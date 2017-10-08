# Startup function for 1 dof advection equation
"""
  This function invokes the solver for the advection equation, using the
  specified input file

  Inputs:
    input_file: a string containing the path to an input file, or just a file
                name.  If it is just a file name, it is taken to be in the
                users pwd.

  Outputs:
    mesh: an AbstractMesh object
    sbp: the SBP operator used in solving the equation
    eqn: the AbstractSolutionData object used to solve the equation
    opts: the options dictonary

"""
function run_advection(input_file::AbstractString)

  mesh, sbp, eqn, opts, pmesh = createObjects(input_file)
  solve_advection(mesh, sbp, eqn, opts, pmesh)

  return mesh, sbp, eqn, opts
end


"""
  This function creates and initializes the mesh, sbp, eqn, and opts objects

  Inputs:
    file_name: input file name

  Outputs:
    mesh: an AbstractMesh.  The concrete type is determined by the options
          dictionary
    sbp: an AbstractSBP.  The concrete type is determined by the options
         dictionary
    eqn: an EulerData object
    opts: the options dictionary
    pmesh: mesh used for preconditioning, can be same object as mesh
"""
function createObjects(input_file::AbstractString)

  if !MPI.Initialized()
    MPI.Init()
  end

  #function runtest(flag::Int)
  opts = read_input(input_file)  # read input file and gets default values
  checkOptions(opts)  # physics specific options checking
  # timestepping parameters
  dim = opts["dimensions"]
  # flag determines whether to calculate u, dR/du, or dR/dx (1, 2, or 3)

  dofpernode = 1

  sbp, mesh, pmesh, Tsol, Tres, Tmsh, mesh_time = createMeshAndOperator(opts, dofpernode)

  # Create advection equation object
  Tdim = dim
  eqn = AdvectionData_{Tsol, Tres, Tdim, Tmsh}(mesh, sbp, opts)
  eqn.params.time.t_meshinit = mesh_time  # save the mesh init time

  # Initialize the advection equation
  init(mesh, sbp, eqn, opts)

  return mesh, sbp, eqn, opts, pmesh
end


"""
  Given fully initialized mesh, sbp, eqn, opts, this function solves
  the advection equations.  The 4 object should be obtained from 
  createObjects().
  

  Specifically, it applies an initial condition and invokes a nonlinear
  solver according to the options dictionary.

  Inputs:
    mesh: an AbstractMesh
    sbp: an AbstractSBP
    eqn: an AdvectionData
    opts: the options dictionary.  This must be the options dictionary returned
          by createObjects().  Changing values in the options dictionary after
          calling createObjects() results in undefined behavior.
    pmesh: mesh used for preconditioning, can be same object as mesh.
           default value of mesh

  Options Keys:
    calc_error
    calc_trunc_error
    calc_havg
    perturb_ic
    calc_dt
    finalize_mpi
"""
function solve_advection(mesh::AbstractMesh, sbp, eqn::AdvectionData, opts, pmesh=mesh)

  myrank = mesh.myrank

  fill!(eqn.res, 0.0)
  fill!(eqn.res_vec, 0.0)

  q_vec = eqn.q_vec


  # Populate with initial conditions
  @mpi_master println("\nEvaluating initial condition")
  ICfunc_name = opts["IC_name"]
  ICfunc = ICDict[ICfunc_name]
  @mpi_master println("ICfunc = ", ICfunc)
  ICfunc(mesh, sbp, eqn, opts, q_vec) 

  writedlm("solution_ic.dat", eqn.q_vec)

  if opts["calc_error"]
    @mpi_master println("\ncalculating error of file ", opts["calc_error_infname"], 
            " compared to initial condition")
    # read in this processors portion of the solution
    vals = readdlm(get_parallel_fname(opts["calc_error_infname"], mesy.myrank))
    @assert length(vals) == mesh.numDof

    err_vec = vals - eqn.q_vec
    err_local = calcNorm(eqn, err_vec)
    if myrank == 0
      outname = opts["calc_error_outfname"]
      @mpi_master println("printed err = ", err, " to file ", outname)
      f = open(outname, "w")
      println(f, err)
      close(f)
    end
  end

  if opts["calc_trunc_error"]  # calculate truncation error
    @mpi_master println("\nCalculating residual for truncation error")
    tmp = calcResidual(mesh, sbp, eqn, opts, evalResidual)
    if myrank == 0
      f = open("error_trunc.dat", "w")
      println(f, tmp)
      close(f)
    end
  end

  if opts["calc_havg"]
    h_avg = calcMeshH(mesh, sbp, eqn, opts)
    if myrank == 0  # TODO: make this globally accurate
      rmfile("havg.dat")
      f = open("havg.dat", "w")
      println(f, h_avg)
      close(f)
    end
  end



  if opts["perturb_ic"]
    @mpi_master println("\nPerturbing initial condition")
    perturb_mag = opts["perturb_mag"]
    for i=1:mesh.numDof
      q_vec[i] += perturb_mag*rand()
    end
  end

  res_vec_exact = deepcopy(q_vec)

  #rmfile("IC_$myrank.dat")
  #writedlm("IC_$myrank.dat", real(q_vec))
  saveSolutionToMesh(mesh, q_vec)
  writeVisFiles(mesh, "solution_ic")
  global int_advec = 1  # ???

  if opts["calc_dt"]
    alpha_net = sqrt(eqn.params.alpha_x^2 + eqn.params.alpha_y^2)
    opts["delta_t"] = opts["CFL"]*mesh.min_el_size/alpha_net
  end


  MPI.Barrier( mesh.comm)

  call_nlsolver(mesh, sbp, eqn, opts, pmesh)
  postproc(mesh, sbp, eqn, opts)

  MPI.Barrier(mesh.comm)
  if opts["finalize_mpi"]
    MPI.Finalize()
  end

  return mesh, sbp, eqn, opts
end  # end function

"""
  This function does post processing, if requested by the input options.
  Typical post processing includes calculation of errors, norms of important
  quantities, writing of files. etc.

  Inputs:
    mesh
    sbp
    eqn
    opts
"""
function postproc(mesh, sbp, eqn, opts)

  ##### Do postprocessing ######
  myrank = mesh.myrank
  @mpi_master println("\nDoing postprocessing")

  Tsol = eltype(eqn.q)
  Tres = eltype(eqn.res)
  Tmsh = eltype(mesh.dxidx)

  if opts["do_postproc"] && opts["solve"]
    @mpi_master println("final time = ", eqn.t)
    exfname = opts["exact_soln_func"]
    if haskey(ICDict, exfname)
      exfunc = ICDict[exfname]
      q_exact = zeros(Tsol, mesh.numDof)
      exfunc(mesh, sbp, eqn, opts, q_exact)
      q_diff = eqn.q_vec - q_exact
      diff_norm = calcNorm(eqn, q_diff)
      discrete_norm = norm(q_diff/length(q_diff))
      discrete_norm = MPI.Allreduce(discrete_norm*discrete_norm, MPI.SUM, mesh.comm)
      discrete_norm = sqrt(discrete_norm)

      if myrank == 0
        println("solution error norm = ", diff_norm)
        println("solution discrete L2 norm = ", discrete_norm)
      end

      saveSolutionToMesh(mesh, real(q_diff))
      writeVisFiles(mesh, "solution_error")

      sol_norm = calcNorm(eqn, eqn.q_vec)
      exact_norm = calcNorm(eqn, q_exact)
      @mpi_master println("numerical solution norm = ", sol_norm)
      @mpi_master println("exact solution norm = ", exact_norm)

      # calculate the average mesh size
      h_avg = calcMeshH(mesh, sbp, eqn, opts)
      if myrank == 0
#        println("mesh.min_node_distance = ", mesh.min_node_dist)
        # print to file
        outname = opts["calc_error_outfname"]
        @mpi_master println("printed err = ", diff_norm, " to file ", outname)
        f = open(outname, "w")
        println(f, diff_norm, " ", h_avg)
        close(f)
      end
      #=
      #----  Calculate functional on a boundary  -----#

      if opts["calc_functional"]
        num_functionals = opts["num_functionals"]
        for j = 1:num_functionals
          functional = OptimizationData{Tsol}(mesh, sbp, opts)
          evalFunctional(mesh, sbp, eqn, opts, functional, functional_number=j)
        end  # End for j = 1:num_functionals
        # evalFunctional(mesh, sbp, eqn, opts) # Legacy
      end    # End if opts["calc_functional"]


      #-----  Calculate adjoint vector for a functional  -----#

      if opts["calc_adjoint"]
        eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
        if mesh.isDG
          boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
        end

        # TODO: Presently adjoint computation only for 1 functional. Figure out
        # API based on future use.
        j = 1
        key = string("geom_edges_functional", j)
        functional_edges = opts[key]
        functional_number = j
        functional_name = getFunctionalName(opts, j)

        adjoint_vec = zeros(Tsol, mesh.numDof)
        calcAdjoint(mesh, sbp, eqn, opts, functional_name, functional_number, adjoint_vec)


        # Write adjoint vector to file and mesh
        file_object = open("adjoint_vector.dat", "w")
        for iter = 1:length(adjoint_vec)
          println(file_object, real(adjoint_vec[iter]))
        end
        close(file_object)
        saveSolutionToMesh(mesh, real(adjoint_vec))
        writeVisFiles(mesh, "adjoint_field")
      end  # end if opts["calc_adjoint"]
      =#
#=
      outname = opts["calc_error_outfname"]
      f = open(outname, "w")
      println(f, mesh.numEl, " ", diff_norm, " ", discrete_norm)
      close(f)
=#
    end  # end if haskey(ICname)
  end  # end if do_postproc

  return nothing
end


