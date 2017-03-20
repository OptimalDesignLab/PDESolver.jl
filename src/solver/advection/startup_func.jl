# Startup function for 1 dof advection equation

function run_advection(input_file::AbstractString)
    
  if !MPI.Initialized()
    MPI.Init()
  end

  #function runtest(flag::Int)
  opts = read_input(input_file)  # read input file and gets default values
  checkOptions(opts)  # physics specific options checking
  # timestepping parameters
  delta_t = opts["delta_t"]
  t_max = opts["t_max"]
  dim = opts["dimensions"]
  # flag determines whether to calculate u, dR/du, or dR/dx (1, 2, or 3)
  flag = opts["run_type"]

  dofpernode = 1

  sbp, mesh, pmesh, Tsol, Tres, Tmsh, mesh_time = createMeshAndOperator(opts, dofpernode)

  #mesh.numColors = 4
  #mesh.maxColors = 4

  if opts["write_timing"]
    MPI.Barrier(mesh.comm)
    if mesh.myrank == 0
      f = open("timing.dat", "a+")
      println(f, mesh_time)
      close(f)
    end
  end

  myrank = mesh.myrank

  # Create advection equation object
  Tdim = dim
  opts["Tdim"] = Tdim
  eqn = AdvectionData_{Tsol, Tres, Tdim, Tmsh}(mesh, sbp, opts)

  q_vec = eqn.q_vec


  # Initialize the advection equation
  init(mesh, sbp, eqn, opts)

  # Populate with initial conditions
  println("\nEvaluating initial condition")
  ICfunc_name = opts["IC_name"]
  ICfunc = ICDict[ICfunc_name]
  println("ICfunc = ", ICfunc)
  ICfunc(mesh, sbp, eqn, opts, q_vec) 

  writedlm("solution_ic.dat", eqn.q_vec)

  if opts["calc_error"]
    println("\ncalculating error of file ", opts["calc_error_infname"], 
            " compared to initial condition")
    # read in this processors portion of the solution
    vals = readdlm(get_parallel_fname(opts["calc_error_infname"], mesy.myrank))
    @assert length(vals) == mesh.numDof

    err_vec = vals - eqn.q_vec
    err_local = calcNorm(eqn, err_vec)
    if myrank == 0
      outname = opts["calc_error_outfname"]
      println("printed err = ", err, " to file ", outname)
      f = open(outname, "w")
      println(f, err)
      close(f)
    end
  end

  if opts["calc_trunc_error"]  # calculate truncation error
    println("\nCalculating residual for truncation error")
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
    println("\nPerturbing initial condition")
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
  global int_advec = 1

  if opts["calc_dt"]
    alpha_net = sqrt(eqn.params.alpha_x^2 + eqn.params.alpha_y^2)
    opts["delta_t"] = opts["CFL"]*mesh.min_el_size/alpha_net
  end


  MPI.Barrier( mesh.comm)
  # evalResidual(mesh, sbp, eqn, opts, t)

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
  println("\nDoing postprocessing")

  Tsol = eltype(eqn.q)
  Tres = eltype(eqn.res)
  Tmsh = eltype(mesh.dxidx)
  myrank = mesh.myrank

  if opts["do_postproc"] && opts["solve"]
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
        println("printed err = ", diff_norm, " to file ", outname)
        f = open(outname, "w")
        println(f, diff_norm, " ", h_avg)
        close(f)
      end

      #----  Calculate functional on a boundary  -----#
      
      if opts["calc_functional"]
        eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
        if mesh.isDG
          boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
        end
        
        # Calculate functional over edges
        num_functionals = opts["num_functionals"]
        for j = 1:num_functionals
          # Geometric edge at which the functional needs to be integrated
          key_j = string("geom_edges_functional", j)
          functional_edges = opts[key_j]
          functional_name = getFunctionalName(opts, j)

          functional_val = zero(Tsol)
          functional_val = calcBndryfunctional(mesh, sbp, eqn, opts, 
                           functional_name, functional_edges)
          println("\nNumerical functional value on geometric edges ", 
                  functional_edges, " = ", functional_val)
          
          analytical_functional_val = opts["analytical_functional_val"]
          println("analytical_functional_val = ", analytical_functional_val)
          
          absolute_functional_error = norm((functional_val - 
                                           analytical_functional_val), 2)
          relative_functional_error = absolute_functional_error/
                                      norm(analytical_functional_val, 2)
          
          mesh_metric = 1/sqrt(mesh.numEl/2)  # Only for a square domain with triangular elements
          
          # write functional error to file
          outname = string(opts["functional_error_outfname"], j, ".dat")
          println("printed relative functional error = ", 
                  relative_functional_error, " to file ", outname, '\n')
          f = open(outname, "w")
          println(f, relative_functional_error, " ", mesh_metric)
          close(f)
        end  # End for j = 1:num_functional
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


