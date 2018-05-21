import PDESolver.solvePDE

#=
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
=#

"""
  This function creates and initializes the mesh, sbp, eqn, and opts objects

  **Inputs**

   * opts: options dictionary

  **Outputs**

   * mesh: an AbstractMesh.  The concrete type is determined by the options
          dictionary
   * sbp: an AbstractSBP.  The concrete type is determined by the options
         dictionary
   * eqn: an AdvectionData object
   * opts: the options dictionary
   * pmesh: mesh used for preconditioning, can be same object as mesh
"""
function createObjects(opts::Dict)

  dim = opts["dimensions"]
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

"""
function solvePDE(mesh::AbstractMesh, sbp::AbstractSBP, eqn::AdvectionData,
                  opts::Dict, pmesh::AbstractMesh=mesh)

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
    tmp = physicsRhs(mesh, sbp, eqn, opts, eqn.res_vec, (evalResidual,))
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

  # this must be last, because all previous calculations need to use the
  # original initial condition
  if opts["is_restart"]
    loadRestartState(mesh, sbp, eqn, opts, pmesh)
  end


  MPI.Barrier( mesh.comm)

  call_nlsolver(mesh, sbp, eqn, opts, pmesh)
  postproc(mesh, sbp, eqn, opts)

  MPI.Barrier(mesh.comm)

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
      outname = opts["calc_error_outfname"]
      f = open(outname, "w")
      println(f, mesh.numEl, " ", diff_norm, " ", discrete_norm)
      close(f)
=#
    end  # end if haskey(ICname)
  end  # end if do_postproc

  return nothing
end
