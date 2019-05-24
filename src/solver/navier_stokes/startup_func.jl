# functions required to start a solve

import PDESolver.solvePDE


"""
  This function creates and initializes the mesh, sbp, eqn, and opts objects

  Inputs:
    opts: options dictionary

  Outputs:
    mesh: an AbstractMesh.  The concrete type is determined by the options
          dictionary
    sbp: an AbstractOperator.  The concrete type is determined by the options
         dictionary
    eqn: an NSData object
    opts: the options dictionary
    pmesh: mesh used for preconditioning, can be same object as mesh
"""
function createObjects(opts::Dict)

  # make sure Euler options were set
  EulerEquationMod.checkOptions(opts)

  Tdim = opts["dimensions"]
  dofpernode = Tdim + 2

  sbp, mesh, pmesh, Tsol, Tres, Tmsh, mesh_time = createMeshAndOperator(opts, dofpernode)

  # create euler equation
  eqn = NSData_{Tsol, Tres, Tdim, Tmsh}(mesh, sbp, opts)

  # initialize physics module and populate any fields in mesh and eqn that
  # depend on the physics module
  init(mesh, sbp, eqn, opts, pmesh)

  return mesh, sbp, eqn, opts, pmesh
end


"""
  Given fully initialized mesh, sbp, eqn, opts, this function solves
  the Euler equations.  The 4 object should be obtained from createObjects().


  Specifically, it applies an initial condition and invokes a nonlinear
  solver according to the options dictionary.

  Inputs:
    mesh: an AbstractMesh
    sbp: an AbstractOperator
    eqn: an NSData
    opts: the options dictionary.  This must be the options dictionary returned
          by createObjects().  Changing values in the options dictionary after
          calling createObjects() results in undefined behavior.
    pmesh: mesh used for preconditioning, can be same object as mesh.
           default value of mesh

"""
function solvePDE(mesh::AbstractMesh, sbp::AbstractOperator, eqn::NSData, opts::Dict, pmesh::AbstractMesh=mesh)
  #delta_t = opts["delta_t"]   # delta_t: timestep for RK

  myrank = mesh.myrank

  # zero some arrays out in case the user forgot

  fill!(eqn.res, 0.0)
  fill!(eqn.res_vec, 0.0)
  res_vec = eqn.res_vec
  q_vec = eqn.q_vec       # solution at current timestep

  #TODO: move as many of these if statements into functions as possible

  # calculate residual of some other function for res_reltol0
  # TODO: add a boolean options here?
  Relfunc_name = opts["Relfunc_name"]
  if haskey(ICDict, Relfunc_name)
    @mpi_master println("\ncalculating residual for relative residual tolerance")
    Relfunc = ICDict[Relfunc_name]
    @mpi_master println("Relfunc = ", Relfunc)
    Relfunc(mesh, sbp, eqn, opts, q_vec)
    tmp = physicsRhs(mesh, sbp, eqn, opts, eqn.res_vec, (evalResidual,))
    res_real = real(eqn.res_vec)
    opts["res_reltol0"] = tmp
    @mpi_master println("res_reltol0 = ", tmp)

    saveSolutionToMesh(mesh, res_real)
    writeVisFiles(mesh, "solution_relfunc")
  end

  # populate u0 with initial condition
  @mpi_master println("\nEvaluating initial condition")
  ICfunc_name = opts["IC_name"]
  ICfunc = ICDict[ICfunc_name]
  @mpi_master println("ICfunc = ", ICfunc)
  ICfunc(mesh, sbp, eqn, opts, q_vec)

  if opts["calc_error"]
    @mpi_master println("\ncalculating error of file ",
                       opts["calc_error_infname"],
                      " compared to initial condition")

    # read in this processors portion of the solution
    vals = copy(eqn.q_vec)
    readSolutionFiles(mesh, sbp, eqn, opts, opts["calc_error_infname"], vals)

    err_vec = abs.(vals - eqn.q_vec)
    err = calcNorm(eqn, err_vec)

    # calculate avg mesh size
    h_avg = calcMeshH(mesh, sbp, eqn, opts)

    @mpi_master begin
      outname = opts["calc_error_outfname"]
      @mpi_master println("printint err = ", err, " to file ", outname)
      f = open(outname, "w")
      println(f, err, " ", h_avg)
      close(f)
    end

    # write visualization
    saveSolutionToMesh(mesh, vec(err_vec))
    writeVisFiles(mesh, "solution_error")
  end

  if opts["calc_trunc_error"]  # calculate truncation error
    @mpi_master println("\nCalculating residual for truncation error")
    tmp = physicsRhs(mesh, sbp, eqn, opts, eqn.res_vec, (evalResidual,))

    @mpi_master begin
      f = open("error_trunc.dat", "w")
      println(f, tmp)
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

  #writeSolutionFiles(mesh, sbp, eqn, opts, "IC")
  saveSolutionToMesh(mesh, q_vec)

  writeVisFiles(mesh, "solution_ic")
  if opts["calc_dt"]
    wave_speed = EulerEquationMod.calcMaxWaveSpeed(mesh, sbp, eqn, opts)
    @mpi_master println("max wave speed = ", wave_speed)
    @mpi_master println("min element size = ", mesh.min_el_size)
    delta_t = opts["CFL"]*mesh.min_el_size/wave_speed
    @mpi_master println("for a CFL of ", opts["CFL"], " delta_t = ", delta_t)
    opts["delta_t"] = delta_t
  end

  # this must be last, because all previous calculations need to use the
  # original initial condition
  if opts["is_restart"]
    loadRestartState(mesh, sbp, eqn, opts, pmesh)
  end

  call_nlsolver(mesh, sbp, eqn, opts, pmesh)
  postproc(mesh, sbp, eqn, opts)

  EulerEquationMod.cleanup(eqn.euler_eqn, opts)

  MPI.Barrier(mesh.comm)

  return mesh, sbp, eqn, opts
end  # end function


"""
  This function does post processing, if requested by the input options.
  Typical post processing includes calculation of errors, norms of important
  quantities, writing of files. etc.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
"""
function postproc(mesh, sbp, eqn, opts)

  ##### Do postprocessing ######
  myrank = mesh.myrank
  @mpi_master println("\nDoing postprocessing")


  if opts["do_postproc"] && opts["solve"]
    exfname = opts["exact_soln_func"]
    if haskey(ICDict, exfname)
      exfunc = ICDict[exfname]
      q_exact = zeros(eqn.q_vec)
      exfunc(mesh, sbp, eqn, opts, q_exact)
      
      myrank = mesh.myrank
      q_diff = eqn.q_vec - q_exact
      saveSolutionToMesh(mesh, abs.(real(q_diff)))
      writeVisFiles(mesh, "solution_error")


      diff_norm = calcNorm(eqn, q_diff)
      @mpi_master println("solution error norm = ", diff_norm)
      # TODO: make this mesh.min_el_size?
      h_avg = calcMeshH(mesh, sbp, eqn, opts)

      # print to file
      @mpi_master begin
        outname = opts["calc_error_outfname"]
        f = open(outname, "w")
        println(f, diff_norm, " ", h_avg)
        close(f)
      end
    end  # end if haskey(ICname)
  end  # end if do_postproc

  return nothing
end
