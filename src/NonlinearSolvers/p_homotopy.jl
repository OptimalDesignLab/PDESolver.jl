# do homotopy over the operator degree p

mutable struct pHomotopyData
  solve_homotopy::Bool
  p_min::Int
  p_max::Int
  directory_names::Vector{String}
  directory_idx::Int
  directory_orig::String  # starting directory
  comm::MPI.Comm
  myrank::Int
  commsize::Int


  function pHomotopyData(opts, comm=MPI.Comm)

    myrank = MPI.Comm_rank(comm)
    commsize = MPI.Comm_size(comm)
    solve_homotopy = opts["phomotopy_solve_homotopy"]
    p_min = opts["order"]
    p_max = opts["phomotopy_p_max"]
    directory_names = makeDirectories(p_min, p_max, myrank)
    directory_idx = 1
    directory_orig = pwd()

    @assert length(opts["phomotopy_euler_taus"]) == p_max

    MPI.Barrier(comm)  # make sure all filesystem operations have finished

    return new(solve_homotopy, p_min, p_max, directory_names, directory_idx,
               directory_orig, comm, myrank, commsize)
  end
end


"""
  Creates the directories for the "easy" solves.

  **Inputs**

   * p_min: degree of the first operator to solve
   * p_max: degree of the last operator to solve
   * myrank: the MPI rank of this process
"""
function makeDirectories(p_min::Integer, p_max::Integer, myrank::Int)

  names = Vector{String}(p_max)
  # make directories
  for i=p_min:p_max
    dname = "p$(i)_easy"
    dname_full = joinpath(pwd(), dname)

    if myrank == 0  # only have root process create directories
      if isdir(dname_full)
        println(BSTDOUT, "removing directory ", dname); flush(BSTDOUT)
        rm(dname_full, recursive=true)
      end
      mkdir(dname_full)
    end

    names[i] = dname_full
  end

  return names
end

function nextDirectory(pdata::pHomotopyData)

  dname = pdata.directory_names[pdata.directory_idx]
  pdata.directory_idx += 1

  return dname
end


"""
  Update opts for doing the first solve with implicit Euler
"""
function updateOptsFirst(opts)

  degree = opts["order"]
  opts["Flux_name"] = opts["phomotopy_flux_regular"]
  opts["newton_globalize_euler"] = true
  opts["setup_globalize_euler"] = true
  opts["euler_tau"] = opts["phomotopy_euler_taus"][degree]

  return nothing
end



"""
  Update opts for next mesh while solving on p=2:p_max
"""
function updateOptsIntermediate(opts)

  degree = opts["order"]
  opts["Flux_name"] = opts["phomotopy_flux_regular"]
  opts["newton_globalize_euler"] = true
  opts["setup_globalize_euler"] = true
  opts["euler_tau"] = opts["phomotopy_euler_taus"][degree + 1]

  return nothing
end

"""
  Update opts and the eqn object for the final solve
"""
function updateOptsFinal(mesh, sbp, eqn, opts)

  opts["Flux_name"] = opts["phomotopy_flux_regular"]
  opts["shock_sensor_name"] = opts["phomotopy_shock_sensor_hard"]
  opts["newton_globalize_euler"] = true
  opts["setup_globalize_euler"] = true
  opts["euler_tau"] = opts["phomotopy_euler_tau_final"]

  # update the eqn object with the flux function and shock sensor in the
  # updated options dictionary
  sensor = createShockSensor(mesh, sbp, eqn, opts)
  setShockSensor(eqn, sensor)
  setFluxFunction(mesh, sbp, eqn, opts)

  # open new log files is current directory
  closeLoggingFiles(eqn, opts); closeLoggingFiles(eqn, opts)

  return nothing
end


"""
  Function to create objects of one degree higher and write solution files
  of both the current solution and the solution interpolated to the new
  objects.
"""
function createNewObjects(pdata::pHomotopyData, mesh, sbp, eqn, opts)

  myrank = pdata.myrank
  updateOptsIntermediate(opts)

  writeSolutionFiles(mesh, sbp, eqn, opts, "solution_p$(sbp.degree)_easy.dat")

  cd(nextDirectory(pdata))
  @mpi_master println(BSTDOUT, "\nworking directory is now ", pwd())

  mesh_h, sbp_h, eqn_h, opts_h = createEnrichedObjects(mesh, sbp, eqn, opts)
  writeSolutionFiles(mesh_h, sbp_h, eqn_h, opts_h, "solution_p$(sbp_h.degree)_interp.dat")

  return mesh_h, sbp_h, eqn_h, opts_h
end

"""
  This function solves a given problem in increasing degree operators,
  using the previous solution as the initial guess for the next degree
  operator.  The problem will be solved twice with the final degree operator,
  once with the initial shock sensor, and once with the final shock sensor.
  This can speed convergence for shock problems by using an "easy" to solve
  shock sensor for all except the final solve.
  
  This function has two modes of operation:

   1. The standard way this function operates is to solve the first
      problem (typically p = 1) with the predictor-corrector homotopy, and then
      solve all remaining problems (p > 1) with implicit Euler, 

   2. It is also possible to use this function starting from p != 1 by
      providing an initial condition that is a good initial guess for implicit
      Euler.  For this mode the option `phomotopy_solve_homotopy` should be
      set to false, so implicit Euler will be used for the first solve rather
      than homotopy.

  The input file should be set for the conditions that will be used for the
  homotopy solve.  There are options below for changing certain options for
  the subsequent solves. 

  **Options Keys**

  Some of the options listed below are general options that take on specific
  meaning for this function.  Option names with the prefix `phomotopy` are
  specific to this function.

   * phomotopy_solve_homotopy: if false, do mode 1 described above, if true,
                               do mode 2.
   * order: the degree of the operator used for the *first* solve.
   * phomotopy_p_max: the maximum degree operator to solve
   * Flux_name: the flux function used by DG methods for the *first* solve.
   * phomotopy_flux_regular: the flux function that will be used for all
                             except the first solve.  This enables using
                             a different flux function for p=1, which can
                             be advantageous.
   * shock_sensor_name: the shock sensor that will be used for all except the
                        final solve
   * phomotopy_shock_sensor_hard: the shock sensor that will be used for the
                                  final solve
   * phomotopy_euler_taus: the initial time-step for implicit Euler for
                           each degree operator (length `phomotopy_p_max`).
                           (See below)
   * phomotopy_euler_tau_final: the implicit Euler timestep for the final solve.
                                (`phomotopy_euler_taus[end]` is used for the 
                                solve with p=`pmax` with the original shock
                                sensor, `phomotopy_euler_tau_final` is used
                                with `phomotopy_shock_sensor hard`.

  For mode 2, some of the keys are used differently or ignored.  Specifically:

   * `phomotopy_flux_regular` is used for all solves, including the first
   * `phomotopy_euler_taus[1]` is used for the first solve (in mode 1 the first
                               solve uses homotopy, so no implicit Euler
                               timestep is required.

  **Additional Details**

   All solves except the final one are done in newly-created subdirectories
   "p\$i_easy", where `i` is the operator degree.  The final solve is done in
   the pwd of the solver when this function is entered.  Existing directories
   with the same names will be deleted.  Only directories between opts["order"]
   and opts["phomotopy_p_max"] will be deleted.

   In each directory, a solution file will be saved with the converged solution.
   for p > opts["order"], there will also be a file containing the solution
   used as the initial guess (the solution interpolated from the previous
   solve).  These files, combined with mode 2, provide a way to resume a
   p homotopy sequence if one of the solves failed.
"""
function pHomotopy(mesh::AbstractMesh, sbp::AbstractSBP,
                   eqn::AbstractSolutionData, opts)

# the eqn object should initially be configured with the easy flux function
# and the easy shock sensor

  pdata = pHomotopyData(opts, eqn.comm)
  myrank = pdata.myrank

  cd(nextDirectory(pdata))
  @mpi_master println(BSTDOUT, "\nworking directory is now ", pwd())
  # open new log files is current directory
  closeLoggingFiles(eqn, opts); closeLoggingFiles(eqn, opts)

  if pdata.solve_homotopy
    # solve p=1 with homotopy
    predictorCorrectorHomotopy(evalResidual, evalHomotopy, mesh, sbp, eqn, opts)
  else
    # solve with implicit Euler
    updateOptsFirst(opts)
    setFluxFunction(mesh, sbp, eqn, opts)
    newton(evalResidual, mesh, sbp, eqn, opts)
  end

  mesh_old = mesh; sbp_old = sbp; eqn_old = eqn; opts_old = opts

  # solve with increasing degree p operators
  for p=(pdata.p_min + 1):pdata.p_max

    mesh, sbp, eqn, opts = createNewObjects(pdata, mesh_old, sbp_old, eqn_old,
                                            opts_old)
    @mpi_master println(BSTDOUT, "about to solve p = $p")
    flush(BSTDOUT)
    newton(evalResidual, mesh, sbp, eqn, opts)
    mesh_old = mesh; sbp_old = sbp; eqn_old = eqn; opts_old = opts
  end

  # solve on the final degree operator, with the final shock sensor

  # the final solve is a little different because we already have objects of
  # the right degree, we just need to update them to the final configuration
  cd(pdata.directory_orig)
  writeSolutionFiles(mesh, sbp, eqn, opts, "solution_p$(sbp.degree)_easy.dat")
  updateOptsFinal(mesh, sbp, eqn, opts)

  newton(evalResidual, mesh, sbp, eqn, opts)

  return mesh, sbp, eqn, opts
end
