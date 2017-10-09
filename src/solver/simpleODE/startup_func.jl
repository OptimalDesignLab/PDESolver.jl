# startup file for simple ODE test problem:
#   q = t^2 + x^2

"""
  This function invokes the solver for the simpleODE equation, using the
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
function run_simpleode(input_file::AbstractString)

  mesh, sbp, eqn, opts, pmesh = createObjects(input_file)
  solve_simpleODE(mesh, sbp, eqn, opts, pmesh)

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

  opts = read_input(input_file)

  dim = opts["dimensions"]
  dofpernode = 1

  sbp, mesh, pmesh, Tsol, Tres, Tmsh, mesh_time = createMeshAndOperator(opts, dofpernode)

  println("\ntypeof(mesh) = ", typeof(mesh))
  println("is subtype of DG mesh = ", typeof(mesh) <: AbstractDGMesh)
  println("mesh.isDG = ", mesh.isDG)

  # Create a simpleODE equation object
  Tdim = dim
  eqn = SimpleODEData_{Tsol, Tres, Tdim, Tmsh}(mesh, sbp, opts)
  eqn.params.time.t_meshinit = mesh_time


  init(mesh, sbp, eqn, opts)

  return mesh, sbp, eqn, opts, pmesh
end


"""
  Given fully initialized mesh, sbp, eqn, opts, this function solves
  the simpleODE equation.  The 4 object should be obtained from 
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
function solve_simpleODE(mesh::AbstractMesh, sbp, eqn::SimpleODEData, opts, pmesh=mesh)

  # timestepping parameters
  delta_t = opts["delta_t"]
  t_max = opts["t_max"]
  # flag determines whether to calculate u, dR/du, or dR/dx (1, 2, or 3)
  flag = opts["run_type"]

  myrank = mesh.myrank

  fill!(eqn.res, 0.0)
  fill!(eqn.res_vec, 0.0)

  q_vec = eqn.q_vec
  # Apply IC's
  println("\nEvaluating initial condition")
  ICfunc_name = opts["IC_name"]
  ICfunc = ICDict[ICfunc_name]
  println("ICfunc = ", ICfunc)
  ICfunc(mesh, sbp, eqn, opts, q_vec) 
  println("finished initializing q")

  saveSolutionToMesh(mesh, real(eqn.q_vec))
  writeVisFiles(mesh, "solution_ic")

  writedlm("solution_ic.dat", real(eqn.q_vec))
  writedlm("residual_ic_$myrank.dat", real(eqn.res_vec))

  # because simpleODE is different in how it calls rk4, it does not use
  # call_nlsolver
#  call_nlsolver(mesh, sbp, eqn, opts, pmesh)

  if opts["solve"]

    # handle flags

    if flag == 1        # RK4

      # NOTE: needs to have a custom pde_pre_func & pde_post_func suited to ODEs, named ode_pre_func & ode_post_func

      rk4(evalResidual, opts["delta_t"], t_max, eqn.q_vec, eqn.res_vec, ode_pre_func, ode_post_func, 
          (mesh, sbp, eqn), opts; 
          majorIterationCallback=eqn.majorIterationCallback, res_tol=opts["res_abstol"], real_time=opts["real_time"])

    elseif flag == 20   # Crank-Nicolson
      @time crank_nicolson(evalResidual, opts["delta_t"], t_max, mesh, sbp, eqn, 
                           opts, opts["res_abstol"], opts["real_time"])

    else

      throw(ErrorException("only RK4 and CN are implemented now; user selected different solve flag"))
      return nothing

    end

    eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    evalResidual(mesh, sbp, eqn, opts, eqn.params.t)

    eqn.res_vec[:] = 0.0
    eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)

    # printing solution and residual to file
    writedlm("solution_final.dat", real(eqn.q_vec))
    writedlm("residual_final_$myrank.dat", real(eqn.res_vec))
   

    saveSolutionToMesh(mesh, real(eqn.q_vec))
    writeVisFiles(mesh, "solution_done")

    postproc(mesh, sbp, eqn, opts)
    # TODO: comparison with exact solution

  end   # end of if opts["solve"]

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
  if opts["do_postproc"] && opts["solve"]
    exfname = opts["exact_soln_func"]
    if haskey(ICDict, exfname)
      exfunc = ICDict[exfname]
      q_exact = zeros(Tsol, mesh.numDof)
      exfunc(mesh, sbp, eqn, opts, q_exact)

    end
  end     # end of if opts["do_postproc"]

  return nothing
end


