# startup file for simple ODE test problem:
#   q = t^2 + x^2

function run_simpleode(input_file::AbstractString)

  if !MPI.Initialized()
    MPI.Init()
  end

  opts = read_input(input_file)

  # timestepping parameters
  delta_t = opts["delta_t"]
  t_max = opts["t_max"]
  dim = opts["dimensions"]
  # flag determines whether to calculate u, dR/du, or dR/dx (1, 2, or 3)
  flag = opts["run_type"]

  dofpernode = 1

  sbp, mesh, pmesh, Tsol, Tres, Tmsh, mesh_time = createMeshAndOperator(opts, dofpernode)

  if opts["write_timing"]
    MPI.Barrier(mesh.comm)
    if mesh.myrank == 0
      f = open("timing.dat", "a+")
      println(f, mesh_time)
      close(f)
    end
  end

  myrank = mesh.myrank

  println("\ntypeof(mesh) = ", typeof(mesh))
  println("is subtype of DG mesh = ", typeof(mesh) <: AbstractDGMesh)
  println("mesh.isDG = ", mesh.isDG)

  # Create a simpleODE equation object
  Tdim = dim
  eqn = SimpleODEData_{Tsol, Tres, Tdim, Tmsh}(mesh, sbp, opts)

  q_vec = eqn.q_vec

  init(mesh, sbp, eqn, opts)

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

  if opts["solve"]

    # handle flags

    if flag == 1        # RK4

      # NOTE: needs to have a custom pde_pre_func & pde_post_func suited to ODEs, named ode_pre_func & ode_post_func

      rk4(evalSimpleODE, opts["delta_t"], t_max, eqn.q_vec, eqn.res_vec, ode_pre_func, ode_post_func, 
          (mesh, sbp, eqn), opts; 
          majorIterationCallback=eqn.majorIterationCallback, res_tol=opts["res_abstol"], real_time=opts["real_time"])

    elseif flag == 20   # Crank-Nicolson
      @time crank_nicolson(evalSimpleODE, opts["delta_t"], t_max, mesh, sbp, eqn, 
                           opts, opts["res_abstol"], opts["real_time"])

    else

      throw(ErrorException("only RK4 and CN are implemented now; user selected different solve flag"))
      return nothing

    end

    eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    evalSimpleODE(mesh, sbp, eqn, opts, eqn.params.t)

    eqn.res_vec[:] = 0.0
    eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)

    # printing solution and residual to file
    writedlm("solution_final.dat", real(eqn.q_vec))
    writedlm("residual_final_$myrank.dat", real(eqn.res_vec))
   

    saveSolutionToMesh(mesh, real(eqn.q_vec))
    writeVisFiles(mesh, "solution_done")

    # TODO: comparison with exact solution

    if opts["do_postproc"]
      exfname = opts["exact_soln_func"]
      if haskey(ICDict, exfname)
        exfunc = ICDict[exfname]
        q_exact = zeros(Tsol, mesh.numDof)
        exfunc(mesh, sbp, eqn, opts, q_exact)

      end
    end     # end of if opts["do_postproc"]

  end   # end of if opts["solve"]

  return mesh, sbp, eqn, opts
end  # end function
