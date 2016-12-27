# function used by all physics' startup functions to invoke a nonlinear
# solver

"""
  This function takes in the 4 principle object, fully initialized, and calls
  a nonlinear solver on them, according to the options in the dictionary.
  The evalResidual function is passed to the nonlinear solver

  Inputs:
    mesh: a mesh object
    sbp: an SBP operator
    eqn: an equation object
    opts: options dictionary, used to determine which nonlinear solver to call
    pmesh: mesh used for calculating preconditioning jacobian in Newton's
           method, default to using mesh if not specified

  Outputs:
    none

  Aliasing restrictions: none (specificaly, mesh and pmesh *can* be the same
                         object)
"""
function call_nlsolver(mesh::AbstractMesh, sbp::AbstractSBP, 
                       eqn::AbstractSolutionData, opts::Dict, 
                       pmesh::AbstractMesh=mesh)
# from Advection---------------------------------------------------------------
  flag = opts["run_type"]::Int
  if opts["solve"]
    
    solve_time = @elapsed if flag == 1 # normal run
      # RK4 solver
      delta_t = opts["delta_t"]
      t_max = opts["t_max"]
      @time rk4(evalResidual, delta_t, t_max, mesh, sbp, eqn, opts, 
                res_tol=opts["res_abstol"], real_time=opts["real_time"])
      println("finish rk4")
  #    printSolution("rk4_solution.dat", eqn.res_vec)
    
    elseif flag == 2 # forward diff dR/du
    #= 
      # define nested function
      function dRdu_rk4_wrapper(u_vals::AbstractVector, res_vec::AbstractVector)
        eqn.q_vec = u_vals
        eqn.q_vec = res_vec
        rk4(evalResidual, delta_t, t_max, mesh, sbp, eqn)
        return nothing
      end

      # use ForwardDiff package to generate function that calculate jacobian
      calcdRdu! = forwarddiff_jacobian!(dRdu_rk4_wrapper, Float64, 
                  fadtype=:dual, n = mesh.numDof, m = mesh.numDof)

      jac = zeros(Float64, mesh.numDof, mesh.numDof)  # array to be populated
      calcdRdu!(eqn.q_vec, jac)
    =#
    elseif flag == 3 # calculate dRdx

      # dRdx here

    elseif flag == 4 || flag == 5
      @time newton(evalResidual, mesh, sbp, eqn, opts, pmesh, itermax=opts["itermax"], 
                   step_tol=opts["step_tol"], res_abstol=opts["res_abstol"], 
                   res_reltol=opts["res_reltol"], res_reltol0=opts["res_reltol0"])

      printSolution("newton_solution.dat", eqn.res_vec)

    elseif flag == 9
      # to non-pde rk4 run
      function pre_func(mesh, sbp, eqn,  opts)
        println("pre_func was called")
        return nothing
      end

      function post_func(mesh, sbp, eqn, opts)
  #      for i=1:mesh.numDof
  #        eqn.res_vec[i] *= eqn.Minv[i]
  #      end
        nrm = norm(eqn.res_vec)
        println("post_func returning residual norm = ", nrm)
        return nrm
      end
      delta_t = opts["delta_t"]
      t_max = opts["t_max"]

      rk4(evalResidual, delta_t, t_max, eqn.q_vec, eqn.res_vec, pre_func, post_func, (mesh, sbp, eqn), opts, majorIterationCallback=eqn.majorIterationCallback, real_time=opts["real_time"])

    elseif flag == 10
      function test_pre_func(mesh, sbp, eqn, opts)
        
        eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
      end

      function test_post_func(mesh, sbp, eqn, opts, calc_norm=true)
        return calcNorm(eqn, eqn.res_vec)
      end
      delta_t = opts["delta_t"]
      t_max = opts["t_max"]

      rk4(evalResidual, delta_t, t_max, eqn.q_vec, eqn.res_vec, test_pre_func,
          test_post_func, (mesh, sbp, eqn), opts, 
          majorIterationCallback=eqn.majorIterationCallback, real_time=opts["real_time"])

    elseif flag == 20

      @time t = crank_nicolson(evalResidual, opts["delta_t"], opts["t_max"], 
                               mesh, sbp, eqn, opts, opts["res_abstol"], 
                               opts["real_time"])

      eqn.t = t

  #   else
  #     throw(ErrorException("No flag specified: no solve will take place"))
  #     return nothing

    end       # end of if/elseif blocks checking flag

    println("total solution time printed above")
    params = eqn.params
    myrank = mesh.myrank

    if opts["write_timing"]
      MPI.Barrier(mesh.comm)
      if mesh.myrank == 0
        f = open("timing.dat", "a+")
        println(f, solve_time)
        close(f)
      end
    end

    # evaluate residual at final q value
    need_res = false
    if need_res
      eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
      evalResidual(mesh, sbp, eqn, opts, eqn.t)

      eqn.res_vec[:] = 0.0
      eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
    end

    if opts["write_finalsolution"]
      println("writing final solution")
      writedlm("solution_final_$myrank.dat", real(eqn.q_vec))
    end

    if opts["write_finalresidual"]
      writedlm("residual_final_$myrank.dat", real(eqn.res_vec))
    end

    myrank = mesh.myrank
  #  f = open("profile_$myrank.dat", "a+")
  #  Profile.print(f, format=:flat, C=true)
  #  close(f)

    saveSolutionToMesh(mesh, real(eqn.q_vec))
  #  printSolution(mesh, real(eqn.q_vec))
  #  printCoordinates(mesh)
    writeVisFiles(mesh, "solution_done")

    # write timings
  #  timings = [params.t_volume, params.t_face, params.t_source, params.t_sharedface, params.t_bndry, params.t_send, params.t_wait, params.t_allreduce, params.t_jacobian, params.t_solve, params.t_barrier, params.t_barrier2, params.t_barrier3]
  #  writedlm("timing_breakdown_$myrank.dat", vcat(timings, params.t_barriers))
    fname = "timing_breakdown_$myrank"
    write_timings(params.time, fname)

    MPI.Barrier(mesh.comm)
    if opts["finalize_mpi"]
      MPI.Finalize()
    end
  end  # end if (opts[solve])

  return nothing
end
