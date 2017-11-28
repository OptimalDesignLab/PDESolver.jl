# Dissipation based homotopy predictor-corrector globalization for solving
# steady problems using Newton's method
# based on Brown and Zingg, "A monolithic Homotopy Continuation Algorithm 
# with applications to Computational Fluid Dynamics"
# Journal of Computational Physics 321, (2016), 55-75
# specifically, Algorithm 2

#global evalPhysicsResidual  # evalResidual
#global evalHomotopyResidual # evalHomotopy
"""
  This function solves steady problems using a dissipation-based
  predictor-correcor homotopy globalization for Newtons method.

  Inputs:
    physics_func: the function to solve, ie. func(q) = 0  mathematically.
                  func must have the signature func(mesh, sbp, eqn, opts)

    g_func: the function that evalutes G(q), the dissipation
    sbp: an SBP operator
    eqn: a AbstractSolutionData.  On entry, eqn.q_vec must be the
         initial condition.  On exit, eqn.q_vec will be the solution to
         func(q) = 0
    opts: options dictionary

  Keyword Arguments:
    pmesh: mesh used for calculating preconditioner


  This function uses Newtons method internally, and supports all the
  different jacobian types and jacobian calculation methods that 
  Newton does.

  On entry, eqn.q_vec should contain the initial guess for q.  On exit
  it will contain the solution for func(q) = 0.

  This function is reentrant.
"""
function predictorCorrectorHomotopy{Tsol, Tres, Tmsh}(physics_func::Function,
                                    g_func::Function,
                                    mesh::AbstractMesh{Tmsh}, 
                                    sbp::AbstractSBP, 
                                    eqn::AbstractSolutionData{Tsol, Tres}, 
                                    opts; pmesh=mesh)

#  global evalPhysicsResidual = physics_func
#  global evalHomotopyResidual = g_func
  #----------------------------------------------------------------------------
  # define the homotopy function H and dH/dLambda
  # defines these as nested functions so predictorCorrectorHomotopy is
  # re-entrant
  res_homotopy = zeros(eqn.res)  # used my homotopyPhysics
  """
    This function makes it appear as though the combined homotopy function H
    is a physics.  This works because an elementwise combinations of physics
    is still a physics.

    physics_func is used for the physcs function R and g_func is used for 
    the homotopy function G.  The combined homotopy function is

      (1 - lambda)R(q) + lambda*G(q)

    Inputs: 
      mesh
      sbp
      eqn
      opts
      t
  """
  function homotopyPhysics(mesh, sbp, eqn, opts, t)

    # this function is only for use with Newton's method, where parallel
    # communication is started outside the physics
    # q_vec -> q

#    res_homotopy = zeros(eqn.res)
    fill!(eqn.res, 0.0)
    fill!(res_homotopy, 0.0)


    # calculate physics residual
    # call this function before g_func, to receive parallel communication
    physics_func(mesh, sbp, eqn, opts, t)

    # calculate homotopy function
    g_func(mesh, sbp, eqn, opts, res_homotopy)

    # combine (use lambda from outer function)
    lambda_c = 1 - lambda # complement of lambda
    for i=1:length(eqn.res)
      eqn.res[i] =  lambda_c*eqn.res[i] + lambda*res_homotopy[i]
    end

  #  println("homotopy physics exiting with residual norm ", norm(vec(eqn.res)))
    return nothing
  end


  #----------------------------------------------------------------------------
  # setup

  Tjac = real(Tres)

 

  time = eqn.params.time
  lambda = 1.0  # homotopy parameter
  myrank = mesh.myrank

  # some parameters
  lambda_min = 0.0
  itermax = opts["itermax"]::Int
  res_reltol=opts["res_reltol"]::Float64
  res_abstol=opts["res_abstol"]::Float64


  # counters/loop variables
  iter = 1
  homotopy_tol = 1e-2
  delta_max = 1.0  # step size limit, set to 1 initially,
  psi_max = 10*pi/180  # angle between tangent limiter
  psi = 0.0  # angle between tangent vectors
  tan_norm = 0.0  # current tangent vector norm
  tan_norm_1 = 0.0  # previous tangent vector norm
  res_norm = 0.0  # norm of residual (not homotopy)
  res_norm_0 = 0.0  # residual norm of initial guess
  h = 0.05  # step size
  lambda -= h  # the jacobian is ill-conditioned at lambda=1, so skip it
  # log file
  @mpi_master fconv = BufferedIO("convergence.dat", "a+")

  # needed arrays
  q_vec0 = zeros(eqn.q_vec)
  delta_q = zeros(eqn.q_vec)
  tan_vec = zeros(Tjac, length(eqn.q_vec))  # tangent vector
  tan_vec_1 = zeros(tan_vec)  # previous tangent vector
  dHdLambda_real = zeros(Tjac, length(eqn.q_vec))  


  # stuff for newtonInner
  # because the homotopy function is a pseudo-physics, we can reuse the
  # Newton PC and LO stuff, supplying homotopyPhysics in ctx_residual
  rhs_func = physicsRhs
  ctx_residual = (homotopyPhysics,)
  pc, lo = getNewtonPCandLO(mesh, sbp, eqn, opts)
  ls = StandardLinearSolver(pc, lo, eqn.comm, opts)


  # configure NewtonData
  newton_data, rhs_vec = setupNewton(mesh, pmesh, sbp, eqn, opts, ls)
  newton_data.itermax = 30
 
  # calculate physics residual
  res_norm = real(physicsRhs(mesh, sbp, eqn, opts, eqn.res_vec, (physics_func,)))
  res_norm_0 = res_norm

  # print to log file
  @mpi_master println(fconv, 0, " ", res_norm, " ", 0.0)
  @mpi_master flush(fconv)

  eqn.majorIterationCallback(0, mesh, sbp, eqn, opts, STDOUT)

  #----------------------------------------------------------------------------
  # main loop
  while res_norm > res_norm_0*res_reltol && res_norm > res_abstol && iter < itermax  # predictor loop

    @mpi_master begin
      println(BSTDOUT, "\npredictor iteration ", iter, ", lambda = ", lambda)
      println(BSTDOUT, "res_norm = ", res_norm)
      println(BSTDOUT, "res_norm_0 = ", res_norm_0)
      println(BSTDOUT, "res_norm/res_norm_0 = ", res_norm/res_norm_0)
      println(BSTDOUT, "res_norm = ", res_norm)
      println(BSTDOUT, "h = ", h)
    end

    # calculate homotopy residual
    homotopy_norm = physicsRhs(mesh, sbp, eqn, opts, eqn.res_vec, (homotopyPhysics,))

    homotopy_norm0 = homotopy_norm
    copy!(q_vec0, eqn.q_vec)  # save initial q to calculate delta_q later

    # if we have finished traversing the homotopy path, solve the 
    # homotopy problem (= the physics problem because lambda is zero)
    # tightly
    if abs(lambda - lambda_min) <= eps()
      @mpi_master println(BSTDOUT, "tightening homotopy tolerance")
      homotopy_tol = res_reltol
      reltol = res_reltol*1e-3  # smaller than newton tolerance
      abstol = res_abstol*1e-3  # smaller than newton tolerance
      setTolerances(newton_data.ls, reltol, abstol, -1, -1)

      @mpi_master begin
        println(BSTDOUT, "setting homotopy tolerance to ", homotopy_tol)
        println(BSTDOUT, "ksp reltol = ", reltol)
        println(BSTDOUT, "ksp abstol = ", abstol)
      end
    end

    # do corrector steps
    newtonInner(newton_data, mesh, sbp, eqn, opts, rhs_func, ls, 
                rhs_vec, ctx_residual)

    # compute delta_q
    for i=1:length(eqn.q)
      delta_q[i] = eqn.q_vec[i] - q_vec0[i]
    end
    delta = calcEuclidianNorm(mesh.comm, delta_q)
    if iter == 2
      delta_max = delta
    end

    # predictor step calculation
    if abs(lambda - lambda_min) > eps()
      # calculate dHdLambda at new q value
      calcdHdLambda(mesh, sbp, eqn, opts, lambda, physics_func, g_func, rhs_vec)
      for i=1:length(rhs_vec)
        dHdLambda_real[i] = real(rhs_vec[i])
      end

      # calculate tangent vector dH/dq * t = dH/dLambda
      @mpi_master println(BSTDOUT, "solving for tangent vector")
      calcPCandLO(ls, mesh, sbp, eqn, opts, ctx_residual, 0.0)

#      jac_func(newton_data, mesh, sbp, eqn, opts, jac, ctx_residual)
      linearSolve(ls, dHdLambda_real, tan_vec)
#      matrixSolve(newton_data, eqn, mesh, opts, jac, tan_vec, dHdLambda_real, STDOUT)

      # normalize tangent vector
      tan_norm = calcEuclidianNorm(mesh.comm, tan_vec)
      tan_norm = sqrt(tan_norm*tan_norm + 1)
#      tan_norm = sqrt(dot(tan_vec, tan_vec) + 1)  # + 1 ?
      scale!(tan_vec, 1/tan_norm)

      psi = psi_max
      if iter > 1
        tan_term = dot(tan_vec_1, tan_vec)
        time.t_allreduce += @elapsed tan_term = MPI.Allreduce(tan_term, MPI.SUM, eqn.comm)
        tan_norm_term = (1/tan_norm)*(1/tan_norm_1) 
        arg = tan_term + tan_norm_term
        arg = clamp(arg, -1.0, 1.0)
        psi = acos( arg )
      end

      # save the tangent vector
      copy!(tan_vec_1, tan_vec)
      tan_norm_1 = tan_norm

      # calculate step size
      fac = max(psi/psi_max, delta/delta_max)
      h /= fac

      # take predictor step
      scale!(tan_vec, h)
      for i=1:length(eqn.q_vec)
        eqn.q_vec[i] += tan_vec[i]
      end

      lambda = max(lambda_min, lambda - h)
    end  # end if lambda too large

    # calculate physics residual at new state q
    res_norm = real(physicsRhs(mesh, sbp, eqn, opts, eqn.res_vec, (physics_func,),))

    # print to log file
    @mpi_master println(fconv, iter, " ", res_norm, " ", h )
    @mpi_master flush(fconv)

    eqn.majorIterationCallback(iter, mesh, sbp, eqn, opts, STDOUT)

    iter += 1
  end  # end while loop

  print(BSTDOUT, "\n")

  # inform user of final status
  @mpi_master if iter >= itermax
    println(BSTDERR, "Warning: predictor-corrector did not converge in $iter iterations")
  
  elseif res_norm <= res_abstol
    println(BSTDOUT, "predictor-corrector converged with absolute residual norm $res_norm")
  elseif res_norm/res_norm_0 <= res_reltol
    tmp = res_norm/res_norm_0
    println(BSTDOUT, "predictor-corrector converged with relative residual norm $tmp")
  end

  free(newton_data)
#  cleanupNewton(newton_data, mesh, mesh, sbp, eqn, opts)

  flush(BSTDOUT)
  flush(BSTDERR)

  return nothing
end



"""
  This function calculates dH/dLambda, where H is the homotopy function
  calculated by homotopyPhysics.  The differentiation is done analytically

  Inputs:
    mesh
    sbp
    eqn: eqn.res and eqn.res_vec are overwritten
    opts
    lambda: homotopy parameter lambda
    physics_func: function that evalutes the physics residual
    g_func: function that evalutes g(q)

  Inputs/Outputs
    res_vec: vector to store dH/dLambda in

  Aliasing restrictions: res_vec and eqn.res_vec may not alias
"""
function calcdHdLambda(mesh, sbp, eqn, opts, lambda, physics_func, g_func, res_vec)

  # it appears this only gets called after parallel communication is done
  # so no need to start communication here

#  lambda = eqn.params.homotopy_lambda
  res_homotopy = zeros(eqn.res)


  # calculate physics residual
  physics_func(mesh, sbp, eqn, opts)
  assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)


  # calculate homotopy function
  g_func(mesh, sbp, eqn, opts, res_homotopy)
  assembleSolution(mesh, sbp, eqn, opts, res_homotopy, res_vec)

  # combine them
  for i=1:length(res_vec)
    res_vec[i] -= eqn.res_vec[i]
  end

  return nothing
end



