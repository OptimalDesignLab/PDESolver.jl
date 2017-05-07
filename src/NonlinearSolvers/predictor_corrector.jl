# Dissipation based homotopy predictor-corrector globalization for solving
# steady problems using Newton's method
# based on Brown and Zingg, "A monolithic Homotopy Continuation Algorithm 
# with applications to Computational Fluid Dynamics"
# Journal of Computational Physics 321, (2016), 55-75
# specifically, Algorithm 2

global evalPhysicsResidual  # evalResidual
global evalHomotopyResidual # evalHomotopy
"""
  This function solves steady problems using a dissipation-based
  predictor-correcor homotopy globalization for Newtons method.

  Inputs:
    func: the function to solve, ie. func(q) = 0  mathematically.
          func must have the signature func(mesh, sbp, eqn, opts)

    homotopy_func: the function that evalutes G(q), the dissipation
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

"""
function predictorCorrectorHomotopy{Tsol, Tres, Tmsh}(physics_func::Function,
                                    homotopy_func::Function,
                                    mesh::AbstractMesh{Tmsh}, 
                                    sbp::AbstractSBP, 
                                    eqn::AbstractSolutionData{Tsol, Tres}, 
                                    opts; pmesh=mesh)

  global evalPhysicsResidual = physics_func
  global evalHomotopyResidual = homotopy_func

  Tjac = real(Tres)

  newton_data, jac, rhs_vec = setupNewton(mesh, pmesh, sbp, eqn, opts, 
                                          homotopyPhysics)
  

  lambda = 1.0  # homotopy parameter
  eqn.params.homotopy_lambda = lambda

  # some parameters
  lambda_min = 0.0
  itermax = 100
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
  h = 0.01  # step size

  # log file
  fconv = BufferedIO("convergence.dat", "a+")

  # needed arrays
  q_vec0 = zeros(eqn.q_vec)
  delta_q = zeros(eqn.q_vec)
  tan_vec = zeros(Tjac, length(eqn.q_vec))  # tangent vector
  tan_vec_1 = zeros(tan_vec)  # previous tangent vector
  dHdLambda_real = zeros(Tjac, length(eqn.q_vec))  
#  dHdLambda = zeros(eqn.q_vec)


  # stuff for newtonInner
  rhs_func = physicsRhs
  jac_func = physicsJac
  ctx_residual = (homotopyPhysics,)

  # calculate physics residual
  res_norm = real(calcResidual(mesh, sbp, eqn, opts, evalPhysicsResidual))
  res_norm_0 = res_norm

  # print to log file
  println(fconv, 0, " ", res_norm)
  flush(fconv)

  eqn.majorIterationCallback(0, mesh, sbp, eqn, opts, STDOUT)

  while res_norm > res_norm_0*res_reltol && res_norm > res_abstol && iter < itermax  # predictor loop

    println("\npredictor iteration ", iter, ", lambda = ", lambda)
    println("res_norm = ", res_norm)
    println("res_norm_0 = ", res_norm_0)
    println("res_norm/res_norm_0 = ", res_norm/res_norm_0)
    println("res_norm = ", res_norm)

    # calculate homotopy residual
    homotopy_norm = calcResidual(mesh, sbp, eqn, opts, homotopyPhysics)

    homotopy_norm0 = homotopy_norm
    copy!(q_vec0, eqn.q_vec)  # save initial q to calculate delta_q later

    # if we have finished traversing the homotopy path, solve the 
    # homotopy problem (= the physics problem because lambda is zero)
    # tightly
    if abs(lambda - lambda_min) <= eps()
      println("tightening homotopy tolerance")
      homotopy_tol = res_reltol
    end

    # do corrector steps
    # TODO: pass tolerances to newtonInnter
    newtonInner(newton_data, mesh, sbp, eqn, opts, rhs_func, jac_func, jac, 
                rhs_vec, ctx_residual, res_reltol=homotopy_tol, 
                res_abstol=res_abstol, itermax=30)

    # compute delta_q
    for i=1:length(eqn.q)
      delta_q[i] = eqn.q[i] - q_vec0[i]
    end

    println("delta_q norm = ", norm(delta_q))

    delta = real(norm(delta_q))
    if iter == 2
      delta_max = delta
    end

    # predictor step calculation
    if abs(lambda - lambda_min) > eps()
      # calculate dHdLambda at new q value
      calcdHdLambda(mesh, sbp, eqn, opts, rhs_vec)
      for i=1:length(rhs_vec)
        dHdLambda_real[i] = real(rhs_vec[i])
      end

      # calculate tangent vector dH/dq * t = dH/dLambda
      jac_func(newton_data, mesh, sbp, eqn, opts, jac, ctx_residual)

      matrixSolve(newton_data, eqn, mesh, opts, jac, tan_vec, dHdLambda_real, STDOUT)
      # normalize tangent vector
      tan_norm = sqrt(dot(tan_vec, tan_vec) + 1)  # + 1 ?
      scale!(tan_vec, 1/tan_norm)

      psi = psi_max
      if iter > 1
        tan_term = dot(tan_vec_1, tan_vec) 
        tan_norm_term = (1/tan_norm)*(1/tan_norm_1) 
        println("tan_term = ", tan_term)
        println("tan_norm_term = ", tan_norm_term)
        println("acos argument = ", tan_term + tan_norm_term)
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

      println("iteration ", iter, ", psi = ", psi, ", delta/delta_max = ", delta/delta_max, ", step size = ", h)

      # take predictor step
      scale!(tan_vec, h)
      for i=1:length(eqn.q_vec)
        eqn.q_vec[i] += tan_vec[i]
      end

      lambda = max(lambda_min, lambda - h)
      eqn.params.homotopy_lambda = lambda


    end  # end if lambda too large

    # calculate physics residual at new state q
    println("old res_norm = ", res_norm)
    res_norm = real(calcResidual(mesh, sbp, eqn, opts, evalPhysicsResidual))
    println("new_res_norm = ", res_norm)

    # print to log file
    println(fconv, iter, " ", res_norm, " ", h )
    flush(fconv)

    eqn.majorIterationCallback(iter, mesh, sbp, eqn, opts, STDOUT)

    iter += 1
  end  # end while loop

  print("\n")

  # inform user of final status
  if iter >= itermax
    println(STDERR, "Warning: predictor-corrector did not converge in $iter iterations")
  
  elseif res_norm <= res_abstol
    println("predictor=corrector converged with absolute residual norm $res_norm")
  elseif res_norm/res_norm_0 <= res_reltol
    tmp = res_norm/res_norm_0
    println("predictor=corrector converged with relative residual norm $tmp")
  end

  return nothing
end



"""
  This function makes it appear as though the combined homotopy function H
  is a physics.  This works because an elementwise combinations of physics
  is still a physics.

  evalPhysicsResidual is used for the physcs function R and evalHomotopyResidual is used for 
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

  global evalPhysicsResidual

  lambda = eqn.params.homotopy_lambda
  res_homotopy = zeros(eqn.res)


  # calculate physics residual
  evalPhysicsResidual(mesh, sbp, eqn, opts, t)


  # calculate homotopy function
  evalHomotopyResidual(mesh, sbp, eqn, opts, res_homotopy)

  # combine 
  lambda_c = 1 - lambda # complement of lambda
  for i=1:length(eqn.res)
    eqn.res[i] =  lambda_c*eqn.res[i] + lambda*res_homotopy[i]
  end

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

  Inputs/Outputs
    res_vec: vector to store dH/dLambda in

  Aliasing restrictions: res_vec and eqn.res_vec may not alias
"""
function calcdHdLambda(mesh, sbp, eqn, opts, res_vec)

  # it appears this only gets called after parallel communication is done
  # so no need to start communication here

  global evalPhysicsResidual

  lambda = eqn.params.homotopy_lambda
  res_homotopy = zeros(eqn.res)


  # calculate physics residual
  evalPhysicsResidual(mesh, sbp, eqn, opts)
  assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)


  # calculate homotopy function
  evalHomotopyResidual(mesh, sbp, eqn, opts, res_homotopy)
  assembleSolution(mesh, sbp, eqn, opts, res_homotopy, res_vec)

  # combine them
  for i=1:length(res_vec)
    res_vec[i] -= eqn.res_vec[i]
  end

  return nothing
end

