# Dissipation based homotopy predictor-corrector globalization for solving
# steady problems using Newton's method
# based on Brown and Zingg, "A monolithic Homotopy Continuation Algorithm 
# with applications to Computational Fluid Dynamics"
# Journal of Computational Physics 321, (2016), 55-75
# specifically, Algorithm 2

"""
  This function solves steady problems using a dissipation-based
  predictor-correcor homotopy globalization for Newtons method.

  Inputs:
    func: the function to solve, ie. func(q) = 0  mathematically.
          func must have the signature func(mesh, sbp, eqn, opts)

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
function predictorCorrectorHomotopy{Tsol, Tres, Tmsh}(func::Function, 
                                    mesh::AbstractMesh{Tmsh}, 
                                    sbp::AbstractSBP, 
                                    eqn::AbstractSolutionData{Tsol, Tres}, 
                                    opts, pmesh=mesh)

  homotopy_func =  # TODO: the homotopy function
  newton_data, jac, rhs_vec = setupNewton(mesh, pmesh, sbp, eqn, opts, 
                                          homotopy_func)
  

  lambda = 1.0  # homotopy parameter

  # some parameters
  lambda_min = 0.0
  itermax = 100

  # counters/loop variables
  iter = 1

  # needed arrays


  while iter < itermax  # predictor loop

    # calculate homotopy residual

    if # norm of homotopy residual too large
      newtonInner(newton_data, mesh, sbp, eqn, opts, rhs_func, jac_func, jac, 
                   rhs_vec)
    end


    # calculate some quantities of interest
    # tangent vector, delta_q, delta, phi

    # calculate step length h

    # compute predictor step q += h*delta_q

    iter += 1
  end

  return nothing
end
