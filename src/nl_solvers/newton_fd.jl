using ArrayViews
using Debug


@debug function newton_fd(func, mesh, sbp, eqn; itermax=10, step_tol=1e-6)
  # this function drives the non-linear residual to some specified tolerance
  # using Newton's Method
  # the jacobian is formed using finite differences
  # the initial condition is stored in eqn.SL0
  # itermax is the maximum number of iterations

  m = length(eqn.SL)
  Tsol = typeof(eqn.SL[1])
  jac = zeros(Tsol, m, m)  # storage of the jacobian matrix
  res_0 = zeros(eqn.SL0)  # function evaluated at u0
  delta_SL = zeros(eqn.SL0)  # newton update
  step_norm = zero(Tsol)  # norm of newton update

  for i=1:itermax
  println("Newton iteration: ", i)
  # compute residual at initial condition
  func(mesh, sbp, eqn, eqn.SL0, eqn.SL)
  res_0 = deepcopy(eqn.SL) 
#  println("SL0 = ", res_0)

    epsilon = 1e-6  # finite difference perturbation
    # calculate jacobian
    for j=1:m
      if j==1
	eqn.SL0[j] +=  epsilon
      else
	eqn.SL0[j-1] -= epsilon # undo previous iteration pertubation
	eqn.SL0[j] += epsilon
      end

      # evaluate residual
      fill!(eqn.SL, 0.0)
      func(mesh, sbp, eqn, eqn.SL0, eqn.SL)
 #     println("column ", j, " of jacobian, SL = ", eqn.SL)
      calcJacRow(unsafe_view(jac, :, j), res_0, eqn.SL, epsilon)
      println("SL norm = ", norm(SL)/m)
      
    end

#    @bp
    # now jac is complete
    delta_SL[:] = jac\(-res_0)  #  calculate Newton update
    eqn.SL0[:] += delta_SL  # update SL0

    step_norm = norm(delta_SL)/m

#    println("delta_SL = ", delta_SL)
    println("delta_sl = ")
    for i=1:m
      println(i, "   ", delta_SL[i])
    end
    println("step_norm = ", step_norm)
#    println("jac = ", jac)

    print("\n")

    # check stopping criteria
    if (step_norm < step_tol)
      println("Newton iteration converged with step_norm = ", step_norm)

      # compute residual at final SL0
      func(mesh, sbp, eqn, eqn.SL0, eqn.SL)
      return nothing
    end
  end  # end loop over newton iterations

  println("Warning: Newton iteration did not converge in ", itermax, " iterations")
  println("  Final step size: ", step_norm)
  return nothing
end


function calcJacRow(jac_row, res_0, res, epsilon)
# calculate a row of the jacobian from res_0, the function evaluated 
# at the original point, and res, the function evaluated at a perturbed point

m = length(res_0)

for i=1:m
  jac_row[i] = (res[i] - res_0[i])/epsilon
end

return nothing

end




