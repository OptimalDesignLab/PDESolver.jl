export newton, newton_check, newton_check_fd

@doc """
  This function uses Newton's method to reduce the residual.  The Jacobian
  is calculated using one of several available methods.

  Arguments:
    * func  : function that evalutes the residual
    * mesh : mesh to use in evaluating the residual
    * sbp : sbp operator to be used to evaluate the residual
    * eqn : EulerData to use to evaluate the residual
    * opts : options dictonary

    Optional Arguments
    * itermax : maximum number of Newton iterations
    * step_tol : step size stopping tolerance
    * res_tol : residual stopping tolerance

    func must have the signature func(mesh, sbp, eqn, opts, eqn.SL0, eqn.SL) 

"""->
function newton(func, mesh, sbp, eqn, opts; itermax=200, step_tol=1e-6, res_tol=1e-6)
  # this function drives the non-linear residual to some specified tolerance
  # using Newton's Method
  # the jacobian is formed using finite differences
  # the initial condition is stored in eqn.SL0
  # itermax is the maximum number of iterations


  # options
  write_rhs = opts["write_rhs"]::Bool
  write_jac = opts["write_jac"]::Bool
  print_cond = opts["print_cond"]::Bool
  write_sol = opts["write_sol"]::Bool
  write_vis = opts["write_vis"]::Bool
  jac_method = opts["jac_method"]::Int

  step_fac = 1.0 # step size limiter
  m = length(eqn.SL)
  Tsol = typeof(eqn.SL[1])
  Tjac = typeof(real(eqn.SL[1]))  # type of jacobian, residual
  jac = zeros(Tjac, m, m)  # storage of the jacobian matrix
  res_0 = zeros(Tjac, m)  # function evaluated at u0
  res_0_norm = 0.0  # norm of res_0
  delta_SL = zeros(Tjac, m)  # newton update
  step_norm = zero(Tjac)  # norm of newton update
  step_norm_1 = zero(Tjac) # norm of previous newton update


  # open file to write convergence data to
  # append to be on the safe side
  fconv = open("convergence.dat", "a+")


  # evaluating residual at initial condition
  println("evaluating residual at initial condition")
  res_0_norm = calcResidual(mesh, sbp, eqn, opts, func, res_0)

  # write rhs to file
  if write_rhs
    writedlm("rhs1.dat", res_0)
  end

  if res_0_norm < res_tol
   println("Initial condition satisfies res_tol with residual norm ", res_0_norm)
   println("writing to convergence.dat")
   println(fconv, i, " ", res_0_norm, " ", 0.0)

   close(fconv)
   return nothing
 end


  # do Newton's method if not converged
  print("\n")


  for i=1:itermax
    println("Newton iteration: ", i)
    println("step_fac = ", step_fac)

    # calculate jacobian using selected method
    if jac_method == 1
#      println("calculating finite difference jacobian")
      calcJacFD(mesh, sbp, eqn, opts, func, res_0, jac)

    elseif jac_method == 2
#      println("calculating complex step jacobian")
      calcJacobianComplex(mesh, sbp, eqn, opts, func, jac)
    end

    # print as determined by options
    if write_jac
#      fname = string("jacobian", i, ".dat")
#      printMatrix(fname, jac)
      writedlm("jacobian$i.dat", jac)
      println("finished printing jacobian")
    end

    # calculate Jacobian condition number
    if print_cond
      cond_j = cond(jac)
      println("Condition number of jacobian = ", cond_j)
    end

    # calculate Newton step
    delta_SL[:] = jac\(-res_0)  #  calculate Newton update
    step_norm = norm(delta_SL)/m
    println("step_norm = ", step_norm)

    # perform Newton update
    eqn.SL0[:] += step_fac*delta_SL  # update SL0

    # write starting values for next iteration to file
    if write_sol
      writedlm("SL0$i.dat", eqn.SL0)
    end

    # write paraview files
    if write_vis
      vals = abs(real(eqn.SL0))  # remove unneded imaginary part
      saveSolutionToMesh(mesh, vals)
      fname = string("solution_newton", i)
      writeVisFiles(mesh, fname)
    end
 

    # write to convergence file
    println(fconv, i, " ", res_0_norm, " ", step_norm)
    flush(fconv)

    # calculate residual at updated location, used for next iteration rhs
    res_0_norm = calcResidual(mesh, sbp, eqn, opts, func, res_0)

    # write rhs to file
    if write_rhs
      tmp = i+1
      writedlm("rhs$tmp.dat", res_0)
    end



   if res_0_norm < res_tol
     println("Newton iteration converged with residual norm ", res_0_norm)
     close(fconv)

     return nothing
   end

    if (step_norm < step_tol)
      println("Newton iteration converged with step_norm = ", step_norm)
      println("Final residual = ", res_0_norm)
      close(fconv)

      return nothing
    end

#=
    # adjust step size limiter
    if (step_norm < step_norm_1)  # decreasing step size
      step_fac *= 1.2

      if step_fac > 1.0
	step_fac = 1.0
      end
    end
=#
#    if (step_norm > step_norm_1)
#      step_fac /= 1.1
#    end


    print("\n")
    step_norm_1 = step_norm
  end  # end loop over newton iterations

  println("Warning: Newton iteration did not converge in ", itermax, " iterations")
  println("  Final step size: ", step_norm)
  println("  Final residual: ", res_0_norm)
  close(fconv)
  return nothing
end


function calcResidual(mesh, sbp, eqn, opts, func, res_0)
# calculate the residual and its norm

  m = length(res_0)

  fill!(eqn.SL, 0.0)
  func(mesh, sbp, eqn, opts, eqn.SL0, eqn.SL)
#  res_0[:] = real(eqn.SL)  # is there an unnecessary copy here?

  for j=1:m
    res_0[j] = real(eqn.SL[j])
  end

  res_0_norm = norm(eqn.SL)/m
  println("residual norm = ", res_0_norm)

 return res_0_norm
end


function calcJacFD(mesh, sbp, eqn, opts, func, res_0, jac)
# calculate the jacobian using finite difference

  (m,n) = size(jac)
  entry_orig = zero(eltype(eqn.SL0))
  epsilon = 1e-6  # finite difference perturbation
  # calculate jacobian
  for j=1:m
#      println("  jacobian iteration ", j)
    if j==1
      entry_orig = eqn.SL0[j]
      eqn.SL0[j] +=  epsilon
    else
      eqn.SL0[j-1] = entry_orig # undo previous iteration pertubation
      entry_orig = eqn.SL0[j]
      eqn.SL0[j] += epsilon
    end

    # evaluate residual
    fill!(eqn.SL, 0.0)
    func(mesh, sbp, eqn, opts, eqn.SL0, eqn.SL)
#     println("column ", j, " of jacobian, SL = ", eqn.SL)
    calcJacRow(unsafe_view(jac, :, j), res_0, eqn.SL, epsilon)
#      println("SL norm = ", norm(SL)/m)
    
  end

  # undo final perturbation
  eqn.SL0[m] = entry_orig

  return nothing
end



function calcJacRow{T <: Real}(jac_row, res_0, res::AbstractArray{T,1}, epsilon)
# calculate a row of the jacobian from res_0, the function evaluated 
# at the original point, and res, the function evaluated at a perturbed point

m = length(res_0)

for i=1:m
  jac_row[i] = (res[i] - res_0[i])/epsilon
end

return nothing

end





function calcJacobianComplex(mesh, sbp, eqn, opts, func, jac)

  epsilon = 1e-20  # complex step perturbation
  entry_orig = zero(eltype(eqn.SL0))
  (m,n) = size(jac)
  # calculate jacobian
  for j=1:m
    if j==1
      entry_orig = eqn.SL0[j]
      eqn.SL0[j] +=  complex(0, epsilon)
    else
      eqn.SL0[j-1] = entry_orig # undo previous iteration pertubation
      entry_orig = eqn.SL0[j]
      eqn.SL0[j] += complex(0, epsilon)
    end

    # evaluate residual
    fill!(eqn.SL, 0.0)
    func(mesh, sbp, eqn, opts, eqn.SL0, eqn.SL)
#     println("column ", j, " of jacobian, SL = ", eqn.SL)
    calcJacRow(unsafe_view(jac, :, j), eqn.SL, epsilon)
#      println("SL norm = ", norm(SL)/m)
    
  end  # end loop over rows of jacobian


  # undo final perturbation
  eqn.SL0[m] = entry_orig
#

  return nothing
end


function calcJacRow{T <: Complex}(jac_row, res::AbstractArray{T, 1}, epsilon)
# calculate a row of the jacobian from res_0, the function evaluated 
# at the original point, and res, the function evaluated at a perturbed point

m = length(res)

for i=1:m
  jac_row[i] = imag(res[i])/epsilon
end

return nothing

end

@doc """
### newton_check

  Uses complex step to compare jacobian vector product to directional derivative.

"""->
function newton_check(func, mesh, sbp, eqn, opts)
  # this function drives the non-linear residual to some specified tolerance
  # using Newton's Method
  # the jacobian is formed using finite differences
  # the initial condition is stored in eqn.SL0
  # itermax is the maximum number of iterations

  step_fac = 0.5  # step size limiter
  m = length(eqn.SL)
  Tsol = typeof(eqn.SL[1])
  Tjac = typeof(real(eqn.SL[1]))  # type of jacobian, residual
  jac = zeros(Tjac, m, m)  # storage of the jacobian matrix
  direction_der = zeros(mesh.numDof)
#  v = rand(mesh.numDof)
   v = readdlm("randvec.txt")

  epsilon = 1e-20  # complex step perturbation
  fill!(eqn.SL, 0.0)  # zero out SL
  # compute directional derivative
  for i=1:mesh.numDof
    eqn.SL0[i] += complex(0, epsilon*v[i])  # apply perturbation
  end

  func(mesh, sbp, eqn, opts, eqn.SL0, eqn.SL)
  println("evaluated directional derivative")

  # calculate derivative
  for i=1:mesh.numDof
    direction_der[i] = imag(eqn.SL[i])/epsilon
    eqn.SL0[i] -= complex(0, epsilon*v[i])  # undo perturbation
  end



    println("Calculating Jacobian")

    # calculate jacobian
    for j=1:m
      println("\ncalculating column ", j, " of the jacobian")
      if j==1
	eqn.SL0[j] +=  complex(0, epsilon)
      else
	eqn.SL0[j-1] -= complex(0, epsilon) # undo previous iteration pertubation
	eqn.SL0[j] += complex(0, epsilon)
      end

      # evaluate residual
      fill!(eqn.SL, 0.0)
      func(mesh, sbp, eqn, opts, eqn.SL0, eqn.SL)
 #     println("column ", j, " of jacobian, SL = ", eqn.SL)
      calcJacRow(unsafe_view(jac, :, j), eqn.SL, epsilon)
#      println("SL norm = ", norm(SL)/m)
      
    end  # end loop over rows of jacobian

    # undo final perturbation
    eqn.SL0[m] -= complex(0, epsilon)

    # now jac is complete

    fname = string("jacobian", ".dat")
    printMatrix(fname, jac)
    println("finished printing jacobian")

    cond_j = cond(jac)
    println("Condition number of jacobian = ", cond_j)
    svals = svdvals(jac)
    println("svdvals = \n", svals)

    jac_mult = jac*v

    # copy difference between directional derivative and
    # jacobian multiplication into SL for return

    for i=1:mesh.numDof
      eqn.SL[i] = direction_der[i] - jac_mult[i]
    end

    err_norm = norm(eqn.SL)/mesh.numDof
    println("step_norm = ", err_norm)
#    println("jac = ", jac)

    print("\n")

    println("finished newton_check")
  return nothing
end



@doc """
### newton_check

  This method calculates a single column of the jacobian with the complex step method.
"""->
function newton_check(func, mesh, sbp, eqn, opts, j)
# calculate a single column of hte jacobian
    
      jac_col = zeros(Float64, mesh.numDof)
      println("\ncalculating column ", j, " of the jacobian")

      epsilon = 1e-20

      eqn.SL0[j] += complex(0, epsilon)

      # evaluate residual
      fill!(eqn.SL, 0.0)
      func(mesh, sbp, eqn, opts, eqn.SL0, eqn.SL)
 #     println("column ", j, " of jacobian, SL = ", eqn.SL)
      calcJacRow(jac_col, eqn.SL, epsilon)
#      println("SL norm = ", norm(SL)/m)

      return jac_col
end 

@doc """
### newton_check_fd

  This method calcualtes a single column of the jacobian with the finite difference method.

"""->
function newton_check_fd(func, mesh, sbp, eqn, opts, j)
# calculate a single column of hte jacobian
    
      jac_col = zeros(Float64, mesh.numDof)
      println("\ncalculating column ", j, " of the jacobian")

     func(mesh, sbp, eqn, opts, eqn.SL0, eqn.SL)
     res_0 = copy(eqn.SL)

      epsilon = 1e-6

      eqn.SL0[j] += epsilon

      # evaluate residual
      fill!(eqn.SL, 0.0)
      func(mesh, sbp, eqn, opts, eqn.SL0, eqn.SL)
 #     println("column ", j, " of jacobian, SL = ", eqn.SL)

      calcJacRow(jac_col, res_0, eqn.SL, epsilon)
#      println("SL norm = ", norm(SL)/m)

      return jac_col
end 
 
