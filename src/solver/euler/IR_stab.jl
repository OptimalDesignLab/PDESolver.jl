#  IRStab.jl: functions for stabilizing the IR flux

"""
  Computes dq/dv, where q are the conservative variables and v are the
  IR entropy variables.  This is equiavlent to calcA0 scaled by gamma_1,
  but computed from the conservative variables, which is much less expensive.

  Methods are available for 2 and 3 dimensions.
  A0 is overwritten with the result
"""
@inline function getIRA0{Tsol}(params::ParamType{2}, 
                           q::AbstractArray{Tsol,1}, 
                           A0::AbstractArray{Tsol, 2})


  gamma = params.gamma
  gamma_1 = params.gamma_1
  p = calcPressure(params, q)

  rho = q[1]
  rhou = q[2]
  rhov = q[3]
  rhoe = q[4]

  rhoinv = 1/rho

  h = (rhoe + p)*rhoinv
  a2 = gamma*p*rhoinv  # speed of sound

  A0[1,1] = rho
  A0[2,1] = rhou
  A0[3,1] = rhov
  A0[4,1] = rhoe

  A0[1,2] = rhou
  A0[2,2] = rhou*rhou*rhoinv + p
  A0[3,2] = rhou*rhov*rhoinv
  A0[4,2] = rhou*h

  A0[1,3] = rhov
  A0[2,3] = rhou*rhov/rho
  A0[3,3] = rhov*rhov*rhoinv + p
  A0[4,3] = rhov*h

  A0[1,4] = rhoe
  A0[2,4] = h*rhou
  A0[3,4] = h*rhov
  A0[4,4] = rho*h*h - a2*p/gamma_1

  return nothing
end

@inline function getIRA0{Tsol}(params::ParamType{3}, 
                           q::AbstractArray{Tsol,1}, 
                           A0::AbstractArray{Tsol, 2})


  gamma = params.gamma
  gamma_1 = params.gamma_1
  p = calcPressure(params, q)

  rho = q[1]
  rhou = q[2]
  rhov = q[3]
  rhow = q[4]
  rhoe = q[5]

  rhoinv = 1/rho

  h = (rhoe + p)*rhoinv
  a2 = gamma*p*rhoinv  # speed of sound

  A0[1,1] = rho
  A0[2,1] = rhou
  A0[3,1] = rhov
  A0[4,1] = rhow
  A0[5,1] = rhoe

  A0[1,2] = rhou
  A0[2,2] = rhou*rhou*rhoinv + p
  A0[3,2] = rhou*rhov*rhoinv
  A0[4,2] = rhou*rhow*rhoinv
  A0[5,2] = rhou*h

  A0[1,3] = rhov
  A0[2,3] = rhou*rhov/rho
  A0[3,3] = rhov*rhov*rhoinv + p
  A0[4,3] = rhov*rhow*rhoinv
  A0[5,3] = rhov*h


  A0[1,4] = rhow
  A0[2,4] = rhow*rhou*rhoinv
  A0[3,4] = rhow*rhov*rhoinv
  A0[4,4] = rhow*rhow*rhoinv + p
  A0[5,4] = rhow*h

  A0[1,5] = rhoe
  A0[2,5] = h*rhou
  A0[3,5] = h*rhov
  A0[4,5] = h*rhow
  A0[5,5] = rho*h*h - a2*p/gamma_1

  return nothing
end



"""
  This function computes the IR stabilization using simple averaging of
  qL and qR.  The term is subtracted off of F.

  This function is dimension agnostic, but only works for conservative
  variables.

  Aliasing restrictions: params.q_vals3, see also getIRStab_inner
"""
function getIRStab1{Tmsh, Tsol, Tres, Tdim}(
                      params::ParamType{Tdim, :conservative}, 
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tres},
                      dir::AbstractArray{Tmsh},  F::AbstractArray{Tres,1})

  q_avg = params.q_vals3
  for i=1:length(q_avg)
    q_avg[i] = 0.5*(qL[i] + qR[i])
  end
  getIRStab_inner(params, qL, qR, q_avg, aux_vars, dir, F)

  return nothing
end


"""
  Updates the vector F with the stabilization term from Carpenter, Fisher,
  Nielsen, Frankel, Entrpoy stable spectral collocation schemes for the 
  Navier-Stokes equatiosn: Discontinuous interfaces.  The term is subtracted
  off from F.

  The q_avg vector should some average of qL and qR, but the type of 
  averaging is left up to the user.

  This function is agnostic to dimension, but only works for conservative
  variables.

  Aliasing: from params the following arrays are used: A0, v_vals
              v_vals2.

"""
function getIRStab_inner{Tmsh, Tsol, Tres, Tdim}(
                      params::ParamType{Tdim, :conservative}, 
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      q_avg::AbstractArray{Tsol}, aux_vars::AbstractArray{Tres},
                      dir::AbstractArray{Tmsh},  F::AbstractArray{Tres,1})
#  println("entered getIRStab_inner")

  A0 = params.A0
  vL = params.v_vals
  vR = params.v_vals2
  gamma = params.gamma
  gamma_1inv = 1/params.gamma_1
  p = calcPressure(params, q_avg)

  convertToEntropy(params, qL, vL)
  convertToEntropy(params, qR, vR)

  for i=1:length(vL)
    vL[i] = gamma_1inv*(vR[i] - vL[i]) # scale by 1/gamma_1 to make IR entropy
                                       # variables, now vL has vL - vR
  end

#  println("delta v = \n", vL)
  # common-subexpression-elimination has a strong influence on the weak minded
  getIRA0(params, q_avg, A0)

  # multiply into vR
  smallmatvec!(A0, vL, vR)

  # calculate lambda_max at average state
  rhoinv = 1/q_avg[1]
  a = sqrt(gamma*p*rhoinv)  # speed of sound

  Un = zero(Tres)
  dA = zero(Tmsh)
#  Un = nx*q_avg[2]*rhoinv + ny*q_avg[3]*rhoinv + nz*q_avg[4]*rhoinv
#  dA = sqrt(nx*nx + ny*ny + nz*nz)
  for i=1:Tdim
    Un += dir[i]*q_avg[i+1]*rhoinv
    dA += dir[i]*dir[i]
  end
  dA = sqrt(dA)
#=
  println("qL = \n", qL)
  println("qR = \n", qR)
  println("q_avg = \n", q_avg)
  println("dir = \n", dir)
=#
  # the two eigenvalues are Un + dA*a and Un - dA*a, so depending on
  # the sign of Un, the maximum is abs(Un) + dA*a
#  lambda_max1 = absvalue(Un) + dA*a
#  println("lambda_max1 = ", lambda_max1)

  # use the Roe solver code
  rhoA = absvalue(Un) + dA*a
  lambda1 = Un + dA*a
  lambda2 = Un - dA*a
  lambda3 = Un
  sat_Vn = convert(Tsol, 0.025)
  sat_Vl = convert(Tsol, 0.025)
  tau = 1
  lambda1 = (tau*max(absvalue(lambda1),sat_Vn *rhoA) - lambda1)
  lambda2 = (tau*max(absvalue(lambda2),sat_Vn *rhoA) - lambda2)
  lambda3 = (tau*max(absvalue(lambda3),sat_Vl *rhoA) - lambda3)
  lambda_max1 = max(lambda1, lambda2)
  lambda_max1 = max(lambda_max1, lambda3)


  # DEBUGGING: try definition from Carpenters paper

#  println("Un = ", Un, ", a = ", a, ", dA = ", dA)
#  lambda_max2 = getLambdaMax(params, qL, qR, dir)
#  println("lambda_max1 = ", lambda_max1, ", lambda_max2 = ", lambda_max2)
  lambda_max = lambda_max1

#  println("lambda_max = ", lambda_max)
  fac = 1
  for i=1:length(vR)
    F[i] -= fac* 0.5*lambda_max*vR[i]
  end

  return nothing
end

function getLambdaMax{Tsol, Tmsh, Tdim}(params::ParamType{Tdim}, qL::AbstractVector{Tsol}, 
                      qR::AbstractVector{Tsol}, dir::AbstractVector{Tmsh})
# compute lambda_max approximation from Carpenter's Entropy Stable Collocation
# Schemes paper

   gamma = params.gamma
  Tres = promote_type(Tsol, Tmsh)
  UnL = zero(Tres)
  UnR = zero(Tres)
  rhoLinv = 1/qL[1]
  rhoRinv = 1/qR[1]

  pL = calcPressure(params, qL)
  pR = calcPressure(params, qR)

  aL = sqrt(gamma*pL*rhoLinv)  # speed of sound
  aR = sqrt(gamma*pR*rhoRinv)  # speed of sound
#  Un = nx*q_avg[2]*rhoinv + ny*q_avg[3]*rhoinv + nz*q_avg[4]*rhoinv
#  dA = sqrt(nx*nx + ny*ny + nz*nz)
  for i=1:Tdim
    UnL += dir[i]*qL[i+1]*rhoLinv
    UnR += dir[i]*qR[i+1]*rhoRinv
  end

#  println("UnL = ", UnL, ", UnR = ", UnR, ", aL = ", aL, ", aR = ", aR)

  lambda_max = 0.5*(UnL^4 + aL^4 + UnR^4 + aR^4)
  lambda_max = lambda_max^(1/4)

  return lambda_max
end
