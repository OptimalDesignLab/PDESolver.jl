#  IRStab.jl: functions for stabilizing the IR flux
#             These functions are called from bc_solvers.jl to produce
#             fluxes that include dissipation terms.

include("IR_stab_diff.jl")

"""
  Computes dq/dv, where q are the conservative variables and v are the
  IR entropy variables.  This is equiavlent to calcA0 scaled by gamma_1,
  but computed from the conservative variables, which is much less expensive.

  Methods are available for 2 and 3 dimensions.
  A0 is overwritten with the result
"""
@inline function getIRA0(params::ParamType{2}, 
                     q::AbstractArray{Tsol,1}, 
                     A0::AbstractArray{Tsol, 2}) where Tsol


  gamma = params.gamma
  gamma_1 = params.gamma_1
  p = calcPressure(params, q)

  rho = q[1]
  rhou = q[2]
  rhov = q[3]
  rhoe = q[4]

  rhoinv = 1/rho

  h = (rhoe + p)*rhoinv  # this isn't really h, but including the factor of
                         # 1/rho is convenient
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

@inline function getIRA0(params::ParamType{3}, 
                     q::AbstractArray{Tsol,1}, 
                     A0::AbstractArray{Tsol, 2}) where Tsol


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
  Gets inverse matrix of [`getIRA0`](@ref)
"""
function getIRA0inv(params::ParamType{2}, q::AbstractArray{Tsol, 1},
                    A0inv::AbstractMatrix{Tsol}) where {Tsol}

  gamma = params.gamma
  gamma_1 = params.gamma_1
  gamma_1i = 1/gamma_1

  k1 = 0.5*(q[2]^2 + q[3]^2)/q[1]
  k1_dot1 = -k1/q[1]
  k1_dot2 = q[2]/q[1]
  k1_dot3 = q[3]/q[1]
#  k1_dot4 = 0

  rho_int = q[4] - k1
  rho_int_dot1 = -k1_dot1
  rho_int_dot2 = -k1_dot2
  rho_int_dot3 = -k1_dot3
  rho_int_dot4 =  1

  s = log(gamma_1*rho_int/(q[1]^gamma))
  s_dot1 = rho_int_dot1/rho_int - gamma/q[1]
  s_dot2 = rho_int_dot2/rho_int
  s_dot3 = rho_int_dot3/rho_int
  s_dot4 = rho_int_dot4/rho_int

  fac = 1.0/rho_int
  t1 = -fac*fac
  fac_dot1 = t1*rho_int_dot1
  fac_dot2 = t1*rho_int_dot2
  fac_dot3 = t1*rho_int_dot3
  fac_dot4 = t1*rho_int_dot4

  t2 = (rho_int*(gamma + 1 - s) - q[4])

  A0inv[1, 1] = (rho_int_dot1*(gamma + 1) - (rho_int*s_dot1 + s*rho_int_dot1))*fac*gamma_1i + t2*fac_dot1*gamma_1i
  A0inv[2, 1] = q[2]*fac_dot1*gamma_1i
  A0inv[3, 1] = q[3]*fac_dot1*gamma_1i
  A0inv[4, 1] = -(fac + q[1]*fac_dot1)*gamma_1i

  A0inv[1, 2] = ( (rho_int_dot2*(gamma + 1) - (rho_int*s_dot2 +
                   s*rho_int_dot2))*fac + fac_dot2*t2)*gamma_1i
  A0inv[2, 2] = (fac + q[2]*fac_dot2)*gamma_1i
  A0inv[3, 2] = q[3]*fac_dot2*gamma_1i
  A0inv[4, 2] = -(q[1]*fac_dot2)*gamma_1i

  A0inv[1, 3] = ( (rho_int_dot3*(gamma + 1) - (rho_int*s_dot3 +
                   s*rho_int_dot3))*fac + fac_dot3*t2)*gamma_1i
  A0inv[2, 3] = q[2]*fac_dot3*gamma_1i
  A0inv[3, 3] = (fac + q[3]*fac_dot3)*gamma_1i
  A0inv[4, 3] = -(q[1]*fac_dot3)*gamma_1i

  A0inv[1, 4] = ( (rho_int_dot4*(gamma + 1) - (rho_int*s_dot4 +
                   s*rho_int_dot4) - 1)*fac + fac_dot4*t2)*gamma_1i
  A0inv[2, 4] = q[2]*fac_dot4*gamma_1i
  A0inv[3, 4] = q[3]*fac_dot4*gamma_1i
  A0inv[4, 4] = -(q[1]*fac_dot4)*gamma_1i

  return nothing
end

function getIRA0inv(params::ParamType{3}, q::AbstractArray{Tsol, 1},
                    A0inv::AbstractMatrix{Tsol}) where {Tsol}

  gamma = params.gamma
  gamma_1 = params.gamma_1
  gamma_1i = 1/gamma_1

  k1 = 0.5*(q[2]^2 + q[3]^2 + q[4]^2)/q[1]
  k1_dot1 = -k1/q[1]
  k1_dot2 = q[2]/q[1]
  k1_dot3 = q[3]/q[1]
  k1_dot4 = q[4]/q[1]
#  k1_dot4 = 0

  rho_int = q[5] - k1
  rho_int_dot1 = -k1_dot1
  rho_int_dot2 = -k1_dot2
  rho_int_dot3 = -k1_dot3
  rho_int_dot4 = -k1_dot4
  rho_int_dot5 =  1

  s = log(gamma_1*rho_int/(q[1]^gamma))
  s_dot1 = rho_int_dot1/rho_int - gamma/q[1]
  s_dot2 = rho_int_dot2/rho_int
  s_dot3 = rho_int_dot3/rho_int
  s_dot4 = rho_int_dot4/rho_int
  s_dot5 = rho_int_dot5/rho_int

  fac = 1.0/rho_int
  t1 = -fac*fac
  fac_dot1 = t1*rho_int_dot1
  fac_dot2 = t1*rho_int_dot2
  fac_dot3 = t1*rho_int_dot3
  fac_dot4 = t1*rho_int_dot4
  fac_dot5 = t1*rho_int_dot5

  t2 = (rho_int*(gamma + 1 - s) - q[5])

  A0inv[1, 1] = (rho_int_dot1*(gamma + 1) - (rho_int*s_dot1 + s*rho_int_dot1))*fac*gamma_1i + t2*fac_dot1*gamma_1i
  A0inv[2, 1] = q[2]*fac_dot1*gamma_1i
  A0inv[3, 1] = q[3]*fac_dot1*gamma_1i
  A0inv[4, 1] = q[4]*fac_dot1*gamma_1i
  A0inv[5, 1] = -(fac + q[1]*fac_dot1)*gamma_1i

  A0inv[1, 2] = ( (rho_int_dot2*(gamma + 1) - (rho_int*s_dot2 +
                   s*rho_int_dot2))*fac + fac_dot2*t2)*gamma_1i
  A0inv[2, 2] = (fac + q[2]*fac_dot2)*gamma_1i
  A0inv[3, 2] = q[3]*fac_dot2*gamma_1i
  A0inv[4, 2] = q[4]*fac_dot2*gamma_1i
  A0inv[5, 2] = -(q[1]*fac_dot2)*gamma_1i

  A0inv[1, 3] = ( (rho_int_dot3*(gamma + 1) - (rho_int*s_dot3 +
                   s*rho_int_dot3))*fac + fac_dot3*t2)*gamma_1i
  A0inv[2, 3] = q[2]*fac_dot3*gamma_1i
  A0inv[3, 3] = (fac + q[3]*fac_dot3)*gamma_1i
  A0inv[4, 3] = q[4]*fac_dot3*gamma_1i
  A0inv[5, 3] = -(q[1]*fac_dot3)*gamma_1i

  A0inv[1, 4] = (rho_int_dot4*(gamma + 1) - (rho_int*s_dot4 + s*rho_int_dot4))*fac*gamma_1i + t2*fac_dot4*gamma_1i
  A0inv[2, 4] = q[2]*fac_dot4*gamma_1i
  A0inv[3, 4] = q[3]*fac_dot4*gamma_1i
  A0inv[4, 4] = (fac + q[4]*fac_dot4)*gamma_1i
  A0inv[5, 4] = -(q[1]*fac_dot4)*gamma_1i


  A0inv[1, 5] = ( (rho_int_dot5*(gamma + 1) - (rho_int*s_dot5 +
                   s*rho_int_dot5) - 1)*fac + fac_dot5*t2)*gamma_1i
  A0inv[2, 5] = q[2]*fac_dot5*gamma_1i
  A0inv[3, 5] = q[3]*fac_dot5*gamma_1i
  A0inv[4, 5] = q[4]*fac_dot5*gamma_1i
  A0inv[5, 5] = -(q[1]*fac_dot5)*gamma_1i

  return nothing
end




"""
  This function applies any [`AbstractEntropyKernel`](@ref) when defining
  a type 1 face integral (the normal type) for an entropy-stable scheme
  using a diagonal E operator.

  **Inputs**

   * params: ParamType
   * kernel: the `AbstractEntropyKernel` to apply
   * qL: solution at left state
   * qR: solution at right state
   * aux_vars: auxiliary varialbes
   * dir: normal vector

  **Inputs/Outputs**

   * F: flux vector to have the entropy kernel contribution added to (well,
        subtracted because the contribution is negative).
"""
#function getEntropyLFStab(
function applyEntropyKernel_diagE(
                      params::ParamType{Tdim, :conservative},
                      kernel::AbstractEntropyKernel,
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tres},
                      dir::AbstractArray{Tmsh},  F::AbstractArray{Tres,1}) where {Tmsh, Tsol, Tres, Tdim}

  q_avg = params.apply_entropy_kernel_diagE_data.q_avg
  for i=1:length(q_avg)
    q_avg[i] = 0.5*(qL[i] + qR[i])
  end

  applyEntropyKernel_diagE_inner(params, kernel, qL, qR, q_avg, aux_vars, dir, F)

  return nothing
end


"""
  Applies the specified [`AbstractEntropyKernel`](@ref)

  **Inputs**

   * params
   * kernel: the kernel to apply
   * qL: left state
   * qR: right state
   * q_avg: the state at which to evaluate the kernel
   * aux_vars
   * dir: normal vector
 
  **Inputs/Outputs**

   * F: flux vector to update with contribtion
"""
#function getEntropyLFStab_inner(
function applyEntropyKernel_diagE_inner(
                      params::ParamType{Tdim, :conservative}, 
                      kernel::AbstractEntropyKernel,
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      q_avg::AbstractArray{Tsol}, aux_vars::AbstractArray{Tres},
                      dir::AbstractArray{Tmsh},  F::AbstractArray{Tres,1}) where {Tmsh, Tsol, Tres, Tdim}
#  println("entered getEntropyLFStab_inner")

  @unpack params.apply_entropy_kernel_diagE_data vL vR F_tmp
  gamma = params.gamma
  gamma_1inv = 1/params.gamma_1
  p = calcPressure(params, q_avg)

  convertToIR(params, qL, vL)
  convertToIR(params, qR, vR)

  for i=1:length(vL)
    vL[i] = vR[i] - vL[i]
  end

#  F_tmp = zeros(Tres, length(F))
  applyEntropyKernel(kernel, params, q_avg, vL, dir, F_tmp)

  for i=1:length(F_tmp)
    F[i] -= F_tmp[i]
  end

  return nothing
end

function getLambdaMax(params::ParamType{Tdim}, 
    qL::AbstractVector{Tsol}, qR::AbstractVector{Tsol}, 
    dir::AbstractVector{Tmsh}) where {Tsol, Tmsh, Tdim}
# compute lambda_max approximation from Carpenter's Entropy Stable Collocation
# Schemes paper

  gamma = params.gamma
  Tres = promote_type(Tsol, Tmsh)
  UnL = zero(Tres)
  UnR = zero(Tres)
  rhoLinv = 1/qL[1]
  rhoRinv = 1/qR[1]
  dA = zero(Tmsh)

  pL = calcPressure(params, qL)
  pR = calcPressure(params, qR)

  aL = sqrt(gamma*pL*rhoLinv)  # speed of sound
  aR = sqrt(gamma*pR*rhoRinv)  # speed of sound
#  Un = nx*q_avg[2]*rhoinv + ny*q_avg[3]*rhoinv + nz*q_avg[4]*rhoinv
#  dA = sqrt(nx*nx + ny*ny + nz*nz)
  for i=1:Tdim
    UnL += dir[i]*qL[i+1]*rhoLinv
    UnR += dir[i]*qR[i+1]*rhoRinv
    dA += dir[i]*dir[i]
  end

  dA = sqrt(dA)
  aL *= dA
  aR *= dA

#  println("UnL = ", UnL, ", UnR = ", UnR, ", aL = ", aL, ", aR = ", aR)

  lambda_max = 0.5*(UnL^4 + aL^4 + UnR^4 + aR^4)
  lambda_max = lambda_max^(1/4)

  return lambda_max
end

"""
  Calculates the maximum magnitude eigenvalue of the Euler flux 
  jacobian at the arithmatic average of two states.

  This functions works in both 2D and 3D
  Inputs:
    params:  ParamType, conservative variable
    qL: left state
    qR: right state
    dir: direction vector (does *not* have to be unit vector)

  Outputs:
    lambda_max: eigenvalue of maximum magnitude

"""
function getLambdaMaxSimple(params::ParamType{Tdim}, 
                      qL::AbstractVector{Tsol}, qR::AbstractVector{Tsol}, 
                      dir::AbstractVector{Tmsh}) where {Tsol, Tmsh, Tdim}
# calculate maximum eigenvalue at simple average state

  gamma = params.gamma
  Tres = promote_type(Tsol, Tmsh)

  q_avg = params.get_lambda_max_simple_data.q_avg

  for i=1:length(q_avg)
    q_avg[i] = 0.5*(qL[i] + qR[i])
  end

  Un = zero(Tres)
  dA = zero(Tres)
  rhoinv = 1/q_avg[1]
  p = calcPressure(params, q_avg)
  a = sqrt(gamma*p*rhoinv)  # speed of sound

  for i=1:Tdim
    Un += dir[i]*q_avg[i+1]*rhoinv
    dA += dir[i]*dir[i]
  end

  dA = sqrt(dA)

  lambda_max = absvalue(Un) + dA*a

  return lambda_max
end


function getLambdaMaxRoe(params::ParamType{Tdim}, 
                      qL::AbstractVector{Tsol}, qR::AbstractVector{Tsol}, 
                      dir::AbstractVector{Tmsh}) where {Tsol, Tmsh, Tdim}
# compute lambda_max approximation from Carpenter's Entropy Stable Collocation
# Schemes paper

  gamma = params.gamma
  Tres = promote_type(Tsol, Tmsh)
  Un = zero(Tres)
  dA = zero(Tmsh)
  rhoLinv = 1/qL[1]
  rhoRinv = 1/qR[1]

  pL = calcPressure(params, qL)
  pR = calcPressure(params, qR)

  aL = sqrt(gamma*pL*rhoLinv)  # speed of sound
  aR = sqrt(gamma*pR*rhoRinv)  # speed of sound
  a = 0.5*(aL + aR)
#  Un = nx*q_avg[2]*rhoinv + ny*q_avg[3]*rhoinv + nz*q_avg[4]*rhoinv
#  dA = sqrt(nx*nx + ny*ny + nz*nz)
  for i=1:Tdim
    Un += dir[i]*0.5*(qL[i+1]*rhoLinv + qR[i+1]*rhoLinv)
    dA += dir[i]*dir[i]
  end

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
  lambda_max1 = max(absvalue(lambda1), absvalue(lambda2))
  lambda_max1 = max(absvalue(lambda_max1), absvalue(lambda3))



  return lambda_max1
end
