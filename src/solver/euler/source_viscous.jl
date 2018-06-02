
#
# Given PRIMITIVE variable and its 1st order and 
# 2nd order derivatives, calculate the source terms 
#
function calcMmsSource(params::ParamType{3, :conservative},
                       q::AbstractArray{Tsrc, 1},
                       q_x::AbstractArray{Tsrc, 2},
                       q_xx::AbstractArray{Tsrc, 3},
                       src::AbstractArray{Tsrc, 1}) where Tsrc
  gamma = params.gamma
  gamma_1 = params.gamma - 1

  p   = q[1]*q[5]/gamma
  E   = q[5]/(gamma*gamma_1) + 0.5*(q[2]*q[2] + q[3]*q[3] + q[4]*q[4])

  p_x   = 1.0/gamma*(q_x[1,1]*q[5] + q[1]*q_x[1,5])
  p_y   = 1.0/gamma*(q_x[2,1]*q[5] + q[1]*q_x[2,5])
  p_z   = 1.0/gamma*(q_x[3,1]*q[5] + q[1]*q_x[3,5])
  E_x   = q_x[1,5]/(gamma*gamma_1) + (q[2]*q_x[1,2] + q[3]*q_x[1,3] + q[4]*q_x[1,4])	
  E_y   = q_x[2,5]/(gamma*gamma_1) + (q[2]*q_x[2,2] + q[3]*q_x[2,3] + q[4]*q_x[2,4])
  E_z   = q_x[3,5]/(gamma*gamma_1) + (q[2]*q_x[3,2] + q[3]*q_x[3,3] + q[4]*q_x[3,4])

  #
  # contribution from inviscid terms
  #
  src[1]  = q_x[1,1]*q[2] + q[1]*q_x[1,2] 
  src[1] += q_x[2,1]*q[3] + q[1]*q_x[2,3] 
  src[1] += q_x[3,1]*q[4] + q[1]*q_x[3,4]

  src[2]  = q_x[1,1]*q[2]*q[2] + q[1]*q_x[1,2]*q[2] + q[1]*q[2]*q_x[1,2] + p_x
  src[2] += q_x[2,1]*q[2]*q[3] + q[1]*q_x[2,2]*q[3] + q[1]*q[2]*q_x[2,3]
  src[2] += q_x[3,1]*q[2]*q[4] + q[1]*q_x[3,2]*q[4] + q[1]*q[2]*q_x[3,4]

  src[3]  = q_x[1,1]*q[3]*q[2] + q[1]*q_x[1,3]*q[2] + q[1]*q[3]*q_x[1,2]
  src[3] += q_x[2,1]*q[3]*q[3] + q[1]*q_x[2,3]*q[3] + q[1]*q[3]*q_x[2,3] + p_y
  src[3] += q_x[3,1]*q[3]*q[4] + q[1]*q_x[3,3]*q[4] + q[1]*q[3]*q_x[3,4]

  src[4]  = q_x[1,1]*q[4]*q[2] + q[1]*q_x[1,4]*q[2] + q[1]*q[4]*q_x[1,2]
  src[4] += q_x[2,1]*q[4]*q[3] + q[1]*q_x[2,4]*q[3] + q[1]*q[4]*q_x[2,3]
  src[4] += q_x[3,1]*q[4]*q[4] + q[1]*q_x[3,4]*q[4] + q[1]*q[4]*q_x[3,4] + p_z

  src[5]  = q_x[1,1]*E*q[2] + q[1]*E_x*q[2] + q[1]*E*q_x[1,2] + p_x*q[2] + p*q_x[1,2]
  src[5] += q_x[2,1]*E*q[3] + q[1]*E_y*q[3] + q[1]*E*q_x[2,3] + p_y*q[3] + p*q_x[2,3]
  src[5] += q_x[3,1]*E*q[4] + q[1]*E_z*q[4] + q[1]*E*q_x[3,4] + p_z*q[4] + p*q_x[3,4]

  if !params.isViscous 
    return nothing
  end

  muK = Array(typeof(q[5]), 2)
  getMuK(q[5], muK)
  rmu = muK[1]
  rK = muK[2]
  two3rd = 2.0/3.0

  #
  # in case someone forgets the lower diagonal of q_xx
  #
  q_xx[2,1,2] = q_xx[1,2,2]
  q_xx[3,1,2] = q_xx[1,3,2]
  q_xx[3,2,2] = q_xx[2,3,2]
  q_xx[2,1,3] = q_xx[1,2,3]
  q_xx[3,1,3] = q_xx[1,3,3]
  q_xx[3,2,3] = q_xx[2,3,3]
  q_xx[2,1,4] = q_xx[1,2,4]
  q_xx[3,1,4] = q_xx[1,3,4]
  q_xx[3,2,4] = q_xx[2,3,4]

  # tau
  txx = rmu * two3rd * (2*q_x[1,2] - q_x[2,3] - q_x[3,4])
  tyy = rmu * two3rd * (2*q_x[2,3] - q_x[1,2] - q_x[3,4])
  tzz = rmu * two3rd * (2*q_x[3,4] - q_x[1,2] - q_x[2,3])
  txy = rmu * (q_x[2,2] + q_x[1,3]) 
  txz = rmu * (q_x[3,2] + q_x[1,4]) 
  tyz = rmu * (q_x[2,4] + q_x[3,3]) 
  tyx = txy
  tzx = txz
  tzy = tyz

  # gradient of tau
  txx_x = rmu*two3rd*(2*q_xx[1,1,2] - q_xx[2,1,3] - q_xx[3,1,4])
  txx_y = rmu*two3rd*(2*q_xx[1,2,2] - q_xx[2,2,3] - q_xx[3,2,4])
  txx_z = rmu*two3rd*(2*q_xx[1,3,2] - q_xx[2,3,3] - q_xx[3,3,4])

  tyy_x = rmu*two3rd*(2*q_xx[2,1,3] - q_xx[1,1,2] - q_xx[3,1,4])
  tyy_y = rmu*two3rd*(2*q_xx[2,2,3] - q_xx[1,2,2] - q_xx[3,2,4])
  tyy_z = rmu*two3rd*(2*q_xx[2,3,3] - q_xx[1,3,2] - q_xx[3,3,4])

  tzz_x = rmu*two3rd*(2*q_xx[3,1,4] - q_xx[1,1,2] - q_xx[2,1,3])
  tzz_y = rmu*two3rd*(2*q_xx[3,2,4] - q_xx[1,2,2] - q_xx[2,2,3])
  tzz_z = rmu*two3rd*(2*q_xx[3,3,4] - q_xx[1,3,2] - q_xx[2,3,3])

  txy_x = rmu*(q_xx[2,1,2] + q_xx[1,1,3])
  txy_y = rmu*(q_xx[2,2,2] + q_xx[1,2,3])
  txy_z = rmu*(q_xx[2,3,2] + q_xx[1,3,3])
  txz_x = rmu*(q_xx[3,1,2] + q_xx[1,1,4])
  txz_y = rmu*(q_xx[3,2,2] + q_xx[1,2,4])
  txz_z = rmu*(q_xx[3,3,2] + q_xx[1,3,4])
  tyz_x = rmu*(q_xx[2,1,4] + q_xx[3,1,3])
  tyz_y = rmu*(q_xx[2,2,4] + q_xx[3,2,3])
  tyz_z = rmu*(q_xx[2,3,4] + q_xx[3,3,3])
  tyx_x = txy_x
  tyx_y = txy_y
  tyx_z = txy_z
  tzx_x = txz_x
  tzx_y = txz_y
  tzx_z = txz_z
  tzy_x = tyz_x
  tzy_y = tyz_y
  tzy_z = tyz_z

  # these coefficients are from nondimensionalization.
  Pr = 0.72
  c1 = params.Ma/params.Re
  c2 = c1/(Pr*gamma_1)
  src[2] -= c1*(txx_x + txy_y + txz_z)
  src[3] -= c1*(tyx_x + tyy_y + tyz_z)
  src[4] -= c1*(tzx_x + tzy_y + tzz_z)
  src[5] -= c1*(txx_x*q[2] + txx*q_x[1,2] + txy_x*q[3] + txy*q_x[1,3] + txz_x*q[4] + txz*q_x[1,4]) 
  src[5] -= c1*(tyx_y*q[2] + tyx*q_x[2,2] + tyy_y*q[3] + tyy*q_x[2,3] + tyz_y*q[4] + tyz*q_x[2,4]) 
  src[5] -= c1*(tzx_z*q[2] + tzx*q_x[3,2] + tzy_z*q[3] + tzy*q_x[3,3] + tzz_z*q[4] + tzz*q_x[3,4]) 
  src[5] -= c2*rK*(q_xx[1,1,5] + q_xx[2,2,5] + q_xx[3,3,5])

  return nothing
end

function calcMmsSource(params::ParamType{2, :conservative},
                       q::AbstractArray{Tsrc, 1},
                       q_x::AbstractArray{Tsrc, 2},
                       q_xx::AbstractArray{Tsrc, 3},
                       src::AbstractArray{Tsrc, 1}) where Tsrc
  gamma = params.gamma
  gamma_1 = params.gamma - 1

  p   = q[1]*q[4]/gamma
  E   = q[4]/(gamma*gamma_1) + 0.5*(q[2]*q[2] + q[3]*q[3])

  p_x   = 1.0/gamma*(q_x[1,1]*q[4] + q[1]*q_x[1,4])
  p_y   = 1.0/gamma*(q_x[2,1]*q[4] + q[1]*q_x[2,4])
  E_x   = q_x[1,4]/(gamma*gamma_1) + (q[2]*q_x[1,2] + q[3]*q_x[1,3])	
  E_y   = q_x[2,4]/(gamma*gamma_1) + (q[2]*q_x[2,2] + q[3]*q_x[2,3])

  #
  # contribution from inviscid terms
  #
  src[1]  = q_x[1,1]*q[2] + q[1]*q_x[1,2] 
  src[1] += q_x[2,1]*q[3] + q[1]*q_x[2,3] 

  src[2]  = q_x[1,1]*q[2]*q[2] + q[1]*q_x[1,2]*q[2] + q[1]*q[2]*q_x[1,2] + p_x
  src[2] += q_x[2,1]*q[2]*q[3] + q[1]*q_x[2,2]*q[3] + q[1]*q[2]*q_x[2,3]

  src[3]  = q_x[1,1]*q[3]*q[2] + q[1]*q_x[1,3]*q[2] + q[1]*q[3]*q_x[1,2]
  src[3] += q_x[2,1]*q[3]*q[3] + q[1]*q_x[2,3]*q[3] + q[1]*q[3]*q_x[2,3] + p_y

  src[4]  = q_x[1,1]*E*q[2] + q[1]*E_x*q[2] + q[1]*E*q_x[1,2] + p_x*q[2] + p*q_x[1,2]
  src[4] += q_x[2,1]*E*q[3] + q[1]*E_y*q[3] + q[1]*E*q_x[2,3] + p_y*q[3] + p*q_x[2,3]

  if !params.isViscous 
    return nothing
  end

  #
  # in case someone forgets the lower diagonal of q_xx
  #
  q_xx[2,1,2] = q_xx[1,2,2]
  q_xx[2,1,3] = q_xx[1,2,3]

  muK = Array(typeof(q[4]), 2)
  getMuK(q[4], muK)
  rmu = muK[1]
  rK = muK[2]
  two3rd = 2.0/3.0

  txx = rmu * two3rd * (2*q_x[1,2] - q_x[2,3])
  tyy = rmu * two3rd * (2*q_x[2,3] - q_x[1,2])
  txy = rmu * (q_x[2,2] + q_x[1,3]) 
  tyx = txy


  txx_x = rmu*two3rd*(2*q_xx[1,1,2] - q_xx[2,1,3])
  txx_y = rmu*two3rd*(2*q_xx[1,2,2] - q_xx[2,2,3])
  tyy_x = rmu*two3rd*(2*q_xx[2,1,3] - q_xx[1,1,2])
  tyy_y = rmu*two3rd*(2*q_xx[2,2,3] - q_xx[1,2,2])

  txy_x = rmu*(q_xx[2,1,2] + q_xx[1,1,3])
  txy_y = rmu*(q_xx[2,2,2] + q_xx[1,2,3])
  tyx_x = txy_x
  tyx_y = txy_y

  Pr = 0.72
  c1 = params.Ma/params.Re
  c2 = c1/(Pr*gamma_1)
  src[2] -= c1*(txx_x + txy_y)
  src[3] -= c1*(tyx_x + tyy_y)
  src[4] -= c1*(txx_x*q[2] + txx*q_x[1,2] + txy_x*q[3] + txy*q_x[1,3]) 
  src[4] -= c1*(tyx_y*q[2] + tyx*q_x[2,2] + tyy_y*q[3] + tyy*q_x[2,3]) 
  src[4] -= c2*rK*(q_xx[1,1,4] + q_xx[2,2,4])

  return nothing
end
#
#TODO: As far as we know `q`, `q_x`, `q_xx`, we can move
# the rest into a function, which is less error prone, and 
# will save a lot of space
#
mutable struct SRCPolynomial <: SRCType
end
function (obj::SRCPolynomial)(
              src::AbstractVector,
              xyz::AbstractVector, 
              params::ParamType{3}, 
              t)

  Tdim = 3
  sigma = 0.01
  gamma = params.gamma
  gamma_1 = params.gamma_1
  aoa = params.aoa
  beta = params.sideslip_angle
  q = zeros(typeof(src[1]), Tdim+2)
  q_x = zeros(typeof(src[1]), Tdim, Tdim+2)
  q_xx = zeros(typeof(src[1]), Tdim, Tdim, Tdim+2)
  qRef = zeros(typeof(src[1]), Tdim+2)
  qRef[1] = 1.0
  qRef[2] = params.Ma * cos(beta) * cos(aoa)
  qRef[3] = params.Ma * sin(beta) * -1
  qRef[4] = params.Ma * cos(beta) * sin(aoa)
  qRef[5] = 1.0

  x = xyz[1]
  y = xyz[2]
  z = xyz[3]

  q[1] = (x-x*x)*(y-y*y)* (z - z*z) 
  q[2] = (x-x*x)*(y-y*y)* (z - z*z) 
  q[3] = (x-x*x)*(y-y*y)* (z - z*z) 
  q[4] = (x-x*x)*(y-y*y)* (z - z*z) 
  q[5] = (x-x*x)*(y-y*y)* (z - z*z) 

  q_x[1,1] = (1.0 - 2.0*x) * (y - y*y) * (z - z*z)
  q_x[1,2] = (1.0 - 2.0*x) * (y - y*y) * (z - z*z)
  q_x[1,3] = (1.0 - 2.0*x) * (y - y*y) * (z - z*z)
  q_x[1,4] = (1.0 - 2.0*x) * (y - y*y) * (z - z*z)
  q_x[1,5] = (1.0 - 2.0*x) * (y - y*y) * (z - z*z)

  q_x[2,1] = (x - x*x) * (1.0 - 2.0*y) * (z - z*z)
  q_x[2,2] = (x - x*x) * (1.0 - 2.0*y) * (z - z*z)
  q_x[2,3] = (x - x*x) * (1.0 - 2.0*y) * (z - z*z)
  q_x[2,4] = (x - x*x) * (1.0 - 2.0*y) * (z - z*z)
  q_x[2,5] = (x - x*x) * (1.0 - 2.0*y) * (z - z*z)

  q_x[3,1] = (x - x*x) * (y - y*y) * (1.0 - 2.0*z) 
  q_x[3,2] = (x - x*x) * (y - y*y) * (1.0 - 2.0*z) 
  q_x[3,3] = (x - x*x) * (y - y*y) * (1.0 - 2.0*z) 
  q_x[3,4] = (x - x*x) * (y - y*y) * (1.0 - 2.0*z) 
  q_x[3,5] = (x - x*x) * (y - y*y) * (1.0 - 2.0*z) 

  q[1] = (sigma*q[1] + 1.0)*qRef[1]
  q[2] = (sigma*q[2] + 1.0)*qRef[2] 
  q[3] = (sigma*q[3] + 1.0)*qRef[3] 
  q[4] = (sigma*q[4] + 1.0)*qRef[4] 
  q[5] = (sigma*q[5] + 1.0)*qRef[5] 

  q_x[1,1] = q_x[1,1]*qRef[1]*sigma
  q_x[2,1] = q_x[2,1]*qRef[1]*sigma
  q_x[3,1] = q_x[3,1]*qRef[1]*sigma

  q_x[1,2] = q_x[1,2]*qRef[2]*sigma
  q_x[2,2] = q_x[2,2]*qRef[2]*sigma
  q_x[3,2] = q_x[3,2]*qRef[2]*sigma

  q_x[1,3] = q_x[1,3]*qRef[3]*sigma
  q_x[2,3] = q_x[2,3]*qRef[3]*sigma
  q_x[3,3] = q_x[3,3]*qRef[3]*sigma

  q_x[1,4] = q_x[1,4]*qRef[4]*sigma
  q_x[2,4] = q_x[2,4]*qRef[4]*sigma
  q_x[3,4] = q_x[3,4]*qRef[4]*sigma

  q_x[1,5] = q_x[1,5]*qRef[5]*sigma
  q_x[2,5] = q_x[2,5]*qRef[5]*sigma
  q_x[3,5] = q_x[3,5]*qRef[5]*sigma

  if !params.isViscous 
    calcMmsSource(params, q, q_x, q_xx, src)
    return nothing
  end
  #
  # contribution from viscous terms
  #
  q_xx[1,1,2] = -2.0 * (y - y*y) * (z - z*z) 
  q_xx[1,1,3] = -2.0 * (y - y*y) * (z - z*z) 
  q_xx[1,1,4] = -2.0 * (y - y*y) * (z - z*z) 
  q_xx[1,1,5] = -2.0 * (y - y*y) * (z - z*z) 

  q_xx[2,2,2] = -2.0 * (x - x*x) * (z - z*z)
  q_xx[2,2,3] = -2.0 * (x - x*x) * (z - z*z)
  q_xx[2,2,4] = -2.0 * (x - x*x) * (z - z*z)
  q_xx[2,2,5] = -2.0 * (x - x*x) * (z - z*z)

  q_xx[3,3,2] = -2.0 * (x - x*x) * (y - y*y)
  q_xx[3,3,3] = -2.0 * (x - x*x) * (y - y*y)
  q_xx[3,3,4] = -2.0 * (x - x*x) * (y - y*y)
  q_xx[3,3,5] = -2.0 * (x - x*x) * (y - y*y)

  q_xx[1,2,2] = (1.0 - 2.0*x) * (1.0 - 2.0*y) * (z - z*z)
  q_xx[1,2,3] = (1.0 - 2.0*x) * (1.0 - 2.0*y) * (z - z*z)
  q_xx[1,2,4] = (1.0 - 2.0*x) * (1.0 - 2.0*y) * (z - z*z)

  q_xx[1,3,2] = (1.0 - 2.0*x) * (1.0 - 2.0*z) * (y - y*y)
  q_xx[1,3,3] = (1.0 - 2.0*x) * (1.0 - 2.0*z) * (y - y*y)
  q_xx[1,3,4] = (1.0 - 2.0*x) * (1.0 - 2.0*z) * (y - y*y)

  q_xx[2,3,2] = (1.0 - 2.0*y) * (1.0 - 2.0*z) * (x - x*x)
  q_xx[2,3,3] = (1.0 - 2.0*y) * (1.0 - 2.0*z) * (x - x*x)
  q_xx[2,3,4] = (1.0 - 2.0*y) * (1.0 - 2.0*z) * (x - x*x)

  q_xx[1,1,2] *= sigma*qRef[2]
  q_xx[2,2,2] *= sigma*qRef[2]
  q_xx[3,3,2] *= sigma*qRef[2]
  q_xx[1,2,2] *= sigma*qRef[2]
  q_xx[1,3,2] *= sigma*qRef[2]
  q_xx[2,3,2] *= sigma*qRef[2]
  q_xx[2,1,2] = q_xx[1,2,2]
  q_xx[3,1,2] = q_xx[1,3,2]
  q_xx[3,2,2] = q_xx[2,3,2]

  q_xx[1,1,3] *= sigma*qRef[3]
  q_xx[2,2,3] *= sigma*qRef[3]
  q_xx[3,3,3] *= sigma*qRef[3]
  q_xx[1,2,3] *= sigma*qRef[3]
  q_xx[1,3,3] *= sigma*qRef[3]
  q_xx[2,3,3] *= sigma*qRef[3]
  q_xx[2,1,3] = q_xx[1,2,3]
  q_xx[3,1,3] = q_xx[1,3,3]
  q_xx[3,2,3] = q_xx[2,3,3]

  q_xx[1,1,4] *= sigma*qRef[4]
  q_xx[2,2,4] *= sigma*qRef[4]
  q_xx[3,3,4] *= sigma*qRef[4]
  q_xx[1,2,4] *= sigma*qRef[4]
  q_xx[1,3,4] *= sigma*qRef[4]
  q_xx[2,3,4] *= sigma*qRef[4]
  q_xx[2,1,4] = q_xx[1,2,4]
  q_xx[3,1,4] = q_xx[1,3,4]
  q_xx[3,2,4] = q_xx[2,3,4]

  q_xx[1,1,5] *= sigma*qRef[5]
  q_xx[2,2,5] *= sigma*qRef[5]
  q_xx[3,3,5] *= sigma*qRef[5]

  calcMmsSource(params, q, q_x, q_xx, src)

  return nothing
end

function (obj::SRCPolynomial)(
              src::AbstractVector,
              coords::AbstractVector, 
              params::ParamType{2}, 
              t)
  Tdim = 2
  sigma = 0.01
  gamma = params.gamma
  gamma_1 = params.gamma_1
  aoa = params.aoa
  q = zeros(typeof(src[1]), Tdim+2)
  q_x = zeros(typeof(src[1]), Tdim, Tdim+2)
  q_xx = zeros(typeof(src[1]), Tdim, Tdim, Tdim+2)
  qRef = zeros(typeof(src[1]), Tdim+2)
  qRef[1] = 1.0
  qRef[2] = params.Ma*cos(aoa)
  qRef[3] = params.Ma*sin(aoa)
  qRef[4] = 1.0
  x = coords[1]
  y = coords[2]

  q[1] = (x-x*x)*(y-y*y) 
  q[2] = (x-x*x)*(y-y*y)
  q[3] = (x-x*x)*(y-y*y)
  q[4] = (x-x*x)*(y-y*y)
  q_x[1,1] = (1.0 - 2.0*x) * (y - y*y) 
  q_x[1,2] = (1.0 - 2.0*x) * (y - y*y)
  q_x[1,3] = (1.0 - 2.0*x) * (y - y*y)
  q_x[1,4] = (1.0 - 2.0*x) * (y - y*y)
  q_x[2,1] = (x - x*x) * (1.0 - 2.0*y)
  q_x[2,2] = (x - x*x) * (1.0 - 2.0*y)
  q_x[2,3] = (x - x*x) * (1.0 - 2.0*y)
  q_x[2,4] = (x - x*x) * (1.0 - 2.0*y)

  q[1] = (sigma*q[1] + 1.0)*qRef[1] 
  q[2] = (sigma*q[2] + 1.0)*qRef[2]
  q[3] = (sigma*q[3] + 1.0)*qRef[3]
  q[4] = (sigma*q[4] + 1.0)*qRef[4]
  q_x[1,1] = q_x[1,1]*qRef[1]*sigma
  q_x[2,1] = q_x[2,1]*qRef[1]*sigma
  q_x[1,2] = q_x[1,2]*qRef[2]*sigma
  q_x[2,2] = q_x[2,2]*qRef[2]*sigma
  q_x[1,3] = q_x[1,3]*qRef[3]*sigma
  q_x[2,3] = q_x[2,3]*qRef[3]*sigma
  q_x[1,4] = q_x[1,4]*qRef[4]*sigma
  q_x[2,4] = q_x[2,4]*qRef[4]*sigma

  if !params.isViscous 
    calcMmsSource(params, q, q_x, q_xx, src)
    return nothing
  end
  #
  # contribution from viscous terms
  #
  q_xx[1,1,2] = -2.0 * (y - y*y)
  q_xx[2,2,2] = -2.0 * (x - x*x)
  q_xx[1,2,2] = (1.0 - 2.0*x) * (1.0 - 2.0*y) 

  q_xx[1,1,3] = -2.0 * (y - y*y) 
  q_xx[2,2,3] = -2.0 * (x - x*x)
  q_xx[1,2,3] = (1.0 - 2.0*x) * (1.0 - 2.0*y)

  q_xx[1,1,4] = -2.0 * (y - y*y) 
  q_xx[2,2,4] = -2.0 * (x - x*x)

  q_xx[1,1,2] *= sigma*qRef[2]
  q_xx[1,2,2] *= sigma*qRef[2]
  q_xx[2,2,2] *= sigma*qRef[2]

  q_xx[1,1,3] *= sigma*qRef[3]
  q_xx[1,2,3] *= sigma*qRef[3]
  q_xx[2,2,3] *= sigma*qRef[3]

  q_xx[1,1,4] *= sigma*qRef[4]
  q_xx[2,2,4] *= sigma*qRef[4]

  q_xx[2,1,2] = q_xx[1,2,2]
  q_xx[2,1,3] = q_xx[1,2,3]

  calcMmsSource(params, q, q_x, q_xx, src)

  return nothing
end

mutable struct SRCChannel <: SRCType
end
function (obj::SRCChannel)(
              src::AbstractVector,
              xyz::AbstractVector, 
              params::ParamType{3}, 
              t)
  Tdim = 3
  sigma = 0.01
  gamma = params.gamma
  gamma_1 = params.gamma_1
  aoa = params.aoa
  beta = params.sideslip_angle
  qRef = zeros(typeof(src[1]), Tdim+2)
  q = zeros(typeof(src[1]), Tdim+2)
  q_x = zeros(typeof(src[1]), Tdim, Tdim+2)
  q_xx = zeros(typeof(src[1]), Tdim, Tdim, Tdim+2)
  qRef[1] = 1.0
  qRef[2] = params.Ma * cos(beta) * cos(aoa)
  qRef[3] = params.Ma * sin(beta) * -1
  qRef[4] = params.Ma * cos(beta) * sin(aoa)
  qRef[5] = 1.0

  x = xyz[1]
  y = xyz[2]
  z = xyz[3]

  q[1] = (1 + sigma *x*y*z) * qRef[1]  
  ux = sin(pi*x) + 1
  uy = sin(pi*y) + 1
  uz = sin(pi*z) + 1
  q[2]  = (1 + sigma * ux * uy * uz) * qRef[2]
  vx = sin(pi*x) + 1
  vy = sin(pi*y) + 1
  vz = sin(pi*z) + 1
  q[3]  = (1 + sigma * ux * uy * uz) * qRef[3]
  wx = sin(pi*x) + 1
  wy = sin(pi*y) + 1
  wz = sin(pi*z) + 1
  q[4]  = wx * wy * wz * qRef[4]
  q[4]  = (1 + sigma * ux * uy * uz) * qRef[4]
  q[5]  = qRef[5] 


  q_x[1,1] = qRef[1] * y*z*sigma 
  q_x[2,1] = qRef[1] * x*z*sigma 
  q_x[3,1] = qRef[1] * x*y*sigma 

  q_x[1,2] = qRef[2] * pi * cos(pi*x) * uy * uz * sigma
  q_x[2,2] = qRef[2] * pi * cos(pi*y) * ux * uz * sigma
  q_x[3,2] = qRef[2] * pi * cos(pi*z) * ux * uy * sigma

  q_x[1,3] = qRef[3] * pi * cos(pi*x) * uy * uz * sigma
  q_x[2,3] = qRef[3] * pi * cos(pi*y) * ux * uz * sigma
  q_x[3,3] = qRef[3] * pi * cos(pi*z) * ux * uy * sigma

  q_x[1,4] = qRef[4] * pi * cos(pi*x) * uy * uz * sigma
  q_x[2,4] = qRef[4] * pi * cos(pi*y) * ux * uz * sigma
  q_x[3,4] = qRef[4] * pi * cos(pi*z) * ux * uy * sigma

  q_x[1,5] = 0
  q_x[2,5] = 0
  q_x[3,5] = 0

  if !params.isViscous 
    q[2] += 0.2 * qRef[2] 
    q[3] += 0.2 * qRef[3] 
    q[4] += 0.2 * qRef[4] 
    calcMmsSource(params, q, q_x, q_xx, src)
    return nothing
  end
  #
  # contribution from viscous terms
  #
  q_xx[1,1,2] = -sigma*qRef[2]*pi*pi * sin(pi*x) * uy * uz
  q_xx[2,2,2] = -sigma*qRef[2]*pi*pi * sin(pi*y) * ux * uz
  q_xx[3,3,2] = -sigma*qRef[2]*pi*pi * sin(pi*z) * ux * uy
  q_xx[1,2,2] =  sigma*qRef[2]*pi*pi * cos(pi*x) * cos(pi*y) * uz
  q_xx[1,3,2] =  sigma*qRef[2]*pi*pi * cos(pi*x) * cos(pi*z) * uy
  q_xx[2,3,2] =  sigma*qRef[2]*pi*pi * cos(pi*y) * cos(pi*z) * ux
  q_xx[2,1,2] =  q_xx[1,2,2]
  q_xx[3,1,2] =  q_xx[1,3,2]
  q_xx[3,2,2] =  q_xx[2,3,2]

  q_xx[1,1,3] = -sigma*qRef[3]*pi*pi * sin(pi*x) * uy * uz
  q_xx[2,2,3] = -sigma*qRef[3]*pi*pi * sin(pi*y) * ux * uz
  q_xx[3,3,3] = -sigma*qRef[3]*pi*pi * sin(pi*z) * ux * uy
  q_xx[1,2,3] =  sigma*qRef[3]*pi*pi * cos(pi*x) * cos(pi*y) * uz
  q_xx[1,3,3] =  sigma*qRef[3]*pi*pi * cos(pi*x) * cos(pi*z) * uy
  q_xx[2,3,3] =  sigma*qRef[3]*pi*pi * cos(pi*y) * cos(pi*z) * ux
  q_xx[2,1,3] =  q_xx[1,2,3]
  q_xx[3,1,3] =  q_xx[1,3,3]
  q_xx[3,2,3] =  q_xx[2,3,3]

  q_xx[1,1,4] = -sigma*qRef[4]*pi*pi * sin(pi*x) * uy * uz
  q_xx[2,2,4] = -sigma*qRef[4]*pi*pi * sin(pi*y) * ux * uz
  q_xx[3,3,4] = -sigma*qRef[4]*pi*pi * sin(pi*z) * ux * uy
  q_xx[1,2,4] =  sigma*qRef[4]*pi*pi * cos(pi*x) * cos(pi*y) * uz
  q_xx[1,3,4] =  sigma*qRef[4]*pi*pi * cos(pi*x) * cos(pi*z) * uy
  q_xx[2,3,4] =  sigma*qRef[4]*pi*pi * cos(pi*y) * cos(pi*z) * ux
  q_xx[2,1,4] =  q_xx[1,2,4]
  q_xx[3,1,4] =  q_xx[1,3,4]
  q_xx[3,2,4] =  q_xx[2,3,4]

  q_xx[1,1,5] = 0
  q_xx[2,2,5] = 0
  q_xx[3,3,5] = 0

  calcMmsSource(params, q, q_x, q_xx, src)
  return nothing
end

function (obj::SRCChannel)(
              src::AbstractVector,
              coords::AbstractVector, 
              params::ParamType{2}, 
              t)
  Tdim = 2
  pi = 3.14159265358979323846264338
  sigma = 0.1
  gamma = params.gamma
  gamma_1 = params.gamma_1
  aoa = params.aoa
  q = zeros(typeof(src[1]), Tdim+2)
  q_x = zeros(typeof(src[1]), Tdim, Tdim+2)
  q_xx = zeros(typeof(src[1]), Tdim, Tdim, Tdim+2)
  qRef = zeros(typeof(src[1]), Tdim+2)
  qRef[1] = 1.0
  qRef[2] = params.Ma*cos(aoa)
  qRef[3] = params.Ma*sin(aoa)
  qRef[4] = 1.0
  x = coords[1]
  y = coords[2]

  # Exact solution in form of primitive variables

  q[1] = qRef[1] * (exp(sin(0.5*pi*(x+y))) *sigma +  1.0)
  ux  = (exp(x) * sin(pi*x) * sigma + 1) * qRef[2]
  uy   = exp(y) * sin(pi*y) 
  q[2] = ux * uy
  vx  = (exp(x) * sin(pi*x) * sigma + 1) * qRef[3]
  vy   = exp(y) * sin(pi*y) 
  q[3] = vx * vy
  q[4] = (1 + sigma* exp(0.1*x+0.1*y)) * qRef[4]

  #
  # contribution from inviscid terms
  #
  q_x[1,1] = qRef[1] * exp(sin(0.5*pi*(x+y))) * 0.5*pi*cos(0.5*pi*(x+y)) * sigma
  q_x[2,1] = qRef[1] * exp(sin(0.5*pi*(x+y))) * 0.5*pi*cos(0.5*pi*(x+y)) * sigma
  ux_x = qRef[2] * exp(x) * (sin(pi*x) + pi* cos(pi*x)) * sigma
  uy_y = exp(y) * (sin(pi*y) + pi * cos(pi*y))
  q_x[1,2] = ux_x * uy 
  q_x[2,2] = ux * uy_y
  vx_x = qRef[3] * exp(x) * (sin(pi*x) + pi* cos(pi*x)) * sigma
  vy_y = exp(y) * (sin(pi*y) + pi * cos(pi*y))
  q_x[1,3] = vx_x * vy 
  q_x[2,3] = vx * vy_y
  q_x[1,4] = exp(0.1*x+0.1*y) * qRef[4] * sigma * 0.1
  q_x[2,4] = exp(0.1*x+0.1*y) * qRef[4] * sigma * 0.1
  # q_x[1,4] = 0.0
  # q_x[2,4] = 0.0

  if !params.isViscous 
    q[2] += qRef[2] * 0.2
    calcMmsSource(params, q, q_x, q_xx, src)
    return nothing
  end

  ux_xx  = exp(x) * (sin(pi*x) + pi*cos(pi*x))
  ux_xx += exp(x) * (cos(pi*x) - pi*sin(pi*x)) * pi 
  ux_xx *= sigma * qRef[2]
  uy_yy  = exp(y) * (sin(pi*y) + pi*cos(pi*y))
  uy_yy += exp(y) * (cos(pi*y) - pi*sin(pi*y)) * pi
  q_xx[1,1,2] = ux_xx * uy
  q_xx[1,2,2] = ux_x * uy_y
  q_xx[2,2,2] = ux * uy_yy
  q_xx[2,1,2] = q_xx[1,2,2]

  vx_xx  = exp(x) * (sin(pi*x) + pi*cos(pi*x))
  vx_xx += exp(x) * (cos(pi*x) - pi*sin(pi*x)) * pi
  vx_xx *= sigma * qRef[3]
  vy_yy  = exp(y) * (sin(pi*y) + pi*cos(pi*y))
  vy_yy += exp(y) * (cos(pi*y) - pi*sin(pi*y)) * pi
  q_xx[1,1,3] = vx_xx * vy
  q_xx[1,2,3] = vx_x * vy_y
  q_xx[2,2,3] = vx * vy_yy
  q_xx[2,1,3] = q_xx[1,2,3]

  q_xx[1,1,4] = exp(0.1*x+0.1*y) * qRef[4] * sigma * 0.01
  q_xx[2,2,4] = exp(0.1*x+0.1*y) * qRef[4] * sigma * 0.01
  # q_xx[1,1,4] = 0.0
  # q_xx[2,2,4] = 0.0

  calcMmsSource(params, q, q_x, q_xx, src)

  return nothing
end

#######################
#  ___________
# |           |
# |    ---    |
# |   |   |   |
# |    ---    |
# |           |
#  -----------
#
#######################
# inner block = nonslip wall
# outer block = freestream

# This MMS is not working very well. 
# (I don't remember clearly, but it's not converging).
mutable struct SRCDoubleSquare <: SRCType
end
function (obj::SRCDoubleSquare)(
              src::AbstractVector,
              coords::AbstractVector, 
              params::ParamType{2}, 
              t)
  Tdim = 2
  pi = 3.14159265358979323846264338
  sigma = 0.01
  gamma = params.gamma
  gamma_1 = params.gamma - 1
  aoa = params.aoa
  q = zeros(typeof(src[1]), Tdim+2)
  q_x = zeros(typeof(src[1]), Tdim, Tdim+2)
  q_xx = zeros(typeof(src[1]), Tdim, Tdim, Tdim+2)
  qRef = zeros(typeof(src[1]), Tdim+2)
  qRef[1] = 1.0
  qRef[2] = params.Ma*cos(aoa)
  qRef[3] = params.Ma*sin(aoa)
  qRef[4] = 1.0
  x = coords[1]
  y = coords[2]
  #
  # Exact solution in form of primitive variables
  #
  si = [0.5, 1.5]
  a = [-2.375, 16.875, -45.0, 55.0, -30.0, 6.0]
  gx    = 0.0
  gy    = 0.0
  gx_x  = 0.0
  gy_y  = 0.0
  gx_xx = 0.0
  gy_yy = 0.0

  if x >= si[1] && x < si[2]
    gx = a[1] + a[2]*x + a[3]*x*x + a[4]*x^3 + a[5]*x^4 + a[6]*x^5
    gx_x = a[2] + 2*a[3]*x + 3*a[4]*x*x + 4*a[5]*x*x*x + 5*a[6]*x^4
    gx_xx = 2*a[3] + 6*a[4]*x + 12*a[5]*x*x + 20*a[6]*x*x*x 
  elseif x >= si[2] 
    gx = 1.0
  elseif x <= -si[1] && x > -si[2]
    gx = a[1] - a[2]*x + a[3]*x*x - a[4]*x^3 + a[5]*x^4 - a[6]*x^5
    gx_x = -a[2] + 2*a[3]*x - 3*a[4]*x*x + 4*a[5]*x*x*x - 5*a[6]*x^4
    gx_xx = 2*a[3] - 6*a[4]*x + 12*a[5]*x*x - 20*a[6]*x*x*x 
  elseif x <= -si[2]
    gx = 1.0
  end  
  if y >= si[1] && y < si[2]
    gy = a[1] + a[2]*y + a[3]*y*y + a[4]*y^3 + a[5]*y^4 + a[6]*y^5
    gy_y = a[2] + 2*a[3]*y + 3*a[4]*y*y + 4*a[5]*y*y*y + 5*a[6]*y^4
    gy_yy = 2*a[3] + 6*a[4]*y + 12*a[5]*y*y + 20*a[6]*y*y*y 
  elseif y >= si[2] 
    gy = 1.0
  elseif y <= -si[1] && y > -si[2]
    gy = a[1] - a[2]*y + a[3]*y*y - a[4]*y^3 + a[5]*y^4 - a[6]*y^5
    gy_y = -a[2] + 2*a[3]*y - 3*a[4]*y*y + 4*a[5]*y*y*y - 5*a[6]*y^4
    gy_yy = 2*a[3] - 6*a[4]*y + 12*a[5]*y*y - 20*a[6]*y*y*y 
  elseif y <= -si[2]
    gy = 1.0
  end  

  q[1] = qRef[1]
  q[2] = qRef[2] * (gx + gy - gx*gy) 
  q[3] = qRef[3] * (gx + gy - gx*gy) 
  q[4] = qRef[4]

  q_x[1,1] = 0.0
  q_x[2,1] = 0.0
  q_x[1,2] = qRef[2] * (gx_x - gx_x * gy) 
  q_x[2,2] = qRef[2] * (gy_y - gx * gy_y) 
  q_x[1,3] = qRef[3] * (gx_x - gx_x * gy) 
  q_x[2,3] = qRef[3] * (gy_y - gx * gy_y) 
  q_x[1,4] = 0.0 
  q_x[2,4] = 0.0 

  if !params.isViscous 
    calcMmsSource(params, q, q_x, q_xx, src)
    return nothing
  end

  q_xx[1,1,2] = qRef[2] * (gx_xx - gx_xx * gy) 
  q_xx[1,2,2] = qRef[2] * (-gx_x * gy_y) 
  q_xx[2,2,2] = qRef[2] * (gy_yy - gx * gy_yy)
  q_xx[2,1,2] = q_xx[1,2,2]
  q_xx[1,1,3] = qRef[3] * (gx_xx - gx_xx * gy) 
  q_xx[1,2,3] = qRef[3] * (-gx_x * gy_y) 
  q_xx[2,2,3] = qRef[3] * (gy_yy - gx * gy_yy)
  q_xx[2,1,3] = q_xx[1,2,3]
  q_xx[1,1,4] = 0.0
  q_xx[2,2,4] = 0.0

  calcMmsSource(params, q, q_x, q_xx, src)

  return nothing
end


mutable struct SRCTrigonometric <: SRCType
end
function (obj::SRCTrigonometric)(
              src::AbstractVector,
              xyz::AbstractVector, 
              params::ParamType{3}, 
              t)
  Tdim = 3
  sigma = 0.0001
  gamma = params.gamma
  gamma_1 = params.gamma - 1
  aoa = params.aoa
  beta = params.sideslip_angle
  qRef = zeros(typeof(src[1]), Tdim+2)
  q    = zeros(typeof(src[1]), Tdim+2)
  q_x  = zeros(typeof(src[1]), Tdim, Tdim+2)
  q_xx = zeros(typeof(src[1]), Tdim, Tdim, Tdim+2)
  qRef[1] = 1.0
  qRef[2] = params.Ma * cos(beta) * cos(aoa)
  qRef[3] = params.Ma * sin(beta) * -1
  qRef[4] = params.Ma * cos(beta) * sin(aoa)
  qRef[5] = 1.0

  xyz1 = 1 * pi * xyz
  xyz2 = 2 * pi * xyz
  xyz4 = 4 * pi * xyz
  sin_val_1 = sin(xyz1)
  cos_val_1 = cos(xyz1)
  sin_val_2 = sin(xyz2)
  cos_val_2 = cos(xyz2)
  sin_val_4 = sin(xyz4)
  cos_val_4 = cos(xyz4)

  q[1] = sin_val_2[1] * sin_val_2[2] * sin_val_2[3] 
  q[2] = sin_val_4[1] * sin_val_4[2] * sin_val_4[3]
  q[3] = sin_val_2[1] * sin_val_2[2] * sin_val_2[3]
  q[4] = sin_val_1[1] * sin_val_1[2] * sin_val_1[3] 
  q[5] = (1.0 - cos_val_4[1]) * (1.0 - cos_val_4[2]) * (1.0 - cos_val_4[3])

  q_x[1,1] = 2*pi * cos_val_2[1] * sin_val_2[2] * sin_val_2[3] 
  q_x[1,2] = 4*pi * cos_val_4[1] * sin_val_4[2] * sin_val_4[3] 
  q_x[1,3] = 2*pi * cos_val_2[1] * sin_val_2[2] * sin_val_2[3] 
  q_x[1,4] =   pi * cos_val_1[1] * sin_val_1[2] * sin_val_1[3] 
  q_x[1,5] = 4*pi * sin_val_4[1] * (1.0 - cos_val_4[2]) * (1.0 - cos_val_4[3])  

  q_x[2,1] = 2*pi * cos_val_2[2] * sin_val_2[1] * sin_val_2[3] 
  q_x[2,2] = 4*pi * cos_val_4[2] * sin_val_4[1] * sin_val_4[3] 
  q_x[2,3] = 2*pi * cos_val_2[2] * sin_val_2[1] * sin_val_2[3] 
  q_x[2,4] =   pi * cos_val_1[2] * sin_val_1[1] * sin_val_1[3] 
  q_x[2,5] = 4*pi * sin_val_4[2] * (1.0 - cos_val_4[1]) * (1.0 - cos_val_4[3])

  q_x[3,1] = 2*pi * cos_val_2[3] * sin_val_2[1] * sin_val_2[2] 
  q_x[3,2] = 4*pi * cos_val_4[3] * sin_val_4[1] * sin_val_4[2] 
  q_x[3,3] = 2*pi * cos_val_2[3] * sin_val_2[1] * sin_val_2[2] 
  q_x[3,4] =    pi* cos_val_1[3] * sin_val_1[1] * sin_val_1[2] 
  q_x[3,5] = 4*pi * sin_val_4[3] * (1.0 - cos_val_4[1]) * (1.0 - cos_val_4[2])


  q[1] = (sigma*q[1] + 1.0)*qRef[1] 
  q[2] = (sigma*q[2] + 1.0)*qRef[2]
  q[3] = (sigma*q[3] + 1.0)*qRef[3]
  q[4] = (sigma*q[4] + 1.0)*qRef[4]
  q[5] = (sigma*q[5] + 1.0)*qRef[5]

  q_x[1,1] = q_x[1,1]*qRef[1]*sigma
  q_x[2,1] = q_x[2,1]*qRef[1]*sigma
  q_x[3,1] = q_x[3,1]*qRef[1]*sigma

  q_x[1,2] = q_x[1,2]*qRef[2]*sigma
  q_x[2,2] = q_x[2,2]*qRef[2]*sigma
  q_x[3,2] = q_x[3,2]*qRef[2]*sigma

  q_x[1,3] = q_x[1,3]*qRef[3]*sigma
  q_x[2,3] = q_x[2,3]*qRef[3]*sigma
  q_x[3,3] = q_x[3,3]*qRef[3]*sigma

  q_x[1,4] = q_x[1,4]*qRef[4]*sigma
  q_x[2,4] = q_x[2,4]*qRef[4]*sigma
  q_x[3,4] = q_x[3,4]*qRef[4]*sigma

  q_x[1,5] = q_x[1,5]*qRef[5]*sigma
  q_x[2,5] = q_x[2,5]*qRef[5]*sigma
  q_x[3,5] = q_x[3,5]*qRef[5]*sigma

  if !params.isViscous 
    calcMmsSource(params, q, q_x, q_xx, src)
    return nothing
  end

  q_xx[1,1,2] = -16*pi*pi * sin_val_4[1] * sin_val_4[2] * sin_val_4[3]
  q_xx[1,1,3] = - 4*pi*pi * sin_val_2[1] * sin_val_2[2] * sin_val_2[3] 
  q_xx[1,1,4] =    -pi*pi * sin_val_1[1] * sin_val_1[2] * sin_val_1[3] 
  q_xx[1,1,5] =  16*pi*pi * cos_val_4[1] * (1.0 - cos_val_4[2]) * (1.0 - cos_val_4[3])

  q_xx[2,2,2] = -16*pi*pi * sin_val_4[2] * sin_val_4[1] * sin_val_4[3]
  q_xx[2,2,3] =  -4*pi*pi * sin_val_2[2] * sin_val_2[1] * sin_val_2[3] 
  q_xx[2,2,4] =    -pi*pi * sin_val_1[2] * sin_val_1[1] * sin_val_1[3] 
  q_xx[2,2,5] =  16*pi*pi * cos_val_4[2] * (1.0 - cos_val_4[1]) * (1.0 - cos_val_4[3])

  q_xx[3,3,2] = -16*pi*pi * sin_val_4[3] * sin_val_4[1] * sin_val_4[2]
  q_xx[3,3,3] =  -4*pi*pi * sin_val_2[3] * sin_val_2[1] * sin_val_2[2] 
  q_xx[3,3,4] =    -pi*pi * sin_val_1[3] * sin_val_1[1] * sin_val_1[2] 
  q_xx[3,3,5] =  16*pi*pi * cos_val_4[3] * (1.0 - cos_val_4[1]) * (1.0 - cos_val_4[2])

  q_xx[1,2,2] = 16*pi*pi * cos_val_4[1] * cos_val_4[2] * sin_val_4[3] 
  q_xx[1,2,3] =  4*pi*pi * cos_val_2[1] * cos_val_2[2] * sin_val_2[3]
  q_xx[1,2,4] =    pi*pi * cos_val_1[1] * cos_val_1[2] * sin_val_1[3]

  q_xx[1,3,2] = 16*pi*pi * cos_val_4[1] * cos_val_4[3] * sin_val_4[2] 
  q_xx[1,3,3] =  4*pi*pi * cos_val_2[1] * cos_val_2[3] * sin_val_2[2]
  q_xx[1,3,4] =    pi*pi * cos_val_1[1] * cos_val_1[3] * sin_val_1[2]

  q_xx[2,3,2] = 16*pi*pi * cos_val_4[2] * cos_val_4[3] * sin_val_4[1] 
  q_xx[2,3,3] =  4*pi*pi * cos_val_2[2] * cos_val_2[3] * sin_val_2[1]
  q_xx[2,3,4] =    pi*pi * cos_val_1[2] * cos_val_1[3] * sin_val_1[1]

  q_xx[1,1,2] *= sigma*qRef[2]
  q_xx[2,2,2] *= sigma*qRef[2]
  q_xx[3,3,2] *= sigma*qRef[2]
  q_xx[1,2,2] *= sigma*qRef[2]
  q_xx[1,3,2] *= sigma*qRef[2]
  q_xx[2,3,2] *= sigma*qRef[2]
  q_xx[2,1,2] = q_xx[1,2,2]
  q_xx[3,1,2] = q_xx[1,3,2]
  q_xx[3,2,2] = q_xx[2,3,2]

  q_xx[1,1,3] *= sigma*qRef[3]
  q_xx[2,2,3] *= sigma*qRef[3]
  q_xx[3,3,3] *= sigma*qRef[3]
  q_xx[1,2,3] *= sigma*qRef[3]
  q_xx[1,3,3] *= sigma*qRef[3]
  q_xx[2,3,3] *= sigma*qRef[3]
  q_xx[2,1,3] = q_xx[1,2,3]
  q_xx[3,1,3] = q_xx[1,3,3]
  q_xx[3,2,3] = q_xx[2,3,3]

  q_xx[1,1,4] *= sigma*qRef[4]
  q_xx[2,2,4] *= sigma*qRef[4]
  q_xx[3,3,4] *= sigma*qRef[4]
  q_xx[1,2,4] *= sigma*qRef[4]
  q_xx[1,3,4] *= sigma*qRef[4]
  q_xx[2,3,4] *= sigma*qRef[4]
  q_xx[2,1,4] = q_xx[1,2,4]
  q_xx[3,1,4] = q_xx[1,3,4]
  q_xx[3,2,4] = q_xx[2,3,4]

  q_xx[1,1,5] *= sigma*qRef[5]
  q_xx[2,2,5] *= sigma*qRef[5]
  q_xx[3,3,5] *= sigma*qRef[5]

  calcMmsSource(params, q, q_x, q_xx, src)

  return nothing
end

function (obj::SRCTrigonometric)(
              src::AbstractVector,
              coords::AbstractVector, 
              params::ParamType{2}, 
              t)
  Tdim = 2
  sigma = 0.01
  gamma = params.gamma
  gamma_1 = gamma - 1.0
  aoa = params.aoa
  q    = zeros(typeof(src[1]), Tdim+2)
  q_x  = zeros(typeof(src[1]), Tdim, Tdim+2)
  q_xx = zeros(typeof(src[1]), Tdim, Tdim, Tdim+2)
  qRef = zeros(typeof(src[1]), Tdim+2)
  qRef[1] = 1.0
  qRef[2] = params.Ma*cos(aoa)
  qRef[3] = params.Ma*sin(aoa)
  qRef[4] = 1.0
  x = coords[1]
  y = coords[2]
  x2 = 2*x*pi
  y2 = 2*y*pi
  x3 = 3*x*pi
  y3 = 3*y*pi
  x4 = 4*x*pi
  y4 = 4*y*pi
  sx2 = sin(x2)
  sy2 = sin(y2)
  sx3 = sin(x3)
  sy3 = sin(y3)
  sx4 = sin(x4)
  sy4 = sin(y4)
  cx2 = cos(x2)
  cx3 = cos(x3)
  cx4 = cos(x4)
  cy2 = cos(y2)
  cy3 = cos(y3)
  cy4 = cos(y4)
  #
  # Exact solution in form of primitive variables
  #
  q[1] = sx2 * sy2
  q[2] = sx4 * sy4
  q[3] = sx3 * sy3
  q[4] = (1.0 - cx2) * (1.0 - cy2)
  q_x[1,1] = 2*pi * cx2 * sy2
  q_x[2,1] = 2*pi * sx2 * cy2
  q_x[1,2] = 4*pi* cx4 * sy4
  q_x[2,2] = 4*pi* sx4 * cy4
  q_x[1,3] = 3*pi* cx3 * sy3
  q_x[2,3] = 3*pi* cy3 * sx3
  q_x[1,4] = 2*pi* sx2 * (1.0 - cy2)
  q_x[2,4] = 2*pi* sy2 * (1.0 - cx2) 

  q[1] = (sigma*q[1] + 1.0)*qRef[1] 
  q[2] = (sigma*q[2] + 1.0)*qRef[2]
  q[3] = (sigma*q[3] + 1.0)*qRef[3]
  q[4] = (sigma*q[4] + 1.0)*qRef[4]

  q_x[1,1] = q_x[1,1]*qRef[1]*sigma
  q_x[2,1] = q_x[2,1]*qRef[1]*sigma
  q_x[1,2] = q_x[1,2]*qRef[2]*sigma
  q_x[2,2] = q_x[2,2]*qRef[2]*sigma
  q_x[1,3] = q_x[1,3]*qRef[3]*sigma
  q_x[2,3] = q_x[2,3]*qRef[3]*sigma
  q_x[1,4] = q_x[1,4]*qRef[4]*sigma
  q_x[2,4] = q_x[2,4]*qRef[4]*sigma

  if !params.isViscous 
    calcMmsSource(params, q, q_x, q_xx, src)
    return nothing
  end

  q_xx[1,1,2] = -16*pi*pi * sx4 * sy4
  q_xx[2,2,2] = -16*pi*pi * sx4 * sy4
  q_xx[1,2,2] =  16*pi*pi * cx4 * cy4
  q_xx[1,1,3] = -9*pi*pi * sx3 * sy3
  q_xx[2,2,3] = -9*pi*pi * sx3 * sy3
  q_xx[1,2,3] =  9*pi*pi * cx3 * cy3
  q_xx[1,1,4] =  4*pi*pi * cx2 * (1.0 - cy2)
  q_xx[2,2,4] =  4*pi*pi * cy2 * (1.0 - cx2)
  q_xx[1,1,2] = q_xx[1,1,2] * qRef[2] * sigma
  q_xx[1,2,2] = q_xx[1,2,2] * qRef[2] * sigma
  q_xx[2,2,2] = q_xx[2,2,2] * qRef[2] * sigma
  q_xx[2,1,2] = q_xx[1,2,2]
  q_xx[1,1,3] = q_xx[1,1,3] * qRef[3] * sigma
  q_xx[1,2,3] = q_xx[1,2,3] * qRef[3] * sigma
  q_xx[2,2,3] = q_xx[2,2,3] * qRef[3] * sigma
  q_xx[2,1,3] = q_xx[1,2,3]
  q_xx[1,1,4] = q_xx[1,1,4] * qRef[4] * sigma
  q_xx[2,2,4] = q_xx[2,2,4] * qRef[4] * sigma

  calcMmsSource(params, q, q_x, q_xx, src)
  return nothing
end
