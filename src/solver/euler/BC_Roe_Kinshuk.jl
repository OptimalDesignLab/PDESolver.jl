# Euler roe solver
function getEulerFlux(u_vals::AbstractArray, nx, ny)
  include("./euler.jl")  # solver functions
  eqn = EulerEquation(sbp)

  pressure = calcPressure(u_vals, eqn)
  f1[1] = u_vals[2]
  f1[2] = (u_vals[2]^2)/u_vals[1] + pressure
  f1[3] = (u_vals[2]*u_vals[3])/u_vals[1]
  f1[4] = (u_vals[4] + pressure)*u_vals[2]/u_vals[1]

  f2[1] = u_vals[3]
  f2[2] = (u_vals[2]*u_vals[3])/u_vals[1]
  f2[3] = (u_vals[3]^2)/u_vals[1] + pressure
  f2[4] = (u_vals[4] + pressure)*u_vals[3]/u_vals[1]

  eulerflux = f1*nx + f2*ny
  return eulerflux
end

function isentropicVortexBC(q, x, dxidx, nrm)
  #sgn::Int, tau::double, dxidx::AbstractArray{T,4}, q::Float64, qg::Float64, sat::AbstractArray{T,1})
  
  # Variable Specification
  # SAT : simultaneous approximation term

  # Call jared calcisentropicvotex get qexact
  # 

  #= nx::Float64; ny::Float64; nz::Float64; dA::Float64; uL::Float64;
  vL::Float64; wL::Float64; HL::Float64; uR::Float64; vR::Float64;
  wR::Float64; HR::Float64; sqL::Float64; sqR::Float64; phi::Float64;
  u::Float64; v::Float64; w::Float64; H::Float64; a::Float64; Un::Float64
  lambda1::Float64; lambda2::Float64; lambda3::Float64; rhoA::Float64
  dq1::Float64; dq2::Float64; dq3::Float64; dq4::Float64; dq5::Float64
  tmp1::Float64; tmp2::Float64; tmp3::Float64; fac::Float64; =#

  E1dq = zeros(Float64, 4)
  E2dq = zeros(Float64, 4)

  # getting qg
  qg = zeros(Float64, 4)
  calcIsentropicVortex(x, eqn, qg)

  # viscFac
  # mu
  d1_0 = 1.0
  d0_0 = 0.0
  gami = 1.4 - 1
  sat_Vn = 0.025
  sat_Vl = 0.025

  # Begin main executuion
  nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
  ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]

  # nx = dxidx[1]; ny = dxidx[2]; # nz = dxidx[3]
  # dA = sqrt(nx*nx + ny*ny + nz*nz)
  dA = sqrt(nx*nx + ny*ny)
  
  fac = d1_0/q[1]
  uL = q[2]*fac; vL = q[3]*fac; # wL = q[4]*fac
  # phi = d0_5*(uL*uL + vL*vL + wL*wL)
  phi = d0_5*(uL*uL + vL*vL)
  HL = gamma*q[4]*fac - gami*phi
  
  fac = d1_0/qg[1]
  uR = qg[2]*fac; vR = qg[3]*fac; # wR = qg[4]*fac
  # phi = d0_5*(uR*uR + vR*vR + wR*wR)
  phi = d0_5*(uR*uR + vR*vR)
  HR = gamma*qg[4]*fac - gami*phi

  sqL = sqrt(q[1]); sqR = sqrt(qg[1])
  fac = d1_0/(sqL + sqR)
  u = (sqL*uL + sqR*uR)*fac
  v = (sqL*vL + sqR*vR)*fac
  # w = (sqL*wL + sqR*wR)*fac
  
  H = (sqL*HL + sqR*HR)*fac
  # phi = d0_5*(u*u + v*v + w*w)
  phi = d0_5*(u*u + v*v)
  
  a = sqrt(gami*(H - phi))
  # Un = u*nx + v*ny + w*nz
  Un = u*nx + v*ny

  lambda1 = Un + dA*a
  lambda2 = Un - dA*a
  lambda3 = Un
  rhoA = abs(Un) + dA*a
  lambda1 = d0_5*(tau*max(abs(lambda1),sat_Vn *rhoA) + sgn*lambda1)
  lambda2 = d0_5*(tau*max(abs(lambda2),sat_Vn *rhoA) + sgn*lambda2)
  lambda3 = d0_5*(tau*max(abs(lambda3),sat_Vl *rhoA) + sgn*lambda3)

  dq1 = q[1] - qg[1] 
  dq2 = q[2] - qg[2]
  dq3 = q[3] - qg[3]
  dq4 = q[4] - qg[4]
  # dq5 = q[5] - qg[5]

  #-- diagonal matrix multiply
  sat[1] = lambda3*dq1
  sat[2] = lambda3*dq2
  sat[3] = lambda3*dq3
  sat[4] = lambda3*dq4
  # sat[5] = lambda3*dq5

  #-- get E1*dq
  # E1dq[1] = phi*dq1 - u*dq2 - v*dq3 - w*dq4 + dq5
  E1dq[1] = phi*dq1 - u*dq2 - v*dq3 + dq4
  E1dq[2] = E1dq[1]*u
  E1dq[3] = E1dq[1]*v
  # E1dq[4] = E1dq[1]*w
  # E1dq[5] = E1dq[1]*H
  E1dq[4] = E1dq[1]*H

  #-- get E2*dq
  E2dq[1] = d0_0
  E2dq[2] = -Un*dq1 + nx*dq2 + ny*dq3 + nz*dq4
  E2dq[3] = E2dq[2]*ny
  # E2dq[4] = E2dq[2]*nz
  E2dq[4] = E2dq[2]*Un
  E2dq[2] = E2dq[2]*nx

  #-- add to sat
  tmp1 = d0_5*(lambda1 + lambda2) - lambda3
  tmp2 = gami/(a*a)
  tmp3 = d1_0/(dA*dA)
  sat[:] = sat[:] + tmp1*(tmp2*E1dq[:] + tmp3*E2dq[:])
  
  #-- get E3*dq
  E1dq[1] = -Un*dq1 + nx*dq2 + ny*dq3 + nz*dq4
  E1dq[2] = E1dq[1]*u
  E1dq[3] = E1dq[1]*v
  # E1dq[4] = E1dq[1]*w
  E1dq[4] = E1dq[1]*H

  #-- get E4*dq
  E2dq[1] = d0_0
  E2dq[2] = phi*dq1 - u*dq2 - v*dq3 - w*dq4 + dq5
  E2dq[3] = E2dq[2]*ny
  # E2dq[4] = E2dq[2]*nz
  E2dq[4] = E2dq[2]*Un
  E2dq[2] = E2dq[2]*nx

  #-- add to sat
  tmp1 = d0_5*(lambda1 - lambda2)/(dA*a)
  sat[:] = sat[:] + tmp1*(E1dq[:] + gami*E2dq[:])
  
  return sat + eulerflux(q, nx, ny)
end # ends the function eulerRoeSAT
