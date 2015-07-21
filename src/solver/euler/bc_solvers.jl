
# this file contains all the flux solvers for weakly imposed boundary conditions

function RoeSolver{Tmsh, Tsol, Tres}( q::AbstractArray{Tsol,1}, qg::AbstractArray{Tsol, 1}, aux_vars::AbstractArray{Tsol, 1}, dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1}, flux::AbstractArray{Tres, 1}, params::ParamType{2})

  E1dq = zeros(Tres, 4)
  E2dq = zeros(Tres, 4)

#  println("entered isentropicOvrtexBC (low level)")
#  println("Tsol = ", Tsol)

  # Declaring constants 
  d1_0 = 1.0
  d0_0 = 0.0
  d0_5 = 0.5
  tau = 1.0
  sgn = -1.0
  gamma = 1.4
  gami = gamma - 1
  sat_Vn = convert(Tsol, 0.025)
  sat_Vl = convert(Tsol, 0.025)

  # Begin main executuion
  nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
  ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]

  dA = sqrt(nx*nx + ny*ny)
  
  fac = d1_0/q[1]
#   println(typeof(fac))
#   println(typeof(q[4]))
  uL = q[2]*fac; vL = q[3]*fac;
  phi = d0_5*(uL*uL + vL*vL)

  HL = gamma*q[4]*fac - gami*phi
  
  fac = d1_0/qg[1]
  uR = qg[2]*fac; vR = qg[3]*fac;
  phi = d0_5*(uR*uR + vR*vR)
  HR = gamma*qg[4]*fac - gami*phi

  sqL = sqrt(q[1]); sqR = sqrt(qg[1])
  fac = d1_0/(sqL + sqR)
  u = (sqL*uL + sqR*uR)*fac
  v = (sqL*vL + sqR*vR)*fac
  
  H = (sqL*HL + sqR*HR)*fac
  phi = d0_5*(u*u + v*v)
  
  a = sqrt(gami*(H - phi))
  Un = u*nx + v*ny

  lambda1 = Un + dA*a
  lambda2 = Un - dA*a
  lambda3 = Un
  rhoA = absvalue(Un) + dA*a

#  println("sat_Vn = ", sat_Vn)
#  println("lambda1 = ", lambda1)
#  println("absvalue(lambda1) = ", absvalue(lambda1))
  lambda1 = d0_5*(tau*max(absvalue(lambda1),sat_Vn *rhoA) + sgn*lambda1)
  lambda2 = d0_5*(tau*max(absvalue(lambda2),sat_Vn *rhoA) + sgn*lambda2)
  lambda3 = d0_5*(tau*max(absvalue(lambda3),sat_Vl *rhoA) + sgn*lambda3)

  dq1 = q[1] - qg[1] 
  dq2 = q[2] - qg[2]
  dq3 = q[3] - qg[3]
  dq4 = q[4] - qg[4]

  #-- diagonal matrix multiply
  sat = zeros(Tres, 4)
  sat[1] = lambda3*dq1
  sat[2] = lambda3*dq2
  sat[3] = lambda3*dq3
  sat[4] = lambda3*dq4

  #-- get E1*dq
  E1dq[1] = phi*dq1 - u*dq2 - v*dq3 + dq4
  E1dq[2] = E1dq[1]*u
  E1dq[3] = E1dq[1]*v
  E1dq[4] = E1dq[1]*H

  #-- get E2*dq
  E2dq[1] = d0_0
  E2dq[2] = -Un*dq1 + nx*dq2 + ny*dq3
  E2dq[3] = E2dq[2]*ny
  E2dq[4] = E2dq[2]*Un
  E2dq[2] = E2dq[2]*nx

  #-- add to sat
  tmp1 = d0_5*(lambda1 + lambda2) - lambda3
  tmp2 = gami/(a*a)
  tmp3 = d1_0/(dA*dA)
  sat[:] = sat[:] + tmp1*(tmp2*E1dq[:] + tmp3*E2dq[:])
  
  #-- get E3*dq
  E1dq[1] = -Un*dq1 + nx*dq2 + ny*dq3
  E1dq[2] = E1dq[1]*u
  E1dq[3] = E1dq[1]*v
  E1dq[4] = E1dq[1]*H

  #-- get E4*dq
  E2dq[1] = d0_0
  E2dq[2] = phi*dq1 - u*dq2 - v*dq3 + dq4
  E2dq[3] = E2dq[2]*ny
  E2dq[4] = E2dq[2]*Un
  E2dq[2] = E2dq[2]*nx

  #-- add to sat
  tmp1 = d0_5*(lambda1 - lambda2)/(dA*a)
  sat[:] = sat[:] + tmp1*(E1dq[:] + gami*E2dq[:])

  euler_flux = zeros(Tsol, 4)
  calcEulerFlux(params, q, aux_vars, [nx, ny], euler_flux)

#  flux[:] = sat + getEulerFlux(q, nx, ny, eqn)
#  flux[:] = -(sat + euler_flux)
  for i=1:4  # ArrayViews does not support flux[:] = .
#=
    println("i = ", i)
    println("flux = ", flux)
    println("typeof(flux) = ", typeof(flux))

    arr = flux.arr
    println("size(arr) = ", size(arr))
    println("offset = ", flux.offset)
    println("len = ", flux.len)
    println("shp = ", flux.shp)
    println("flux[i] = ", flux[i])
    println("sat[i] = ", sat[i])
    println("euler_flux[i] = ", euler_flux[i])
=#
    flux[i] = -(sat[i] + euler_flux[i])
#     flux[i] = euler_flux[i]
  end

#  println("sat = ", sat)
#  println("euler = ", euler_flux)
#  println("Roe flux = ", flux)
 
#  return sat + getEulerFlux(q, nx, ny)
   return nothing

end # ends the function eulerRoeSAT


