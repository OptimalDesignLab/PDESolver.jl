
# this file contains all the flux solvers for weakly imposed boundary conditions
@doc """
### EulerEquationMod.RoeSolver
  This calculates the Roe flux for boundary conditions at a node

  Inputs:
  q  : conservative variables of the fluid a
  qg : conservative variables of the boundary
  flux_parametric : (scaled) Euler flux in the xi and eta directions
  aux_vars : vector of all auxiliary variables at this node
  dxidx : dxidx matrix at the node
  nrm : sbp face normal vector
  params :: parameter structure

  Outputs:
    flux : vector to populate with solution

  flux_parametric is accessed using *linear* indexing only.  The first 4 entries must be
  the xi direction flux, the next 4 must be the eta direction flux.  This
  makes it possible to pass view(flux_parametric, :, j, i :) and have it work correctly


"""->
function RoeSolver{Tmsh, Tsol, Tres}(q::AbstractArray{Tsol,1}, 
                                     qg::AbstractArray{Tsol, 1}, 
                                     flux_parametric::AbstractArray{Tres}, 
                                     aux_vars::AbstractArray{Tres, 1}, 
                                     dxidx::AbstractArray{Tmsh,2}, 
                                     nrm::AbstractArray{Tmsh,1}, 
                                     flux::AbstractArray{Tres, 1}, 
                                     params::ParamType{2})

  # SAT terms are used for ensuring consistency with the physical problem. Its 
  # similar to upwinding which adds dissipation to the problem. SATs on the 
  # boundary can be thought of as having two overlapping nodes and because of
  # the discontinuous nature of SBP adds some dissipation.

  E1dq = zeros(Tres, 4)
  E2dq = zeros(Tres, 4)

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
  sat_fac = 1  # multiplier for SAT term
  # Begin main executuion
  nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
  ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]

  # println("nx, ny = ", nx, ", ", ny)
  dA = sqrt(nx*nx + ny*ny)
  # println("dA = ",  (dA))

  fac = d1_0/q[1]
  uL = q[2]*fac; vL = q[3]*fac;
  phi = d0_5*(uL*uL + vL*vL)

  HL = gamma*q[4]*fac - gami*phi
  
  fac = d1_0/qg[1]
  uR = qg[2]*fac; vR = qg[3]*fac;
  phi = d0_5*(uR*uR + vR*vR)
  HR = gamma*qg[4]*fac - gami*phi

  sqL = sqrt(q[1]) 
  sqR = sqrt(qg[1])
  fac = d1_0/(sqL + sqR)
  u = (sqL*uL + sqR*uR)*fac
  v = (sqL*vL + sqR*vR)*fac
  
  H = (sqL*HL + sqR*HR)*fac
  phi = d0_5*(u*u + v*v)
  # println("H = ", H, " phi = ", phi)
  
  a = sqrt(gami*(H - phi))
  Un = u*nx + v*ny

  lambda1 = Un + dA*a
  lambda2 = Un - dA*a
  lambda3 = Un
  rhoA = absvalue(Un) + dA*a

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

  # println("E1dq = ",  (E1dq))
  # println("E2dq = ",  (E2dq))

  #-- add to sat
  tmp1 = d0_5*(lambda1 - lambda2)/(dA*a)
  sat[:] = sat[:] + tmp1*(E1dq[:] + gami*E2dq[:])

  euler_flux = zeros(Tsol, 4)
  euler_flux2 = zeros(Tsol, 4)


  calcEulerFlux(params, q, aux_vars, [nx, ny], euler_flux)
  # println("euler_flux = ", euler_flux)

  # calculate Euler flux in wall normal directiona
#  for i=1:4
#    euler_flux[i] = flux_parametric[i]*nrm[1] + flux_parametric[i+4]*nrm[2]
#  end


 #    println("euler_flux = ",  (euler_flux))
  
#  println("euler_flux = ", euler_flux)
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
#    flux[i] = -(sat_fac*sat[i] + 0.5*euler_flux[i] + 0.5*euler_flux2[i])
    flux[i] = -sat_fac*sat[i]  # when weak differentiate has transpose = false
    # flux[i] = -(sat_fac*sat[i] + euler_flux[i]) # when weak differentiate has transpose = true
    # flux[i] =  -euler_flux[i]
#=
if nx < 0.0  # inlet
       flux[1] = qg[1]
       flux[2] = qg[2]
       flux[3] = qg[3]
       flux[4] = q[4]
     else
       flux[1] = q[1]
       flux[2] = q[2] 
       flux[3] = q[3]
       flux[4] = qg[4]
     end
=#
#     flux[i] = euler_flux[i]
  end

#  println("sat = ", sat)
#  println("euler = ", euler_flux)
#  println("Roe flux = ", flux)

#  println("flux = ", flux)
#  print("\n")

#  return sat + getEulerFlux(q, nx, ny)
   return nothing

end # ends the function eulerRoeSAT


