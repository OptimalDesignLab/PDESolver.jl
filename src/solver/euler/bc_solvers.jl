
# this file contains all the flux solvers for weakly imposed boundary conditions
@doc """
### EulerEquationMod.RoeSolver
  This calculates the Roe flux for boundary conditions at a node. The inputs 
  must be in *conservative* variables.

  Inputs:
  q  : conservative variables of the fluid
  qg : conservative variables of the boundary
  aux_vars : vector of all auxiliary variables at this node
  dxidx : dxidx matrix at the node
  nrm : sbp face normal vector
  params : ParamType

  Outputs:
    flux : vector to populate with solution

  Aliasing restrictions:  none of the inputs can alias params.res_vals1,
                          params.res_vals2, params.q_vals, params.flux_vals1, or                          params.sat, or params.nrm


"""->
function RoeSolver{Tmsh, Tsol, Tres}(q::AbstractArray{Tsol,1}, 
                                     qg::AbstractArray{Tsol, 1}, 
                                     aux_vars::AbstractArray{Tres, 1}, 
                                     dxidx::AbstractArray{Tmsh,2}, 
                                     nrm::AbstractArray{Tmsh,1}, 
                                     flux::AbstractArray{Tres, 1}, 
                                     params::ParamType{2})

  # SAT terms are used for ensuring consistency with the physical problem. Its 
  # similar to upwinding which adds dissipation to the problem. SATs on the 
  # boundary can be thought of as having two overlapping nodes and because of
  # the discontinuous nature of SBP adds some dissipation.


  E1dq = params.res_vals1
  E2dq = params.res_vals2

  # Declaring constants 
  d1_0 = 1.0
  d0_0 = 0.0
  d0_5 = 0.5
  tau = 1.0
#  sgn = -1.0
  gamma = params.gamma
  gami = params.gamma_1
  sat_Vn = convert(Tsol, 0.025)
  sat_Vl = convert(Tsol, 0.025)
  sat_fac = 1  # multiplier for SAT term

  # Begin main executuion
  nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
  ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]

  dA = sqrt(nx*nx + ny*ny)

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
 
  a = sqrt(gami*(H - phi))
  Un = u*nx + v*ny

  lambda1 = Un + dA*a
  lambda2 = Un - dA*a
  lambda3 = Un
  rhoA = absvalue(Un) + dA*a

  lambda1 = d0_5*(tau*max(absvalue(lambda1),sat_Vn *rhoA) - lambda1)
  lambda2 = d0_5*(tau*max(absvalue(lambda2),sat_Vn *rhoA) - lambda2)
  lambda3 = d0_5*(tau*max(absvalue(lambda3),sat_Vl *rhoA) - lambda3)


  dq1 = q[1] - qg[1] 
  dq2 = q[2] - qg[2]
  dq3 = q[3] - qg[3]
  dq4 = q[4] - qg[4]

  #-- diagonal matrix multiply
#  sat = zeros(Tres, 4)
  sat = params.sat_vals
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
  for i=1:length(sat)
    sat[i] = sat[i] + tmp1*(tmp2*E1dq[i] + tmp3*E2dq[i])
  end
  
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
  for i=1:length(sat)
    sat[i] = sat[i] + tmp1*(E1dq[i] + gami*E2dq[i])
  end

  euler_flux = params.flux_vals1


  # calculate Euler flux in wall normal directiona
  # because edge numbering is rather arbitary, any memory access is likely to
  # be a cache miss, so we recalculate the Euler flux
  v_vals = params.q_vals
  nrm2 = params.nrm
  nrm2[1] = nx
  nrm2[2] = ny

  convertFromNaturalToWorkingVars(params, q, v_vals)
  calcEulerFlux(params, v_vals, aux_vars, nrm2, euler_flux)

  for i=1:4  # ArrayViews does not support flux[:] = .

    flux[i] = -(sat_fac*sat[i] + euler_flux[i]) 
    # when weak differentiate has transpose = true
  end

  return nothing

end # ends the function eulerRoeSAT

function RoeSolver{Tmsh, Tsol, Tres}(q::AbstractArray{Tsol,1}, 
                                     qg::AbstractArray{Tsol, 1}, 
                                     aux_vars::AbstractArray{Tres, 1}, 
                                     dxidx::AbstractArray{Tmsh,2}, 
                                     nrm::AbstractArray{Tmsh,1}, 
                                     flux::AbstractArray{Tres, 1}, 
                                     params::ParamType{3})

  # SAT terms are used for ensuring consistency with the physical problem. Its 
  # similar to upwinding which adds dissipation to the problem. SATs on the 
  # boundary can be thought of as having two overlapping nodes and because of
  # the discontinuous nature of SBP adds some dissipation.

  E1dq = params.res_vals1
  E2dq = params.res_vals2

  # Declaring constants 
  d1_0 = 1.0
  d0_0 = 0.0
  d0_5 = 0.5
  tau = 1.0
#  sgn = -1.0
  gamma = params.gamma
  gami = params.gamma_1
  sat_Vn = convert(Tsol, 0.025)
  sat_Vl = convert(Tsol, 0.025)
  sat_fac = 1  # multiplier for SAT term

  # Begin main executuion
  nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2] + dxidx[3,1]*nrm[3]
  ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2] + dxidx[3,2]*nrm[3]
  nz = dxidx[1,3]*nrm[1] + dxidx[2,3]*nrm[2] + dxidx[3,3]*nrm[3]

  dA = sqrt(nx*nx + ny*ny + nz*nz)

  fac = d1_0/q[1]
  uL = q[2]*fac; vL = q[3]*fac; wL = q[4]*fac
  phi = d0_5*(uL*uL + vL*vL + wL*wL)

  HL = gamma*q[5]*fac - gami*phi
  
  fac = d1_0/qg[1]
  uR = qg[2]*fac; vR = qg[3]*fac; wR = qg[4]*fac
  phi = d0_5*(uR*uR + vR*vR + wR*wR)
  HR = gamma*qg[5]*fac - gami*phi
  sqL = sqrt(q[1]) 
  sqR = sqrt(qg[1])
  fac = d1_0/(sqL + sqR)
  u = (sqL*uL + sqR*uR)*fac
  v = (sqL*vL + sqR*vR)*fac
  w = (sqL*wL + sqR*wR)*fac
  
  H = (sqL*HL + sqR*HR)*fac
  phi = d0_5*(u*u + v*v + w*w)
#  println("H = ", H)
#  println("phi = ", phi)
  a = sqrt(gami*(H - phi))
  Un = u*nx + v*ny + w*nz

  lambda1 = Un + dA*a
  lambda2 = Un - dA*a
  lambda3 = Un
  rhoA = absvalue(Un) + dA*a
#=
  println("before entropy fix:")
  println("rhoA = ", rhoA)
  println("lambda1 = ", lambda1)
  println("lambda2 = ", lambda2)
  println("lambda3 = ", lambda3)
=#

  lambda1 = d0_5*(tau*max(absvalue(lambda1),sat_Vn *rhoA) - lambda1)
  lambda2 = d0_5*(tau*max(absvalue(lambda2),sat_Vn *rhoA) - lambda2)
  lambda3 = d0_5*(tau*max(absvalue(lambda3),sat_Vl *rhoA) - lambda3)
#=
  println("after entropy fix:")
  println("lambda1 = ", lambda1)
  println("lambda2 = ", lambda2)
  println("lambda3 = ", lambda3)
=#
  dq1 = q[1] - qg[1] 
  dq2 = q[2] - qg[2]
  dq3 = q[3] - qg[3]
  dq4 = q[4] - qg[4]
  dq5 = q[5] - qg[5]

  #-- diagonal matrix multiply
#  sat = zeros(Tres, 4)
  sat = params.sat_vals
  sat[1] = lambda3*dq1
  sat[2] = lambda3*dq2
  sat[3] = lambda3*dq3
  sat[4] = lambda3*dq4
  sat[5] = lambda3*dq5

  #-- get E1*dq
  E1dq[1] = phi*dq1 - u*dq2 - v*dq3 - w*dq4 + dq5
  E1dq[2] = E1dq[1]*u
  E1dq[3] = E1dq[1]*v
  E1dq[4] = E1dq[1]*w
  E1dq[5] = E1dq[1]*H

  #-- get E2*dq
  E2dq[1] = d0_0
  E2dq[2] = -Un*dq1 + nx*dq2 + ny*dq3 + nz*dq4
  E2dq[3] = E2dq[2]*ny
  E2dq[4] = E2dq[2]*nz
  E2dq[5] = E2dq[2]*Un
  E2dq[2] = E2dq[2]*nx

  #-- add to sat
  tmp1 = d0_5*(lambda1 + lambda2) - lambda3
  tmp2 = gami/(a*a)
  tmp3 = d1_0/(dA*dA)
#=
  println("tmp1 = ", tmp1)
  println("tmp2 = ", tmp2)
  println("tmp3 = ", tmp3)
=#
  for i=1:5
    sat[i] = sat[i] + tmp1*(tmp2*E1dq[i] + tmp3*E2dq[i])
  end
  
  #-- get E3*dq
  E1dq[1] = -Un*dq1 + nx*dq2 + ny*dq3 + nz*dq4
  E1dq[2] = E1dq[1]*u
  E1dq[3] = E1dq[1]*v
  E1dq[4] = E1dq[1]*w
  E1dq[5] = E1dq[1]*H

  #-- get E4*dq
  E2dq[1] = d0_0
  E2dq[2] = phi*dq1 - u*dq2 - v*dq3 - w*dq4 + dq5
  E2dq[3] = E2dq[2]*ny
  E1dq[4] = E2dq[2]*nz
  E2dq[5] = E2dq[2]*Un
  E2dq[2] = E2dq[2]*nx

  #-- add to sat
  tmp1 = d0_5*(lambda1 - lambda2)/(dA*a)
  for i=1:5
    sat[i] = sat[i] + tmp1*(E1dq[i] + gami*E2dq[i])
  end

  euler_flux = params.flux_vals1


  # calculate Euler flux in wall normal directiona
  # because edge numbering is rather arbitary, any memory access is likely to
  # be a cache miss, so we recalculate the Euler flux
  v_vals = params.q_vals
  nrm2 = params.nrm
  nrm2[1] = nx
  nrm2[2] = ny
  nrm2[3] = nz

  convertFromNaturalToWorkingVars(params, q, v_vals)
  calcEulerFlux(params, v_vals, aux_vars, nrm2, euler_flux)

  for i=1:5  # ArrayViews does not support flux[:] = .

    flux[i] = -(sat_fac*sat[i] + euler_flux[i]) 
    # when weak differentiate has transpose = true
  end

  return nothing

end # ends the function eulerRoeSAT



function LFSolver(qL, qR, aux_vars, dxidx, nrm, flux, params)

  nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
  ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]

  # determine if the left state is entering or existing the element
  l_to_r = nx*qL[2] + ny*qL[3]

  if l_to_r > 0
    calcEulerFlux(params, qL, aux_vars, [nx, ny], flux)
  else
    aux_vars[1] = calcPressure(params, qR)
    calcEulerFlux(params, qR, aux_vars, [nx, ny], flux)
  end

  # negate it
  for i=1:length(flux)
    flux[i] = -flux[i]
  end

  return nothing
end


function AvgSolver(qL, qR, aux_vars, dxidx, nrm, flux, params)

  nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
  ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]


  for i=1:length(qL)
    params.q_vals[i] = 0.5*(qL[i] + qR[i])
  end

  aux_vars[1] = calcPressure(params, params.q_vals)
  calcEulerFlux(params, params.q_vals, aux_vars, [nx, ny], flux)

  # negate it
  for i=1:length(flux)
    flux[i] = -flux[i]
  end



  return nothing
end



