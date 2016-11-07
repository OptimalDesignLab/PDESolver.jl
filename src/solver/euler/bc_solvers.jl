# this file contains all the flux solvers for weakly imposed boundary conditions

include("IR_stab.jl")  # stabilization for the IR flux

"""
  Calculate the face integrals in an entropy stable manner for a given
  interface.  Unlike standard face integrals, this requires data from
  the entirety of both elements, not just data interpolated to the face

  resL and resR are updated with the results of the computation for the 
  left and right elements, respectively.

  Note that dxidx_face must contains the scaled metric terms interpolated 
  to the face nodes.

  Aliasing restrictions: none, although its unclear what the meaning of this
                         function would be if resL and resR alias

  Performance note: the version in the tests is the same speed as this one
                    for p=1 Omega elements and about 10% faster for 
                    p=4 elements, but would not be able to take advantage of t
                    he sparsity of R for SBP Gamma elements
"""
                  
function calcESFaceIntegral{Tdim, Tsol, Tres, Tmsh}(
                             params::AbstractParamType{Tdim}, 
                             sbpface::AbstractFace, 
                             iface::Interface,
                             qL::AbstractMatrix{Tsol}, 
                             qR::AbstractMatrix{Tsol}, 
                             aux_vars::AbstractMatrix{Tres}, 
                             dxidx_face::Abstract3DArray{Tmsh},
                             functor::FluxType, 
                             resL::AbstractMatrix{Tres}, 
                             resR::AbstractMatrix{Tres})


  Flux_tmp = params.flux_vals1
  numDofPerNode = length(Flux_tmp)
  nrm = params.nrm
  for dim = 1:Tdim
    fill!(nrm, 0.0)
    nrm[dim] = 1

    # loop over the nodes of "left" element that are in the stencil of interp
    for i = 1:sbpface.stencilsize
      p_i = sbpface.perm[i, iface.faceL]
      qi = sview(qL, :, p_i)
      aux_vars_i = sview(aux_vars, :, p_i)  # !!!! why no aux_vars_j???

      # loop over the nodes of "right" element that are in the stencil of interp
      for j = 1:sbpface.stencilsize
        p_j = sbpface.perm[j, iface.faceR]
        qj = sview(qR, :, p_j)

        # accumulate entry p_i, p_j of E
        Eij = zero(Tres)  # should be Tres
        for k = 1:sbpface.numnodes
          # the computation of nrm_k could be moved outside i,j loops and saved
          # in an array of size [3, sbp.numnodes]
          nrm_k = zero(Tmsh)
          for d = 1:Tdim
            nrm_k += sbpface.normal[d, iface.faceL]*dxidx_face[d, dim, k]
          end
          kR = sbpface.nbrperm[k, iface.orient]
          Eij += sbpface.interp[i,k]*sbpface.interp[j,kR]*sbpface.wface[k]*nrm_k
        end  # end loop k
        
        # compute flux and add contribution to left and right elements
        functor(params, qi, qj, aux_vars_i, nrm, Flux_tmp)
        for p=1:numDofPerNode
          resL[p, p_i] -= Eij*Flux_tmp[p]
          resR[p, p_j] += Eij*Flux_tmp[p]
        end

      end
    end
  end  # end loop Tdim


  return nothing
end


"""
  A wrapper for the Roe Solver that computes the scaled normal vector
  in parametric coordinates from the the face normal and the scaled 
  mapping jacobian.

  Useful for boundary conditions.

""" 
function RoeSolver{Tmsh, Tsol, Tres}(params::ParamType,
                                     q::AbstractArray{Tsol,1}, 
                                     qg::AbstractArray{Tsol, 1}, 
                                     aux_vars::AbstractArray{Tres, 1}, 
                                     dxidx::AbstractArray{Tmsh,2}, 
                                     nrm::AbstractArray{Tmsh,1}, 
                                     flux::AbstractArray{Tres, 1})
                                     
  nrm2 = params.nrm2
  calcBCNormal(params, dxidx, nrm, nrm2)
  RoeSolver(params, q, qg, aux_vars, nrm2, flux)

  return nothing
end




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
                          params.res_vals2, params.q_vals, params.flux_vals1, or                          params.sat


"""->
function RoeSolver{Tmsh, Tsol, Tres}(params::ParamType{2},
                                     q::AbstractArray{Tsol,1}, 
                                     qg::AbstractArray{Tsol, 1}, 
                                     aux_vars::AbstractArray{Tres, 1}, 
                                     nrm::AbstractArray{Tmsh,1}, 
                                     flux::AbstractArray{Tres, 1})

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
  nx = nrm[1]
  ny = nrm[2]

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

    flux[i] = (sat_fac*sat[i] + euler_flux[i]) 
    # when weak differentiate has transpose = true
  end

  return nothing

end # ends the function eulerRoeSAT


"""
  The main Roe solver.  Populates `flux` with the computed flux.
"""
function RoeSolver{Tmsh, Tsol, Tres}(params::ParamType{3},
                                     q::AbstractArray{Tsol,1}, 
                                     qg::AbstractArray{Tsol, 1}, 
                                     aux_vars::AbstractArray{Tres, 1}, 
                                     nrm::AbstractArray{Tmsh,1}, 
                                     flux::AbstractArray{Tres, 1})
                                     

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
  nx = nrm[1]
  ny = nrm[2]
  nz = nrm[3]

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

  lambda1 = d0_5*(tau*max(absvalue(lambda1),sat_Vn *rhoA) - lambda1)
  lambda2 = d0_5*(tau*max(absvalue(lambda2),sat_Vn *rhoA) - lambda2)
  lambda3 = d0_5*(tau*max(absvalue(lambda3),sat_Vl *rhoA) - lambda3)
  
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
  v_vals2 = zeros(v_vals)

  convertFromNaturalToWorkingVars(params, q, v_vals)
  convertFromNaturalToWorkingVars(params, qg, v_vals2)
  calcEulerFlux(params, v_vals, aux_vars, nrm, euler_flux)
#  calcEulerFlux_Ducros(params, v_vals, v_vals2, nrm, euler_flux)

  for i=1:5  # ArrayViews does not support flux[:] = .

    flux[i] = (sat_fac*sat[i] + euler_flux[i]) 
    # when weak differentiate has transpose = true
  end

  return nothing

end # ends the function eulerRoeSAT


function calcEulerFlux_standard{Tmsh, Tsol, Tres}(params::ParamType,
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tres}, 
                      dxidx::AbstractMatrix{Tmsh},
                      nrm::AbstractArray{Tmsh},  F::AbstractArray{Tres,1})

  nrm2 = params.nrm2
  calcBCNormal(params, dxidx, nrm, nrm2)
  calcEulerFlux_standard(params, qL, qR, aux_vars, nrm2, F)
  return nothing
end



function calcEulerFlux_standard{Tmsh, Tsol, Tres}(
                      params::ParamType{2, :conservative}, 
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tsol, 1},
                      dir::AbstractArray{Tmsh},  F::AbstractArray{Tres,1})
# calculate the split form numerical flux function corresponding to the standard DG flux

  pL = calcPressure(params, qL); pR = calcPressure(params, qR)
  rho_avg = 0.5*(qL[1] + qR[1])
  rhou_avg = 0.5*(qL[2] + qR[2])
  rhov_avg = 0.5*(qL[3] + qR[3])
  p_avg = 0.5*(pL + pR)
  rhoLinv = 1/qL[1]; rhoRinv = 1/qR[1]



  F[1] = dir[1]*(rhou_avg) + dir[2]*rhov_avg

  tmp1 = 0.5*(qL[2]*qL[2]*rhoLinv + qR[2]*qR[2]*rhoRinv)
  tmp2 = 0.5*(qL[2]*qL[3]*rhoLinv + qR[2]*qR[3]*rhoRinv)
  F[2] = dir[1]*(tmp1 + p_avg) + dir[2]*tmp2

  tmp1 = 0.5*(qL[2]*qL[3]*rhoLinv + qR[2]*qR[3]*rhoRinv)
  tmp2 = 0.5*(qL[3]*qL[3]*rhoLinv + qR[3]*qR[3]*rhoRinv)
  F[3] = dir[1]*tmp1 + dir[2]*(tmp2 + p_avg)


  tmp1 = 0.5*( (qL[4] + pL)*qL[2]*rhoLinv + (qR[4] + pR)*qR[2]*rhoRinv)
  tmp2 = 0.5*( (qL[4] + pL)*qL[3]*rhoLinv + (qR[4] + pR)*qR[3]*rhoRinv)
  F[4] = dir[1]*tmp1 + dir[2]*tmp2

  return nothing
end

function calcEulerFlux_standard{Tmsh, Tsol, Tres}(
                      params::ParamType{3, :conservative}, 
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tres, 1},
                      dir::AbstractArray{Tmsh},  F::AbstractArray{Tres,1})
# calculate the split form numerical flux function corresponding to the standard DG flux
#TODO: pre-calculate 1/qL[1], 1/qR[1]

  pL = calcPressure(params, qL); pR = calcPressure(params, qR)
  rho_avg = 0.5*(qL[1] + qR[1])
  rhou_avg = 0.5*(qL[2] + qR[2])
  rhov_avg = 0.5*(qL[3] + qR[3])
  rhow_avg = 0.5*(qL[4] + qR[4])
  p_avg = 0.5*(pL + pR)
  rhoLinv = 1/qL[1]; rhoRinv = 1/qR[1]

  F[1] = dir[1]*(rhou_avg) + dir[2]*rhov_avg + dir[3]*rhow_avg

  tmp1 = 0.5*(qL[2]*qL[2]*rhoLinv + qR[2]*qR[2]*rhoRinv)
  tmp2 = 0.5*(qL[2]*qL[3]*rhoLinv + qR[2]*qR[3]*rhoRinv)
  tmp3 = 0.5*(qL[2]*qL[4]*rhoLinv + qR[2]*qR[4]*rhoRinv)
  F[2] = dir[1]*(tmp1 + p_avg) + dir[2]*tmp2 + dir[3]*tmp3

  tmp1 = 0.5*(qL[2]*qL[3]*rhoLinv + qR[2]*qR[3]*rhoRinv)
  tmp2 = 0.5*(qL[3]*qL[3]*rhoLinv + qR[3]*qR[3]*rhoRinv)
  tmp3 = 0.5*(qL[4]*qL[3]*rhoLinv + qR[4]*qR[3]*rhoRinv)
  F[3] = dir[1]*tmp1 + dir[2]*(tmp2 + p_avg) + dir[3]*tmp3

  tmp1 = 0.5*(qL[2]*qL[4]*rhoLinv + qR[2]*qR[4]*rhoRinv)
  tmp2 = 0.5*(qL[3]*qL[4]*rhoLinv + qR[3]*qR[4]*rhoRinv)
  tmp3 = 0.5*(qL[4]*qL[4]*rhoLinv + qR[4]*qR[4]*rhoRinv)
  F[4] = dir[1]*tmp1 + dir[2]*tmp2 + dir[3]*(tmp3 + p_avg)


  tmp1 = 0.5*( (qL[5] + pL)*qL[2]*rhoLinv + (qR[5] + pR)*qR[2]*rhoRinv)
  tmp2 = 0.5*( (qL[5] + pL)*qL[3]*rhoLinv + (qR[5] + pR)*qR[3]*rhoRinv)
  tmp3 = 0.5*( (qL[5] + pL)*qL[4]*rhoLinv + (qR[5] + pR)*qR[4]*rhoRinv)
  F[5] = dir[1]*tmp1 + dir[2]*tmp2 + dir[3]*tmp3

  return nothing
end



function calcEulerFlux_Ducros{Tmsh, Tsol, Tres}(
                      params::ParamType, 
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tres}, 
                      dxidx::AbstractMatrix{Tmsh},
                      nrm::AbstractArray{Tmsh},  F::AbstractArray{Tres,1})

  nrm2 = params.nrm2
  calcBCNormal(params, dxidx, nrm, nrm2)
  calcEulerFlux_Ducros(params, qL, qR, aux_vars, nrm2, F)
  return nothing
end



function calcEulerFlux_Ducros{Tmsh, Tsol, Tres}(params::ParamType{2, :conservative}, 
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tres},
                      dir::AbstractArray{Tmsh},  F::AbstractArray{Tres,1})
# calculate the split form numerical flux function proposed by Ducros et al.

  pL = calcPressure(params, qL); pR = calcPressure(params, qR)
  uL = qL[2]/qL[1]; uR = qR[2]/qR[1]
  vL = qL[3]/qL[1]; vR = qR[3]/qR[1]
  
  u_avg = 0.5*(uL + uR)
  v_avg = 0.5*(vL + vR)

  rho_avg = 0.5*(qL[1] + qR[1])
  rhou_avg = 0.5*(qL[2] + qR[2])
  rhov_avg = 0.5*(qL[3] + qR[3])
  E_avg = 0.5*(qL[4] + qR[4])
  p_avg = 0.5*(pL + pR)

  F[1] = dir[1]*rho_avg*u_avg + dir[2]*rho_avg*v_avg
  F[2] = dir[1]*(rhou_avg*u_avg + p_avg) + dir[2]*(rhou_avg*v_avg)
  F[3] = dir[1]*(rhov_avg*u_avg) + dir[2]*(rhov_avg*v_avg + p_avg)
  F[4] = dir[1]*(E_avg + p_avg)*u_avg + dir[2]*(E_avg + p_avg)*v_avg

  return nothing
end

function calcEulerFlux_Ducros{Tmsh, Tsol, Tres}(params::ParamType{3, :conservative}, 
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tres},
                      dir::AbstractArray{Tmsh},  F::AbstractArray{Tres,1})
# calculate the split form numerical flux function proposed by Ducros et al.

  pL = calcPressure(params, qL); pR = calcPressure(params, qR)
  uL = qL[2]/qL[1]; uR = qR[2]/qR[1]
  vL = qL[3]/qL[1]; vR = qR[3]/qR[1]
  wL = qL[4]/qL[1]; wR = qR[4]/qR[1]

  u_avg = 0.5*(uL + uR)
  v_avg = 0.5*(vL + vR)
  w_avg = 0.5*(wL + wR)

  rho_avg = 0.5*(qL[1] + qR[1])
  rhou_avg = 0.5*(qL[2] + qR[2])
  rhov_avg = 0.5*(qL[3] + qR[3])
  rhow_avg = 0.5*(qL[4] + qR[4])
  E_avg = 0.5*(qL[5] + qR[5])
  p_avg = 0.5*(pL + pR)

  F[1] = dir[1]*rho_avg*u_avg + dir[2]*rho_avg*v_avg + dir[3]*rho_avg*w_avg
  F[2] = dir[1]*(rhou_avg*u_avg + p_avg) + dir[2]*(rhou_avg*v_avg) + 
         dir[3]rhou_avg*w_avg
  F[3] = dir[1]*(rhov_avg*u_avg) + dir[2]*(rhov_avg*v_avg + p_avg) + 
         dir[3]*rhov_avg*w_avg
  F[4] = dir[1]*rhow_avg*u_avg + dir[2]*rhow_avg*v_avg + 
         dir[3]*(rhow_avg*w_avg + p_avg)
  F[5] = dir[1]*(E_avg + p_avg)*u_avg + dir[2]*(E_avg + p_avg)*v_avg + 
         dir[3]*( (E_avg + p_avg)*w_avg)

  return nothing
end


function calcEulerFlux_IR{Tmsh, Tsol, Tres}(params::ParamType,
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tres}, 
                      dxidx::AbstractMatrix{Tmsh},
                      nrm::AbstractArray{Tmsh},  F::AbstractArray{Tres,1})

  nrm2 = params.nrm2
  calcBCNormal(params, dxidx, nrm, nrm2)
  calcEulerFlux_IR(params, qL, qR, aux_vars, nrm2, F)
  return nothing
end


# IR flux
function calcEulerFlux_IR{Tmsh, Tsol, Tres}(params::ParamType{2, :conservative}, 
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tres},
                      dir::AbstractArray{Tmsh},  F::AbstractArray{Tres,1})

  gamma = params.gamma
  gamma_1 = params.gamma_1
  pL = calcPressure(params, qL); pR = calcPressure(params, qR)
  z1L = sqrt(qL[1]/pL); z1R = sqrt(qR[1]/pR)
  z2L = z1L*qL[2]/qL[1]; z2R = z1R*qR[2]/qR[1]
  z3L = z1L*qL[3]/qL[1]; z3R = z1R*qR[3]/qR[1]
  z4L = sqrt(qL[1]*pL); z4R = sqrt(qR[1]*pR)

  rho_hat = 0.5*(z1L + z1R)*logavg(z4L, z4R)
  u_hat = (z2L + z2R)/(z1L + z1R)
  v_hat = (z3L + z3R)/(z1L + z1R)
  p1_hat = (z4L + z4R)/(z1L + z1R)
  p2_hat = ((gamma + 1)/(2*gamma) )*logavg(z4L, z4R)/logavg(z1L, z1R) + ( gamma_1/(2*gamma) )*(z4L + z4R)/(z1L + z1R)
  h_hat = gamma*p2_hat/(rho_hat*gamma_1) + 0.5*(u_hat*u_hat + v_hat*v_hat)

  
  Un = dir[1]*u_hat + dir[2]*v_hat
  F[1] = rho_hat*Un
  F[2] = rho_hat*u_hat*Un + dir[1]*p1_hat
  F[3] = rho_hat*v_hat*Un + dir[2]*p1_hat
  F[4] = rho_hat*h_hat*Un
  #=
  F[1] = dir[1]*rho_hat*u_hat + dir[2]*rho_hat*v_hat
  F[2] = dir[1]*(rho_hat*u_hat*u_hat + p1_hat) + dir[2]*rho_hat*u_hat*v_hat
  F[3] = dir[1]*rho_hat*u_hat*v_hat + dir[2]*(rho_hat*v_hat*v_hat + p1_hat)
  F[4] = dir[1]*rho_hat*u_hat*h_hat + dir[2]*rho_hat*v_hat*h_hat
  =#
  return nothing
end

function calcEulerFlux_IR{Tmsh, Tsol, Tres}(params::ParamType{3, :conservative}, 
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tres},
                      dir::AbstractArray{Tmsh},  F::AbstractArray{Tres,1})

  gamma = params.gamma
  gamma_1 = params.gamma_1
  pL = calcPressure(params, qL); pR = calcPressure(params, qR)
  z1L = sqrt(qL[1]/pL); z1R = sqrt(qR[1]/pR)
  z2L = z1L*qL[2]/qL[1]; z2R = z1R*qR[2]/qR[1]
  z3L = z1L*qL[3]/qL[1]; z3R = z1R*qR[3]/qR[1]
  z4L = z1L*qL[4]/qL[1]; z4R = z1R*qR[4]/qR[1]
  z5L = sqrt(qL[1]*pL); z5R = sqrt(qR[1]*pR)

  rho_hat = 0.5*(z1L + z1R)*logavg(z5L, z5R)
  u_hat = (z2L + z2R)/(z1L + z1R)
  v_hat = (z3L + z3R)/(z1L + z1R)
  w_hat = (z4L + z4R)/(z1L + z1R)
  p1_hat = (z5L + z5R)/(z1L + z1R)
  p2_hat = ((gamma + 1)/(2*gamma) )*logavg(z5L, z5R)/logavg(z1L, z1R) + ( gamma_1/(2*gamma) )*(z5L + z5R)/(z1L + z1R)
  h_hat = gamma*p2_hat/(rho_hat*gamma_1) + 0.5*(u_hat*u_hat + v_hat*v_hat + w_hat*w_hat)

  F[1] = dir[1]*rho_hat*u_hat + dir[2]*rho_hat*v_hat + dir[3]*rho_hat*w_hat
  F[2] = dir[1]*(rho_hat*u_hat*u_hat + p1_hat) + dir[2]*rho_hat*u_hat*v_hat + 
         dir[3]*rho_hat*u_hat*w_hat
  F[3] = dir[1]*rho_hat*u_hat*v_hat + dir[2]*(rho_hat*v_hat*v_hat + p1_hat) + 
         dir[3]*rho_hat*v_hat*w_hat
  F[4] = dir[1]*rho_hat*u_hat*w_hat + dir[2]*rho_hat*v_hat*w_hat + 
         dir[3]*(rho_hat*w_hat*w_hat + p1_hat)
  F[5] = dir[1]*rho_hat*u_hat*h_hat + dir[2]*rho_hat*v_hat*h_hat + dir[3]*rho_hat*w_hat*h_hat

  return nothing
end

# stabilized IR flux

function calcEulerFlux_IRStable{Tmsh, Tsol, Tres}(params::ParamType,
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tres}, 
                      dxidx::AbstractMatrix{Tmsh},
                      nrm::AbstractArray{Tmsh},  F::AbstractArray{Tres,1})

  println("in calcEulerFlux_IRStable matrix version")
  nrm2 = params.nrm2
  calcBCNormal(params, dxidx, nrm, nrm2)
  calcEulerFlux_IRStable(params, qL, qR, aux_vars, nrm2, F)
  return nothing
end


function calcEulerFlux_IRStable{Tmsh, Tsol, Tres, Tdim}(
                      params::ParamType{Tdim, :conservative}, 
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractVector{Tres},
                      dir::AbstractVector{Tmsh},  F::AbstractArray{Tres,1})

  calcEulerFlux_IR(params, qL, qR, aux_vars, dir, F)
  getIRStab1(params, qL, qR, aux_vars, dir, F)

  return nothing
end



function logavg(aL, aR)
# calculate the logarithmic average needed by the IR flux
  xi = aL/aR
  f = (xi - 1)/(xi + 1)
  u = f*f
  eps = 1e-2
  if u < eps
    F = @evalpoly( u, 1, 1/3, 1/5, 1/7, 1/9)
#    F = 1.0 + u/3.0 + u*u/5.0 + u*u*u/7.0 + u*u*u*u/9.0
  else
    F = (log(xi)/2.0)/f
  end

  return (aL + aR)/(2*F)
end



