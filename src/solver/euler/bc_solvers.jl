# this file contains all the flux solvers for weakly imposed boundary conditions

# Some of these flux functions are used by the FaceElementIntegrals.
# These flux functions should implement 2 forms the function: a version that
# computes the flux in a single direction and a version that computes the flux
# in d dimensions simultaneously (where d is the dimensionality of the system)
# In some cases the second form is significantly more computationally efficient
# See calcEulerFlux_IR() for an example

@doc """
### EulerEquationMod.RoeSolver
  This calculates the Roe flux for boundary conditions at a node. The inputs
  must be in *conservative* variables.

  **Inputs**
   * params : ParamType
   * q  : conservative variables of the fluid
   * qg : conservative variables of the boundary
   * aux_vars : vector of all auxiliary variables at this node
   * nrm : scaled face normal vector (x-y space, outward normal of the face q
        lives on)

  **Outputs**
   * flux : vector to populate with solution

  Aliasing restrictions:  none of the inputs can alias params.res_vals1,
                          params.res_vals2, params.q_vals, params.flux_vals1, or
                          params.sat


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

  # Declaring constants
  d1_0 = 1.0
  d0_0 = 0.0
  d0_5 = 0.5
  tau = 1.0
  gamma = params.gamma
  gami = params.gamma_1
  sat_fac = 1  # multiplier for SAT term

  # Begin main executuion
  nx = nrm[1]
  ny = nrm[2]

  # Compute the Roe Averaged states
  # The left state of Roe are the actual solution variables
  fac = d1_0/q[1]
  uL = q[2]*fac; vL = q[3]*fac;
  phi = d0_5*(uL*uL + vL*vL)
  HL = gamma*q[4]*fac - gami*phi # Total enthalpy, H = e + 0.5*(u^2 + v^2) + p/rho,
                                 # where e is the internal energy per unit mass

  # The right side of the Roe solver comprises the boundary conditions
  fac = d1_0/qg[1]
  uR = qg[2]*fac; vR = qg[3]*fac;
  phi = d0_5*(uR*uR + vR*vR)
  HR = gamma*qg[4]*fac - gami*phi # Total Enthalpy

  # Averaged states
  sqL = sqrt(q[1])
  sqR = sqrt(qg[1])
  fac = d1_0/(sqL + sqR)
  u = (sqL*uL + sqR*uR)*fac
  v = (sqL*vL + sqR*vR)*fac
  H = (sqL*HL + sqR*HR)*fac

  dq = params.v_vals2 # zeros(Tsol, 4)
  for i=1:length(dq)
    dq[i] = q[i] - qg[i]
  end
  sat = params.sat_vals
  roe_vars = params.roe_vars
  roe_vars[1] = u
  roe_vars[2] = v
  roe_vars[3] = H
  calcSAT(params, roe_vars, dq, nrm, sat)

  euler_flux = params.flux_vals1
  # calculate Euler flux in wall normal directiona
  # because edge numbering is rather arbitary, any memory access is likely to
  # be a cache miss, so we recalculate the Euler flux
  v_vals = params.q_vals
  nrm2 = params.nrm
  nrm2[1] = nx   # why are we assigning to nrm2?
  nrm2[2] = ny

  convertFromNaturalToWorkingVars(params, q, v_vals)
  calcEulerFlux(params, v_vals, aux_vars, nrm2, euler_flux)

  for i=1:4  # ArrayViews does not support flux[:] = .
    flux[i] = (sat_fac*sat[i] + euler_flux[i])
  end

  return nothing

end # ends the function eulerRoeSAT


@doc """
###EulerEquationMod.RoeSolver_revm

Reverse mode of `EulerEquationMod.RoeSolver`. This function computes the
reverse mode of the Roe flux w.r.t the mesh metrics

**Inputs**

* `params` : Parameter object
* `q`  : Conservative variable of the fluid
* `qg` : Conservative variable of the boundary or the adjacent element
* `aux_vars` : Auxiliary variables
* `nrm` : Element face normal vector in the physical space
* `flux_bar` : Flux value which needs to get differentiated

**Output**

* `nrm_bar` : derivaitve of the flux_bar w.r.t the mesh metrics

Aliasing Restrictions: Same as the forward function

"""->

function RoeSolver_revm{Tmsh, Tsol, Tres}(params::ParamType{2},
                                     q::AbstractArray{Tsol,1},
                                     qg::AbstractArray{Tsol, 1},
                                     aux_vars::AbstractArray{Tres, 1},
                                     nrm::AbstractArray{Tmsh,1},
                                     flux_bar::AbstractArray{Tres, 1},
                                     nrm_bar::AbstractArray{Tmsh, 1})

  # Forward sweep
  tau = 1.0
  gamma = params.gamma
  gami = params.gamma_1
  sat_fac = 1  # multiplier for SAT term

  # Begin main executuion
  nx = nrm[1]
  ny = nrm[2]

  # Compute the Roe Averaged states
  fac = 1.0/q[1]
  uL = q[2]*fac; vL = q[3]*fac;
  phi = 0.5*(uL*uL + vL*vL)
  HL = gamma*q[4]*fac - gami*phi
  fac = 1.0/qg[1]
  uR = qg[2]*fac
  vR = qg[3]*fac
  phi = 0.5*(uR*uR + vR*vR)
  HR = gamma*qg[4]*fac - gami*phi # Total Enthalpy
  sqL = sqrt(q[1])
  sqR = sqrt(qg[1])
  fac = 1.0/(sqL + sqR)
  u = (sqL*uL + sqR*uR)*fac
  v = (sqL*vL + sqR*vR)*fac
  H = (sqL*HL + sqR*HR)*fac
  dq = params.v_vals2
  dq[:] = q[:] - qg[:]
  v_vals = params.q_vals
  nrm2 = params.nrm
  nrm2[1] = nx
  nrm2[2] = ny
  convertFromNaturalToWorkingVars(params, q, v_vals)

  # Reverse Sweep
  # for i=1:4  # ArrayViews does not support flux[:] = .
  #   flux[i] = (sat_fac*sat[i] + euler_flux[i])
  # end
  euler_flux_bar = params.flux_vals1 # zeros(Tsol, 4)
  sat_bar = params.sat_vals
  fill!(euler_flux_bar, 0.0)
  fill!(sat_bar, 0.0)
  for i = 4:-1:1
    euler_flux_bar[i] += flux_bar[i]
    sat_bar[i] += sat_fac*flux_bar[i]
  end

  # calcEulerFlux(params, v_vals, aux_vars, nrm2, euler_flux)
  nrm2_bar = zeros(Tmsh, 2)
  calcEulerFlux_revm(params, v_vals, aux_vars, nrm2, euler_flux_bar, nrm2_bar)

  # nrm2[2] = ny
  ny_bar = nrm2_bar[2]
  # nrm2[1] = nx
  nx_bar = nrm2_bar[1]

  #  calcSAT(params, nrm, dq, sat, [u, v], H)
  calcSAT_revm(params, nrm, dq, [u, v], H, sat_bar, nrm_bar)
  nrm_bar[1] += nx_bar
  nrm_bar[2] += ny_bar

  return nothing
end



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

  dq = params.v_vals2 # zeros(Tsol, 4)
  for i=1:length(dq)
    dq[i] = q[i] - qg[i]
  end

  roe_vars = params.roe_vars
  roe_vars[1] = u
  roe_vars[2] = v
  roe_vars[3] = w
  roe_vars[4] = H
  sat = params.sat_vals

  calcSAT(params, roe_vars, dq, nrm, sat)

  # calculate Euler flux in wall normal directiona
  # because edge numbering is rather arbitary, any memory access is likely to
  # be a cache miss, so we recalculate the Euler flux
  euler_flux = params.flux_vals1
  v_vals = params.q_vals
  convertFromNaturalToWorkingVars(params, q, v_vals)
  calcEulerFlux(params, v_vals, aux_vars, nrm, euler_flux)

  for i=1:5  # ArrayViews does not support flux[:] = .
    flux[i] = (sat_fac*sat[i] + euler_flux[i])
  end

  return nothing

end # ends the function eulerRoeSAT


function RoeSolver_revm{Tmsh, Tsol, Tres}(params::ParamType{3},
                        q::AbstractArray{Tsol,1}, qg::AbstractArray{Tsol, 1},
                        aux_vars::AbstractArray{Tres, 1},
                        nrm::AbstractArray{Tmsh,1},
                        flux_bar::AbstractArray{Tres, 1},
                        nrm_bar::AbstractArray{Tmsh,1})

  E1dq = params.res_vals1
  E2dq = params.res_vals2
  # E1dq = zeros(Tres, 5)
  # E2dq = zeros(Tres, 5)

  # Declaring constants
  d1_0 = one(Tres)
  d0_0 = zero(Tres)
  d0_5 = 0.5*one(Tres)
  tau = one(Tres)

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

  # start reverse sweep

  euler_flux_bar = zeros(flux_bar)
  sat_bar = zeros(flux_bar)
  for i=1:5
    # flux[i] = (sat_fac*sat[i] + euler_flux[i])
    euler_flux_bar[i] += flux_bar[i]
    sat_bar[i] += sat_fac*flux_bar[i]
  end

  # calcEulerFlux(q, nrm, euler_flux)
  # calcEulerFlux_revm(params, v_vals, aux_vars, nrm2, euler_flux_bar, nrm2_bar)
  calcEulerFlux_revm(params, q, aux_vars, nrm, euler_flux_bar, nrm_bar)

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
  E2dq[4] = E2dq[2]*nz
  E2dq[5] = E2dq[2]*Un
  E2dq[2] = E2dq[2]*nx

  tmp1 = d0_5*(lambda1 - lambda2)/(dA*a)
  tmp1_bar = zero(tmp1)
  E1dq_bar = zeros(E1dq)
  E2dq_bar = zeros(E2dq)
  for i = 1:5
    # sat[i] = sat[i] + tmp1*(E1dq[i] + gami*E2dq[i])
    tmp1_bar += (E1dq[i] + gami*E2dq[i])*sat_bar[i]
    E1dq_bar[i] += tmp1*sat_bar[i]
    E2dq_bar[i] += tmp1*gami*sat_bar[i]
  end
  # tmp1 = d0_5*(lambda1 - lambda2)/(dA*a)
  lambda1_bar = d0_5*tmp1_bar/(dA*a)
  lambda2_bar = -d0_5*tmp1_bar/(dA*a)
  dA_bar = -d0_5*(lambda1 - lambda2)*tmp1_bar/(dA*dA*a)

  fac = phi*dq1 - u*dq2 - v*dq3 - w*dq4 + dq5
  #E2dq[5] = fac*Un
  Un_bar = E2dq_bar[5]*fac
  # E2dq[4] = fac*nz
  nz_bar = E2dq_bar[4]*fac
  # E2dq[3] = fac*ny
  ny_bar = E2dq_bar[3]*fac
  # E2dq[2] = fac*nx
  nx_bar = E2dq_bar[2]*fac

  fac = -Un*dq1 + nx*dq2 + ny*dq3 + nz*dq4
  fac_bar = zero(fac)
  # E1dq[5] = fac*H
  # E1dq[4] = fac*w
  # E1dq[3] = fac*v
  # E1dq[2] = fac*u
  # E1dq[1] = fac
  fac_bar += (E1dq_bar[5]*H + E1dq_bar[4]*w + E1dq_bar[3]*v + E1dq_bar[2]*u
  + E1dq_bar[1])
  # fac = -Un*dq1 + nx*dq2 + ny*dq3 + nz*dq4
  Un_bar += -dq1*fac_bar
  nx_bar += dq2*fac_bar
  ny_bar += dq3*fac_bar
  nz_bar += dq4*fac_bar

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

  tmp1 = d0_5*(lambda1 + lambda2) - lambda3
  tmp2 = gami/(a*a)
  tmp3 = d1_0/(dA*dA)
  tmp1_bar = zero(tmp1)
  tmp3_bar = zero(tmp3)
  fill!(E1dq_bar, zero(Tres))
  fill!(E2dq_bar, zero(Tres))
  for i=1:5
    # sat[i] = sat[i] + tmp1*(tmp2*E1dq[i] + tmp3*E2dq[i])
    tmp1_bar += (tmp2*E1dq[i] + tmp3*E2dq[i])*sat_bar[i]
    tmp3_bar += tmp1*E2dq[i]*sat_bar[i]
    # E1dq_bar[i] += tmp1*tmp2*sat_bar[i] # E1dq is independent of nrm
    E2dq_bar[i] += tmp1*tmp3*sat_bar[i]
  end
  # tmp1 = d0_5*(lambda1 + lambda2) - lambda3
  lambda1_bar += d0_5*tmp1_bar
  lambda2_bar += d0_5*tmp1_bar
  lambda3_bar = -tmp1_bar
  # tmp3 = d1_0/(dA*dA)
  dA_bar += -2.0*tmp3_bar/(dA*dA*dA)

  fac = -Un*dq1 + nx*dq2 + ny*dq3 + nz*dq4
  fac_bar = zero(fac)
  # E2dq[5] = fac*Un
  Un_bar += fac*E2dq_bar[5]
  fac_bar += Un*E2dq_bar[5]
  # E2dq[4] = fac*nz
  nz_bar += fac*E2dq_bar[4]
  fac_bar += nz*E2dq_bar[4]
  # E2dq[3] = fac*ny
  ny_bar += fac*E2dq_bar[3]
  fac_bar += ny*E2dq_bar[3]
  # E2dq[2] = fac*nx
  nx_bar += fac*E2dq_bar[2]
  fac_bar += nx*E2dq_bar[2]
  # fac = -Un*dq1 + nx*dq2 + ny*dq3 + nz*dq4
  Un_bar += -dq1*fac_bar
  nx_bar += dq2*fac_bar
  ny_bar += dq3*fac_bar
  nz_bar += dq4*fac_bar

  # E1*dq is independent of nrm
  # E1dq[1] = phi*dq1 - u*dq2 - v*dq3 - w*dq4 + dq5
  # E1dq[2] = E1dq[1]*u
  # E1dq[3] = E1dq[1]*v
  # E1dq[4] = E1dq[1]*w
  # E1dq[5] = E1dq[1]*H

  # sat[1] = lambda3*dq1
  # sat[2] = lambda3*dq2
  # sat[3] = lambda3*dq3
  # sat[4] = lambda3*dq4
  # sat[5] = lambda3*dq5
  lambda3_bar += (dq1*sat_bar[1] + dq2*sat_bar[2] + dq3*sat_bar[3] + dq4*sat_bar[4]
  + dq5*sat_bar[5])

  rhoA_bar = zero(rhoA)

  # lambda1 = d0_5*(tau*max_val - lambda1)
  max_val_bar = d0_5*tau*lambda1_bar
  lambda1_bar = -d0_5*lambda1_bar
  # max_val = max(absvalue(lambda1),sat_Vn *rhoA)
  tmp1_bar, tmp2_bar = max_deriv_rev(absvalue(Un + dA*a), sat_Vn*rhoA, max_val_bar)
  lambda1_bar += tmp1_bar*absvalue_deriv(Un + dA*a)
  rhoA_bar += tmp2_bar*sat_Vn

  # lambda2 = d0_5*(tau*max_val - lambda2)
  max_val_bar = d0_5*tau*lambda2_bar
  lambda2_bar = -d0_5*lambda2_bar
  # max_val =  max(absvalue(lambda2),sat_Vn *rhoA)
  tmp1_bar, tmp2_bar = max_deriv_rev(absvalue(Un - dA*a), sat_Vn*rhoA, max_val_bar)
  lambda2_bar += tmp1_bar*absvalue_deriv(Un - dA*a)
  rhoA_bar += tmp2_bar*sat_Vn

  # lambda3 = d0_5*(tau*max_val - lambda3)
  max_val_bar = d0_5*tau*lambda3_bar
  lambda3_bar = -d0_5*lambda3_bar
  # max_val =  max(absvalue(lambda3),sat_Vn *rhoA)
  tmp1_bar, tmp2_bar = max_deriv_rev(absvalue(Un), sat_Vn*rhoA, max_val_bar)
  lambda3_bar += tmp1_bar*absvalue_deriv(Un)
  rhoA_bar += tmp2_bar*sat_Vn

  #rhoA = absvalue(Un) + dA*a
  Un_bar += absvalue_deriv(Un)*rhoA_bar
  dA_bar += a*rhoA_bar
  #lambda3 = Un
  Un_bar += lambda3_bar
  #lambda2 = Un - dA*a
  Un_bar += lambda2_bar
  dA_bar += -a*lambda2_bar
  #lambda1 = Un + dA*a
  Un_bar += lambda1_bar
  dA_bar += a*lambda1_bar

  # Un = u*nx + v*ny + w*nz
  nx_bar += u*Un_bar
  ny_bar += v*Un_bar
  nz_bar += w*Un_bar

  # dA = sqrt(nx*nx + ny*ny + nz*nz)
  nx_bar += nx*dA_bar/dA
  ny_bar += ny*dA_bar/dA
  nz_bar += nz*dA_bar/dA

  # nx = nrm[1]
  # ny = nrm[2]
  # nz = nrm[3]
  nrm_bar[1] += nx_bar
  nrm_bar[2] += ny_bar
  nrm_bar[3] += nz_bar

  return nothing

end # ends the function eulerRoeSAT_revm

@doc """
###EulerEquationMod.calcSAT

Computes the simultaneous approximation term for use in computing the numerical
flux, namely:

  0.5*(|A_hat| - A)(qL - qR)

**Arguments**

 * `params` : Parameter object of type ParamType
 * `roe_vars: [u, v, H] at roe average state
 * `dq`  : difference between left and right states
 * `nrm` : Normal to face in the physical space
 * `sat` : Simultaneous approximation Term

"""->
function calcSAT{Tmsh, Tsol}(params::ParamType{2},
                             roe_vars::AbstractArray{Tsol, 1},
                             dq::AbstractArray{Tsol,1},
                             nrm::AbstractArray{Tmsh,1},
                             sat::AbstractArray{Tsol,1})
# roe_vars = [u, v, H] at Roe average 

  # SAT parameters
  sat_Vn = convert(Tsol, 0.025)
  sat_Vl = convert(Tsol, 0.025)
  tau = 1.0

  u = roe_vars[1]
  v = roe_vars[2]
  H = roe_vars[3]

  gami = params.gamma_1

  # Begin main executuion
  nx = nrm[1]
  ny = nrm[2]

  dA = sqrt(nx*nx + ny*ny)

  Un = u*nx + v*ny # Normal Velocity

  phi = 0.5*(u*u + v*v)

  a = sqrt(gami*(H - phi)) # speed of sound

  lambda1 = Un + dA*a
  lambda2 = Un - dA*a
  lambda3 = Un

  rhoA = absvalue(Un) + dA*a

  # Compute Eigen Values of the Flux Jacobian
  # The eigen values calculated above cannot be used directly. Near stagnation
  # points lambda3 approaches zero while near sonic lines lambda1 and lambda2
  # approach zero. This has a possibility of creating numerical difficulties.
  # As a result, the eigen values are limited by the following expressions.

  lambda1 = 0.5*(tau*max(absvalue(lambda1),sat_Vn *rhoA) - lambda1)
  lambda2 = 0.5*(tau*max(absvalue(lambda2),sat_Vn *rhoA) - lambda2)
  lambda3 = 0.5*(tau*max(absvalue(lambda3),sat_Vl *rhoA) - lambda3)


  dq1 = dq[1]
  dq2 = dq[2]
  dq3 = dq[3]
  dq4 = dq[4]

  sat[1] = lambda3*dq1
  sat[2] = lambda3*dq2
  sat[3] = lambda3*dq3
  sat[4] = lambda3*dq4


  E1dq = params.res_vals1
  E2dq = params.res_vals2

  #-- get E1*dq
  E1dq[1] = phi*dq1 - u*dq2 - v*dq3 + dq4
  E1dq[2] = E1dq[1]*u
  E1dq[3] = E1dq[1]*v
  E1dq[4] = E1dq[1]*H

  #-- get E2*dq
  E2dq[1] = 0.0
  E2dq[2] = -Un*dq1 + nx*dq2 + ny*dq3
  E2dq[3] = E2dq[2]*ny
  E2dq[4] = E2dq[2]*Un
  E2dq[2] = E2dq[2]*nx

  #-- add to sat
  tmp1 = 0.5*(lambda1 + lambda2) - lambda3
  tmp2 = gami/(a*a)
  tmp3 = 1.0/(dA*dA)

  for i=1:length(sat)
    sat[i] = sat[i] + tmp1*(tmp2*E1dq[i] + tmp3*E2dq[i])
  end

  #-- get E3*dq
  E1dq[1] = -Un*dq1 + nx*dq2 + ny*dq3
  E1dq[2] = E1dq[1]*u
  E1dq[3] = E1dq[1]*v
  E1dq[4] = E1dq[1]*H

  #-- get E4*dq
  E2dq[1] = 0.0
  E2dq[2] = phi*dq1 - u*dq2 - v*dq3 + dq4
  E2dq[3] = E2dq[2]*ny
  E2dq[4] = E2dq[2]*Un
  E2dq[2] = E2dq[2]*nx

  #-- add to sat
  tmp1 = 0.5*(lambda1 - lambda2)/(dA*a)
  for i=1:length(sat)
    sat[i] = sat[i] + tmp1*(E1dq[i] + gami*E2dq[i])
  end

  return nothing
end  # End function calcSAT


function calcSAT{Tsol, Tmsh}(params::ParamType{3},
                             roe_vars::AbstractArray{Tsol, 1},
                             dq::AbstractArray{Tsol,1},
                             nrm::AbstractArray{Tmsh,1},
                             sat::AbstractArray{Tsol,1})
  # roe_vars = [u, v, w, H] at Roe average 

  # SAT parameters
  sat_Vn = convert(Tsol, 0.025)
  sat_Vl = convert(Tsol, 0.025)
  tau = 1.0

  u = roe_vars[1]
  v = roe_vars[2]
  w = roe_vars[3]
  H = roe_vars[4]

  gami = params.gamma_1

  # Begin main executuion
  nx = nrm[1]
  ny = nrm[2]
  nz = nrm[3]

  dA = sqrt(nx*nx + ny*ny + nz*nz)

  phi = 0.5*(u*u + v*v + w*w)
  a = sqrt(gami*(H - phi))

  Un = u*nx + v*ny + w*nz


  lambda1 = Un + dA*a
  lambda2 = Un - dA*a
  lambda3 = Un
  rhoA = absvalue(Un) + dA*a

  lambda1 = 0.5*(tau*max(absvalue(lambda1),sat_Vn *rhoA) - lambda1)
  lambda2 = 0.5*(tau*max(absvalue(lambda2),sat_Vn *rhoA) - lambda2)
  lambda3 = 0.5*(tau*max(absvalue(lambda3),sat_Vl *rhoA) - lambda3)

  dq1 = dq[1]
  dq2 = dq[2]
  dq3 = dq[3]
  dq4 = dq[4]
  dq5 = dq[5]

  E1dq = params.res_vals1
  E2dq = params.res_vals2



  #-- diagonal matrix multiply
#  sat = zeros(Tres, 4)
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
  E2dq[1] = 0.0
  E2dq[2] = -Un*dq1 + nx*dq2 + ny*dq3 + nz*dq4
  E2dq[3] = E2dq[2]*ny
  E2dq[4] = E2dq[2]*nz
  E2dq[5] = E2dq[2]*Un
  E2dq[2] = E2dq[2]*nx

  #-- add to sat
  tmp1 = 0.5*(lambda1 + lambda2) - lambda3
  tmp2 = gami/(a*a)
  tmp3 = 1.0/(dA*dA)

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
  E2dq[1] = 0.0
  E2dq[2] = phi*dq1 - u*dq2 - v*dq3 - w*dq4 + dq5
  E2dq[3] = E2dq[2]*ny
  E2dq[4] = E2dq[2]*nz
  E2dq[5] = E2dq[2]*Un
  E2dq[2] = E2dq[2]*nx

  t1 = phi*dq1 - u*dq2 - v*dq3 - w*dq4 + dq5 # debug

  #-- add to sat
  tmp1 = 0.5*(lambda1 - lambda2)/(dA*a)

  for i=1:5
    sat[i] = sat[i] + tmp1*(E1dq[i] + gami*E2dq[i])
  end

  return nothing
end  # end function calcSAT




@doc """
###EulerEquationMod.calcSAT_revm

Reverse mode of calcSAT

**Inputs**
* `params` : Parameter object of type ParamType
* `nrm` : Normal to face in the physical space
* `dq`  : Boundary condition penalty variable
* `vel` : Velocities along X & Y directions in the physical space
* `H`   : Total enthalpy
* `sat_bar` : Inpute seed for sat flux whose derivative needs to be computed

**Output**

* `nrm_bar` : derivative of `sat_bar` w.r.t physical normal vector

"""->

function calcSAT_revm{Tmsh, Tsol}(params::ParamType{2}, nrm::AbstractArray{Tmsh,1},
                      dq::AbstractArray{Tsol,1}, vel::AbstractArray{Tsol, 1},
                      H::Tsol, sat_bar::AbstractArray{Tsol, 1},
                      nrm_bar::AbstractArray{Tmsh,1})

  # Forward Sweep
  sat_Vn = convert(Tsol, 0.025)
  sat_Vl = convert(Tsol, 0.025)
  u = vel[1]
  v = vel[2]
  gami = params.gamma_1

  # Begin main executuion
  nx = nrm[1]
  ny = nrm[2]
  dA = sqrt(nx*nx + ny*ny)
  Un = u*nx + v*ny # Normal Velocity
  phi = 0.5*(u*u + v*v)
  a = sqrt(gami*(H - phi)) # speed of sound
  lambda1 = Un + dA*a
  lambda2 = Un - dA*a
  lambda3 = Un
  rhoA = absvalue(Un) + dA*a

  lambda1 = 0.5*(max(absvalue(lambda1),sat_Vn *rhoA) - lambda1)
  lambda2 = 0.5*(max(absvalue(lambda2),sat_Vn *rhoA) - lambda2)
  lambda3 = 0.5*(max(absvalue(lambda3),sat_Vl *rhoA) - lambda3)

  dq1 = dq[1]
  dq2 = dq[2]
  dq3 = dq[3]
  dq4 = dq[4]

  sat = zeros(Tsol, 4)
  sat[1] = lambda3*dq1
  sat[2] = lambda3*dq2
  sat[3] = lambda3*dq3
  sat[4] = lambda3*dq4

  E1dq = params.res_vals1
  E2dq = params.res_vals2
  E3dq = zeros(Tsol, 4)
  E4dq = zeros(Tsol, 4)

  #-- get E1*dq
  E1dq[1] = phi*dq1 - u*dq2 - v*dq3 + dq4
  E1dq[2] = E1dq[1]*u
  E1dq[3] = E1dq[1]*v
  E1dq[4] = E1dq[1]*H

  #-- get E2*dq
  E2dq[1] = 0.0
  E2dq[2] = -Un*dq1 + nx*dq2 + ny*dq3
  E2dq[3] = E2dq[2]*ny
  E2dq[4] = E2dq[2]*Un
  E2dq[2] = E2dq[2]*nx

  #-- add to sat
  tmp1 = 0.5*(lambda1 + lambda2) - lambda3
  tmp2 = gami/(a*a)
  tmp3 = 1.0/(dA*dA)

  #-- get E3*dq
  E3dq[1] = -Un*dq1 + nx*dq2 + ny*dq3
  E3dq[2] = E3dq[1]*u
  E3dq[3] = E3dq[1]*v
  E3dq[4] = E3dq[1]*H

  #-- get E4*dq
  E4dq[1] = 0.0
  E4dq[2] = phi*dq1 - u*dq2 - v*dq3 + dq4
  E4dq[3] = E4dq[2]*ny
  E4dq[4] = E4dq[2]*Un
  E4dq[2] = E4dq[2]*nx

  #-- add to sat
  tmp1 = 0.5*(lambda1 - lambda2)/(dA*a)

  # Reverse sweep
  E3dq_bar = zeros(Tsol, 4)
  E4dq_bar = zeros(Tsol, 4)
  tmp1_bar = zero(Tsol)
  for i = 1:length(sat_bar)
    E3dq_bar[i] += sat_bar[i]*tmp1
    E4dq_bar[i] = sat_bar[i]*tmp1*gami
    tmp1_bar += sat_bar[i]*(E3dq[i] + gami*E4dq[i])
  end

  # tmp1 = 0.5*(lambda1 - lambda2)/(dA*a)
  lambda1_bar = (0.5/(dA*a))*tmp1_bar
  lambda2_bar = -(0.5/(dA*a))*tmp1_bar
  dA_bar = -(0.5*(lambda1 - lambda2)/(dA*dA*a))*tmp1_bar


  #  E4dq[2] = phi*dq1 - u*dq2 - v*dq3 + dq4
  #  E4dq[3] = E4dq[2]*ny
  #  E4dq[4] = E4dq[2]*Un
  #  E4dq[2] = E4dq[2]*nx
  nx_bar = E4dq_bar[2]*(phi*dq1 - u*dq2 - v*dq3 + dq4)
  Un_bar = E4dq_bar[4]*(phi*dq1 - u*dq2 - v*dq3 + dq4)
  ny_bar = E4dq_bar[3]*(phi*dq1 - u*dq2 - v*dq3 + dq4)

  # E3dq[1] = -Un*dq1 + nx*dq2 + ny*dq3
  # E3dq[2] = E3dq[1]*u
  # E3dq[3] = E3dq[1]*v
  # E3dq[4] = E3dq[1]*H
  E3dq_bar[1] += E3dq_bar[4]*H + E3dq_bar[3]*v + E3dq_bar[2]*u
  Un_bar += -E3dq_bar[1]*dq1
  nx_bar += E3dq_bar[1]*dq2
  ny_bar += E3dq_bar[1]*dq3

  #  for i=1:length(sat)
  #    sat[i] = sat[i] + tmp1*(tmp2*E1dq[i] + tmp3*E2dq[i])
  #  end
  tmp1 = 0.5*(lambda1 + lambda2) - lambda3 # For E1dq & E2dq matrices
  tmp1_bar = zero(Tsol) # Reset tmp1_bar to 0 since the variable is being reused
  E1dq_bar = zeros(Tsol, 4)
  E2dq_bar = zeros(Tsol, 4)
  tmp3_bar = zero(Tsol)
  for i = 1:length(sat_bar)
    E1dq_bar[i] += sat_bar[i]*tmp1*tmp2
    E2dq_bar[i] += sat_bar[i]*tmp1*tmp3
    tmp1_bar += sat_bar[i]*(tmp2*E1dq[i] + tmp3*E2dq[i])
    tmp3_bar += sat_bar[i]*tmp1*E2dq[i]
  end

  # tmp1 = 0.5*(lambda1 + lambda2) - lambda3
  # tmp2 = gami/(a*a)
  # tmp3 = 1.0/(dA*dA)
  dA_bar += -2*tmp3_bar/(dA^3)
  lambda1_bar += tmp1_bar*0.5
  lambda2_bar += tmp1_bar*0.5
  lambda3_bar = -tmp1_bar

  #  E2dq[1] = 0.0
  #  intvar = -Un*dq1 + nx*dq2 + ny*dq3
  #  E2dq[2] = intvar*nx
  #  E2dq[3] = intvar*ny
  #  E2dq[4] = intvar*Un
  intvar = -Un*dq1 + nx*dq2 + ny*dq3
  Un_bar += E2dq_bar[4]*intvar
  ny_bar += E2dq_bar[3]*intvar
  nx_bar += E2dq_bar[2]*intvar
  intvar_bar = E2dq_bar[4]*Un + E2dq_bar[3]*ny + E2dq_bar[2]*nx
  Un_bar -= intvar_bar*dq1
  nx_bar += intvar_bar*dq2
  ny_bar += intvar_bar*dq3

  # None of the following need to be differentiated as they dont depen on nrm
  # E1dq[1] = phi*dq1 - u*dq2 - v*dq3 + dq4
  # E1dq[2] = E1dq[1]*u
  # E1dq[3] = E1dq[1]*v
  # E1dq[4] = E1dq[1]*H

  # sat[1] = lambda3*dq1
  # sat[2] = lambda3*dq2
  # sat[3] = lambda3*dq3
  # sat[4] = lambda3*dq4
  lambda3_bar += dq1*sat_bar[1] + dq2*sat_bar[2] + dq3*sat_bar[3] + dq4*sat_bar[4]

  rhoA_bar = zero(Tsol)
  # lambda3 = 0.5*(max(absvalue(lambda3),sat_Vl *rhoA) - lambda3)
  # Breaking the above down. lambda3 on RHS = Un
  # intVar1 = = absvalue(Un)
  # intVar2 = sat_Vl*rhoA
  # intVar3 = max(intVar1, intVar2)
  # lambda3 = 0.5*(intVar3 - Un)
  intVar3_bar = 0.5*lambda3_bar
  Un_bar -= 0.5*lambda3_bar
  intVar1_bar, intVar2_bar = max_deriv_rev(absvalue(Un), sat_Vl*rhoA, intVar3_bar)
  rhoA_bar += sat_Vl*intVar2_bar
  Un_bar += intVar1_bar*absvalue_deriv(Un)

  L2 = Un - dA*a
  # lambda2 = 0.5*(max(absvalue(L2),sat_Vn *rhoA) - L2)
  # intVar1 = absvalue(L2)
  # intVar2 = sat_Vn*rhoA
  # intVar3 = max(intVar1, intVar2)
  # lambda2 = 0.5*(intVar3 - L2)
  intVar3_bar = 0.5*lambda2_bar
  L2_bar = -0.5*lambda2_bar
  intVar1_bar, intVar2_bar = max_deriv_rev(absvalue(L2), sat_Vn*rhoA, intVar3_bar)
  rhoA_bar += sat_Vn*intVar2_bar
  L2_bar += intVar1_bar*absvalue_deriv(L2)

  L1 = Un + dA*a
  # lambda1 = 0.5*(max(absvalue(L1),sat_Vn *rhoA) - L1)
  # intVar1 = absvalue(L1)
  # intVar2 = sat_Vn*rhoA
  # intVar3 = max(intVar1, intVar2)
  # lambda1 = 0.5*(intVar3 - L1)
  intVar3_bar = 0.5*lambda1_bar
  L1_bar = -0.5*lambda1_bar
  intVar1_bar, intVar2_bar = max_deriv_rev(absvalue(Un + dA*a), sat_Vn*rhoA, intVar3_bar)
  rhoA_bar += intVar2_bar*sat_Vn
  L1_bar += intVar1_bar*absvalue_deriv(L1)

  # rhoA = absvalue(Un) + dA*a
  dA_bar += rhoA_bar*a
  if real(Un) >= 0.0
    Un_bar += rhoA_bar
  else
    Un_bar -= rhoA_bar
  end

  # lambda1 = Un + dA*a
  # lambda2 = Un - dA*a
  # lambda3 = Un
  Un_bar += L1_bar + L2_bar #  + lambda3_bar
  dA_bar += a*L1_bar - a*L2_bar

  # nx = nrm[1]
  # ny = nrm[2]
  # dA = sqrt(nx*nx + ny*ny)
  # Un = u*nx + v*ny
  nx_bar += u*Un_bar
  ny_bar += v*Un_bar
  nx_bar += nx*dA_bar/dA
  ny_bar += ny*dA_bar/dA
  nrm_bar[1] += nx_bar
  nrm_bar[2] += ny_bar

  return nothing
end

"""
  Calculates the Lax-Friedrich flux function on the conservative variables
"""
function calcLFFlux{Tmsh, Tsol, Tres, Tdim}(
                      params::ParamType{Tdim, :conservative},
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tsol, 1},
                      dir::AbstractArray{Tmsh, 1},  F::AbstractArray{Tres,1})

  # compute Euler flux of left and right states
  fluxL = params.flux_vals1
  fluxR = params.flux_vals2
  calcEulerFlux(params, qL, aux_vars, dir, fluxL)
  calcEulerFlux(params, qR, aux_vars, dir, fluxR)
  lambda_max = getLambdaMaxSimple(params, qL, qR, dir)

  for i=1:length(F)
    F[i] = 0.5*(fluxL[i] + fluxR[i] - lambda_max*(qR[i] - qL[i]))
  end

  return nothing
end

#=
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
=#


function calcEulerFlux_standard{Tmsh, Tsol, Tres}(
                      params::ParamType{2, :conservative},
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tsol, 1},
                      dir::AbstractArray{Tmsh, 1},  F::AbstractArray{Tres,1})
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

"""
  This function computes the split form flux corresponding to the standard
  flux in Tdim directions at once for a given state.  This is more efficient
  Than calling the single direction method Tdim times.  Methods are available
  for 2 and 3 dimensions.

  Inputs:
    params: ParamType
    qL: left state vector
    qR: right state vector
    aux_vars: auxiliary variable vector for qL
    dir: a Tdim x Tdim matrix with each column containing a normal vector

  Inputs/Outputs:
    F: a numDofPerNode x Tdim matrix where each column will be populated with
       the flux in the direction specified by the corresponding column of nrm

  Aliasing restrictions: none
"""
function calcEulerFlux_standard{Tmsh, Tsol, Tres}(
                      params::ParamType{2, :conservative},
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tsol, 1},
                      dir::AbstractArray{Tmsh, 2},  F::AbstractArray{Tres,2})
# calculate the split form numerical flux function corresponding to the standard DG flux

  pL = calcPressure(params, qL); pR = calcPressure(params, qR)
  rho_avg = 0.5*(qL[1] + qR[1])
  rhou_avg = 0.5*(qL[2] + qR[2])
  rhov_avg = 0.5*(qL[3] + qR[3])
  p_avg = 0.5*(pL + pR)
  rhoLinv = 1/qL[1]; rhoRinv = 1/qR[1]


  for i=1:2
    F[1, i] = dir[1, i]*(rhou_avg) + dir[2, i]*rhov_avg

    tmp1 = 0.5*(qL[2]*qL[2]*rhoLinv + qR[2]*qR[2]*rhoRinv)
    tmp2 = 0.5*(qL[2]*qL[3]*rhoLinv + qR[2]*qR[3]*rhoRinv)
    F[2, i] = dir[1, i]*(tmp1 + p_avg) + dir[2, i]*tmp2

    tmp1 = 0.5*(qL[2]*qL[3]*rhoLinv + qR[2]*qR[3]*rhoRinv)
    tmp2 = 0.5*(qL[3]*qL[3]*rhoLinv + qR[3]*qR[3]*rhoRinv)
    F[3, i] = dir[1, i]*tmp1 + dir[2, i]*(tmp2 + p_avg)


    tmp1 = 0.5*( (qL[4] + pL)*qL[2]*rhoLinv + (qR[4] + pR)*qR[2]*rhoRinv)
    tmp2 = 0.5*( (qL[4] + pL)*qL[3]*rhoLinv + (qR[4] + pR)*qR[3]*rhoRinv)
    F[4, i] = dir[1, i]*tmp1 + dir[2, i]*tmp2
  end

  return nothing
end


function calcEulerFlux_standard{Tmsh, Tsol, Tres}(
                      params::ParamType{3, :conservative},
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tres, 1},
                      dir::AbstractArray{Tmsh, 1},  F::AbstractArray{Tres,1})
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

function calcEulerFlux_standard{Tmsh, Tsol, Tres}(
                      params::ParamType{3, :conservative},
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tres, 1},
                      dir::AbstractArray{Tmsh, 2},  F::AbstractArray{Tres,2})
# calculate the split form numerical flux function corresponding to the standard DG flux
#TODO: pre-calculate 1/qL[1], 1/qR[1]

  pL = calcPressure(params, qL); pR = calcPressure(params, qR)
  rho_avg = 0.5*(qL[1] + qR[1])
  rhou_avg = 0.5*(qL[2] + qR[2])
  rhov_avg = 0.5*(qL[3] + qR[3])
  rhow_avg = 0.5*(qL[4] + qR[4])
  p_avg = 0.5*(pL + pR)
  rhoLinv = 1/qL[1]; rhoRinv = 1/qR[1]

  for i=1:3
    F[1, i] = dir[1, i]*(rhou_avg) + dir[2, i]*rhov_avg + dir[3, i]*rhow_avg

    tmp1 = 0.5*(qL[2]*qL[2]*rhoLinv + qR[2]*qR[2]*rhoRinv)
    tmp2 = 0.5*(qL[2]*qL[3]*rhoLinv + qR[2]*qR[3]*rhoRinv)
    tmp3 = 0.5*(qL[2]*qL[4]*rhoLinv + qR[2]*qR[4]*rhoRinv)
    F[2, i] = dir[1, i]*(tmp1 + p_avg) + dir[2, i]*tmp2 + dir[3, i]*tmp3

    tmp1 = 0.5*(qL[2]*qL[3]*rhoLinv + qR[2]*qR[3]*rhoRinv)
    tmp2 = 0.5*(qL[3]*qL[3]*rhoLinv + qR[3]*qR[3]*rhoRinv)
    tmp3 = 0.5*(qL[4]*qL[3]*rhoLinv + qR[4]*qR[3]*rhoRinv)
    F[3, i] = dir[1, i]*tmp1 + dir[2, i]*(tmp2 + p_avg) + dir[3, i]*tmp3

    tmp1 = 0.5*(qL[2]*qL[4]*rhoLinv + qR[2]*qR[4]*rhoRinv)
    tmp2 = 0.5*(qL[3]*qL[4]*rhoLinv + qR[3]*qR[4]*rhoRinv)
    tmp3 = 0.5*(qL[4]*qL[4]*rhoLinv + qR[4]*qR[4]*rhoRinv)
    F[4, i] = dir[1, i]*tmp1 + dir[2, i]*tmp2 + dir[3, i]*(tmp3 + p_avg)


    tmp1 = 0.5*( (qL[5] + pL)*qL[2]*rhoLinv + (qR[5] + pR)*qR[2]*rhoRinv)
    tmp2 = 0.5*( (qL[5] + pL)*qL[3]*rhoLinv + (qR[5] + pR)*qR[3]*rhoRinv)
    tmp3 = 0.5*( (qL[5] + pL)*qL[4]*rhoLinv + (qR[5] + pR)*qR[4]*rhoRinv)
    F[5, i] = dir[1, i]*tmp1 + dir[2, i]*tmp2 + dir[3, i]*tmp3
  end

  return nothing
end



#=
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
=#

"""
  Calculates the numerical flux function associated with the Ducros flux
  splitting.  Methods are available for 2D and 3D.

  **Inputs**:

   * params:
   * qL: the left state
   * qR: the right state
   * aux_vars: the aux vars for the left state
   * dir: the direction vector


  **Inputs/Outputs**:

   * F: vector to be populated with the flux

"""
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

#=
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
=#

"""
  This function calculates the Ismail-Roe numerical flux at a node in a
  specified direction

  **Inputs**:

   * params: ParamType
   * qL: left state vector
   * qR: right state vector
   * aux_vars: auxiliary variable vector for qL
   * dir: a direction vector of length Tdim

  **Inputs/Outputs**:

   * F: a numDofPerNode x Tdim matrix where each column will be populated with
       the flux in the direction specified by the corresponding column of nrm

  Aliasing Restrictions: none
"""
function calcEulerFlux_IR{Tmsh, Tsol, Tres}(params::ParamType{2, :conservative},
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tres},
                      dir::AbstractArray{Tmsh, 1},  F::AbstractArray{Tres,1})

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


"""
  This function computes the Ismail-Roe
  flux in Tdim directions at once for a given state.  This is more efficient
  Than calling the single direction method Tdim times.  Methods are available
  for 2 and 3 dimensions.

  **Inputs**:

   * params: ParamType
   * qL: left state vector
   * qR: right state vector
   * aux_vars: auxiliary variable vector for qL
   * dir: a Tdim x Tdim matrix with each column containing a normal vector

  **Inputs/Outputs**:

   * F: a numDofPerNode x Tdim matrix where each column will be populated with
       the flux in the direction specified by the corresponding column of nrm

  Aliasing restrictions: none
"""
function calcEulerFlux_IR{Tmsh, Tsol, Tres}(params::ParamType{2, :conservative},
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tres},
                      dir::AbstractArray{Tmsh, 2},  F::AbstractArray{Tres,2})

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

  for i=1:2
    mv_n = rho_hat*(dir[1, i]*u_hat + dir[2, i]*v_hat)  # normal momentum
    F[1, i] = mv_n
    F[2, i] = mv_n*u_hat + dir[1, i]*p1_hat
    F[3, i] = mv_n*v_hat + dir[2, i]*p1_hat
    F[4, i] = mv_n*h_hat
  end
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
                      dir::AbstractArray{Tmsh, 1},  F::AbstractArray{Tres,1})

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


# multi-dimension version
function calcEulerFlux_IR{Tmsh, Tsol, Tres}(params::ParamType{3, :conservative},
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tres},
                      dir::AbstractArray{Tmsh, 2},  F::AbstractArray{Tres,2})

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

  @simd for i=1:3
    mv_n = rho_hat*(dir[1,i]*u_hat + dir[2, i]*v_hat + dir[3,i]*w_hat)
    F[1, i] = mv_n
    F[2, i] = mv_n*u_hat + dir[1, i]*p1_hat
    F[3, i] = mv_n*v_hat + dir[2, i]*p1_hat
    F[4, i] = mv_n*w_hat + dir[3, i]*p1_hat
    F[5, i] = mv_n*h_hat
  end

  return nothing
end


#TODO: move documentation to second method
# stabilized IR flux
#=
"""
  This function calculates the flux across an interface using the IR
  numerical flux function and a Lax-Friedrich type of entropy dissipation.

  Currently this is implemented for conservative variables only.

  Methods are available that take in dxidx and a normal vector in parametric
  space and compute and normal vector xy space and that take in a
  normal vector directly.

  Inputs:
    qL, qR: vectors conservative variables at left and right states
    aux_vars: aux_vars for qL
    dxidx: scaled mapping jacobian (2x2 or 3x3 in 3d)
    nrm: normal vector in parametric space

  Inputs/Outputs:
    F: vector to be updated with the result

  Aliasing restrictions:
    nothing may alias params.nrm2.  See also getEntropyLFStab
"""
=#
#=
function calcEulerFlux_IRSLF{Tmsh, Tsol, Tres}(params::ParamType,
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tres},
                      dxidx::AbstractMatrix{Tmsh},
                      nrm::AbstractArray{Tmsh},  F::AbstractArray{Tres,1})

  nrm2 = params.nrm2
  calcBCNormal(params, dxidx, nrm, nrm2)
  calcEulerFlux_IRSLF(params, qL, qR, aux_vars, nrm2, F)
  return nothing
end
=#

"""
  This is the second method that takes in a normal vector directly.
  See the first method for a description of what this function does.

  Inputs
    qL, qR: vectors conservative variables at left and right states
    aux_vars: aux_vars for qL
    nrm: a normal vector in xy space

  Inputs/Outputs
    F: vector to be updated with the result

  Alising restrictions:
    See getEntropyLFStab

"""
function calcEulerFlux_IRSLF{Tmsh, Tsol, Tres, Tdim}(
                      params::ParamType{Tdim, :conservative},
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractVector{Tres},
                      dir::AbstractVector{Tmsh},  F::AbstractArray{Tres,1})

  calcEulerFlux_IR(params, qL, qR, aux_vars, dir, F)
  getEntropyLFStab(params, qL, qR, aux_vars, dir, F)

  return nothing
end

"""
  This function is similar to calcEulerFlux_IRSLF, but uses Lax-Wendroff
  dissipation rather than Lax-Friedrich.

  Aliasing restrictions: see getEntropyLWStab
"""
function calcEulerFlux_IRSLW{Tmsh, Tsol, Tres}(params::ParamType,
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tres},
                      dxidx::AbstractMatrix{Tmsh},
                      nrm::AbstractArray{Tmsh},  F::AbstractArray{Tres,1})
  #TODO: remove this method
  nrm2 = params.nrm2
  calcBCNormal(params, dxidx, nrm, nrm2)
  calcEulerFlux_IRSLW(params, qL, qR, aux_vars, nrm2, F)
  return nothing
end

function calcEulerFlux_IRSWF{Tmsh, Tsol, Tres, Tdim}(
                      params::ParamType{Tdim, :conservative},
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractVector{Tres},
                      dir::AbstractVector{Tmsh},  F::AbstractArray{Tres,1})

  calcEulerFlux_IR(params, qL, qR, aux_vars, dir, F)
  getEntropyLWStab(params, qL, qR, aux_vars, dir, F)

  return nothing
end




function logavg(aL, aR)
# calculate the logarithmic average needed by the IR flux
  xi = aL/aR
  f = (xi - 1)/(xi + 1)
  u = f*f
  eps = 1e-3
  if u < eps
    F = @evalpoly( u, 1, 1/3, 1/5, 1/7, 1/9)
#    F = 1.0 + u/3.0 + u*u/5.0 + u*u*u/7.0 + u*u*u*u/9.0
  else
    F = (log(xi)/2.0)/f
  end

  return (aL + aR)/(2*F)
end
