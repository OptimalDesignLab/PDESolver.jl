# differentiated version of functions in bc_solvers.jl


"""
  Differentiated version of the [`RoeSolver`](@ref).  Computes the jacobian
  of the flux with respect to `q` and `qg`.  Methods are available for 2D
  and 3D.

  The caller must zero out the output arrays (if desired) between calls

  **Inputs**

   * params: ParamType, conservative variables only
   * q: vector of conservative variables for the left state
   * qg: vector of conservative variables for the right state
   * aux_vars: auxiliary variables for the left state
   * nrm: scaled normal vector in x-y space (outward wrt the element q lives on

  **Inputs/Outputs**

   * fluxL_dot: flux jacobian wrt `q`, numDofPerNode x numDofPerNode (summed into)
   * fluxR_dot: flux jacobian wrt `qg`, numDofPerNode x numDofPerNode (summed into)

  Aliasing restrictions:

  params: `sat`, `roe_vars`, `roe_vars_dot`, `euler_fluxjac`, `p_dot` must not be used
"""
function RoeSolver_diff(params::ParamType{2, :conservative},
                   q::AbstractArray{Tsol,1},
                   qg::AbstractArray{Tsol, 1},
                   aux_vars::AbstractArray{Tres, 1},
                   nrm::AbstractArray{Tmsh,1},
                   fluxL_dot::AbstractArray{Tres, 2},
                   fluxR_dot::AbstractArray{Tres, 2}) where {Tmsh, Tsol, Tres}

  # SAT terms are used for ensuring consistency with the physical problem. Its
  # similar to upwinding which adds dissipation to the problem. SATs on the
  # boundary can be thought of as having two overlapping nodes and because of
  # the discontinuous nature of SBP adds some dissipation.

  data = params.roefluxdata
  @unpack data dq sat roe_vars euler_flux v_vals nrm2 roe_vars_dot

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
  # All the _dot variables are wrt q
  fac = d1_0/q[1]
  fac_dot1 = -fac*fac
#  fac_dot1 = -d1_0/(q[1]*q[1])

  uL = q[2]*fac
  uL_dot1 = q[2]*fac_dot1
  uL_dot2 = fac

  vL = q[3]*fac
  vL_dot1 = q[3]*fac_dot1
  vL_dot3 = fac

  phi = d0_5*(uL*uL + vL*vL)
  phi_dot1 = uL*uL_dot1 + vL*vL_dot1
  phi_dot2 = uL*uL_dot2
  phi_dot3 = vL*vL_dot3

  HL = gamma*q[4]*fac - gami*phi # Total enthalpy, H = e + 0.5*(u^2 + v^2) + p/rho,
                                 # where e is the internal energy per unit mass
  HL_dot1 = gamma*q[4]*fac_dot1 - gami*phi_dot1
  HL_dot2 = -gami*phi_dot2
  HL_dot3 = -gami*phi_dot3
  HL_dot4 = gamma*fac  # q[4]_dot = 1

  # The right side of the Roe solver comprises the boundary conditions
  # all the _dot variables are wrt qg now
  fac = d1_0/qg[1]
  fac_dot1 = -fac*fac
#  fac_dot1 = -d1_0/(qg[1]*qg[1])

  uR = qg[2]*fac
  uR_dot1 = qg[2]*fac_dot1
  uR_dot2 = fac

  vR = qg[3]*fac
  vR_dot1 = qg[3]*fac_dot1
  vR_dot3 = fac

  phi = d0_5*(uR*uR + vR*vR)
  phi_dot1 = uR*uR_dot1 + vR*vR_dot1
  phi_dot2 = uR*uR_dot2
  phi_dot3 = vR*vR_dot3


  HR = gamma*qg[4]*fac - gami*phi # Total Enthalpy
  HR_dot1 = gamma*qg[4]*fac_dot1 - gami*phi_dot1
  HR_dot2 = -gami*phi_dot2
  HR_dot3 = -gami*phi_dot3
  HR_dot4 = gamma*fac

  # Averaged states
  sqL = sqrt(q[1])
  sqL_dot1 = 0.5/sqL

  sqR = sqrt(qg[1])
  sqR_dot1 = 0.5/sqR

  fac = d1_0/(sqL + sqR)
  t1 = -1/((sqL + sqR)*(sqL + sqR))
  fac_dotL1 = t1*sqL_dot1
  fac_dotR1 = t1*sqR_dot1

  u = (sqL*uL + sqR*uR)*fac
  t2 = sqR*uR
  t3 = sqL*uL
  t4 = sqL*fac
  t5 = sqR*fac
  u_dotL1 = t3*fac_dotL1 + sqL*fac*uL_dot1 + uL*fac*sqL_dot1 + t2*fac_dotL1
  u_dotR1 = t2*fac_dotR1 + sqR*fac*uR_dot1 + uR*fac*sqR_dot1 + t3*fac_dotR1


  u_dotL2 = t4*uL_dot2
  u_dotR2 = t5*uR_dot2

  v = (sqL*vL + sqR*vR)*fac
  t2 = sqL*vL
  t3 = sqR*vR
  v_dotL1 = t2*fac_dotL1 + sqL*fac*vL_dot1 + vL*fac*sqL_dot1 + t3*fac_dotL1
  v_dotR1 = t3*fac_dotR1 + sqR*fac*vR_dot1 + vR*fac*sqR_dot1 + t2*fac_dotR1

  v_dotL3 = t4*vL_dot3
  v_dotR3 = t5*vR_dot3

  H = (sqL*HL + sqR*HR)*fac
  t2 = sqL*HL
  t3 = sqR*HR
  H_dotL1 = t2*fac_dotL1 + sqL*fac*HL_dot1 + HL*fac*sqL_dot1 + t3*fac_dotL1
  H_dotR1 = t3*fac_dotR1 + sqR*fac*HR_dot1 + HR*fac*sqR_dot1 + t2*fac_dotR1
 
  H_dotL2 = t4*HL_dot2 
  H_dotR2 = t5*HR_dot2

  H_dotL3 = t4*HL_dot3
  H_dotR3 = t5*HR_dot3

  H_dotL4 = t4*HL_dot4
  H_dotR4 = t5*HR_dot4


  for i=1:length(dq)
    dq[i] = q[i] - qg[i]
  end

  # dq_dotL* = 1, dq_dotR* = -1, so omit them

  roe_vars[1] = u
  roe_vars[2] = v
  roe_vars[3] = H

  roe_vars_dot[1]  = u_dotL1
  roe_vars_dot[2]  = u_dotR1
  roe_vars_dot[3]  = u_dotL2
  roe_vars_dot[4]  = u_dotR2

  roe_vars_dot[5]  = v_dotL1
  roe_vars_dot[6]  = v_dotR1
  roe_vars_dot[7]  = v_dotL3
  roe_vars_dot[8]  = v_dotR3

  roe_vars_dot[9]  = H_dotL1
  roe_vars_dot[10] = H_dotR1
  roe_vars_dot[11] = H_dotL2
  roe_vars_dot[12] = H_dotR2
  roe_vars_dot[13] = H_dotL3
  roe_vars_dot[14] = H_dotR3
  roe_vars_dot[15] = H_dotL4
  roe_vars_dot[16] = H_dotR4
  

  # pass in fluxL_dot and fluxR_dot here, then add the Euler flux
  # contribution below
  calcSAT_diff(params, roe_vars, roe_vars_dot,  dq, nrm, fluxL_dot, fluxR_dot)

#  calcSAT(params, nrm, dq, sat, u, v, H, use_efix)
  
  nrm2[1] = nx   # why are we assigning to nrm2?
  nrm2[2] = ny

#  convertFromNaturalToWorkingVars(params, q, v_vals)
  calcEulerFlux_diff(params, q, aux_vars, nrm2, fluxL_dot)


  return nothing

end # ends the function RoeSolver


function RoeSolver_revq(params::ParamType{2, :conservative},
                   q::AbstractArray{Tsol,1},
                   q_bar::AbstractArray{Tsol, 1},
                   qg::AbstractArray{Tsol, 1},
                   qg_bar::AbstractArray{Tsol, 1},
                   aux_vars::AbstractArray{Tres, 1},
                   nrm::AbstractArray{Tmsh,1},
                   flux_bar::AbstractArray{Tres, 1}) where {Tmsh, Tsol, Tres}

  # SAT terms are used for ensuring consistency with the physical problem. Its
  # similar to upwinding which adds dissipation to the problem. SATs on the
  # boundary can be thought of as having two overlapping nodes and because of
  # the discontinuous nature of SBP adds some dissipation.

  data = params.roefluxdata
  @unpack data dq sat roe_vars euler_flux v_vals nrm2

  # Declaring constants
  d1_0 = 1.0
  d0_0 = 0.0
  d0_5 = 0.5
  tau = 1.0
  gamma = params.gamma
  gami = params.gamma_1
  sat_fac = 1  # multiplier for SAT term

  # Begin main execution
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

  for i=1:length(dq)
    dq[i] = q[i] - qg[i]
  end

  roe_vars[1] = u
  roe_vars[2] = v
  roe_vars[3] = H

#  calcSAT(params, roe_vars, dq, nrm, sat)

  # calculate Euler flux in wall normal directiona
  # because edge numbering is rather arbitary, any memory access is likely to
  # be a cache miss, so we recalculate the Euler flux
#  nrm2[1] = nx   # why are we assigning to nrm2?
#  nrm2[2] = ny

#  convertFromNaturalToWorkingVars(params, q, v_vals)
#  calcEulerFlux(params, v_vals, aux_vars, nrm2, euler_flux)

#  for i=1:4  # ArrayViews does not support flux[:] = .
#    flux[i] = (sat_fac*sat[i] + euler_flux[i])
#  end

  #----------------------------------------------------------------------------
  # reverse sweep

  numDofPerNode = length(dq)
  @unpack data sat_bar v_vals_bar roe_vars_bar dq_bar

  fill!(sat_bar, 0); fill!(roe_vars_bar, 0); fill!(dq_bar, 0)
  fill!(v_vals_bar, 0)
  for i=1:4
    sat_bar[i] = sat_fac*flux_bar[i]
    v_vals[i] = q[i]
  end
  euler_flux_bar = flux_bar

  calcEulerFlux_revq(params, v_vals, v_vals_bar, aux_vars, nrm, euler_flux_bar)

  #TODO: this is where convertFromNatualToWorkingVars_rev should be
  for i=1:numDofPerNode
    q_bar[i] += v_vals_bar[i]
  end

  calcSAT_revq(params, roe_vars, roe_vars_bar, dq, dq_bar, nrm, sat_bar)

  for i=1:length(dq)
    q_bar[i]  += dq_bar[i]
    qg_bar[i] -= dq_bar[i]
  end

  u_bar = roe_vars_bar[1]
  v_bar = roe_vars_bar[2]
  H_bar = roe_vars_bar[3]

  # averaged states
  sqL_bar = HL*fac*H_bar; sqR_bar = HR*fac*H_bar
  HL_bar  = sqL*fac*H_bar; HR_bar = sqR*fac*H_bar
  fac_bar = H_bar*(sqL*HL + sqR*HR)

  sqL_bar += vL*fac*v_bar; sqR_bar += vR*fac*v_bar
  vL_bar   = sqL*fac*v_bar; vR_bar  = sqR*fac*v_bar
  fac_bar += (sqL*vL + sqR*vR)*v_bar

  sqL_bar += uL*fac*u_bar; sqR_bar += uR*fac*u_bar
  uL_bar   = sqL*fac*u_bar; uR_bar  = sqR*fac*u_bar
  fac_bar += (sqL*uL + sqR*uR)*u_bar

  sqL_bar += -fac*fac*fac_bar; sqR_bar += -fac*fac*fac_bar
  q_bar[1] += 0.5/sqL*sqL_bar; qg_bar[1] += 0.5/sqR*sqR_bar

  # right state
  fac = d1_0/qg[1]
  fac_bar = zero(fac_bar)  # this is a different variable in single assignment
  qg_bar[4] += gamma*fac*HR_bar
  fac_bar += gamma*qg[4]*HR_bar
  phi_bar = -gami*HR_bar

  uR_bar += uR*phi_bar
  vR_bar += vR*phi_bar

  qg_bar[2] += fac*uR_bar; qg_bar[3] += fac*vR_bar
  fac_bar += qg[2]*uR_bar; fac_bar += qg[3]*vR_bar

  qg_bar[1] += -fac*fac*fac_bar

  # left state
  fac = d1_0/q[1]
  fac_bar = zero(fac_bar); phi_bar = zero(phi_bar)
  q_bar[4] += gamma*fac*HL_bar
  fac_bar += gamma*q[4]*HL_bar
  phi_bar += -gami*HL_bar

  uL_bar += uL*phi_bar
  vL_bar += vL*phi_bar

  q_bar[2] += fac*uL_bar; q_bar[3] += fac*vL_bar
  fac_bar += q[2]*uL_bar; fac_bar += q[3]*vL_bar
  q_bar[1] += -fac*fac*fac_bar


  return nothing
end





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

function RoeSolver_revm(params::ParamType{2},
                   q::AbstractArray{Tsol,1},
                   qg::AbstractArray{Tsol, 1},
                   aux_vars::AbstractArray{Tres, 1},
                   nrm::AbstractArray{Tmsh,1},
                   nrm_bar::AbstractArray{Tmsh, 1},
                   flux_bar::AbstractArray{Tres, 1},
                   ) where {Tmsh, Tsol, Tres}

  data = params.roefluxdata
  @unpack data dq sat roe_vars euler_flux v_vals nrm2 euler_flux_bar sat_bar

  # Forward sweep
  tau = 1.0
  gamma = params.gamma
  gami = params.gamma_1
  sat_fac = 1  # multiplier for SAT term

  # Begin main execution
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
  for i=1:length(dq)
    dq[i] = q[i] - qg[i]
  end

  nrm2[1] = nx
  nrm2[2] = ny
  convertFromNaturalToWorkingVars(params, q, v_vals)

  # Reverse Sweep
  # for i=1:4  # ArrayViews does not support flux[:] = .
  #   flux[i] = (sat_fac*sat[i] + euler_flux[i])
  # end
  fill!(euler_flux_bar, 0.0)
  fill!(sat_bar, 0.0)
  for i = 4:-1:1
    euler_flux_bar[i] += flux_bar[i]
    sat_bar[i] += sat_fac*flux_bar[i]
  end

  # calcEulerFlux(params, v_vals, aux_vars, nrm2, euler_flux)
  nrm2_bar = zeros(Tmsh, 2)
  calcEulerFlux_revm(params, v_vals, aux_vars, nrm2, nrm2_bar, euler_flux_bar)

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
function calcSAT_revm(params::ParamType{2}, nrm::AbstractArray{Tmsh,1},
          dq::AbstractArray{Tsol,1}, vel::AbstractArray{Tsol, 1},
          H::Tsol, sat_bar::AbstractArray{Tsol, 1},
          nrm_bar::AbstractArray{Tmsh,1}) where {Tmsh, Tsol}

  data = params.calcsatdata
  @unpack data E1dq E2dq E3dq E4dq

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
  @unpack data E3dq_bar E4dq_bar
  fill!(E3dq_bar, 0); fill!(E4dq_bar, 0)
#  E3dq_bar = zeros(Tsol, 4)
#  E4dq_bar = zeros(Tsol, 4)
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
  Reverse mode wrt. q of [`calcSAT`](@ref)

  **Inputs**

   * params
   * roe_vars
   * dq
   * nrm
   * sat_bar

  **Inputs/Outputs**

   * roe_vars_bar
   *  dq_bar
"""
function calcSAT_revq(params::ParamType{2},
                 roe_vars::AbstractArray{Tsol, 1},
                 roe_vars_bar::AbstractArray{Tres, 1},
                 dq::AbstractArray{Tsol,1},
                 dq_bar::AbstractArray{Tres, 1},
                 nrm::AbstractArray{Tmsh,1},
                 sat_bar::AbstractArray{Tsol,1}) where {Tmsh, Tsol, Tres}
# roe_vars = [u, v, H] at Roe average 

  data = params.calcsatdata
  @unpack data E1dq E2dq E3dq E4dq

  numDofPerNode = length(dq)

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

  lambda1_orig = lambda1  # save these for reverse sweep
  lambda2_orig = lambda2
  lambda3_orig = lambda3

  lambda1 = 0.5*(tau*max(absvalue(lambda1),sat_Vn *rhoA) - lambda1)
  lambda2 = 0.5*(tau*max(absvalue(lambda2),sat_Vn *rhoA) - lambda2)
  lambda3 = 0.5*(tau*max(absvalue(lambda3),sat_Vl *rhoA) - lambda3)


  dq1 = dq[1]
  dq2 = dq[2]
  dq3 = dq[3]
  dq4 = dq[4]

  
  #=
  sat[1] = lambda3*dq1
  sat[2] = lambda3*dq2
  sat[3] = lambda3*dq3
  sat[4] = lambda3*dq4
  =#

  #-- get E1*dq
  t1 = phi*dq1 - u*dq2 - v*dq3 + dq4
  E1dq[1] = t1
  E1dq[2] = t1*u
  E1dq[3] = t1*v
  E1dq[4] = t1*H


  #-- get E2*dq
  t2 = -Un*dq1 + nx*dq2 + ny*dq3
  E2dq[1] = 0.0
  E2dq[2] = t2*nx
  E2dq[3] = t2*ny
  E2dq[4] = t2*Un

  #-- add to sat
  tmp1 = 0.5*(lambda1 + lambda2) - lambda3
  tmp2 = gami/(a*a)
  tmp3 = 1.0/(dA*dA)

#  for i=1:length(sat)
#    sat[i] = sat[i] + tmp1*(tmp2*E1dq[i] + tmp3*E2dq[i])
#  end

  #-- get E3*dq
  t3 = -Un*dq1 + nx*dq2 + ny*dq3

  E3dq[1] = t3
  E3dq[2] = t3*u
  E3dq[3] = t3*v
  E3dq[4] = t3*H

  #-- get E4*dq
  t4 = phi*dq1 - u*dq2 - v*dq3 + dq4
  E4dq[1] = 0.0
  E4dq[2] = t4*nx
  E4dq[3] = t4*ny
  E4dq[4] = t4*Un

  #-- add to sat
  tmp4 = 0.5*(lambda1 - lambda2)/(dA*a)
#  for i=1:length(sat)
#    sat[i] = sat[i] + tmp4*(E3dq[i] + gami*E4dq[i])
#  end

  #----------------------------------------------------------------------------
  # reverse sweep
  @unpack data E1dq_bar E2dq_bar E3dq_bar E4dq_bar
  t1_bar = zero(Tres)
  t2_bar = zero(Tres)
  t3_bar = zero(Tres)
  t4_bar = zero(Tres)
  tmp1_bar = zero(Tres)
  tmp2_bar = zero(Tres)
  tmp3_bar = zero(Tres)
  tmp4_bar = zero(Tres)
  Un_bar = zero(Tres)
  phi_bar = zero(Tres)
  dq1_bar = zero(Tres)
  dq2_bar = zero(Tres)
  dq3_bar = zero(Tres)
  dq4_bar = zero(Tres)
  lambda1_bar = zero(Tres)
  lambda2_bar = zero(Tres)
  lambda3_bar = zero(Tres)
  rhoA_bar = zero(Tres)
  a_bar = zero(Tres)
  u_bar = zero(Tres)
  v_bar = zero(Tres)
  H_bar = zero(Tres)

  fill!(E1dq_bar, 0)
  fill!(E2dq_bar, 0)
  fill!(E3dq_bar, 0)
  fill!(E4dq_bar, 0)

  tmp_bar = zero(Tres)

  for i=1:length(sat_bar)
    tmp4_bar += (E3dq[i] + gami*E4dq[i])*sat_bar[i]
    E3dq_bar[i] += tmp4*sat_bar[i]
    E4dq_bar[i] += tmp4*gami*sat_bar[i]
  end

  lambda1_bar +=  0.5*tmp4_bar/(dA*a)
  lambda2_bar += -0.5*tmp4_bar/(dA*a)
  a_bar       += -tmp4*tmp4_bar/a


  # E4*dq
  
  t4_bar += E4dq_bar[4]*Un
  Un_bar += t4*E4dq_bar[4]

  t4_bar += E4dq_bar[3]*ny

  t4_bar += E4dq_bar[2]*nx

  phi_bar +=  dq1*t4_bar
  dq1_bar +=  phi*t4_bar
  u_bar   += -dq2*t4_bar
  dq2_bar += -u*t4_bar
  v_bar   += -dq3*t4_bar
  dq3_bar += -v*t4_bar
  dq4_bar +=  t4_bar


  # E3*dq
  t3_bar += E3dq_bar[4]*H
  H_bar  += t3*E3dq_bar[4]

  t3_bar += E3dq_bar[3]*v
  v_bar  += E3dq_bar[3]*t3

  t3_bar += E3dq_bar[2]*u
  u_bar  += E3dq_bar[2]*t3
  
  t3_bar += E3dq_bar[1]

  Un_bar  += -dq1*t3_bar
  dq1_bar += -Un*t3_bar
  dq2_bar +=  nx*t3_bar
  dq3_bar +=  ny*t3_bar



  for i=1:length(sat_bar)
    tmp1_bar    += sat_bar[i]*(tmp2*E1dq[i] + tmp3*E2dq[i])
    tmp2_bar    += tmp1*E1dq[i]*sat_bar[i]
    E1dq_bar[i] += tmp1*tmp2*sat_bar[i]
#    tmp3_bar    += tmp1*E2dq[i]*sat_bar[i]
    E2dq_bar[i] += tmp1*tmp3*sat_bar[i]
  end

  a_bar += (-2*tmp2/a)*tmp2_bar
  lambda1_bar += 0.5*tmp1_bar
  lambda2_bar += 0.5*tmp1_bar
  lambda3_bar += -tmp1_bar


  # E2*dq

  t2_bar += Un*E2dq_bar[4]
  Un_bar += t2*E2dq_bar[4]

  t2_bar += ny*E2dq_bar[3]

  t2_bar += nx*E2dq_bar[2]

  Un_bar  += -dq1*t2_bar
  dq1_bar += -Un*t2_bar
  dq2_bar +=  nx*t2_bar
  dq3_bar +=  ny*t2_bar



  # E1*dq
  t1_bar += E1dq_bar[4]*H
  H_bar  += E1dq_bar[4]*t1

  t1_bar += E1dq_bar[3]*v
  v_bar  += E1dq_bar[3]*t1

  t1_bar += E1dq_bar[2]*u
  u_bar  += E1dq_bar[2]*t1

  t1_bar += E1dq_bar[1]

  phi_bar +=  dq1*t1_bar
  dq1_bar +=  phi*t1_bar
  u_bar   += -dq2*t1_bar
  dq2_bar += -u*t1_bar
  v_bar   += -dq3*t1_bar
  dq3_bar += -v*t1_bar
  dq4_bar +=  t1_bar



  # sat diagonal term
  lambda3_bar += dq1*sat_bar[1]
  dq1_bar     += lambda3*sat_bar[1]

  lambda3_bar += dq2*sat_bar[2]
  dq2_bar     += lambda3*sat_bar[2]

  lambda3_bar += dq3*sat_bar[3]
  dq3_bar     += lambda3*sat_bar[3]

  lambda3_bar += dq4*sat_bar[4]
  dq4_bar     += lambda3*sat_bar[4]



  # re-pack dq
  dq_bar[1] += dq1_bar
  dq_bar[2] += dq2_bar
  dq_bar[3] += dq3_bar
  dq_bar[4] += dq4_bar


  # restore lambdas to original values from forward sweep
  # lambda_bar variables will be overwritten below
  lambda1 = lambda1_orig
  lambda2 = lambda2_orig
  lambda3 = lambda3_orig

  # entropy fix
  alambda1 = absvalue(lambda1)
  alambda2 = absvalue(lambda2)
  alambda3 = absvalue(lambda3)
  # compute the alambda_bar variables inside the conditionals because they
  # are not always needed

  if alambda1 > sat_Vn*rhoA
    # overwrite lambda here because in the primal code it gets re-assigned
    alambda1_bar = lambda1 > 0 ? lambda1_bar : -lambda1_bar
    lambda1_bar  = 0.5*(tau*alambda1_bar - lambda1_bar)
  else
    rhoA_bar    += 0.5*tau*sat_Vn*lambda1_bar
    lambda1_bar = -0.5*lambda1_bar
  end

  if alambda2 > sat_Vn*rhoA
    alambda2_bar = lambda2 > 0 ? lambda2_bar : -lambda2_bar
    lambda2_bar  = 0.5*(tau*alambda2_bar - lambda2_bar)
  else
    rhoA_bar    += 0.5*tau*sat_Vn*lambda2_bar
    lambda2_bar = -0.5*lambda2_bar
  end

  if alambda3 > sat_Vl*rhoA
    alambda3_bar = lambda3 > 0 ? lambda3_bar : -lambda3_bar
    lambda3_bar  = 0.5*tau*(alambda3_bar - lambda3_bar)
  else
    rhoA_bar    += 0.5*tau*sat_Vl*lambda3_bar
    lambda3_bar = -0.5*lambda3_bar
  end

  # rhoA
  if Un > 0
    Un_bar += rhoA_bar
  else
    Un_bar += -rhoA_bar
  end
  a_bar += dA*rhoA_bar


  # lambda calculation
  Un_bar += lambda1_bar
  a_bar  += dA*lambda1_bar

  Un_bar += lambda2_bar
  a_bar  += -dA*lambda2_bar

  Un_bar += lambda3_bar

  # a calculation
  H_bar   += (0.5/a)*gami*a_bar
  phi_bar += -(0.5/a)*gami*a_bar

  # phi
  u_bar += u*phi_bar
  v_bar += v*phi_bar

  # Un
  u_bar += nx*Un_bar
  v_bar += ny*Un_bar

  roe_vars_bar[1] += u_bar
  roe_vars_bar[2] += v_bar
  roe_vars_bar[3] += H_bar
  

  return nothing
end  # End function calcSAT




function RoeSolver_diff(params::ParamType{3, :conservative},
                   q::AbstractArray{Tsol,1},
                   qg::AbstractArray{Tsol, 1},
                   aux_vars::AbstractArray{Tres, 1},
                   nrm::AbstractArray{Tmsh,1},
                   fluxL_dot::AbstractArray{Tres, 2},
                   fluxR_dot::AbstractArray{Tres, 2}) where {Tmsh, Tsol, Tres}

  # SAT terms are used for ensuring consistency with the physical problem. Its
  # similar to upwinding which adds dissipation to the problem. SATs on the
  # boundary can be thought of as having two overlapping nodes and because of
  # the discontinuous nature of SBP adds some dissipation.

  data = params.roefluxdata
  @unpack data dq sat roe_vars euler_flux v_vals nrm2 roe_vars_dot


  # Declaring constants
  tau = 1.0
  gamma = params.gamma
  gami = params.gamma_1
  sat_fac = 1  # multiplier for SAT term

  # Begin main executuion
  nx = nrm[1]
  ny = nrm[2]
  nz = nrm[3]

  # Compute the Roe Averaged states
  # The left state of Roe are the actual solution variables
  # All the _dot variables are wrt q
  fac = 1.0/q[1]
  fac_dot1 = -fac*fac

  uL = q[2]*fac
  uL_dot1 = q[2]*fac_dot1
  uL_dot2 = fac

  vL = q[3]*fac
  vL_dot1 = q[3]*fac_dot1
  vL_dot3 = fac

  wL = q[4]*fac
  wL_dot1 = q[4]*fac_dot1
  wL_dot4 = fac

  phi = 0.5*(uL*uL + vL*vL + wL*wL)
  phi_dot1 = uL*uL_dot1 + vL*vL_dot1 + wL*wL_dot1
  phi_dot2 = uL*uL_dot2
  phi_dot3 = vL*vL_dot3
  phi_dot4 = wL*wL_dot4

  HL = gamma*q[5]*fac - gami*phi # Total enthalpy, H = e + 0.5*(u^2 + v^2) + p/rho,
                                 # where e is the internal energy per unit mass
  HL_dot1 = gamma*q[5]*fac_dot1 - gami*phi_dot1
  HL_dot2 = -gami*phi_dot2
  HL_dot3 = -gami*phi_dot3
  HL_dot4 = -gami*phi_dot4
  HL_dot5 = gamma*fac  # q[5]_dot = 1

  # The right side of the Roe solver comprises the boundary conditions
  # all the _dot variables are wrt qg now
  fac = 1.0/qg[1]
  fac_dot1 = -fac*fac
#  fac_dot1 = -1.0/(qg[1]*qg[1])

  uR = qg[2]*fac
  uR_dot1 = qg[2]*fac_dot1
  uR_dot2 = fac

  vR = qg[3]*fac
  vR_dot1 = qg[3]*fac_dot1
  vR_dot3 = fac

  wR = qg[4]*fac
  wR_dot1 = qg[4]*fac_dot1
  wR_dot4 = fac


  phi = 0.5*(uR*uR + vR*vR + wR*wR)
  phi_dot1 = uR*uR_dot1 + vR*vR_dot1 + wR*wR_dot1
  phi_dot2 = uR*uR_dot2
  phi_dot3 = vR*vR_dot3
  phi_dot4 = wR*wR_dot4


  HR = gamma*qg[5]*fac - gami*phi # Total Enthalpy
  HR_dot1 = gamma*qg[5]*fac_dot1 - gami*phi_dot1
  HR_dot2 = -gami*phi_dot2
  HR_dot3 = -gami*phi_dot3
  HR_dot4 = -gami*phi_dot4
  HR_dot5 = gamma*fac

  # Averaged states
  sqL = sqrt(q[1])
  sqL_dot1 = 0.5/sqL

  sqR = sqrt(qg[1])
  sqR_dot1 = 0.5/sqR

  fac = 1.0/(sqL + sqR)
  t1 = -1/((sqL + sqR)*(sqL + sqR))
  fac_dotL1 = t1*sqL_dot1
  fac_dotR1 = t1*sqR_dot1

  u = (sqL*uL + sqR*uR)*fac
  t2 = sqR*uR
  t3 = sqL*uL
  t4 = sqL*fac
  t5 = sqR*fac
  u_dotL1 = t3*fac_dotL1 + t4*uL_dot1 + uL*fac*sqL_dot1 + t2*fac_dotL1
  u_dotR1 = t2*fac_dotR1 + t5*uR_dot1 + uR*fac*sqR_dot1 + t3*fac_dotR1


  u_dotL2 = t4*uL_dot2
  u_dotR2 = t5*uR_dot2

  v = (sqL*vL + sqR*vR)*fac
  t2 = sqL*vL
  t3 = sqR*vR
  v_dotL1 = t2*fac_dotL1 + t4*vL_dot1 + vL*fac*sqL_dot1 + t3*fac_dotL1
  v_dotR1 = t3*fac_dotR1 + t5*vR_dot1 + vR*fac*sqR_dot1 + t2*fac_dotR1

  v_dotL3 = t4*vL_dot3
  v_dotR3 = t5*vR_dot3

  w = (sqL*wL + sqR*wR)*fac
  t2 = sqL*wL
  t3 = sqR*wR
  w_dotL1 = t2*fac_dotL1 + t4*wL_dot1 + wL*fac*sqL_dot1 + t3*fac_dotL1
  w_dotR1 = t3*fac_dotR1 + t5*wR_dot1 + wR*fac*sqR_dot1 + t2*fac_dotR1

  w_dotL4 = t4*wL_dot4
  w_dotR4 = t5*wR_dot4



  H = (sqL*HL + sqR*HR)*fac
  t2 = sqL*HL
  t3 = sqR*HR
  H_dotL1 = t2*fac_dotL1 + t4*HL_dot1 + HL*fac*sqL_dot1 + t3*fac_dotL1
  H_dotR1 = t3*fac_dotR1 + t5*HR_dot1 + HR*fac*sqR_dot1 + t2*fac_dotR1
 
  H_dotL2 = t4*HL_dot2 
  H_dotR2 = t5*HR_dot2

  H_dotL3 = t4*HL_dot3
  H_dotR3 = t5*HR_dot3

  H_dotL4 = t4*HL_dot4
  H_dotR4 = t5*HR_dot4

  H_dotL5 = t4*HL_dot5
  H_dotR5 = t5*HR_dot5


  for i=1:length(dq)
    dq[i] = q[i] - qg[i]
  end

  # dq_dotL* = 1, dq_dotR* = -1, so omit them

  roe_vars[1] = u
  roe_vars[2] = v
  roe_vars[3] = w
  roe_vars[4] = H

  roe_vars_dot[1]  = u_dotL1
  roe_vars_dot[2]  = u_dotR1
  roe_vars_dot[3]  = u_dotL2
  roe_vars_dot[4]  = u_dotR2

  roe_vars_dot[5]  = v_dotL1
  roe_vars_dot[6]  = v_dotR1
  roe_vars_dot[7]  = v_dotL3
  roe_vars_dot[8]  = v_dotR3

  roe_vars_dot[9]  =  w_dotL1
  roe_vars_dot[10]  = w_dotR1
  roe_vars_dot[11]  = w_dotL4
  roe_vars_dot[12]  = w_dotR4


  roe_vars_dot[13] = H_dotL1
  roe_vars_dot[14] = H_dotR1
  roe_vars_dot[15] = H_dotL2
  roe_vars_dot[16] = H_dotR2
  roe_vars_dot[17] = H_dotL3
  roe_vars_dot[18] = H_dotR3
  roe_vars_dot[19] = H_dotL4
  roe_vars_dot[20] = H_dotR4
  roe_vars_dot[21] = H_dotL5
  roe_vars_dot[22] = H_dotR5
  

  # pass in fluxL_dot and fluxR_dot here, then add the Euler flux
  # contribution below
  calcSAT_diff(params, roe_vars, roe_vars_dot,  dq, nrm, fluxL_dot, fluxR_dot)

  calcEulerFlux_diff(params, q, aux_vars, nrm, fluxL_dot)

  return nothing

end # ends the function RoeSolver




function calcSAT_diff(params::ParamType{2},
                 roe_vars::AbstractArray{Tsol, 1},
                 roe_vars_dot::AbstractArray{Tsol, 1},
                 dq::AbstractArray{Tsol,1},
                 nrm::AbstractArray{Tmsh,1},
                 sat_jacL::AbstractArray{Tsol,2},
                 sat_jacR::AbstractArray{Tsol,2}) where {Tmsh, Tsol}
# roe_vars = [u, v, H] at Roe average 
# roe_vars_dot contains all the non-zero derivatives of the roe_vars packed
# into a vector

  data = params.calcsatdata
  @unpack data E1dq E2dq

  # dq_dotL* = 1, dq_dotR* = -1, so don't pass them explicitly

  # SAT parameters
  sat_Vn = convert(Tsol, 0.025)
  sat_Vl = convert(Tsol, 0.025)
  tau = 1.0

  u = roe_vars[1]
  v = roe_vars[2]
  H = roe_vars[3]

  u_dotL1 = roe_vars_dot[1]
  u_dotR1 = roe_vars_dot[2]
  u_dotL2 = roe_vars_dot[3]
  u_dotR2 = roe_vars_dot[4]

  v_dotL1 = roe_vars_dot[5]
  v_dotR1 = roe_vars_dot[6]
  v_dotL3 = roe_vars_dot[7]
  v_dotR3 = roe_vars_dot[8]

  H_dotL1 = roe_vars_dot[9]
  H_dotR1 = roe_vars_dot[10]
  H_dotL2 = roe_vars_dot[11]
  H_dotR2 = roe_vars_dot[12]
  H_dotL3 = roe_vars_dot[13]
  H_dotR3 = roe_vars_dot[14]
  H_dotL4 = roe_vars_dot[15]
  H_dotR4 = roe_vars_dot[16]

  gami = params.gamma_1

  # Begin main execution
  nx = nrm[1]
  ny = nrm[2]

  dA = sqrt(nx*nx + ny*ny)

  Un = u*nx + v*ny # Normal Velocity
  Un_dotL1 = nx*u_dotL1 + ny*v_dotL1
  Un_dotR1 = nx*u_dotR1 + ny*v_dotR1

  Un_dotL2 = nx*u_dotL2
  Un_dotR2 = nx*u_dotR2

  Un_dotL3 = ny*v_dotL3
  Un_dotR3 = ny*v_dotR3

  phi = 0.5*(u*u + v*v)

  phi_dotL1 = u*u_dotL1 + v*v_dotL1
  phi_dotR1 = u*u_dotR1 + v*v_dotR1

  phi_dotL2 = u*u_dotL2
  phi_dotR2 = u*u_dotR2
  
  phi_dotL3 = v*v_dotL3
  phi_dotR3 = v*v_dotR3


  a = sqrt(gami*(H - phi)) # speed of sound
  t1 = gami/(2*a)
  a_dotL1 = t1*(H_dotL1 - phi_dotL1)
  a_dotR1 = t1*(H_dotR1 - phi_dotR1)

  a_dotL2 = t1*(H_dotL2 - phi_dotL2)
  a_dotR2 = t1*(H_dotR2 - phi_dotR2)

  a_dotL3 = t1*(H_dotL3 - phi_dotL3)
  a_dotR3 = t1*(H_dotR3 - phi_dotR3)

  a_dotL4 = t1*H_dotL4
  a_dotR4 = t1*H_dotR4

  lambda1 = Un + dA*a
  lambda1_dotL1 = Un_dotL1 + dA*a_dotL1
  lambda1_dotR1 = Un_dotR1 + dA*a_dotR1

  lambda1_dotL2 = Un_dotL2 + dA*a_dotL2
  lambda1_dotR2 = Un_dotR2 + dA*a_dotR2

  lambda1_dotL3 = Un_dotL3 + dA*a_dotL3
  lambda1_dotR3 = Un_dotR3 + dA*a_dotR3

  lambda1_dotL4 = dA*a_dotL4
  lambda1_dotR4 = dA*a_dotR4

  
  lambda2 = Un - dA*a
  lambda2_dotL1 = Un_dotL1 - dA*a_dotL1
  lambda2_dotR1 = Un_dotR1 - dA*a_dotR1

  lambda2_dotL2 = Un_dotL2 - dA*a_dotL2
  lambda2_dotR2 = Un_dotR2 - dA*a_dotR2

  lambda2_dotL3 = Un_dotL3 - dA*a_dotL3
  lambda2_dotR3 = Un_dotR3 - dA*a_dotR3

  lambda2_dotL4 = -dA*a_dotL4
  lambda2_dotR4 = -dA*a_dotR4


  lambda3 = Un
  lambda3_dotL1 = Un_dotL1
  lambda3_dotR1 = Un_dotR1

  lambda3_dotL2 = Un_dotL2
  lambda3_dotR2 = Un_dotR2

  lambda3_dotL3 = Un_dotL3
  lambda3_dotR3 = Un_dotR3

  rhoA = absvalue(Un) + dA*a
  #TODO: see if there is a better way to do this
  if Un > 0
    fac = 1
  else
    fac = -1
  end

  rhoA_dotL1 = fac*Un_dotL1 + dA*a_dotL1
  rhoA_dotR1 = fac*Un_dotR1 + dA*a_dotR1

  rhoA_dotL2 = fac*Un_dotL2 + dA*a_dotL2
  rhoA_dotR2 = fac*Un_dotR2 + dA*a_dotR2

  rhoA_dotL3 = fac*Un_dotL3 + dA*a_dotL3
  rhoA_dotR3 = fac*Un_dotR3 + dA*a_dotR3

  rhoA_dotL4 = dA*a_dotL4
  rhoA_dotR4 = dA*a_dotR4

  # Compute Eigen Values of the Flux Jacobian
  # The eigen values calculated above cannot be used directly. Near stagnation
  # points lambda3 approaches zero while near sonic lines lambda1 and lambda2
  # approach zero. This has a possibility of creating numerical difficulties.
  # As a result, the eigen values are limited by the following expressions.


  # see lambda1 expression below
  if absvalue(lambda1) > sat_Vn*rhoA
    if lambda1 > 0
      fac = 1
    else
      fac = -1
    end
    
    t1 = tau*fac
    #TODO: lambda1_dotL1 - lambgda1_dotL1 = 0, so simplify this
    lambda1_dotL1 = 0.5 * (t1 * lambda1_dotL1 - lambda1_dotL1) 
    lambda1_dotR1 = 0.5 * (t1 * lambda1_dotR1 - lambda1_dotR1) 

    lambda1_dotL2 = 0.5 * (t1 * lambda1_dotL2 - lambda1_dotL2) 
    lambda1_dotR2 = 0.5 * (t1 * lambda1_dotR2 - lambda1_dotR2) 

    lambda1_dotL3 = 0.5 * (t1 * lambda1_dotL3 - lambda1_dotL3) 
    lambda1_dotR3 = 0.5 * (t1 * lambda1_dotR3 - lambda1_dotR3) 

    lambda1_dotL4 = 0.5 * (t1 * lambda1_dotL4 - lambda1_dotL4) 
    lambda1_dotR4 = 0.5 * (t1 * lambda1_dotR4 - lambda1_dotR4) 


  else
    t1 = sat_Vn*tau
    lambda1_dotL1 =  0.5 * (t1 * rhoA_dotL1 - lambda1_dotL1) 
    lambda1_dotR1 =  0.5 * (t1 * rhoA_dotR1 - lambda1_dotR1) 
 
    lambda1_dotL2 =  0.5 * (t1 * rhoA_dotL2 - lambda1_dotL2) 
    lambda1_dotR2 =  0.5 * (t1 * rhoA_dotR2 - lambda1_dotR2) 
    
    lambda1_dotL3 =  0.5 * (t1 * rhoA_dotL3 - lambda1_dotL3) 
    lambda1_dotR3 =  0.5 * (t1 * rhoA_dotR3 - lambda1_dotR3) 
 
    lambda1_dotL4 =  0.5 * (t1 * rhoA_dotL4 - lambda1_dotL4) 
    lambda1_dotR4 =  0.5 * (t1 * rhoA_dotR4 - lambda1_dotR4) 
 
  end

  lambda1 = 0.5*(tau*max(absvalue(lambda1),sat_Vn *rhoA) - lambda1)

  # see lambda2 expression below
  if absvalue(lambda2) > sat_Vn*rhoA
    if lambda2 > 0
      fac = 1
    else
      fac = -1
    end

    t1 = tau*fac
    lambda2_dotL1 = 0.5 * (t1 * lambda2_dotL1 - lambda2_dotL1)
    lambda2_dotR1 = 0.5 * (t1 * lambda2_dotR1 - lambda2_dotR1)

    lambda2_dotL2 = 0.5 * (t1 * lambda2_dotL2 - lambda2_dotL2)
    lambda2_dotR2 = 0.5 * (t1 * lambda2_dotR2 - lambda2_dotR2)

    lambda2_dotL3 = 0.5 * (t1 * lambda2_dotL3 - lambda2_dotL3)
    lambda2_dotR3 = 0.5 * (t1 * lambda2_dotR3 - lambda2_dotR3)

    lambda2_dotL4 = 0.5 * (t1 * lambda2_dotL4 - lambda2_dotL4)
    lambda2_dotR4 = 0.5 * (t1 * lambda2_dotR4 - lambda2_dotR4)

  else

    t1 = sat_Vn*tau
    lambda2_dotL1 = 0.5 * (t1 * rhoA_dotL1 - lambda2_dotL1)
    lambda2_dotR1 = 0.5 * (t1 * rhoA_dotR1 - lambda2_dotR1)
 
    lambda2_dotL2 = 0.5 * (t1 * rhoA_dotL2 - lambda2_dotL2)
    lambda2_dotR2 = 0.5 * (t1 * rhoA_dotR2 - lambda2_dotR2)
    
    lambda2_dotL3 = 0.5 * (t1 * rhoA_dotL3 - lambda2_dotL3)
    lambda2_dotR3 = 0.5 * (t1 * rhoA_dotR3 - lambda2_dotR3)
 
    lambda2_dotL4 = 0.5 * (t1 * rhoA_dotL4 - lambda2_dotL4)
    lambda2_dotR4 = 0.5 * (t1 * rhoA_dotR4 - lambda2_dotR4)
 
  end
  lambda2 = 0.5*(tau*max(absvalue(lambda2),sat_Vn *rhoA) - lambda2)


  # see lambda3 expression below
  if absvalue(lambda3) > sat_Vn*rhoA
    if lambda3 > 0
      fac = 1
    else
      fac = -1
    end

    t1 = tau*fac
    lambda3_dotL1 = 0.5 * (t1 * lambda3_dotL1 - lambda3_dotL1)
    lambda3_dotR1 = 0.5 * (t1 * lambda3_dotR1 - lambda3_dotR1)

    lambda3_dotL2 = 0.5 * (t1 * lambda3_dotL2 - lambda3_dotL2)
    lambda3_dotR2 = 0.5 * (t1 * lambda3_dotR2 - lambda3_dotR2)

    lambda3_dotL3 = 0.5 * (t1 * lambda3_dotL3 - lambda3_dotL3)
    lambda3_dotR3 = 0.5 * (t1 * lambda3_dotR3 - lambda3_dotR3)

    lambda3_dotL4 = 0.0
    lambda3_dotR4 = 0.0


  else
    t1 = sat_Vn*tau
    lambda3_dotL1 = 0.5 * (t1 * rhoA_dotL1 - lambda3_dotL1)
    lambda3_dotR1 = 0.5 * (t1 * rhoA_dotR1 - lambda3_dotR1)
 
    lambda3_dotL2 = 0.5 * (t1 * rhoA_dotL2 - lambda3_dotL2)
    lambda3_dotR2 = 0.5 * (t1 * rhoA_dotR2 - lambda3_dotR2)
    
    lambda3_dotL3 = 0.5 * (t1 * rhoA_dotL3 - lambda3_dotL3)
    lambda3_dotR3 = 0.5 * (t1 * rhoA_dotR3 - lambda3_dotR3)
 

    lambda3_dotL4 = 0.5 * t1 * rhoA_dotL4
    lambda3_dotR4 = 0.5 * t1 * rhoA_dotR4
 
  end

  lambda3 = 0.5*(tau*max(absvalue(lambda3),sat_Vl *rhoA) - lambda3)

                    
  dq1 = dq[1]
  dq2 = dq[2]
  dq3 = dq[3]
  dq4 = dq[4]

  # sat[1] = lambda3*dq1
  sat_jacL[1, 1] += lambda3 + dq1*lambda3_dotL1
  sat_jacL[1, 2] +=         + dq1*lambda3_dotL2
  sat_jacL[1, 3] +=         + dq1*lambda3_dotL3
  sat_jacL[1, 4] +=         + dq1*lambda3_dotL4

  sat_jacR[1, 1] += -lambda3 + dq1*lambda3_dotR1
  sat_jacR[1, 2] +=          + dq1*lambda3_dotR2
  sat_jacR[1, 3] +=          + dq1*lambda3_dotR3
  sat_jacR[1, 4] +=          + dq1*lambda3_dotR4

  # sat[2] = lambda3*dq2
  sat_jacL[2, 1] +=         + dq2*lambda3_dotL1
  sat_jacL[2, 2] += lambda3 + dq2*lambda3_dotL2
  sat_jacL[2, 3] +=         + dq2*lambda3_dotL3
  sat_jacL[2, 4] +=         + dq2*lambda3_dotL4

  sat_jacR[2, 1] +=          + dq2*lambda3_dotR1
  sat_jacR[2, 2] += -lambda3 + dq2*lambda3_dotR2
  sat_jacR[2, 3] +=          + dq2*lambda3_dotR3
  sat_jacR[2, 4] +=          + dq2*lambda3_dotR4

  # sat[3] = lambda3*dq3
  sat_jacL[3, 1] +=         + dq3*lambda3_dotL1
  sat_jacL[3, 2] +=         + dq3*lambda3_dotL2
  sat_jacL[3, 3] += lambda3 + dq3*lambda3_dotL3
  sat_jacL[3, 4] +=         + dq3*lambda3_dotL4

  sat_jacR[3, 1] +=          + dq3*lambda3_dotR1
  sat_jacR[3, 2] +=          + dq3*lambda3_dotR2
  sat_jacR[3, 3] += -lambda3 + dq3*lambda3_dotR3
  sat_jacR[3, 4] +=          + dq3*lambda3_dotR4

  # sat[4] = lambda3*dq4
  sat_jacL[4, 1] +=           dq4*lambda3_dotL1
  sat_jacL[4, 2] +=           dq4*lambda3_dotL2
  sat_jacL[4, 3] +=           dq4*lambda3_dotL3
  sat_jacL[4, 4] += lambda3 + dq4*lambda3_dotL4

  sat_jacR[4, 1] +=            dq4*lambda3_dotR1
  sat_jacR[4, 2] +=            dq4*lambda3_dotR2
  sat_jacR[4, 3] +=            dq4*lambda3_dotR3
  sat_jacR[4, 4] += -lambda3 + dq4*lambda3_dotR4

  #-- get E1*dq
  E1dq[1] = phi*dq1 - u*dq2 - v*dq3 + dq4
  E1dq1_dotL1 = phi + dq1*phi_dotL1     - dq2*u_dotL1     - dq3*v_dotL1
  E1dq1_dotL2 =       dq1*phi_dotL2 - u - dq2*u_dotL2
  E1dq1_dotL3 =       dq1*phi_dotL3                   - v - dq3*v_dotL3
  E1dq1_dotL4 = 1

  E1dq1_dotR1 = -phi + dq1*phi_dotR1     - dq2*u_dotR1     - dq3*v_dotR1
  E1dq1_dotR2 =        dq1*phi_dotR2 + u - dq2*u_dotR2
  E1dq1_dotR3 =        dq1*phi_dotR3                   + v - dq3*v_dotR3
  E1dq1_dotR4 = -1


  E1dq[2] = E1dq[1]*u
  E1dq2_dotL1 = u*E1dq1_dotL1 + E1dq[1]*u_dotL1
  E1dq2_dotL2 = u*E1dq1_dotL2 + E1dq[1]*u_dotL2
  E1dq2_dotL3 = u*E1dq1_dotL3
  E1dq2_dotL4 = u*E1dq1_dotL4

  E1dq2_dotR1 = u*E1dq1_dotR1 + E1dq[1]*u_dotR1
  E1dq2_dotR2 = u*E1dq1_dotR2 + E1dq[1]*u_dotR2
  E1dq2_dotR3 = u*E1dq1_dotR3
  E1dq2_dotR4 = u*E1dq1_dotR4

  E1dq[3] = E1dq[1]*v
  E1dq3_dotL1 = v*E1dq1_dotL1 + E1dq[1]*v_dotL1
  E1dq3_dotL2 = v*E1dq1_dotL2
  E1dq3_dotL3 = v*E1dq1_dotL3 + E1dq[1]*v_dotL3
  E1dq3_dotL4 = v*E1dq1_dotL4

  E1dq3_dotR1 = v*E1dq1_dotR1 + E1dq[1]*v_dotR1
  E1dq3_dotR2 = v*E1dq1_dotR2
  E1dq3_dotR3 = v*E1dq1_dotR3 + E1dq[1]*v_dotR3
  E1dq3_dotR4 = v*E1dq1_dotR4

  E1dq[4] = E1dq[1]*H
  E1dq4_dotL1 = H*E1dq1_dotL1 + E1dq[1]*H_dotL1
  E1dq4_dotL2 = H*E1dq1_dotL2 + E1dq[1]*H_dotL2
  E1dq4_dotL3 = H*E1dq1_dotL3 + E1dq[1]*H_dotL3
  E1dq4_dotL4 = H*E1dq1_dotL4 + E1dq[1]*H_dotL4

  E1dq4_dotR1 = H*E1dq1_dotR1 + E1dq[1]*H_dotR1
  E1dq4_dotR2 = H*E1dq1_dotR2 + E1dq[1]*H_dotR2
  E1dq4_dotR3 = H*E1dq1_dotR3 + E1dq[1]*H_dotR3
  E1dq4_dotR4 = H*E1dq1_dotR4 + E1dq[1]*H_dotR4



  #-- get E2*dq
  E2dq[1] = 0.0
  E2dq[2] = -Un*dq1 + nx*dq2 + ny*dq3
  E2dq2_dotL1 = -Un + -Un_dotL1*dq1
  E2dq2_dotL2 =     + -Un_dotL2*dq1 + nx
  E2dq2_dotL3 =     + -Un_dotL3*dq1      + ny

  E2dq2_dotR1 = Un + -Un_dotR1*dq1
  E2dq2_dotR2 =    + -Un_dotR2*dq1 - nx
  E2dq2_dotR3 =    + -Un_dotR3*dq1      - ny


  E2dq[3] = E2dq[2]*ny
  E2dq3_dotL1 = ny*E2dq2_dotL1
  E2dq3_dotL2 = ny*E2dq2_dotL2
  E2dq3_dotL3 = ny*E2dq2_dotL3

  E2dq3_dotR1 = ny*E2dq2_dotR1
  E2dq3_dotR2 = ny*E2dq2_dotR2
  E2dq3_dotR3 = ny*E2dq2_dotR3

  E2dq[4] = E2dq[2]*Un
  E2dq4_dotL1 = Un*E2dq2_dotL1 + E2dq[2]*Un_dotL1
  E2dq4_dotL2 = Un*E2dq2_dotL2 + E2dq[2]*Un_dotL2
  E2dq4_dotL3 = Un*E2dq2_dotL3 + E2dq[2]*Un_dotL3

  E2dq4_dotR1 = Un*E2dq2_dotR1 + E2dq[2]*Un_dotR1
  E2dq4_dotR2 = Un*E2dq2_dotR2 + E2dq[2]*Un_dotR2
  E2dq4_dotR3 = Un*E2dq2_dotR3 + E2dq[2]*Un_dotR3

  E2dq[2] = E2dq[2]*nx
  E2dq2_dotL1 = nx*E2dq2_dotL1
  E2dq2_dotL2 = nx*E2dq2_dotL2
  E2dq2_dotL3 = nx*E2dq2_dotL3

  E2dq2_dotR1 = nx*E2dq2_dotR1
  E2dq2_dotR2 = nx*E2dq2_dotR2
  E2dq2_dotR3 = nx*E2dq2_dotR3

  #-- add to sat
  tmp1 = 0.5*(lambda1 + lambda2) - lambda3
  tmp1_dotL1 = 0.5*(lambda1_dotL1 + lambda2_dotL1) - lambda3_dotL1
  tmp1_dotL2 = 0.5*(lambda1_dotL2 + lambda2_dotL2) - lambda3_dotL2
  tmp1_dotL3 = 0.5*(lambda1_dotL3 + lambda2_dotL3) - lambda3_dotL3
  tmp1_dotL4 = 0.5*(lambda1_dotL4 + lambda2_dotL4) - lambda3_dotL4

  tmp1_dotR1 = 0.5*(lambda1_dotR1 + lambda2_dotR1) - lambda3_dotR1
  tmp1_dotR2 = 0.5*(lambda1_dotR2 + lambda2_dotR2) - lambda3_dotR2
  tmp1_dotR3 = 0.5*(lambda1_dotR3 + lambda2_dotR3) - lambda3_dotR3
  tmp1_dotR4 = 0.5*(lambda1_dotR4 + lambda2_dotR4) - lambda3_dotR4


  tmp2 = gami/(a*a)
  t1 = -2*tmp2/a
#  t1 = -2*gami/(a*a*a) # = -2*tmp2/a
  tmp2_dotL1 = t1*a_dotL1
  tmp2_dotL2 = t1*a_dotL2
  tmp2_dotL3 = t1*a_dotL3
  tmp2_dotL4 = t1*a_dotL4

  tmp2_dotR1 = t1*a_dotR1
  tmp2_dotR2 = t1*a_dotR2
  tmp2_dotR3 = t1*a_dotR3
  tmp2_dotR4 = t1*a_dotR4

  tmp3 = 1.0/(dA*dA)
  #for i=1:length(sat)
  #  sat[i] = sat[i] + tmp1*(tmp2*E1dq[i] + tmp3*E2dq[i])
  #end

  #TODO: align + signs
  sat_jacL[1,1] += tmp1*tmp2*E1dq1_dotL1   + tmp1*E1dq[1]*tmp2_dotL1 + 
                   tmp2*E1dq[1]*tmp1_dotL1 + 
                                           + tmp3*E2dq[1]*tmp1_dotL1

  sat_jacL[1,2] += tmp1*tmp2*E1dq1_dotL2   + tmp1*E1dq[1]*tmp2_dotL2 + 
                   tmp2*E1dq[1]*tmp1_dotL2 + 
                                           + tmp3*E2dq[1]*tmp1_dotL2

  sat_jacL[1,3] += tmp1*tmp2*E1dq1_dotL3   + tmp1*E1dq[1]*tmp2_dotL3 + 
                   tmp2*E1dq[1]*tmp1_dotL3 + 
                                           + tmp3*E2dq[1]*tmp1_dotL3

  sat_jacL[1,4] += tmp1*tmp2*E1dq1_dotL4   + tmp1*E1dq[1]*tmp2_dotL4 + 
                   tmp2*E1dq[1]*tmp1_dotL4 + 
                                           + tmp3*E2dq[1]*tmp1_dotL4

  sat_jacR[1,1] += tmp1*tmp2*E1dq1_dotR1   + tmp1*E1dq[1]*tmp2_dotR1 + 
                   tmp2*E1dq[1]*tmp1_dotR1 + 
                                           + tmp3*E2dq[1]*tmp1_dotR1

  sat_jacR[1,2] += tmp1*tmp2*E1dq1_dotR2   + tmp1*E1dq[1]*tmp2_dotR2 + 
                   tmp2*E1dq[1]*tmp1_dotR2 + 
                                           + tmp3*E2dq[1]*tmp1_dotR2

  sat_jacR[1,3] += tmp1*tmp2*E1dq1_dotR3   + tmp1*E1dq[1]*tmp2_dotR3 + 
                   tmp2*E1dq[1]*tmp1_dotR3 + 
                                           + tmp3*E2dq[1]*tmp1_dotR3

  sat_jacR[1,4] += tmp1*tmp2*E1dq1_dotR4   + tmp1*E1dq[1]*tmp2_dotR4 + 
                   tmp2*E1dq[1]*tmp1_dotR4 + 
                                           + tmp3*E2dq[1]*tmp1_dotR4


  sat_jacL[2,1] += tmp1*tmp2*E1dq2_dotL1   + tmp1*E1dq[2]*tmp2_dotL1 + 
                   tmp2*E1dq[2]*tmp1_dotL1 + 
                   tmp1*tmp3*E2dq2_dotL1   + tmp3*E2dq[2]*tmp1_dotL1

  sat_jacL[2,2] += tmp1*tmp2*E1dq2_dotL2   + tmp1*E1dq[2]*tmp2_dotL2 + 
                   tmp2*E1dq[2]*tmp1_dotL2 + 
                   tmp1*tmp3*E2dq2_dotL2   + tmp3*E2dq[2]*tmp1_dotL2

  sat_jacL[2,3] += tmp1*tmp2*E1dq2_dotL3   + tmp1*E1dq[2]*tmp2_dotL3 + 
                   tmp2*E1dq[2]*tmp1_dotL3 + 
                   tmp1*tmp3*E2dq2_dotL3   + tmp3*E2dq[2]*tmp1_dotL3

  sat_jacL[2,4] += tmp1*tmp2*E1dq2_dotL4   + tmp1*E1dq[2]*tmp2_dotL4 + 
                   tmp2*E1dq[2]*tmp1_dotL4 + 
                                           + tmp3*E2dq[2]*tmp1_dotL4

  sat_jacR[2,1] += tmp1*tmp2*E1dq2_dotR1   + tmp1*E1dq[2]*tmp2_dotR1 + 
                   tmp2*E1dq[2]*tmp1_dotR1 + 
                   tmp1*tmp3*E2dq2_dotR1   + tmp3*E2dq[2]*tmp1_dotR1

  sat_jacR[2,2] += tmp1*tmp2*E1dq2_dotR2   + tmp1*E1dq[2]*tmp2_dotR2 + 
                   tmp2*E1dq[2]*tmp1_dotR2 + 
                   tmp1*tmp3*E2dq2_dotR2   + tmp3*E2dq[2]*tmp1_dotR2

  sat_jacR[2,3] += tmp1*tmp2*E1dq2_dotR3   + tmp1*E1dq[2]*tmp2_dotR3 + 
                   tmp2*E1dq[2]*tmp1_dotR3 + 
                   tmp1*tmp3*E2dq2_dotR3   + tmp3*E2dq[2]*tmp1_dotR3

  sat_jacR[2,4] += tmp1*tmp2*E1dq2_dotR4   + tmp1*E1dq[2]*tmp2_dotR4 + 
                   tmp2*E1dq[2]*tmp1_dotR4 + 
                                           + tmp3*E2dq[2]*tmp1_dotR4


  sat_jacL[3,1] += tmp1*tmp2*E1dq3_dotL1   + tmp1*E1dq[3]*tmp2_dotL1 + 
                   tmp2*E1dq[3]*tmp1_dotL1 + 
                   tmp1*tmp3*E2dq3_dotL1   + tmp3*E2dq[3]*tmp1_dotL1

  sat_jacL[3,2] += tmp1*tmp2*E1dq3_dotL2   + tmp1*E1dq[3]*tmp2_dotL2 + 
                   tmp2*E1dq[3]*tmp1_dotL2 + 
                   tmp1*tmp3*E2dq3_dotL2   + tmp3*E2dq[3]*tmp1_dotL2

  sat_jacL[3,3] += tmp1*tmp2*E1dq3_dotL3   + tmp1*E1dq[3]*tmp2_dotL3 + 
                   tmp2*E1dq[3]*tmp1_dotL3 + 
                   tmp1*tmp3*E2dq3_dotL3   + tmp3*E2dq[3]*tmp1_dotL3

  sat_jacL[3,4] += tmp1*tmp2*E1dq3_dotL4   + tmp1*E1dq[3]*tmp2_dotL4 + 
                   tmp2*E1dq[3]*tmp1_dotL4 + 
                                           + tmp3*E2dq[3]*tmp1_dotL4

  sat_jacR[3,1] += tmp1*tmp2*E1dq3_dotR1   + tmp1*E1dq[3]*tmp2_dotR1 + 
                   tmp2*E1dq[3]*tmp1_dotR1 + 
                   tmp1*tmp3*E2dq3_dotR1   + tmp3*E2dq[3]*tmp1_dotR1

  sat_jacR[3,2] += tmp1*tmp2*E1dq3_dotR2   + tmp1*E1dq[3]*tmp2_dotR2 + 
                   tmp2*E1dq[3]*tmp1_dotR2 + 
                   tmp1*tmp3*E2dq3_dotR2   + tmp3*E2dq[3]*tmp1_dotR2

  sat_jacR[3,3] += tmp1*tmp2*E1dq3_dotR3   + tmp1*E1dq[3]*tmp2_dotR3 + 
                   tmp2*E1dq[3]*tmp1_dotR3 + 
                   tmp1*tmp3*E2dq3_dotR3   + tmp3*E2dq[3]*tmp1_dotR3

  sat_jacR[3,4] += tmp1*tmp2*E1dq3_dotR4   + tmp1*E1dq[3]*tmp2_dotR4 + 
                   tmp2*E1dq[3]*tmp1_dotR4 + 
                                           + tmp3*E2dq[3]*tmp1_dotR4


  sat_jacL[4,1] += tmp1*tmp2*E1dq4_dotL1   + tmp1*E1dq[4]*tmp2_dotL1 + 
                   tmp2*E1dq[4]*tmp1_dotL1 + 
                   tmp1*tmp3*E2dq4_dotL1   + tmp3*E2dq[4]*tmp1_dotL1

  sat_jacL[4,2] += tmp1*tmp2*E1dq4_dotL2   + tmp1*E1dq[4]*tmp2_dotL2 + 
                   tmp2*E1dq[4]*tmp1_dotL2 + 
                   tmp1*tmp3*E2dq4_dotL2   + tmp3*E2dq[4]*tmp1_dotL2

  sat_jacL[4,3] += tmp1*tmp2*E1dq4_dotL3   + tmp1*E1dq[4]*tmp2_dotL3 + 
                   tmp2*E1dq[4]*tmp1_dotL3 + 
                   tmp1*tmp3*E2dq4_dotL3   + tmp3*E2dq[4]*tmp1_dotL3

  sat_jacL[4,4] += tmp1*tmp2*E1dq4_dotL4   + tmp1*E1dq[4]*tmp2_dotL4 + 
                   tmp2*E1dq[4]*tmp1_dotL4 + 
                                           + tmp3*E2dq[4]*tmp1_dotL4

  sat_jacR[4,1] += tmp1*tmp2*E1dq4_dotR1   + tmp1*E1dq[4]*tmp2_dotR1 + 
                   tmp2*E1dq[4]*tmp1_dotR1 + 
                   tmp1*tmp3*E2dq4_dotR1   + tmp3*E2dq[4]*tmp1_dotR1

  sat_jacR[4,2] += tmp1*tmp2*E1dq4_dotR2   + tmp1*E1dq[4]*tmp2_dotR2 + 
                   tmp2*E1dq[4]*tmp1_dotR2 + 
                   tmp1*tmp3*E2dq4_dotR2   + tmp3*E2dq[4]*tmp1_dotR2

  sat_jacR[4,3] += tmp1*tmp2*E1dq4_dotR3   + tmp1*E1dq[4]*tmp2_dotR3 + 
                   tmp2*E1dq[4]*tmp1_dotR3 + 
                   tmp1*tmp3*E2dq4_dotR3   + tmp3*E2dq[4]*tmp1_dotR3

  sat_jacR[4,4] += tmp1*tmp2*E1dq4_dotR4   + tmp1*E1dq[4]*tmp2_dotR4 + 
                   tmp2*E1dq[4]*tmp1_dotR4 + 
                                           + tmp3*E2dq[4]*tmp1_dotR4

  #-- get E3*dq
  E1dq[1] = -Un*dq1 + nx*dq2 + ny*dq3
  E1dq1_dotL1 = -Un - dq1*Un_dotL1
  E1dq1_dotL2 =     - dq1*Un_dotL2 + nx
  E1dq1_dotL3 =     - dq1*Un_dotL3      + ny

  E1dq1_dotR1 = Un - dq1*Un_dotR1
  E1dq1_dotR2 =    - dq1*Un_dotR2 - nx
  E1dq1_dotR3 =    - dq1*Un_dotR3      - ny


  E1dq[2] = E1dq[1]*u
  E1dq2_dotL1 = u*E1dq1_dotL1 + E1dq[1]*u_dotL1
  E1dq2_dotL2 = u*E1dq1_dotL2 + E1dq[1]*u_dotL2
  E1dq2_dotL3 = u*E1dq1_dotL3

  E1dq2_dotR1 = u*E1dq1_dotR1 + E1dq[1]*u_dotR1
  E1dq2_dotR2 = u*E1dq1_dotR2 + E1dq[1]*u_dotR2
  E1dq2_dotR3 = u*E1dq1_dotR3


  E1dq[3] = E1dq[1]*v
  E1dq3_dotL1 = v*E1dq1_dotL1 + E1dq[1]*v_dotL1
  E1dq3_dotL2 = v*E1dq1_dotL2
  E1dq3_dotL3 = v*E1dq1_dotL3 + E1dq[1]*v_dotL3

  E1dq3_dotR1 = v*E1dq1_dotR1 + E1dq[1]*v_dotR1
  E1dq3_dotR2 = v*E1dq1_dotR2
  E1dq3_dotR3 = v*E1dq1_dotR3 + E1dq[1]*v_dotR3


  E1dq[4] = E1dq[1]*H
  E1dq4_dotL1 = H*E1dq1_dotL1 + E1dq[1]*H_dotL1
  E1dq4_dotL2 = H*E1dq1_dotL2 + E1dq[1]*H_dotL2
  E1dq4_dotL3 = H*E1dq1_dotL3 + E1dq[1]*H_dotL3
  E1dq4_dotL4 =               + E1dq[1]*H_dotL4

  E1dq4_dotR1 = H*E1dq1_dotR1 + E1dq[1]*H_dotR1
  E1dq4_dotR2 = H*E1dq1_dotR2 + E1dq[1]*H_dotR2
  E1dq4_dotR3 = H*E1dq1_dotR3 + E1dq[1]*H_dotR3
  E1dq4_dotR4 =               + E1dq[1]*H_dotR4


  #-- get E4*dq
  t1 = phi*dq1 - u*dq2 - v*dq3 + dq4
  t1_dotL1 = phi + dq1*phi_dotL1     - dq2*u_dotL1     - dq3*v_dotL1
  t1_dotL2 =       dq1*phi_dotL2 - u - dq2*u_dotL2
  t1_dotL3 =       dq1*phi_dotL3                   - v - dq3*v_dotL3
  t1_dotL4 = 1

  t1_dotR1 = -phi + dq1*phi_dotR1     - dq2*u_dotR1     - dq3*v_dotR1
  t1_dotR2 =        dq1*phi_dotR2 + u - dq2*u_dotR2
  t1_dotR3 =        dq1*phi_dotR3                   + v - dq3*v_dotR3
  t1_dotR4 = -1


  E2dq[1] = 0.0
  E2dq[2] = t1*nx
  E2dq2_dotL1 = nx*t1_dotL1
  E2dq2_dotL2 = nx*t1_dotL2
  E2dq2_dotL3 = nx*t1_dotL3
  E2dq2_dotL4 = nx*t1_dotL4

  E2dq2_dotR1 = nx*t1_dotR1
  E2dq2_dotR2 = nx*t1_dotR2
  E2dq2_dotR3 = nx*t1_dotR3
  E2dq2_dotR4 = nx*t1_dotR4


  E2dq[3] = t1*ny
  E2dq3_dotL1 = ny*t1_dotL1
  E2dq3_dotL2 = ny*t1_dotL2
  E2dq3_dotL3 = ny*t1_dotL3
  E2dq3_dotL4 = ny*t1_dotL4

  E2dq3_dotR1 = ny*t1_dotR1
  E2dq3_dotR2 = ny*t1_dotR2
  E2dq3_dotR3 = ny*t1_dotR3
  E2dq3_dotR4 = ny*t1_dotR4


  E2dq[4] = t1*Un
  E2dq4_dotL1 = Un*t1_dotL1 + t1*Un_dotL1
  E2dq4_dotL2 = Un*t1_dotL2 + t1*Un_dotL2
  E2dq4_dotL3 = Un*t1_dotL3 + t1*Un_dotL3
  E2dq4_dotL4 = Un*t1_dotL4

  E2dq4_dotR1 = Un*t1_dotR1 + t1*Un_dotR1
  E2dq4_dotR2 = Un*t1_dotR2 + t1*Un_dotR2
  E2dq4_dotR3 = Un*t1_dotR3 + t1*Un_dotR3
  E2dq4_dotR4 = Un*t1_dotR4

  #-- add to sat
  t1 = 1/(dA*a)
  t2 = 1/a
  tmp1 = 0.5*(lambda1 - lambda2)*t1

  tmp1_dotL1 = 0.5 * t1 * ( (lambda1_dotL1 - lambda2_dotL1) - (lambda1 - lambda2) * t2 * a_dotL1)
  tmp1_dotL2 = 0.5 * t1 * ( (lambda1_dotL2 - lambda2_dotL2) - (lambda1 - lambda2) * t2 * a_dotL2)
  tmp1_dotL3 = 0.5 * t1 * ( (lambda1_dotL3 - lambda2_dotL3) - (lambda1 - lambda2) * t2 * a_dotL3)
  tmp1_dotL4 = 0.5 * t1 * ( (lambda1_dotL4 - lambda2_dotL4) - (lambda1 - lambda2) * t2 * a_dotL4)

  tmp1_dotR1 = 0.5 * t1 * ( (lambda1_dotR1 - lambda2_dotR1) - (lambda1 - lambda2) * t2 * a_dotR1)
  tmp1_dotR2 = 0.5 * t1 * ( (lambda1_dotR2 - lambda2_dotR2) - (lambda1 - lambda2) * t2 * a_dotR2)
  tmp1_dotR3 = 0.5 * t1 * ( (lambda1_dotR3 - lambda2_dotR3) - (lambda1 - lambda2) * t2 * a_dotR3)
  tmp1_dotR4 = 0.5 * t1 * ( (lambda1_dotR4 - lambda2_dotR4) - (lambda1 - lambda2) * t2 * a_dotR4)


  #for i=1:length(sat)
  #  sat[i] = sat[i] + tmp1*(E1dq[i] + gami*E2dq[i])
  #end

  t1 = E1dq[1] + gami*E2dq[1]
  sat_jacL[1, 1] += tmp1*E1dq1_dotL1 + t1*tmp1_dotL1
  sat_jacL[1, 2] += tmp1*E1dq1_dotL2 + t1*tmp1_dotL2
  sat_jacL[1, 3] += tmp1*E1dq1_dotL3 + t1*tmp1_dotL3
  sat_jacL[1, 4] +=                    t1*tmp1_dotL4

  sat_jacR[1, 1] += tmp1*E1dq1_dotR1 + t1*tmp1_dotR1
  sat_jacR[1, 2] += tmp1*E1dq1_dotR2 + t1*tmp1_dotR2
  sat_jacR[1, 3] += tmp1*E1dq1_dotR3 + t1*tmp1_dotR3
  sat_jacR[1, 4] +=                    t1*tmp1_dotR4

  t1 = E1dq[2] + gami*E2dq[2]
  sat_jacL[2, 1] += tmp1*(E1dq2_dotL1 + gami*E2dq2_dotL1) + t1*tmp1_dotL1
  sat_jacL[2, 2] += tmp1*(E1dq2_dotL2 + gami*E2dq2_dotL2) + t1*tmp1_dotL2
  sat_jacL[2, 3] += tmp1*(E1dq2_dotL3 + gami*E2dq2_dotL3) + t1*tmp1_dotL3
  sat_jacL[2, 4] += tmp1*(            + gami*E2dq2_dotL4) + t1*tmp1_dotL4

  sat_jacR[2, 1] += tmp1*(E1dq2_dotR1 + gami*E2dq2_dotR1) + t1*tmp1_dotR1
  sat_jacR[2, 2] += tmp1*(E1dq2_dotR2 + gami*E2dq2_dotR2) + t1*tmp1_dotR2
  sat_jacR[2, 3] += tmp1*(E1dq2_dotR3 + gami*E2dq2_dotR3) + t1*tmp1_dotR3
  sat_jacR[2, 4] += tmp1*(            + gami*E2dq2_dotR4) + t1*tmp1_dotR4

  t1 = E1dq[3] + gami*E2dq[3]
  sat_jacL[3, 1] += tmp1*(E1dq3_dotL1 + gami*E2dq3_dotL1) + t1*tmp1_dotL1
  sat_jacL[3, 2] += tmp1*(E1dq3_dotL2 + gami*E2dq3_dotL2) + t1*tmp1_dotL2
  sat_jacL[3, 3] += tmp1*(E1dq3_dotL3 + gami*E2dq3_dotL3) + t1*tmp1_dotL3
  sat_jacL[3, 4] += tmp1*(            + gami*E2dq3_dotL4) + t1*tmp1_dotL4

  sat_jacR[3, 1] += tmp1*(E1dq3_dotR1 + gami*E2dq3_dotR1) + t1*tmp1_dotR1
  sat_jacR[3, 2] += tmp1*(E1dq3_dotR2 + gami*E2dq3_dotR2) + t1*tmp1_dotR2
  sat_jacR[3, 3] += tmp1*(E1dq3_dotR3 + gami*E2dq3_dotR3) + t1*tmp1_dotR3
  sat_jacR[3, 4] += tmp1*(            + gami*E2dq3_dotR4) + t1*tmp1_dotR4

  t1 = E1dq[4] + gami*E2dq[4]
  sat_jacL[4, 1] += tmp1*(E1dq4_dotL1 + gami*E2dq4_dotL1) + t1*tmp1_dotL1
  sat_jacL[4, 2] += tmp1*(E1dq4_dotL2 + gami*E2dq4_dotL2) + t1*tmp1_dotL2
  sat_jacL[4, 3] += tmp1*(E1dq4_dotL3 + gami*E2dq4_dotL3) + t1*tmp1_dotL3
  sat_jacL[4, 4] += tmp1*(E1dq4_dotL4 + gami*E2dq4_dotL4) + t1*tmp1_dotL4

  sat_jacR[4, 1] += tmp1*(E1dq4_dotR1 + gami*E2dq4_dotR1) + t1*tmp1_dotR1
  sat_jacR[4, 2] += tmp1*(E1dq4_dotR2 + gami*E2dq4_dotR2) + t1*tmp1_dotR2
  sat_jacR[4, 3] += tmp1*(E1dq4_dotR3 + gami*E2dq4_dotR3) + t1*tmp1_dotR3
  sat_jacR[4, 4] += tmp1*(E1dq4_dotR4 + gami*E2dq4_dotR4) + t1*tmp1_dotR4


  return nothing
end  # End function calcSAT


@fastmath function calcSAT_diff(params::ParamType{3},
                             roe_vars::AbstractArray{Tsol, 1},
                             roe_vars_dot::AbstractArray{Tsol, 1},
                             dq::AbstractArray{Tsol,1},
                             nrm::AbstractArray{Tmsh,1},
                             sat_jacL::AbstractArray{Tsol,2},
                             sat_jacR::AbstractArray{Tsol,2}) where {Tmsh, Tsol}
  # roe_vars = [u, v, w, H] at Roe average 
  # roe_vars_dot contains all the non-zero derivatives of the roe_vars packed
  # into a vector

  # dq_dotL* = 1, dq_dotR* = -1, so don't pass them explicitly

  data = params.calcsatdata
  @unpack data E1dq E2dq

  # SAT parameters
  sat_Vn = convert(Tsol, 0.025)
  sat_Vl = convert(Tsol, 0.025)
  tau = 1.0

  u = roe_vars[1]
  v = roe_vars[2]
  w = roe_vars[3]
  H = roe_vars[4]

  u_dotL1  = roe_vars_dot[1] 
  u_dotR1  = roe_vars_dot[2] 
  u_dotL2  = roe_vars_dot[3] 
  u_dotR2  = roe_vars_dot[4] 
                              
  v_dotL1  = roe_vars_dot[5] 
  v_dotR1  = roe_vars_dot[6] 
  v_dotL3  = roe_vars_dot[7] 
  v_dotR3  = roe_vars_dot[8] 
                              
  w_dotL1  = roe_vars_dot[9] 
  w_dotR1  = roe_vars_dot[10]
  w_dotL4  = roe_vars_dot[11]
  w_dotR4  = roe_vars_dot[12]
                              
                              
  H_dotL1  = roe_vars_dot[13]
  H_dotR1  = roe_vars_dot[14]
  H_dotL2  = roe_vars_dot[15]
  H_dotR2  = roe_vars_dot[16]
  H_dotL3  = roe_vars_dot[17]
  H_dotR3  = roe_vars_dot[18]
  H_dotL4  = roe_vars_dot[19]
  H_dotR4  = roe_vars_dot[20]
  H_dotL5  = roe_vars_dot[21]
  H_dotR5  = roe_vars_dot[22]
  
  gami = params.gamma_1

  # Begin main execution
  nx = nrm[1]
  ny = nrm[2]
  nz = nrm[3]

  dA = sqrt(nx*nx + ny*ny + nz*nz)

  Un = u*nx + v*ny + w*nz  # Normal Velocity
  Un_dotL1 = nx*u_dotL1 + ny*v_dotL1 + nz*w_dotL1
  Un_dotR1 = nx*u_dotR1 + ny*v_dotR1 + nz*w_dotR1

  Un_dotL2 = nx*u_dotL2
  Un_dotR2 = nx*u_dotR2

  Un_dotL3 = ny*v_dotL3
  Un_dotR3 = ny*v_dotR3

  Un_dotL4 = nz*w_dotL4
  Un_dotR4 = nz*w_dotR4

  phi = 0.5*(u*u + v*v + w*w)

  phi_dotL1 = u*u_dotL1 + v*v_dotL1 + w*w_dotL1
  phi_dotR1 = u*u_dotR1 + v*v_dotR1 + w*w_dotR1

  phi_dotL2 = u*u_dotL2
  phi_dotR2 = u*u_dotR2
  
  phi_dotL3 = v*v_dotL3
  phi_dotR3 = v*v_dotR3

  phi_dotL4 = w*w_dotL4
  phi_dotR4 = w*w_dotR4


  a = sqrt(gami*(H - phi)) # speed of sound
  t1 = gami/(2*a)
  a_dotL1 = t1*(H_dotL1 - phi_dotL1)
  a_dotR1 = t1*(H_dotR1 - phi_dotR1)

  a_dotL2 = t1*(H_dotL2 - phi_dotL2)
  a_dotR2 = t1*(H_dotR2 - phi_dotR2)

  a_dotL3 = t1*(H_dotL3 - phi_dotL3)
  a_dotR3 = t1*(H_dotR3 - phi_dotR3)

  a_dotL4 = t1*(H_dotL4 - phi_dotL4)
  a_dotR4 = t1*(H_dotR4 - phi_dotR4)

  a_dotL5 = t1*H_dotL5
  a_dotR5 = t1*H_dotR5

  lambda1 = Un + dA*a
  lambda1_dotL1 = Un_dotL1 + dA*a_dotL1
  lambda1_dotR1 = Un_dotR1 + dA*a_dotR1

  lambda1_dotL2 = Un_dotL2 + dA*a_dotL2
  lambda1_dotR2 = Un_dotR2 + dA*a_dotR2

  lambda1_dotL3 = Un_dotL3 + dA*a_dotL3
  lambda1_dotR3 = Un_dotR3 + dA*a_dotR3

  lambda1_dotL4 = Un_dotL4 + dA*a_dotL4
  lambda1_dotR4 = Un_dotR4 + dA*a_dotR4

  lambda1_dotL5 = dA*a_dotL5
  lambda1_dotR5 = dA*a_dotR5

  
  lambda2 = Un - dA*a
  lambda2_dotL1 = Un_dotL1 - dA*a_dotL1
  lambda2_dotR1 = Un_dotR1 - dA*a_dotR1

  lambda2_dotL2 = Un_dotL2 - dA*a_dotL2
  lambda2_dotR2 = Un_dotR2 - dA*a_dotR2

  lambda2_dotL3 = Un_dotL3 - dA*a_dotL3
  lambda2_dotR3 = Un_dotR3 - dA*a_dotR3

  lambda2_dotL4 = Un_dotL4 - dA*a_dotL4
  lambda2_dotR4 = Un_dotR4 - dA*a_dotR4

  lambda2_dotL5 = -dA*a_dotL5
  lambda2_dotR5 = -dA*a_dotR5


  lambda3 = Un
  lambda3_dotL1 = Un_dotL1
  lambda3_dotR1 = Un_dotR1

  lambda3_dotL2 = Un_dotL2
  lambda3_dotR2 = Un_dotR2

  lambda3_dotL3 = Un_dotL3
  lambda3_dotR3 = Un_dotR3

  lambda3_dotL4 = Un_dotL4
  lambda3_dotR4 = Un_dotR4

  rhoA = absvalue(Un) + dA*a
  #TODO: see if there is a better way to do this
  if Un > 0
    fac = 1
  else
    fac = -1
  end

  rhoA_dotL1 = fac*Un_dotL1 + dA*a_dotL1
  rhoA_dotR1 = fac*Un_dotR1 + dA*a_dotR1

  rhoA_dotL2 = fac*Un_dotL2 + dA*a_dotL2
  rhoA_dotR2 = fac*Un_dotR2 + dA*a_dotR2

  rhoA_dotL3 = fac*Un_dotL3 + dA*a_dotL3
  rhoA_dotR3 = fac*Un_dotR3 + dA*a_dotR3

  rhoA_dotL4 = fac*Un_dotL4 + dA*a_dotL4
  rhoA_dotR4 = fac*Un_dotR4 + dA*a_dotR4

  rhoA_dotL5 = dA*a_dotL5
  rhoA_dotR5 = dA*a_dotR5

  # Compute Eigen Values of the Flux Jacobian
  # The eigen values calculated above cannot be used directly. Near stagnation
  # points lambda3 approaches zero while near sonic lines lambda1 and lambda2
  # approach zero. This has a possibility of creating numerical difficulties.
  # As a result, the eigen values are limited by the following expressions.


  # see lambda1 expression below
  if absvalue(lambda1) > sat_Vn*rhoA
    if lambda1 > 0
      fac = 1
    else
      fac = -1
    end

    t1 = tau*fac
    lambda1_dotL1 = 0.5 * (t1 * lambda1_dotL1 - lambda1_dotL1)
    lambda1_dotR1 = 0.5 * (t1 * lambda1_dotR1 - lambda1_dotR1)

    lambda1_dotL2 = 0.5 * (t1 * lambda1_dotL2 - lambda1_dotL2)
    lambda1_dotR2 = 0.5 * (t1 * lambda1_dotR2 - lambda1_dotR2)

    lambda1_dotL3 = 0.5 * (t1 * lambda1_dotL3 - lambda1_dotL3)
    lambda1_dotR3 = 0.5 * (t1 * lambda1_dotR3 - lambda1_dotR3)

    lambda1_dotL4 = 0.5 * (t1 * lambda1_dotL4 - lambda1_dotL4)
    lambda1_dotR4 = 0.5 * (t1 * lambda1_dotR4 - lambda1_dotR4)

    lambda1_dotL5 = 0.5 * (t1 * lambda1_dotL5 - lambda1_dotL5)
    lambda1_dotR5 = 0.5 * (t1 * lambda1_dotR5 - lambda1_dotR5)


  else
    t1 = sat_Vn*tau
    lambda1_dotL1 = 0.5 * (t1 * rhoA_dotL1 - lambda1_dotL1)
    lambda1_dotR1 = 0.5 * (t1 * rhoA_dotR1 - lambda1_dotR1)
 
    lambda1_dotL2 = 0.5 * (t1 * rhoA_dotL2 - lambda1_dotL2)
    lambda1_dotR2 = 0.5 * (t1 * rhoA_dotR2 - lambda1_dotR2)
    
    lambda1_dotL3 = 0.5 * (t1 * rhoA_dotL3 - lambda1_dotL3)
    lambda1_dotR3 = 0.5 * (t1 * rhoA_dotR3 - lambda1_dotR3)
 
    lambda1_dotL4 = 0.5 * (t1 * rhoA_dotL4 - lambda1_dotL4)
    lambda1_dotR4 = 0.5 * (t1 * rhoA_dotR4 - lambda1_dotR4)
 
    lambda1_dotL5 = 0.5 * (t1 * rhoA_dotL5 - lambda1_dotL5)
    lambda1_dotR5 = 0.5 * (t1 * rhoA_dotR5 - lambda1_dotR5)
 
  end

  lambda1 = 0.5*(tau*max(absvalue(lambda1),sat_Vn *rhoA) - lambda1)

  # see lambda2 expression below
  if absvalue(lambda2) > sat_Vn*rhoA
    if lambda2 > 0
      fac = 1
    else
      fac = -1
    end

    t1 = tau*fac
    lambda2_dotL1 = 0.5 * (t1 * lambda2_dotL1 - lambda2_dotL1)
    lambda2_dotR1 = 0.5 * (t1 * lambda2_dotR1 - lambda2_dotR1)

    lambda2_dotL2 = 0.5 * (t1 * lambda2_dotL2 - lambda2_dotL2)
    lambda2_dotR2 = 0.5 * (t1 * lambda2_dotR2 - lambda2_dotR2)

    lambda2_dotL3 = 0.5 * (t1 * lambda2_dotL3 - lambda2_dotL3)
    lambda2_dotR3 = 0.5 * (t1 * lambda2_dotR3 - lambda2_dotR3)

    lambda2_dotL4 = 0.5 * (t1 * lambda2_dotL4 - lambda2_dotL4)
    lambda2_dotR4 = 0.5 * (t1 * lambda2_dotR4 - lambda2_dotR4)

    lambda2_dotL5 = 0.5 * (t1 * lambda2_dotL5 - lambda2_dotL5)
    lambda2_dotR5 = 0.5 * (t1 * lambda2_dotR5 - lambda2_dotR5)

  else

    t1 = sat_Vn*tau
    lambda2_dotL1 = 0.5 * (t1 * rhoA_dotL1 - lambda2_dotL1)
    lambda2_dotR1 = 0.5 * (t1 * rhoA_dotR1 - lambda2_dotR1)
 
    lambda2_dotL2 = 0.5 * (t1 * rhoA_dotL2 - lambda2_dotL2)
    lambda2_dotR2 = 0.5 * (t1 * rhoA_dotR2 - lambda2_dotR2)
    
    lambda2_dotL3 = 0.5 * (t1 * rhoA_dotL3 - lambda2_dotL3)
    lambda2_dotR3 = 0.5 * (t1 * rhoA_dotR3 - lambda2_dotR3)
 
    lambda2_dotL4 = 0.5 * (t1 * rhoA_dotL4 - lambda2_dotL4)
    lambda2_dotR4 = 0.5 * (t1 * rhoA_dotR4 - lambda2_dotR4)
 
    lambda2_dotL5 = 0.5 * (t1 * rhoA_dotL5 - lambda2_dotL5)
    lambda2_dotR5 = 0.5 * (t1 * rhoA_dotR5 - lambda2_dotR5)
 
  end
  lambda2 = 0.5*(tau*max(absvalue(lambda2),sat_Vn *rhoA) - lambda2)


  # see lambda3 expression below
  if absvalue(lambda3) > sat_Vn*rhoA
    if lambda3 > 0
      fac = 1
    else
      fac = -1
    end

    t1 = tau*fac
    lambda3_dotL1 = 0.5 * (t1 * lambda3_dotL1 - lambda3_dotL1)
    lambda3_dotR1 = 0.5 * (t1 * lambda3_dotR1 - lambda3_dotR1)

    lambda3_dotL2 = 0.5 * (t1 * lambda3_dotL2 - lambda3_dotL2)
    lambda3_dotR2 = 0.5 * (t1 * lambda3_dotR2 - lambda3_dotR2)

    lambda3_dotL3 = 0.5 * (t1 * lambda3_dotL3 - lambda3_dotL3)
    lambda3_dotR3 = 0.5 * (t1 * lambda3_dotR3 - lambda3_dotR3)

    lambda3_dotL4 = 0.5 * (t1 * lambda3_dotL4 - lambda3_dotL4)
    lambda3_dotR4 = 0.5 * (t1 * lambda3_dotR4 - lambda3_dotR4)

    lambda3_dotL5 = 0.0
    lambda3_dotR5 = 0.0


  else
    t1 = sat_Vn*tau
    lambda3_dotL1 = 0.5 * (t1 * rhoA_dotL1 - lambda3_dotL1)
    lambda3_dotR1 = 0.5 * (t1 * rhoA_dotR1 - lambda3_dotR1)
 
    lambda3_dotL2 = 0.5 * (t1 * rhoA_dotL2 - lambda3_dotL2)
    lambda3_dotR2 = 0.5 * (t1 * rhoA_dotR2 - lambda3_dotR2)
    
    lambda3_dotL3 = 0.5 * (t1 * rhoA_dotL3 - lambda3_dotL3)
    lambda3_dotR3 = 0.5 * (t1 * rhoA_dotR3 - lambda3_dotR3)
 
    lambda3_dotL4 = 0.5 * (t1 * rhoA_dotL4 - lambda3_dotL4)
    lambda3_dotR4 = 0.5 * (t1 * rhoA_dotR4 - lambda3_dotR4)
 
    lambda3_dotL5 = 0.5 * t1 * rhoA_dotL5
    lambda3_dotR5 = 0.5 * t1 * rhoA_dotR5
 
  end

  lambda3 = 0.5*(tau*max(absvalue(lambda3),sat_Vl *rhoA) - lambda3)

                    
  dq1 = dq[1]
  dq2 = dq[2]
  dq3 = dq[3]
  dq4 = dq[4]
  dq5 = dq[5]

  # sat[1] = lambda3*dq1
  sat_jacL[1, 1] += lambda3 + dq1*lambda3_dotL1
  sat_jacL[1, 2] +=         + dq1*lambda3_dotL2
  sat_jacL[1, 3] +=         + dq1*lambda3_dotL3
  sat_jacL[1, 4] +=         + dq1*lambda3_dotL4
  sat_jacL[1, 5] +=         + dq1*lambda3_dotL5

  sat_jacR[1, 1] += -lambda3 + dq1*lambda3_dotR1
  sat_jacR[1, 2] +=          + dq1*lambda3_dotR2
  sat_jacR[1, 3] +=          + dq1*lambda3_dotR3
  sat_jacR[1, 4] +=          + dq1*lambda3_dotR4
  sat_jacR[1, 5] +=          + dq1*lambda3_dotR5

  # sat[2] = lambda3*dq2
  sat_jacL[2, 1] +=         + dq2*lambda3_dotL1
  sat_jacL[2, 2] += lambda3 + dq2*lambda3_dotL2
  sat_jacL[2, 3] +=         + dq2*lambda3_dotL3
  sat_jacL[2, 4] +=         + dq2*lambda3_dotL4
  sat_jacL[2, 5] +=         + dq2*lambda3_dotL5

  sat_jacR[2, 1] +=          + dq2*lambda3_dotR1
  sat_jacR[2, 2] += -lambda3 + dq2*lambda3_dotR2
  sat_jacR[2, 3] +=          + dq2*lambda3_dotR3
  sat_jacR[2, 4] +=          + dq2*lambda3_dotR4
  sat_jacR[2, 5] +=          + dq2*lambda3_dotR5

  # sat[3] = lambda3*dq3
  sat_jacL[3, 1] +=         + dq3*lambda3_dotL1
  sat_jacL[3, 2] +=         + dq3*lambda3_dotL2
  sat_jacL[3, 3] += lambda3 + dq3*lambda3_dotL3
  sat_jacL[3, 4] +=         + dq3*lambda3_dotL4
  sat_jacL[3, 5] +=         + dq3*lambda3_dotL5

  sat_jacR[3, 1] +=          + dq3*lambda3_dotR1
  sat_jacR[3, 2] +=          + dq3*lambda3_dotR2
  sat_jacR[3, 3] += -lambda3 + dq3*lambda3_dotR3
  sat_jacR[3, 4] +=          + dq3*lambda3_dotR4
  sat_jacR[3, 5] +=          + dq3*lambda3_dotR5

  sat_jacL[4, 1] +=         + dq4*lambda3_dotL1
  sat_jacL[4, 2] +=         + dq4*lambda3_dotL2
  sat_jacL[4, 3] +=         + dq4*lambda3_dotL3
  sat_jacL[4, 4] += lambda3 + dq4*lambda3_dotL4
  sat_jacL[4, 5] +=         + dq4*lambda3_dotL5

  sat_jacR[4, 1] +=          + dq4*lambda3_dotR1
  sat_jacR[4, 2] +=          + dq4*lambda3_dotR2
  sat_jacR[4, 3] +=          + dq4*lambda3_dotR3
  sat_jacR[4, 4] += -lambda3 + dq4*lambda3_dotR4
  sat_jacR[4, 5] +=          + dq4*lambda3_dotR5

  # sat[5] = lambda3*dq5
  sat_jacL[5, 1] +=           dq5*lambda3_dotL1
  sat_jacL[5, 2] +=           dq5*lambda3_dotL2
  sat_jacL[5, 3] +=           dq5*lambda3_dotL3
  sat_jacL[5, 4] +=           dq5*lambda3_dotL4
  sat_jacL[5, 5] += lambda3 + dq5*lambda3_dotL5

  sat_jacR[5, 1] +=            dq5*lambda3_dotR1
  sat_jacR[5, 2] +=            dq5*lambda3_dotR2
  sat_jacR[5, 3] +=            dq5*lambda3_dotR3
  sat_jacR[5, 4] +=            dq5*lambda3_dotR4
  sat_jacR[5, 5] += -lambda3 + dq5*lambda3_dotR5

  #-- get E1*dq
  E1dq[1] = phi*dq1 - u*dq2 - v*dq3  -w*dq4 + dq5
  E1dq1_dotL1 = phi + dq1*phi_dotL1     - dq2*u_dotL1     - dq3*v_dotL1     - dq4*w_dotL1
  E1dq1_dotL2 =       dq1*phi_dotL2 - u - dq2*u_dotL2
  E1dq1_dotL3 =       dq1*phi_dotL3                   - v - dq3*v_dotL3
  E1dq1_dotL4 =       dq1*phi_dotL4                                     - w - dq4*w_dotL4
  E1dq1_dotL5 = 1

  E1dq1_dotR1 = -phi + dq1*phi_dotR1     - dq2*u_dotR1     - dq3*v_dotR1     - dq4*w_dotR1
  E1dq1_dotR2 =        dq1*phi_dotR2 + u - dq2*u_dotR2
  E1dq1_dotR3 =        dq1*phi_dotR3                   + v - dq3*v_dotR3
  E1dq1_dotR4 =        dq1*phi_dotR4                                     + w - dq4*w_dotR4
  E1dq1_dotR5 = -1


  E1dq[2] = E1dq[1]*u
  E1dq2_dotL1 = u*E1dq1_dotL1 + E1dq[1]*u_dotL1
  E1dq2_dotL2 = u*E1dq1_dotL2 + E1dq[1]*u_dotL2
  E1dq2_dotL3 = u*E1dq1_dotL3
  E1dq2_dotL4 = u*E1dq1_dotL4
  E1dq2_dotL5 = u*E1dq1_dotL5

  E1dq2_dotR1 = u*E1dq1_dotR1 + E1dq[1]*u_dotR1
  E1dq2_dotR2 = u*E1dq1_dotR2 + E1dq[1]*u_dotR2
  E1dq2_dotR3 = u*E1dq1_dotR3
  E1dq2_dotR4 = u*E1dq1_dotR4
  E1dq2_dotR5 = u*E1dq1_dotR5

  E1dq[3] = E1dq[1]*v
  E1dq3_dotL1 = v*E1dq1_dotL1 + E1dq[1]*v_dotL1
  E1dq3_dotL2 = v*E1dq1_dotL2
  E1dq3_dotL3 = v*E1dq1_dotL3 + E1dq[1]*v_dotL3
  E1dq3_dotL4 = v*E1dq1_dotL4
  E1dq3_dotL5 = v*E1dq1_dotL5

  E1dq3_dotR1 = v*E1dq1_dotR1 + E1dq[1]*v_dotR1
  E1dq3_dotR2 = v*E1dq1_dotR2
  E1dq3_dotR3 = v*E1dq1_dotR3 + E1dq[1]*v_dotR3
  E1dq3_dotR4 = v*E1dq1_dotR4
  E1dq3_dotR5 = v*E1dq1_dotR5

  E1dq[4] = E1dq[1]*w
  E1dq4_dotL1 = w*E1dq1_dotL1 + E1dq[1]*w_dotL1
  E1dq4_dotL2 = w*E1dq1_dotL2
  E1dq4_dotL3 = w*E1dq1_dotL3
  E1dq4_dotL4 = w*E1dq1_dotL4 + E1dq[1]*w_dotL4
  E1dq4_dotL5 = w*E1dq1_dotL5

  E1dq4_dotR1 = w*E1dq1_dotR1 + E1dq[1]*w_dotR1
  E1dq4_dotR2 = w*E1dq1_dotR2
  E1dq4_dotR3 = w*E1dq1_dotR3
  E1dq4_dotR4 = w*E1dq1_dotR4 + E1dq[1]*w_dotR4
  E1dq4_dotR5 = w*E1dq1_dotR5


  E1dq[5] = E1dq[1]*H
  E1dq5_dotL1 = H*E1dq1_dotL1 + E1dq[1]*H_dotL1
  E1dq5_dotL2 = H*E1dq1_dotL2 + E1dq[1]*H_dotL2
  E1dq5_dotL3 = H*E1dq1_dotL3 + E1dq[1]*H_dotL3
  E1dq5_dotL4 = H*E1dq1_dotL4 + E1dq[1]*H_dotL4
  E1dq5_dotL5 = H*E1dq1_dotL5 + E1dq[1]*H_dotL5

  E1dq5_dotR1 = H*E1dq1_dotR1 + E1dq[1]*H_dotR1
  E1dq5_dotR2 = H*E1dq1_dotR2 + E1dq[1]*H_dotR2
  E1dq5_dotR3 = H*E1dq1_dotR3 + E1dq[1]*H_dotR3
  E1dq5_dotR4 = H*E1dq1_dotR4 + E1dq[1]*H_dotR4
  E1dq5_dotR5 = H*E1dq1_dotR5 + E1dq[1]*H_dotR5



  #-- get E2*dq
  E2dq[1] = 0.0
  E2dq[2] = -Un*dq1 + nx*dq2 + ny*dq3 + nz*dq4
  E2dq2_dotL1 = -Un + -Un_dotL1*dq1
  E2dq2_dotL2 =     + -Un_dotL2*dq1 + nx
  E2dq2_dotL3 =     + -Un_dotL3*dq1      + ny
  E2dq2_dotL4 =     + -Un_dotL4*dq1           + nz

  E2dq2_dotR1 = Un + -Un_dotR1*dq1
  E2dq2_dotR2 =    + -Un_dotR2*dq1 - nx
  E2dq2_dotR3 =    + -Un_dotR3*dq1      - ny
  E2dq2_dotR4 =     + -Un_dotR4*dq1           - nz


  E2dq[3] = E2dq[2]*ny
  E2dq3_dotL1 = ny*E2dq2_dotL1
  E2dq3_dotL2 = ny*E2dq2_dotL2
  E2dq3_dotL3 = ny*E2dq2_dotL3
  E2dq3_dotL4 = ny*E2dq2_dotL4

  E2dq3_dotR1 = ny*E2dq2_dotR1
  E2dq3_dotR2 = ny*E2dq2_dotR2
  E2dq3_dotR3 = ny*E2dq2_dotR3
  E2dq3_dotR4 = ny*E2dq2_dotR4

  E2dq[4] = E2dq[2]*nz
  E2dq4_dotL1 = nz*E2dq2_dotL1
  E2dq4_dotL2 = nz*E2dq2_dotL2
  E2dq4_dotL3 = nz*E2dq2_dotL3
  E2dq4_dotL4 = nz*E2dq2_dotL4

  E2dq4_dotR1 = nz*E2dq2_dotR1
  E2dq4_dotR2 = nz*E2dq2_dotR2
  E2dq4_dotR3 = nz*E2dq2_dotR3
  E2dq4_dotR4 = nz*E2dq2_dotR4


  E2dq[5] = E2dq[2]*Un
  E2dq5_dotL1 = Un*E2dq2_dotL1 + E2dq[2]*Un_dotL1
  E2dq5_dotL2 = Un*E2dq2_dotL2 + E2dq[2]*Un_dotL2
  E2dq5_dotL3 = Un*E2dq2_dotL3 + E2dq[2]*Un_dotL3
  E2dq5_dotL4 = Un*E2dq2_dotL4 + E2dq[2]*Un_dotL4

  E2dq5_dotR1 = Un*E2dq2_dotR1 + E2dq[2]*Un_dotR1
  E2dq5_dotR2 = Un*E2dq2_dotR2 + E2dq[2]*Un_dotR2
  E2dq5_dotR3 = Un*E2dq2_dotR3 + E2dq[2]*Un_dotR3
  E2dq5_dotR4 = Un*E2dq2_dotR4 + E2dq[2]*Un_dotR4

  E2dq[2] = E2dq[2]*nx
  E2dq2_dotL1 = nx*E2dq2_dotL1
  E2dq2_dotL2 = nx*E2dq2_dotL2
  E2dq2_dotL3 = nx*E2dq2_dotL3
  E2dq2_dotL4 = nx*E2dq2_dotL4

  E2dq2_dotR1 = nx*E2dq2_dotR1
  E2dq2_dotR2 = nx*E2dq2_dotR2
  E2dq2_dotR3 = nx*E2dq2_dotR3
  E2dq2_dotR4 = nx*E2dq2_dotR4

  #-- add to sat
  tmp1 = 0.5*(lambda1 + lambda2) - lambda3
  tmp1_dotL1 = 0.5*(lambda1_dotL1 + lambda2_dotL1) - lambda3_dotL1
  tmp1_dotL2 = 0.5*(lambda1_dotL2 + lambda2_dotL2) - lambda3_dotL2
  tmp1_dotL3 = 0.5*(lambda1_dotL3 + lambda2_dotL3) - lambda3_dotL3
  tmp1_dotL4 = 0.5*(lambda1_dotL4 + lambda2_dotL4) - lambda3_dotL4
  tmp1_dotL5 = 0.5*(lambda1_dotL5 + lambda2_dotL5) - lambda3_dotL5

  tmp1_dotR1 = 0.5*(lambda1_dotR1 + lambda2_dotR1) - lambda3_dotR1
  tmp1_dotR2 = 0.5*(lambda1_dotR2 + lambda2_dotR2) - lambda3_dotR2
  tmp1_dotR3 = 0.5*(lambda1_dotR3 + lambda2_dotR3) - lambda3_dotR3
  tmp1_dotR4 = 0.5*(lambda1_dotR4 + lambda2_dotR4) - lambda3_dotR4
  tmp1_dotR5 = 0.5*(lambda1_dotR5 + lambda2_dotR5) - lambda3_dotR5


  tmp2 = gami/(a*a)
  t1 = -2*tmp2/a
#  t1 = -2*gami/(a*a*a) # = -2*tmp2/a
  tmp2_dotL1 = t1*a_dotL1
  tmp2_dotL2 = t1*a_dotL2
  tmp2_dotL3 = t1*a_dotL3
  tmp2_dotL4 = t1*a_dotL4
  tmp2_dotL5 = t1*a_dotL5

  tmp2_dotR1 = t1*a_dotR1
  tmp2_dotR2 = t1*a_dotR2
  tmp2_dotR3 = t1*a_dotR3
  tmp2_dotR4 = t1*a_dotR4
  tmp2_dotR5 = t1*a_dotR5

  tmp3 = 1.0/(dA*dA)
  #for i=1:length(sat)
  #  sat[i] = sat[i] + tmp1*(tmp2*E1dq[i] + tmp3*E2dq[i])
  #end

  sat_jacL[1,1] += tmp1*tmp2*E1dq1_dotL1   + tmp1*E1dq[1]*tmp2_dotL1 + 
                   tmp2*E1dq[1]*tmp1_dotL1 + 
                                           + tmp3*E2dq[1]*tmp1_dotL1

  sat_jacL[1,2] += tmp1*tmp2*E1dq1_dotL2   + tmp1*E1dq[1]*tmp2_dotL2 + 
                   tmp2*E1dq[1]*tmp1_dotL2 + 
                                           + tmp3*E2dq[1]*tmp1_dotL2

  sat_jacL[1,3] += tmp1*tmp2*E1dq1_dotL3   + tmp1*E1dq[1]*tmp2_dotL3 + 
                   tmp2*E1dq[1]*tmp1_dotL3 + 
                                           + tmp3*E2dq[1]*tmp1_dotL3

  sat_jacL[1,4] += tmp1*tmp2*E1dq1_dotL4   + tmp1*E1dq[1]*tmp2_dotL4 + 
                   tmp2*E1dq[1]*tmp1_dotL4 + 
                                           + tmp3*E2dq[1]*tmp1_dotL4

  sat_jacL[1,5] += tmp1*tmp2*E1dq1_dotL5   + tmp1*E1dq[1]*tmp2_dotL5 + 
                   tmp2*E1dq[1]*tmp1_dotL5 + 
                                           + tmp3*E2dq[1]*tmp1_dotL5

  sat_jacR[1,1] += tmp1*tmp2*E1dq1_dotR1   + tmp1*E1dq[1]*tmp2_dotR1 + 
                   tmp2*E1dq[1]*tmp1_dotR1 + 
                                           + tmp3*E2dq[1]*tmp1_dotR1

  sat_jacR[1,2] += tmp1*tmp2*E1dq1_dotR2   + tmp1*E1dq[1]*tmp2_dotR2 + 
                   tmp2*E1dq[1]*tmp1_dotR2 + 
                                           + tmp3*E2dq[1]*tmp1_dotR2

  sat_jacR[1,3] += tmp1*tmp2*E1dq1_dotR3   + tmp1*E1dq[1]*tmp2_dotR3 + 
                   tmp2*E1dq[1]*tmp1_dotR3 + 
                                           + tmp3*E2dq[1]*tmp1_dotR3

  sat_jacR[1,4] += tmp1*tmp2*E1dq1_dotR4   + tmp1*E1dq[1]*tmp2_dotR4 + 
                   tmp2*E1dq[1]*tmp1_dotR4 + 
                                           + tmp3*E2dq[1]*tmp1_dotR4

  sat_jacR[1,5] += tmp1*tmp2*E1dq1_dotR5   + tmp1*E1dq[1]*tmp2_dotR5 + 
                   tmp2*E1dq[1]*tmp1_dotR5 + 
                                           + tmp3*E2dq[1]*tmp1_dotR5


  sat_jacL[2,1] += tmp1*tmp2*E1dq2_dotL1   + tmp1*E1dq[2]*tmp2_dotL1 + 
                   tmp2*E1dq[2]*tmp1_dotL1 + 
                   tmp1*tmp3*E2dq2_dotL1   + tmp3*E2dq[2]*tmp1_dotL1

  sat_jacL[2,2] += tmp1*tmp2*E1dq2_dotL2   + tmp1*E1dq[2]*tmp2_dotL2 + 
                   tmp2*E1dq[2]*tmp1_dotL2 + 
                   tmp1*tmp3*E2dq2_dotL2   + tmp3*E2dq[2]*tmp1_dotL2

  sat_jacL[2,3] += tmp1*tmp2*E1dq2_dotL3   + tmp1*E1dq[2]*tmp2_dotL3 + 
                   tmp2*E1dq[2]*tmp1_dotL3 + 
                   tmp1*tmp3*E2dq2_dotL3   + tmp3*E2dq[2]*tmp1_dotL3

  sat_jacL[2,4] += tmp1*tmp2*E1dq2_dotL4   + tmp1*E1dq[2]*tmp2_dotL4 + 
                   tmp2*E1dq[2]*tmp1_dotL4 + 
                   tmp1*tmp3*E2dq2_dotL4   + tmp3*E2dq[2]*tmp1_dotL4

  sat_jacL[2,5] += tmp1*tmp2*E1dq2_dotL5   + tmp1*E1dq[2]*tmp2_dotL5 + 
                   tmp2*E1dq[2]*tmp1_dotL5 + 
                                           + tmp3*E2dq[2]*tmp1_dotL5

  sat_jacR[2,1] += tmp1*tmp2*E1dq2_dotR1   + tmp1*E1dq[2]*tmp2_dotR1 + 
                   tmp2*E1dq[2]*tmp1_dotR1 + 
                   tmp1*tmp3*E2dq2_dotR1   + tmp3*E2dq[2]*tmp1_dotR1

  sat_jacR[2,2] += tmp1*tmp2*E1dq2_dotR2   + tmp1*E1dq[2]*tmp2_dotR2 + 
                   tmp2*E1dq[2]*tmp1_dotR2 + 
                   tmp1*tmp3*E2dq2_dotR2   + tmp3*E2dq[2]*tmp1_dotR2

  sat_jacR[2,3] += tmp1*tmp2*E1dq2_dotR3   + tmp1*E1dq[2]*tmp2_dotR3 + 
                   tmp2*E1dq[2]*tmp1_dotR3 + 
                   tmp1*tmp3*E2dq2_dotR3   + tmp3*E2dq[2]*tmp1_dotR3

  sat_jacR[2,4] += tmp1*tmp2*E1dq2_dotR4   + tmp1*E1dq[2]*tmp2_dotR4 + 
                   tmp2*E1dq[2]*tmp1_dotR4 + 
                   tmp1*tmp3*E2dq2_dotR4   + tmp3*E2dq[2]*tmp1_dotR4

  sat_jacR[2,5] += tmp1*tmp2*E1dq2_dotR5   + tmp1*E1dq[2]*tmp2_dotR5 + 
                   tmp2*E1dq[2]*tmp1_dotR5 + 
                                           + tmp3*E2dq[2]*tmp1_dotR5


  sat_jacL[3,1] += tmp1*tmp2*E1dq3_dotL1   + tmp1*E1dq[3]*tmp2_dotL1 + 
                   tmp2*E1dq[3]*tmp1_dotL1 + 
                   tmp1*tmp3*E2dq3_dotL1   + tmp3*E2dq[3]*tmp1_dotL1

  sat_jacL[3,2] += tmp1*tmp2*E1dq3_dotL2   + tmp1*E1dq[3]*tmp2_dotL2 + 
                   tmp2*E1dq[3]*tmp1_dotL2 + 
                   tmp1*tmp3*E2dq3_dotL2   + tmp3*E2dq[3]*tmp1_dotL2

  sat_jacL[3,3] += tmp1*tmp2*E1dq3_dotL3   + tmp1*E1dq[3]*tmp2_dotL3 + 
                   tmp2*E1dq[3]*tmp1_dotL3 + 
                   tmp1*tmp3*E2dq3_dotL3   + tmp3*E2dq[3]*tmp1_dotL3

  sat_jacL[3,4] += tmp1*tmp2*E1dq3_dotL4   + tmp1*E1dq[3]*tmp2_dotL4 + 
                   tmp2*E1dq[3]*tmp1_dotL4 + 
                   tmp1*tmp3*E2dq3_dotL4   + tmp3*E2dq[3]*tmp1_dotL4

  sat_jacL[3,5] += tmp1*tmp2*E1dq3_dotL5   + tmp1*E1dq[3]*tmp2_dotL5 + 
                   tmp2*E1dq[3]*tmp1_dotL5 + 
                                           + tmp3*E2dq[3]*tmp1_dotL5

  sat_jacR[3,1] += tmp1*tmp2*E1dq3_dotR1   + tmp1*E1dq[3]*tmp2_dotR1 + 
                   tmp2*E1dq[3]*tmp1_dotR1 + 
                   tmp1*tmp3*E2dq3_dotR1   + tmp3*E2dq[3]*tmp1_dotR1

  sat_jacR[3,2] += tmp1*tmp2*E1dq3_dotR2   + tmp1*E1dq[3]*tmp2_dotR2 + 
                   tmp2*E1dq[3]*tmp1_dotR2 + 
                   tmp1*tmp3*E2dq3_dotR2   + tmp3*E2dq[3]*tmp1_dotR2

  sat_jacR[3,3] += tmp1*tmp2*E1dq3_dotR3   + tmp1*E1dq[3]*tmp2_dotR3 + 
                   tmp2*E1dq[3]*tmp1_dotR3 + 
                   tmp1*tmp3*E2dq3_dotR3   + tmp3*E2dq[3]*tmp1_dotR3

  sat_jacR[3,4] += tmp1*tmp2*E1dq3_dotR4   + tmp1*E1dq[3]*tmp2_dotR4 + 
                   tmp2*E1dq[3]*tmp1_dotR4 + 
                   tmp1*tmp3*E2dq3_dotR4   + tmp3*E2dq[3]*tmp1_dotR4

  sat_jacR[3,5] += tmp1*tmp2*E1dq3_dotR5   + tmp1*E1dq[3]*tmp2_dotR5 + 
                   tmp2*E1dq[3]*tmp1_dotR5 + 
                                           + tmp3*E2dq[3]*tmp1_dotR5


  sat_jacL[4,1] += tmp1*tmp2*E1dq4_dotL1   + tmp1*E1dq[4]*tmp2_dotL1 + 
                   tmp2*E1dq[4]*tmp1_dotL1 + 
                   tmp1*tmp3*E2dq4_dotL1   + tmp3*E2dq[4]*tmp1_dotL1

  sat_jacL[4,2] += tmp1*tmp2*E1dq4_dotL2   + tmp1*E1dq[4]*tmp2_dotL2 + 
                   tmp2*E1dq[4]*tmp1_dotL2 + 
                   tmp1*tmp3*E2dq4_dotL2   + tmp3*E2dq[4]*tmp1_dotL2

  sat_jacL[4,3] += tmp1*tmp2*E1dq4_dotL3   + tmp1*E1dq[4]*tmp2_dotL3 + 
                   tmp2*E1dq[4]*tmp1_dotL3 + 
                   tmp1*tmp3*E2dq4_dotL3   + tmp3*E2dq[4]*tmp1_dotL3

  sat_jacL[4,4] += tmp1*tmp2*E1dq4_dotL4   + tmp1*E1dq[4]*tmp2_dotL4 + 
                   tmp2*E1dq[4]*tmp1_dotL4 + 
                   tmp1*tmp3*E2dq4_dotL4   + tmp3*E2dq[4]*tmp1_dotL4

  sat_jacL[4,5] += tmp1*tmp2*E1dq4_dotL5   + tmp1*E1dq[4]*tmp2_dotL5 + 
                   tmp2*E1dq[4]*tmp1_dotL5 + 
                                           + tmp3*E2dq[4]*tmp1_dotL5

  sat_jacR[4,1] += tmp1*tmp2*E1dq4_dotR1   + tmp1*E1dq[4]*tmp2_dotR1 + 
                   tmp2*E1dq[4]*tmp1_dotR1 + 
                   tmp1*tmp3*E2dq4_dotR1   + tmp3*E2dq[4]*tmp1_dotR1

  sat_jacR[4,2] += tmp1*tmp2*E1dq4_dotR2   + tmp1*E1dq[4]*tmp2_dotR2 + 
                   tmp2*E1dq[4]*tmp1_dotR2 + 
                   tmp1*tmp3*E2dq4_dotR2   + tmp3*E2dq[4]*tmp1_dotR2

  sat_jacR[4,3] += tmp1*tmp2*E1dq4_dotR3   + tmp1*E1dq[4]*tmp2_dotR3 + 
                   tmp2*E1dq[4]*tmp1_dotR3 + 
                   tmp1*tmp3*E2dq4_dotR3   + tmp3*E2dq[4]*tmp1_dotR3

  sat_jacR[4,4] += tmp1*tmp2*E1dq4_dotR4   + tmp1*E1dq[4]*tmp2_dotR4 + 
                   tmp2*E1dq[4]*tmp1_dotR4 + 
                   tmp1*tmp3*E2dq4_dotR4   + tmp3*E2dq[4]*tmp1_dotR4

  sat_jacR[4,5] += tmp1*tmp2*E1dq4_dotR5   + tmp1*E1dq[4]*tmp2_dotR5 + 
                   tmp2*E1dq[4]*tmp1_dotR5 + 
                                           + tmp3*E2dq[4]*tmp1_dotR5


  sat_jacL[5,1] += tmp1*tmp2*E1dq5_dotL1   + tmp1*E1dq[5]*tmp2_dotL1 + 
                   tmp2*E1dq[5]*tmp1_dotL1 + 
                   tmp1*tmp3*E2dq5_dotL1   + tmp3*E2dq[5]*tmp1_dotL1

  sat_jacL[5,2] += tmp1*tmp2*E1dq5_dotL2   + tmp1*E1dq[5]*tmp2_dotL2 + 
                   tmp2*E1dq[5]*tmp1_dotL2 + 
                   tmp1*tmp3*E2dq5_dotL2   + tmp3*E2dq[5]*tmp1_dotL2

  sat_jacL[5,3] += tmp1*tmp2*E1dq5_dotL3   + tmp1*E1dq[5]*tmp2_dotL3 + 
                   tmp2*E1dq[5]*tmp1_dotL3 + 
                   tmp1*tmp3*E2dq5_dotL3   + tmp3*E2dq[5]*tmp1_dotL3

  sat_jacL[5,4] += tmp1*tmp2*E1dq5_dotL4   + tmp1*E1dq[5]*tmp2_dotL4 + 
                   tmp2*E1dq[5]*tmp1_dotL4 + 
                   tmp1*tmp3*E2dq5_dotL4   + tmp3*E2dq[5]*tmp1_dotL4

  sat_jacL[5,5] += tmp1*tmp2*E1dq5_dotL5   + tmp1*E1dq[5]*tmp2_dotL5 + 
                   tmp2*E1dq[5]*tmp1_dotL5 + 
                                           + tmp3*E2dq[5]*tmp1_dotL5

  sat_jacR[5,1] += tmp1*tmp2*E1dq5_dotR1   + tmp1*E1dq[5]*tmp2_dotR1 + 
                   tmp2*E1dq[5]*tmp1_dotR1 + 
                   tmp1*tmp3*E2dq5_dotR1   + tmp3*E2dq[5]*tmp1_dotR1

  sat_jacR[5,2] += tmp1*tmp2*E1dq5_dotR2   + tmp1*E1dq[5]*tmp2_dotR2 + 
                   tmp2*E1dq[5]*tmp1_dotR2 + 
                   tmp1*tmp3*E2dq5_dotR2   + tmp3*E2dq[5]*tmp1_dotR2

  sat_jacR[5,3] += tmp1*tmp2*E1dq5_dotR3   + tmp1*E1dq[5]*tmp2_dotR3 + 
                   tmp2*E1dq[5]*tmp1_dotR3 + 
                   tmp1*tmp3*E2dq5_dotR3   + tmp3*E2dq[5]*tmp1_dotR3

  sat_jacR[5,4] += tmp1*tmp2*E1dq5_dotR4   + tmp1*E1dq[5]*tmp2_dotR4 + 
                   tmp2*E1dq[5]*tmp1_dotR4 + 
                   tmp1*tmp3*E2dq5_dotR4   + tmp3*E2dq[5]*tmp1_dotR4


  sat_jacR[5,5] += tmp1*tmp2*E1dq5_dotR5   + tmp1*E1dq[5]*tmp2_dotR5 + 
                   tmp2*E1dq[5]*tmp1_dotR5 + 
                                           + tmp3*E2dq[5]*tmp1_dotR5

  #-- get E3*dq
  E1dq[1] = -Un*dq1 + nx*dq2 + ny*dq3 + nz*dq4
  E1dq1_dotL1 = -Un - dq1*Un_dotL1
  E1dq1_dotL2 =     - dq1*Un_dotL2 + nx
  E1dq1_dotL3 =     - dq1*Un_dotL3      + ny
  E1dq1_dotL4 =     - dq1*Un_dotL4           + nz

  E1dq1_dotR1 = Un - dq1*Un_dotR1
  E1dq1_dotR2 =    - dq1*Un_dotR2 - nx
  E1dq1_dotR3 =    - dq1*Un_dotR3      - ny
  E1dq1_dotR4 =    - dq1*Un_dotR4           - nz

  E1dq[2] = E1dq[1]*u
  E1dq2_dotL1 = u*E1dq1_dotL1 + E1dq[1]*u_dotL1
  E1dq2_dotL2 = u*E1dq1_dotL2 + E1dq[1]*u_dotL2
  E1dq2_dotL3 = u*E1dq1_dotL3
  E1dq2_dotL4 = u*E1dq1_dotL4

  E1dq2_dotR1 = u*E1dq1_dotR1 + E1dq[1]*u_dotR1
  E1dq2_dotR2 = u*E1dq1_dotR2 + E1dq[1]*u_dotR2
  E1dq2_dotR3 = u*E1dq1_dotR3
  E1dq2_dotR4 = u*E1dq1_dotR4


  E1dq[3] = E1dq[1]*v
  E1dq3_dotL1 = v*E1dq1_dotL1 + E1dq[1]*v_dotL1
  E1dq3_dotL2 = v*E1dq1_dotL2
  E1dq3_dotL3 = v*E1dq1_dotL3 + E1dq[1]*v_dotL3
  E1dq3_dotL4 = v*E1dq1_dotL4

  E1dq3_dotR1 = v*E1dq1_dotR1 + E1dq[1]*v_dotR1
  E1dq3_dotR2 = v*E1dq1_dotR2
  E1dq3_dotR3 = v*E1dq1_dotR3 + E1dq[1]*v_dotR3
  E1dq3_dotR4 = v*E1dq1_dotR4


  E1dq[4] = E1dq[1]*w
  E1dq4_dotL1 = w*E1dq1_dotL1 + E1dq[1]*w_dotL1
  E1dq4_dotL2 = w*E1dq1_dotL2
  E1dq4_dotL3 = w*E1dq1_dotL3
  E1dq4_dotL4 = w*E1dq1_dotL4 + E1dq[1]*w_dotL4

  E1dq4_dotR1 = w*E1dq1_dotR1 + E1dq[1]*w_dotR1
  E1dq4_dotR2 = w*E1dq1_dotR2
  E1dq4_dotR3 = w*E1dq1_dotR3
  E1dq4_dotR4 = w*E1dq1_dotR4 + E1dq[1]*w_dotR4

  E1dq[5] = E1dq[1]*H
  E1dq5_dotL1 = H*E1dq1_dotL1 + E1dq[1]*H_dotL1
  E1dq5_dotL2 = H*E1dq1_dotL2 + E1dq[1]*H_dotL2
  E1dq5_dotL3 = H*E1dq1_dotL3 + E1dq[1]*H_dotL3
  E1dq5_dotL4 = H*E1dq1_dotL4 + E1dq[1]*H_dotL4
  E1dq5_dotL5 =               + E1dq[1]*H_dotL5

  E1dq5_dotR1 = H*E1dq1_dotR1 + E1dq[1]*H_dotR1
  E1dq5_dotR2 = H*E1dq1_dotR2 + E1dq[1]*H_dotR2
  E1dq5_dotR3 = H*E1dq1_dotR3 + E1dq[1]*H_dotR3
  E1dq5_dotR4 = H*E1dq1_dotR4 + E1dq[1]*H_dotR4
  E1dq5_dotR5 =               + E1dq[1]*H_dotR5


  #-- get E4*dq
  t1 = phi*dq1 - u*dq2 - v*dq3 - w*dq4 + dq5
  t1_dotL1 = phi + dq1*phi_dotL1     - dq2*u_dotL1     - dq3*v_dotL1     - dq4*w_dotL1
  t1_dotL2 =       dq1*phi_dotL2 - u - dq2*u_dotL2
  t1_dotL3 =       dq1*phi_dotL3                   - v - dq3*v_dotL3
  t1_dotL4 =       dq1*phi_dotL4                                     - w - dq4*w_dotL4
  t1_dotL5 = 1

  t1_dotR1 = -phi + dq1*phi_dotR1     - dq2*u_dotR1     - dq3*v_dotR1    - dq4*w_dotR1
  t1_dotR2 =        dq1*phi_dotR2 + u - dq2*u_dotR2
  t1_dotR3 =        dq1*phi_dotR3                   + v - dq3*v_dotR3
  t1_dotR4 =        dq1*phi_dotR4                                     + w - dq4*w_dotR4 
  t1_dotR5 = -1


  E2dq[1] = 0.0
  E2dq[2] = t1*nx
  E2dq2_dotL1 = nx*t1_dotL1
  E2dq2_dotL2 = nx*t1_dotL2
  E2dq2_dotL3 = nx*t1_dotL3
  E2dq2_dotL4 = nx*t1_dotL4
  E2dq2_dotL5 = nx*t1_dotL5

  E2dq2_dotR1 = nx*t1_dotR1
  E2dq2_dotR2 = nx*t1_dotR2
  E2dq2_dotR3 = nx*t1_dotR3
  E2dq2_dotR4 = nx*t1_dotR4
  E2dq2_dotR5 = nx*t1_dotR5


  E2dq[3] = t1*ny
  E2dq3_dotL1 = ny*t1_dotL1
  E2dq3_dotL2 = ny*t1_dotL2
  E2dq3_dotL3 = ny*t1_dotL3
  E2dq3_dotL4 = ny*t1_dotL4
  E2dq3_dotL5 = ny*t1_dotL5

  E2dq3_dotR1 = ny*t1_dotR1
  E2dq3_dotR2 = ny*t1_dotR2
  E2dq3_dotR3 = ny*t1_dotR3
  E2dq3_dotR4 = ny*t1_dotR4
  E2dq3_dotR5 = ny*t1_dotR5

  E2dq[4] = t1*nz
  E2dq4_dotL1 = nz*t1_dotL1
  E2dq4_dotL2 = nz*t1_dotL2
  E2dq4_dotL3 = nz*t1_dotL3
  E2dq4_dotL4 = nz*t1_dotL4
  E2dq4_dotL5 = nz*t1_dotL5

  E2dq4_dotR1 = nz*t1_dotR1
  E2dq4_dotR2 = nz*t1_dotR2
  E2dq4_dotR3 = nz*t1_dotR3
  E2dq4_dotR4 = nz*t1_dotR4
  E2dq4_dotR5 = nz*t1_dotR5

  E2dq[5] = t1*Un
  E2dq5_dotL1 = Un*t1_dotL1 + t1*Un_dotL1
  E2dq5_dotL2 = Un*t1_dotL2 + t1*Un_dotL2
  E2dq5_dotL3 = Un*t1_dotL3 + t1*Un_dotL3
  E2dq5_dotL4 = Un*t1_dotL4 + t1*Un_dotL4
  E2dq5_dotL5 = Un*t1_dotL5

  E2dq5_dotR1 = Un*t1_dotR1 + t1*Un_dotR1
  E2dq5_dotR2 = Un*t1_dotR2 + t1*Un_dotR2
  E2dq5_dotR3 = Un*t1_dotR3 + t1*Un_dotR3
  E2dq5_dotR4 = Un*t1_dotR4 + t1*Un_dotR4
  E2dq5_dotR5 = Un*t1_dotR5

  #-- add to sat
  t1 = 1/(dA*a)
  t2 = 1/a
  tmp1 = 0.5*(lambda1 - lambda2)*t1

  tmp1_dotL1 = 0.5 * t1 * ( (lambda1_dotL1 - lambda2_dotL1) - (lambda1 - lambda2) * t2 * a_dotL1)
  tmp1_dotL2 = 0.5 * t1 * ( (lambda1_dotL2 - lambda2_dotL2) - (lambda1 - lambda2) * t2 * a_dotL2)
  tmp1_dotL3 = 0.5 * t1 * ( (lambda1_dotL3 - lambda2_dotL3) - (lambda1 - lambda2) * t2 * a_dotL3)
  tmp1_dotL4 = 0.5 * t1 * ( (lambda1_dotL4 - lambda2_dotL4) - (lambda1 - lambda2) * t2 * a_dotL4)
  tmp1_dotL5 = 0.5 * t1 * ( (lambda1_dotL5 - lambda2_dotL5) - (lambda1 - lambda2) * t2 * a_dotL5)

  tmp1_dotR1 = 0.5 * t1 * ( (lambda1_dotR1 - lambda2_dotR1) - (lambda1 - lambda2) * t2 * a_dotR1)
  tmp1_dotR2 = 0.5 * t1 * ( (lambda1_dotR2 - lambda2_dotR2) - (lambda1 - lambda2) * t2 * a_dotR2)
  tmp1_dotR3 = 0.5 * t1 * ( (lambda1_dotR3 - lambda2_dotR3) - (lambda1 - lambda2) * t2 * a_dotR3)
  tmp1_dotR4 = 0.5 * t1 * ( (lambda1_dotR4 - lambda2_dotR4) - (lambda1 - lambda2) * t2 * a_dotR4)
  tmp1_dotR5 = 0.5 * t1 * ( (lambda1_dotR5 - lambda2_dotR5) - (lambda1 - lambda2) * t2 * a_dotR5)


  #for i=1:length(sat)
  #  sat[i] = sat[i] + tmp1*(E1dq[i] + gami*E2dq[i])
  #end

  t1 = E1dq[1] + gami*E2dq[1]
  sat_jacL[1, 1] += tmp1*E1dq1_dotL1 + t1*tmp1_dotL1
  sat_jacL[1, 2] += tmp1*E1dq1_dotL2 + t1*tmp1_dotL2
  sat_jacL[1, 3] += tmp1*E1dq1_dotL3 + t1*tmp1_dotL3
  sat_jacL[1, 4] += tmp1*E1dq1_dotL4 + t1*tmp1_dotL4
  sat_jacL[1, 5] +=                    t1*tmp1_dotL5

  sat_jacR[1, 1] += tmp1*E1dq1_dotR1 + t1*tmp1_dotR1
  sat_jacR[1, 2] += tmp1*E1dq1_dotR2 + t1*tmp1_dotR2
  sat_jacR[1, 3] += tmp1*E1dq1_dotR3 + t1*tmp1_dotR3
  sat_jacR[1, 4] += tmp1*E1dq1_dotR4 + t1*tmp1_dotR4
  sat_jacR[1, 5] +=                    t1*tmp1_dotR5

  t1 = E1dq[2] + gami*E2dq[2]
  sat_jacL[2, 1] += tmp1*(E1dq2_dotL1 + gami*E2dq2_dotL1) + t1*tmp1_dotL1
  sat_jacL[2, 2] += tmp1*(E1dq2_dotL2 + gami*E2dq2_dotL2) + t1*tmp1_dotL2
  sat_jacL[2, 3] += tmp1*(E1dq2_dotL3 + gami*E2dq2_dotL3) + t1*tmp1_dotL3
  sat_jacL[2, 4] += tmp1*(E1dq2_dotL4 + gami*E2dq2_dotL4) + t1*tmp1_dotL4
  sat_jacL[2, 5] += tmp1*(            + gami*E2dq2_dotL5) + t1*tmp1_dotL5

  sat_jacR[2, 1] += tmp1*(E1dq2_dotR1 + gami*E2dq2_dotR1) + t1*tmp1_dotR1
  sat_jacR[2, 2] += tmp1*(E1dq2_dotR2 + gami*E2dq2_dotR2) + t1*tmp1_dotR2
  sat_jacR[2, 3] += tmp1*(E1dq2_dotR3 + gami*E2dq2_dotR3) + t1*tmp1_dotR3
  sat_jacR[2, 4] += tmp1*(E1dq2_dotR4 + gami*E2dq2_dotR4) + t1*tmp1_dotR4
  sat_jacR[2, 5] += tmp1*(            + gami*E2dq2_dotR5) + t1*tmp1_dotR5

  t1 = E1dq[3] + gami*E2dq[3]
  sat_jacL[3, 1] += tmp1*(E1dq3_dotL1 + gami*E2dq3_dotL1) + t1*tmp1_dotL1
  sat_jacL[3, 2] += tmp1*(E1dq3_dotL2 + gami*E2dq3_dotL2) + t1*tmp1_dotL2
  sat_jacL[3, 3] += tmp1*(E1dq3_dotL3 + gami*E2dq3_dotL3) + t1*tmp1_dotL3
  sat_jacL[3, 4] += tmp1*(E1dq3_dotL4 + gami*E2dq3_dotL4) + t1*tmp1_dotL4
  sat_jacL[3, 5] += tmp1*(            + gami*E2dq3_dotL5) + t1*tmp1_dotL5

  sat_jacR[3, 1] += tmp1*(E1dq3_dotR1 + gami*E2dq3_dotR1) + t1*tmp1_dotR1
  sat_jacR[3, 2] += tmp1*(E1dq3_dotR2 + gami*E2dq3_dotR2) + t1*tmp1_dotR2
  sat_jacR[3, 3] += tmp1*(E1dq3_dotR3 + gami*E2dq3_dotR3) + t1*tmp1_dotR3
  sat_jacR[3, 4] += tmp1*(E1dq3_dotR4 + gami*E2dq3_dotR4) + t1*tmp1_dotR4
  sat_jacR[3, 5] += tmp1*(            + gami*E2dq3_dotR5) + t1*tmp1_dotR5

  t1 = E1dq[4] + gami*E2dq[4]
  sat_jacL[4, 1] += tmp1*(E1dq4_dotL1 + gami*E2dq4_dotL1) + t1*tmp1_dotL1
  sat_jacL[4, 2] += tmp1*(E1dq4_dotL2 + gami*E2dq4_dotL2) + t1*tmp1_dotL2
  sat_jacL[4, 3] += tmp1*(E1dq4_dotL3 + gami*E2dq4_dotL3) + t1*tmp1_dotL3
  sat_jacL[4, 4] += tmp1*(E1dq4_dotL4 + gami*E2dq4_dotL4) + t1*tmp1_dotL4
  sat_jacL[4, 5] += tmp1*(            + gami*E2dq4_dotL5) + t1*tmp1_dotL5

  sat_jacR[4, 1] += tmp1*(E1dq4_dotR1 + gami*E2dq4_dotR1) + t1*tmp1_dotR1
  sat_jacR[4, 2] += tmp1*(E1dq4_dotR2 + gami*E2dq4_dotR2) + t1*tmp1_dotR2
  sat_jacR[4, 3] += tmp1*(E1dq4_dotR3 + gami*E2dq4_dotR3) + t1*tmp1_dotR3
  sat_jacR[4, 4] += tmp1*(E1dq4_dotR4 + gami*E2dq4_dotR4) + t1*tmp1_dotR4
  sat_jacR[4, 5] += tmp1*(            + gami*E2dq4_dotR5) + t1*tmp1_dotR5


  t1 = E1dq[5] + gami*E2dq[5]
  sat_jacL[5, 1] += tmp1*(E1dq5_dotL1 + gami*E2dq5_dotL1) + t1*tmp1_dotL1
  sat_jacL[5, 2] += tmp1*(E1dq5_dotL2 + gami*E2dq5_dotL2) + t1*tmp1_dotL2
  sat_jacL[5, 3] += tmp1*(E1dq5_dotL3 + gami*E2dq5_dotL3) + t1*tmp1_dotL3
  sat_jacL[5, 4] += tmp1*(E1dq5_dotL4 + gami*E2dq5_dotL4) + t1*tmp1_dotL4
  sat_jacL[5, 5] += tmp1*(E1dq5_dotL5 + gami*E2dq5_dotL5) + t1*tmp1_dotL5

  sat_jacR[5, 1] += tmp1*(E1dq5_dotR1 + gami*E2dq5_dotR1) + t1*tmp1_dotR1
  sat_jacR[5, 2] += tmp1*(E1dq5_dotR2 + gami*E2dq5_dotR2) + t1*tmp1_dotR2
  sat_jacR[5, 3] += tmp1*(E1dq5_dotR3 + gami*E2dq5_dotR3) + t1*tmp1_dotR3
  sat_jacR[5, 4] += tmp1*(E1dq5_dotR4 + gami*E2dq5_dotR4) + t1*tmp1_dotR4
  sat_jacR[5, 5] += tmp1*(E1dq5_dotR5 + gami*E2dq5_dotR5) + t1*tmp1_dotR5


  return nothing
end  # End function calcSAT


"""
  Reverse mode wrt q and qg of [`RoeSolver`](@ref)
"""
function RoeSolver_revq(params::ParamType{3, :conservative},
                   q::AbstractArray{Tsol,1},
                   q_bar::AbstractArray{Tsol, 1},
                   qg::AbstractArray{Tsol, 1},
                   qg_bar::AbstractArray{Tsol, 1},
                   aux_vars::AbstractArray{Tres, 1},
                   nrm::AbstractArray{Tmsh,1},
                   flux_bar::AbstractArray{Tres, 1}) where {Tmsh, Tsol, Tres}

  # SAT terms are used for ensuring consistency with the physical problem. Its
  # similar to upwinding which adds dissipation to the problem. SATs on the
  # boundary can be thought of as having two overlapping nodes and because of
  # the discontinuous nature of SBP adds some dissipation.

  data = params.roefluxdata
  @unpack data dq sat roe_vars euler_flux v_vals nrm2

  # Declaring constants
  d1_0 = 1.0
  d0_0 = 0.0
  d0_5 = 0.5
  tau = 1.0
  gamma = params.gamma
  gami = params.gamma_1
  sat_fac = 1  # multiplier for SAT term

  # Begin main execution
  nx = nrm[1]
  ny = nrm[2]

  # Compute the Roe Averaged states
  # The left state of Roe are the actual solution variables
  fac = d1_0/q[1]
  uL = q[2]*fac; vL = q[3]*fac; wL = q[4]*fac
  phi = d0_5*(uL*uL + vL*vL + wL*wL)
  HL = gamma*q[5]*fac - gami*phi # Total enthalpy, H = e + 0.5*(u^2 + v^2) + p/rho,
                                 # where e is the internal energy per unit mass

  # The right side of the Roe solver comprises the boundary conditions
  fac = d1_0/qg[1]
  uR = qg[2]*fac; vR = qg[3]*fac; wR = qg[4]*fac
  phi = d0_5*(uR*uR + vR*vR + wR*wR)
  HR = gamma*qg[5]*fac - gami*phi # Total Enthalpy

  # Averaged states
  sqL = sqrt(q[1])
  sqR = sqrt(qg[1])
  fac = d1_0/(sqL + sqR)
  u = (sqL*uL + sqR*uR)*fac
  v = (sqL*vL + sqR*vR)*fac
  w = (sqL*wL + sqR*wR)*fac
  H = (sqL*HL + sqR*HR)*fac

  for i=1:length(dq)
    dq[i] = q[i] - qg[i]
  end

  roe_vars[1] = u
  roe_vars[2] = v
  roe_vars[3] = w
  roe_vars[4] = H

#  calcSAT(params, roe_vars, dq, nrm, sat)

  # calculate Euler flux in wall normal directiona
  # because edge numbering is rather arbitary, any memory access is likely to
  # be a cache miss, so we recalculate the Euler flux
#  nrm2[1] = nx   # why are we assigning to nrm2?
#  nrm2[2] = ny

#  convertFromNaturalToWorkingVars(params, q, v_vals)
#  calcEulerFlux(params, v_vals, aux_vars, nrm2, euler_flux)

#  for i=1:5  # ArrayViews does not support flux[:] = .
#    flux[i] = (sat_fac*sat[i] + euler_flux[i])
#  end

  #----------------------------------------------------------------------------
  # reverse sweep

  numDofPerNode = length(q)

  @unpack data sat_bar v_vals_bar roe_vars_bar dq_bar

  fill!(sat_bar, 0); fill!(roe_vars_bar, 0); fill!(dq_bar, 0)
  fill!(v_vals_bar, 0)
  for i=1:5
    sat_bar[i] = sat_fac*flux_bar[i]
    v_vals[i] = q[i] 
  end
  euler_flux_bar = flux_bar

  calcEulerFlux_revq(params, v_vals, v_vals_bar, aux_vars, nrm, euler_flux_bar)


  #TODO: this is where convertFromNatualToWorkingVars_rev should be
  for i=1:numDofPerNode
    q_bar[i] += v_vals_bar[i]
  end

  calcSAT_revq(params, roe_vars, roe_vars_bar, dq, dq_bar, nrm, sat_bar)


  for i=1:length(dq)
    q_bar[i]  += dq_bar[i]
    qg_bar[i] -= dq_bar[i]
  end

  u_bar = roe_vars_bar[1]
  v_bar = roe_vars_bar[2]
  w_bar = roe_vars_bar[3]
  H_bar = roe_vars_bar[4]

  # averaged states
  sqL_bar = HL*fac*H_bar; sqR_bar = HR*fac*H_bar
  HL_bar  = sqL*fac*H_bar; HR_bar = sqR*fac*H_bar
  fac_bar = H_bar*(sqL*HL + sqR*HR)

  sqL_bar += wL*fac*w_bar; sqR_bar += wR*fac*w_bar
  wL_bar   = sqL*fac*w_bar; wR_bar  = sqR*fac*w_bar
  fac_bar += (sqL*wL + sqR*wR)*w_bar

  sqL_bar += vL*fac*v_bar; sqR_bar += vR*fac*v_bar
  vL_bar   = sqL*fac*v_bar; vR_bar  = sqR*fac*v_bar
  fac_bar += (sqL*vL + sqR*vR)*v_bar

  sqL_bar += uL*fac*u_bar; sqR_bar += uR*fac*u_bar
  uL_bar   = sqL*fac*u_bar; uR_bar  = sqR*fac*u_bar
  fac_bar += (sqL*uL + sqR*uR)*u_bar

  sqL_bar += -fac*fac*fac_bar; sqR_bar += -fac*fac*fac_bar
  q_bar[1] += 0.5/sqL*sqL_bar; qg_bar[1] += 0.5/sqR*sqR_bar


  # right state
  fac = d1_0/qg[1]
  fac_bar = zero(fac_bar)  # this is a different variable in single assignment
  qg_bar[5] += gamma*fac*HR_bar
  fac_bar += gamma*qg[5]*HR_bar
  phi_bar = -gami*HR_bar

  uR_bar += uR*phi_bar
  vR_bar += vR*phi_bar
  wR_bar += wR*phi_bar

  qg_bar[2] += fac*uR_bar; qg_bar[3] += fac*vR_bar; qg_bar[4] += fac*wR_bar
  fac_bar += qg[2]*uR_bar; fac_bar += qg[3]*vR_bar; fac_bar += qg[4]*wR_bar

  qg_bar[1] += -fac*fac*fac_bar

  # left state
  fac = d1_0/q[1]
  fac_bar = zero(fac_bar); phi_bar = zero(phi_bar)
  q_bar[5] += gamma*fac*HL_bar
  fac_bar += gamma*q[5]*HL_bar
  phi_bar += -gami*HL_bar

  uL_bar += uL*phi_bar
  vL_bar += vL*phi_bar
  wL_bar += wL*phi_bar

  q_bar[2] += fac*uL_bar; q_bar[3] += fac*vL_bar; q_bar[4] += fac*wL_bar
  fac_bar += q[2]*uL_bar; fac_bar += q[3]*vL_bar; fac_bar += q[4]*wL_bar
  q_bar[1] += -fac*fac*fac_bar


  return nothing
end

"""
  3D method
"""
function calcSAT_revq(params::ParamType{3},
                 roe_vars::AbstractArray{Tsol, 1},
                 roe_vars_bar::AbstractArray{Tres, 1},
                 dq::AbstractArray{Tsol,1},
                 dq_bar::AbstractArray{Tres, 1},
                 nrm::AbstractArray{Tmsh,1},
                 sat_bar::AbstractArray{Tsol,1}) where {Tmsh, Tsol, Tres}
# roe_vars = [u, v, H] at Roe average 

  data = params.calcsatdata
  @unpack data E1dq E2dq E3dq E4dq

  numDofPerNode = length(dq)

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

  Un = u*nx + v*ny + w*nz# Normal Velocity

  phi = 0.5*(u*u + v*v + w*w)

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

  lambda1_orig = lambda1  # save these for reverse sweep
  lambda2_orig = lambda2
  lambda3_orig = lambda3

  lambda1 = 0.5*(tau*max(absvalue(lambda1),sat_Vn *rhoA) - lambda1)
  lambda2 = 0.5*(tau*max(absvalue(lambda2),sat_Vn *rhoA) - lambda2)
  lambda3 = 0.5*(tau*max(absvalue(lambda3),sat_Vl *rhoA) - lambda3)


  dq1 = dq[1]
  dq2 = dq[2]
  dq3 = dq[3]
  dq4 = dq[4]
  dq5 = dq[5]

  
  #=
  sat[1] = lambda3*dq1
  sat[2] = lambda3*dq2
  sat[3] = lambda3*dq3
  sat[4] = lambda3*dq3
  sat[5] = lambda3*dq5
  =#

  #-- get E1*dq
  t1 = phi*dq1 - u*dq2 - v*dq3 - w*dq4 + dq5
  E1dq[1] = t1
  E1dq[2] = t1*u
  E1dq[3] = t1*v
  E1dq[4] = t1*w
  E1dq[5] = t1*H


  #-- get E2*dq
  t2 = -Un*dq1 + nx*dq2 + ny*dq3 + nz*dq4
  E2dq[1] = 0.0
  E2dq[2] = t2*nx
  E2dq[3] = t2*ny
  E2dq[4] = t2*nz
  E2dq[5] = t2*Un

  #-- add to sat
  tmp1 = 0.5*(lambda1 + lambda2) - lambda3
  tmp2 = gami/(a*a)
  tmp3 = 1.0/(dA*dA)

#  for i=1:length(sat)
#    sat[i] = sat[i] + tmp1*(tmp2*E1dq[i] + tmp3*E2dq[i])
#  end

  #-- get E3*dq
  t3 = -Un*dq1 + nx*dq2 + ny*dq3 + nz*dq4

  E3dq[1] = t3
  E3dq[2] = t3*u
  E3dq[3] = t3*v
  E3dq[4] = t3*w
  E3dq[5] = t3*H

  #-- get E4*dq
  t4 = phi*dq1 - u*dq2 - v*dq3 + - w*dq4 + dq5
  E4dq[1] = 0.0
  E4dq[2] = t4*nx
  E4dq[3] = t4*ny
  E4dq[4] = t4*nz
  E4dq[5] = t4*Un

  #-- add to sat
  tmp4 = 0.5*(lambda1 - lambda2)/(dA*a)
#  for i=1:length(sat)
#    sat[i] = sat[i] + tmp4*(E3dq[i] + gami*E4dq[i])
#  end

  #----------------------------------------------------------------------------
  # reverse sweep
  @unpack data E1dq_bar E2dq_bar E3dq_bar E4dq_bar
  t1_bar = zero(Tres)
  t2_bar = zero(Tres)
  t3_bar = zero(Tres)
  t4_bar = zero(Tres)
  tmp1_bar = zero(Tres)
  tmp2_bar = zero(Tres)
  tmp3_bar = zero(Tres)
  tmp4_bar = zero(Tres)
  Un_bar = zero(Tres)
  phi_bar = zero(Tres)
  dq1_bar = zero(Tres)
  dq2_bar = zero(Tres)
  dq3_bar = zero(Tres)
  dq4_bar = zero(Tres)
  dq5_bar = zero(Tres)
  lambda1_bar = zero(Tres)
  lambda2_bar = zero(Tres)
  lambda3_bar = zero(Tres)
  rhoA_bar = zero(Tres)
  a_bar = zero(Tres)
  u_bar = zero(Tres)
  v_bar = zero(Tres)
  w_bar = zero(Tres)
  H_bar = zero(Tres)

  fill!(E1dq_bar, 0)
  fill!(E2dq_bar, 0)
  fill!(E3dq_bar, 0)
  fill!(E4dq_bar, 0)

  tmp_bar = zero(Tres)

  for i=1:length(sat_bar)
    tmp4_bar += (E3dq[i] + gami*E4dq[i])*sat_bar[i]
    E3dq_bar[i] += tmp4*sat_bar[i]
    E4dq_bar[i] += tmp4*gami*sat_bar[i]
  end

  lambda1_bar +=  0.5*tmp4_bar/(dA*a)
  lambda2_bar += -0.5*tmp4_bar/(dA*a)
  a_bar       += -tmp4*tmp4_bar/a


  # E4*dq
  
  t4_bar += E4dq_bar[5]*Un
  Un_bar += t4*E4dq_bar[5]

  t4_bar += E4dq_bar[4]*nz

  t4_bar += E4dq_bar[3]*ny

  t4_bar += E4dq_bar[2]*nx

  phi_bar +=  dq1*t4_bar
  dq1_bar +=  phi*t4_bar
  u_bar   += -dq2*t4_bar
  dq2_bar += -u*t4_bar
  v_bar   += -dq3*t4_bar
  dq3_bar += -v*t4_bar
  w_bar   += -dq4*t4_bar
  dq4_bar += -w*t4_bar
  dq5_bar +=  t4_bar


  # E3*dq
  t3_bar += E3dq_bar[5]*H
  H_bar  += t3*E3dq_bar[5]

  t3_bar += E3dq_bar[4]*w
  w_bar  += E3dq_bar[4]*t3

  t3_bar += E3dq_bar[3]*v
  v_bar  += E3dq_bar[3]*t3

  t3_bar += E3dq_bar[2]*u
  u_bar  += E3dq_bar[2]*t3
  
  t3_bar += E3dq_bar[1]

  Un_bar  += -dq1*t3_bar
  dq1_bar += -Un*t3_bar
  dq2_bar +=  nx*t3_bar
  dq3_bar +=  ny*t3_bar
  dq4_bar +=  nz*t3_bar

  for i=1:length(sat_bar)
    tmp1_bar    += sat_bar[i]*(tmp2*E1dq[i] + tmp3*E2dq[i])
    tmp2_bar    += tmp1*E1dq[i]*sat_bar[i]
    E1dq_bar[i] += tmp1*tmp2*sat_bar[i]
#    tmp3_bar    += tmp1*E2dq[i]*sat_bar[i]
    E2dq_bar[i] += tmp1*tmp3*sat_bar[i]
  end

  a_bar += (-2*tmp2/a)*tmp2_bar
  lambda1_bar += 0.5*tmp1_bar
  lambda2_bar += 0.5*tmp1_bar
  lambda3_bar += -tmp1_bar


  # E2*dq

  t2_bar += Un*E2dq_bar[5]
  Un_bar += t2*E2dq_bar[5]

  t2_bar += nz*E2dq_bar[4]

  t2_bar += ny*E2dq_bar[3]

  t2_bar += nx*E2dq_bar[2]

  Un_bar  += -dq1*t2_bar
  dq1_bar += -Un*t2_bar
  dq2_bar +=  nx*t2_bar
  dq3_bar +=  ny*t2_bar
  dq4_bar += nz*t2_bar



  # E1*dq
  t1_bar += E1dq_bar[5]*H
  H_bar  += E1dq_bar[5]*t1

  t1_bar += E1dq_bar[4]*w
  w_bar  += E1dq_bar[4]*t1

  t1_bar += E1dq_bar[3]*v
  v_bar  += E1dq_bar[3]*t1

  t1_bar += E1dq_bar[2]*u
  u_bar  += E1dq_bar[2]*t1

  t1_bar += E1dq_bar[1]

  phi_bar +=  dq1*t1_bar
  dq1_bar +=  phi*t1_bar
  u_bar   += -dq2*t1_bar
  dq2_bar += -u*t1_bar
  v_bar   += -dq3*t1_bar
  dq3_bar += -v*t1_bar
  w_bar   += -dq4*t1_bar
  dq4_bar += -w*t1_bar
  dq5_bar +=  t1_bar



  # sat diagonal term
  lambda3_bar += dq1*sat_bar[1]
  dq1_bar     += lambda3*sat_bar[1]

  lambda3_bar += dq2*sat_bar[2]
  dq2_bar     += lambda3*sat_bar[2]

  lambda3_bar += dq3*sat_bar[3]
  dq3_bar     += lambda3*sat_bar[3]

  lambda3_bar += dq4*sat_bar[4]
  dq4_bar     += lambda3*sat_bar[4]

  lambda3_bar += dq5*sat_bar[5]
  dq5_bar     += lambda3*sat_bar[5]



  # re-pack dq
  dq_bar[1] += dq1_bar
  dq_bar[2] += dq2_bar
  dq_bar[3] += dq3_bar
  dq_bar[4] += dq4_bar
  dq_bar[5] += dq5_bar


  # restore lambdas to original values from forward sweep
  # lambda_bar variables will be overwritten below
  lambda1 = lambda1_orig
  lambda2 = lambda2_orig
  lambda3 = lambda3_orig

  # entropy fix
  alambda1 = absvalue(lambda1)
  alambda2 = absvalue(lambda2)
  alambda3 = absvalue(lambda3)
  # compute the alambda_bar variables inside the conditionals because they
  # are not always needed

  if alambda1 > sat_Vn*rhoA
    # overwrite lambda here because in the primal code it gets re-assigned
    alambda1_bar = lambda1 > 0 ? lambda1_bar : -lambda1_bar
    lambda1_bar  = 0.5*(tau*alambda1_bar - lambda1_bar)
  else
    rhoA_bar    += 0.5*tau*sat_Vn*lambda1_bar
    lambda1_bar = -0.5*lambda1_bar
  end

  if alambda2 > sat_Vn*rhoA
    alambda2_bar = lambda2 > 0 ? lambda2_bar : -lambda2_bar
    lambda2_bar  = 0.5*(tau*alambda2_bar - lambda2_bar)
  else
    rhoA_bar    += 0.5*tau*sat_Vn*lambda2_bar
    lambda2_bar = -0.5*lambda2_bar
  end

  if alambda3 > sat_Vl*rhoA
    alambda3_bar = lambda3 > 0 ? lambda3_bar : -lambda3_bar
    lambda3_bar  = 0.5*tau*(alambda3_bar - lambda3_bar)
  else
    rhoA_bar    += 0.5*tau*sat_Vl*lambda3_bar
    lambda3_bar = -0.5*lambda3_bar
  end

  # rhoA
  if Un > 0
    Un_bar += rhoA_bar
  else
    Un_bar += -rhoA_bar
  end
  a_bar += dA*rhoA_bar


  # lambda calculation
  Un_bar += lambda1_bar
  a_bar  += dA*lambda1_bar

  Un_bar += lambda2_bar
  a_bar  += -dA*lambda2_bar

  Un_bar += lambda3_bar

  # a calculation
  H_bar   += (0.5/a)*gami*a_bar
  phi_bar += -(0.5/a)*gami*a_bar

  # phi
  u_bar += u*phi_bar
  v_bar += v*phi_bar
  w_bar += w*phi_bar

  # Un
  u_bar += nx*Un_bar
  v_bar += ny*Un_bar
  w_bar += nz*Un_bar

  roe_vars_bar[1] += u_bar
  roe_vars_bar[2] += v_bar
  roe_vars_bar[3] += w_bar
  roe_vars_bar[4] += H_bar
  

  return nothing
end  # End function calcSAT




function RoeSolver_revm(params::ParamType{3},
      q::AbstractArray{Tsol,1}, qg::AbstractArray{Tsol, 1},
      aux_vars::AbstractArray{Tres, 1},
      nrm::AbstractArray{Tmsh,1},
      nrm_bar::AbstractArray{Tmsh, 1},
      flux_bar::AbstractArray{Tres, 1},
      ) where {Tmsh, Tsol, Tres}

  @unpack params.calcsatdata E1dq E2dq

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
  # calcEulerFlux_revm(params, v_vals, aux_vars, nrm2, nrm2_bar, euler_flux_bar)
  calcEulerFlux_revm(params, q, aux_vars, nrm, nrm_bar, euler_flux_bar)

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




"""
  Computes the Jacobian of the LF flux with respect to qL and qR
"""
function calcLFFlux_diff(
                      params::ParamType{Tdim, :conservative},
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tsol, 1},
                      dir::AbstractArray{Tmsh, 1},
                      _F_dotL::AbstractMatrix{Tres}, _F_dotR::AbstractMatrix{Tres}) where {Tmsh, Tsol, Tres, Tdim}

  lffluxdata = params.lffluxdata
  @unpack lffluxdata F_dotL F_dotR lambda_dotL lambda_dotR
  F_dotL = lffluxdata.F_dotL; F_dotR = lffluxdata.F_dotR
  fill!(F_dotL, 0); fill!(F_dotR, 0)

  numDofPerNode = length(qL)
#  lambda_dotL = zeros(Tres, numDofPerNode)
#  lambda_dotR = zeros(Tres, numDofPerNode)

#  lambda_dotL = params.lambda_dotL
#  lambda_dotR = params.lambda_dotR

  calcEulerFlux_diff(params, qL, aux_vars, dir, F_dotL)
  calcEulerFlux_diff(params, qR, aux_vars, dir, F_dotR)

  lambda_max = getLambdaMaxSimple_diff(params, qL, qR, dir, lambda_dotL, lambda_dotR)

  for j=1:numDofPerNode
    F_dotL[j, j] += lambda_max
    F_dotR[j, j] -= lambda_max

    for i=1:numDofPerNode
      F_dotL[i, j] -=  qR[i]*lambda_dotL[j] - qL[i]*lambda_dotL[j]
      F_dotR[i, j] -= -qL[i]*lambda_dotR[j] + qR[i]*lambda_dotR[j]
      F_dotL[i, j] *= 0.5
      F_dotR[i, j] *= 0.5

      # update output array
      _F_dotL[i, j] += F_dotL[i, j]
      _F_dotR[i, j] += F_dotR[i, j]
    end
  end

  return nothing
end

#------------------------------------------------------------------------------
# IR flux differentiated
# there are a total of 12 versions: (2d vs 3d) x (single direction vs multi
# direction) x (forward vector mode vs reverse mode for q vs reverse mode for
# reverse mode for nrm).



#------------------------------------------------------------------------------
# 2D, single direction

"""
  Forward vector mode of [`calcIRFlux`](@ref), 2D, single direction version

  **Inputs**

   * params: ParamType
   * qL: solution vector at left state
   * qR: solutoin vector at right state
   * aux_vars
   * nrm: normal vector at the current node, length `dim`

  **Inputs/Outputs**
  
   * FL_dot: jacobian of the flux with respect to `qL` (overwritten)
   * FR_dot: jacobian of the flux with respect to `qR`.(overwritten)
"""
function calcEulerFlux_IR_diff(params::ParamType{2, :conservative},
                   qL::AbstractArray{Tsol,1},
                   qR::AbstractArray{Tsol, 1},
                   aux_vars::AbstractArray{Tres, 1},
                   nrm::AbstractArray{Tmsh, 1},
                   FL_dot::AbstractArray{Tres, 2},
                   FR_dot::AbstractArray{Tres, 2}) where {Tmsh, Tsol, Tres}

  @debug1 begin
    @assert size(nrm, 1) == 2
    @assert length(qL) == length(qR)
    @assert length(qL) == size(FL_dot, 1)
    @assert length(qL) == size(FL_dot, 2)
    @assert size(FL_dot, 1) == size(FR_dot, 1)
    @assert size(FL_dot, 2) == size(FR_dot, 2)
  end


  data = params.irfluxdata
  @unpack data pL_dot pR_dot
  @unpack data z1L_dot z2L_dot z3L_dot z4L_dot z1R_dot z2R_dot z3R_dot z4R_dot

  gamma = params.gamma
  gamma_1 = params.gamma_1

  pL = calcPressure_diff(params, qL, pL_dot)
  pR = calcPressure_diff(params, qR, pR_dot)
  uL = qL[2]/qL[1]; uR = qR[2]/qR[1]
  vL = qL[3]/qL[1]; vR = qR[3]/qR[1]
  z1L = sqrt(qL[1]/pL); z1R = sqrt(qR[1]/pR)
  z2L = z1L*uL; z2R = z1R*uR
  z3L = z1L*vL; z3R = z1R*vR
  z4L = sqrt(qL[1]*pL); z4R = sqrt(qR[1]*pR)

  fastzero!(z1L_dot); fastzero!(z1R_dot)
  fastzero!(z2L_dot); fastzero!(z2R_dot)
  fastzero!(z3L_dot); fastzero!(z3R_dot)
  fastzero!(z4L_dot); fastzero!(z4R_dot)

  # differentiate with respect to q (not including chain rule terms for p)
  z1L_dot[1] = (0.5/z1L)*1/pL; z1R_dot[1] = (0.5/z1R)*1/pR

  z2L_dot[1] = -z2L/qL[1]; z2R_dot[1] = -z2R/qR[1]
  z2L_dot[2] =  z1L/qL[1]; z2R_dot[2] =  z1R/qR[1]

  z3L_dot[1] = -z3L/qL[1]; z3R_dot[1] = -z3R/qR[1]
  z3L_dot[3] =  z1L/qL[1]; z3R_dot[3] =  z1R/qR[1]

  z4L_dot[1] =  (0.5/z4L)*pL; z4R_dot[1] = (0.5/z4R)*pR

  # do the pressure/z1L related terms
  @simd for i=1:4
    z1L_dot[i] += (0.5/z1L)*(-qL[1]/(pL*pL))*pL_dot[i]
    z1R_dot[i] += (0.5/z1R)*(-qR[1]/(pR*pR))*pR_dot[i]

    z2L_dot[i] += uL*z1L_dot[i]
    z2R_dot[i] += uR*z1R_dot[i]

    z3L_dot[i] += vL*z1L_dot[i]
    z3R_dot[i] += vR*z1R_dot[i]

   
    z4L_dot[i] += (0.5/z4L)*qL[1]*pL_dot[i]
    z4R_dot[i] += (0.5/z4R)*qR[1]*pR_dot[i]
  end

  @unpack data avgdata z4avg_dotL z4avg_dotR z1avg_dotL z1avg_dotR
  @unpack data rho_hat_dotL rho_hat_dotR u_hat_dotL u_hat_dotR
  @unpack data v_hat_dotL v_hat_dotR p1_hat_dotL p1_hat_dotR
  @unpack data h_hat_dotL h_hat_dotR

  # z4avg_dotL/r, z1avg_dotL/r, rho_hat, u_hat, v_hat, p1_hat, p2_hat
  z4avg = logavg_diff(avgdata, z4L, z4L_dot, z4R, z4R_dot, z4avg_dotL, z4avg_dotR)
  z1avg = logavg_diff(avgdata, z1L, z1L_dot, z1R, z1R_dot, z1avg_dotL, z1avg_dotR)

  z1 = z1L + z1R
  z1inv = 1/z1
  rho_hat = 0.5*(z1L + z1R)*z4avg
  u_hat = (z2L + z2R)*z1inv
  v_hat = (z3L + z3R)*z1inv
  p1_hat = (z4L + z4R)*z1inv
  p2_hat = ((gamma + 1)/(2*gamma) )*z4avg/z1avg + ( gamma_1/(2*gamma) )*p1_hat
  h_hat = gamma*p2_hat/(rho_hat*gamma_1) + 0.5*(u_hat*u_hat + v_hat*v_hat)

  @simd for i=1:4
    rho_hat_dotL[i] = 0.5*(z4avg*z1L_dot[i] + z1*z4avg_dotL[i])
    rho_hat_dotR[i] = 0.5*(z4avg*z1R_dot[i] + z1*z4avg_dotR[i])

    u_hat_dotL[i] = z2L_dot[i]*z1inv - u_hat*z1inv*z1L_dot[i]
    u_hat_dotR[i] = z2R_dot[i]*z1inv - u_hat*z1inv*z1R_dot[i]

    v_hat_dotL[i] = z3L_dot[i]*z1inv - v_hat*z1inv*z1L_dot[i]
    v_hat_dotR[i] = z3R_dot[i]*z1inv - v_hat*z1inv*z1R_dot[i]

    p1_hat_dotL[i] = z4L_dot[i]*z1inv - p1_hat*z1inv*z1L_dot[i]
    p1_hat_dotR[i] = z4R_dot[i]*z1inv - p1_hat*z1inv*z1R_dot[i]

    # p2_hat is an intermediate variable for h_hat below
    p2_hat_dotL = ((gamma + 1)/(2*gamma))*(z4avg_dotL[i]/z1avg +
                      -z4avg/(z1avg*z1avg)*z1avg_dotL[i]) + 
                      ( gamma_1/(2*gamma))*p1_hat_dotL[i]
    p2_hat_dotR = ((gamma + 1)/(2*gamma))*(z4avg_dotR[i]/z1avg +
                      -z4avg/(z1avg*z1avg)*z1avg_dotR[i]) +
                      ( gamma_1/(2*gamma))*p1_hat_dotR[i]

    h_hat_dotL[i] = (gamma/gamma_1)*(p2_hat_dotL/rho_hat +
                     -p2_hat/(rho_hat*rho_hat)*rho_hat_dotL[i]) +
                     u_hat*u_hat_dotL[i] + v_hat*v_hat_dotL[i]
    h_hat_dotR[i] = (gamma/gamma_1)*(p2_hat_dotR/rho_hat +
                     -p2_hat/(rho_hat*rho_hat)*rho_hat_dotR[i]) +
                      u_hat*u_hat_dotR[i] + v_hat*v_hat_dotR[i]
  end


  mv_n = rho_hat*(nrm[1]*u_hat + nrm[2]*v_hat)  # normal momentum
  #F[1, j] = mv_n
  #F[2, j] = mv_n*u_hat + nrm[1, j]*p1_hat
  #F[3, j] = mv_n*v_hat + nrm[2, j]*p1_hat
  #F[4, j] = mv_n*h_hat

  @simd for i=1:4
    mv_n_dotL = (nrm[1]*u_hat + nrm[2]*v_hat)*rho_hat_dotL[i] + 
                rho_hat*(nrm[1]*u_hat_dotL[i] + nrm[2]*v_hat_dotL[i])
    mv_n_dotR = (nrm[1]*u_hat + nrm[2]*v_hat)*rho_hat_dotR[i] +
                rho_hat*(nrm[1]*u_hat_dotR[i] + nrm[2]*v_hat_dotR[i])

    FL_dot[1, i] += mv_n_dotL
    FL_dot[2, i] += u_hat*mv_n_dotL + mv_n*u_hat_dotL[i] + 
                      nrm[1]*p1_hat_dotL[i]
    FL_dot[3, i] += v_hat*mv_n_dotL + mv_n*v_hat_dotL[i] +
                      nrm[2]*p1_hat_dotL[i]
    FL_dot[4, i] += h_hat*mv_n_dotL + mv_n*h_hat_dotL[i]
    
    FR_dot[1, i] += mv_n_dotR
    FR_dot[2, i] += u_hat*mv_n_dotR + mv_n*u_hat_dotR[i] +
                      nrm[1]*p1_hat_dotR[i]
    FR_dot[3, i] += v_hat*mv_n_dotR + mv_n*v_hat_dotR[i] +
                      nrm[2]*p1_hat_dotR[i]
    FR_dot[4, i] += h_hat*mv_n_dotR + mv_n*h_hat_dotR[i]    
  end

  return nothing
end


function calcEulerFlux_IR_revq(params::ParamType{2, :conservative},
                      qL::AbstractArray{Tsol,1}, qL_bar::AbstractArray{Tsol, 1},
                      qR::AbstractArray{Tsol, 1}, qR_bar::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tres}, dir::AbstractArray{Tmsh, 1},  
                      F_bar::AbstractArray{Tres, 1}) where {Tmsh, Tsol, Tres}


  @debug1 begin
    @assert length(qL) == length(qL_bar)
    @assert length(qR) == length(qR_bar)
    @assert size(F_bar, 1) == length(qL_bar)
    @assert size(dir, 1) == 2
    @assert size(F_bar, 2) == size(dir, 2)
  end

  gamma = params.gamma
  gamma_1 = params.gamma_1
  pL = calcPressure(params, qL); pR = calcPressure(params, qR)
  z1L = sqrt(qL[1]/pL); z1R = sqrt(qR[1]/pR)
  z2L = z1L*qL[2]/qL[1]; z2R = z1R*qR[2]/qR[1]
  z3L = z1L*qL[3]/qL[1]; z3R = z1R*qR[3]/qR[1]
  z4L = sqrt(qL[1]*pL); z4R = sqrt(qR[1]*pR)

  z4avg = logavg(z4L, z4R)
  z1avg = logavg(z1L, z1R)
  rho_hat = 0.5*(z1L + z1R)*z4avg
  u_hat = (z2L + z2R)/(z1L + z1R)
  v_hat = (z3L + z3R)/(z1L + z1R)
  p1_hat = (z4L + z4R)/(z1L + z1R)
  p2_hat = ((gamma + 1)/(2*gamma) )*z4avg/z1avg + ( gamma_1/(2*gamma) )*(z4L + z4R)/(z1L + z1R)
  h_hat = gamma*p2_hat/(rho_hat*gamma_1) + 0.5*(u_hat*u_hat + v_hat*v_hat)

#=
  mv_n = rho_hat*(dir[1]*u_hat + dir[2]*v_hat)  # normal momentum
  F[1] = mv_n
  F[2] = mv_n*u_hat + dir[1]*p1_hat
  F[3] = mv_n*v_hat + dir[2]*p1_hat
  F[4] = mv_n*h_hat
=#
  # reverse sweep
  rho_hat_bar = zero(Tsol)
  u_hat_bar = zero(Tsol); v_hat_bar = zero(Tsol); p1_hat_bar = zero(Tsol)
  p2_hat_bar = zero(Tsol); h_hat_bar = zero(Tsol);
  z4avg_bar = zero(Tsol); z1avg_bar = zero(Tsol)

  z1L_bar = zero(Tsol); z2L_bar = zero(Tsol); z3L_bar = zero(Tsol); z4L_bar = zero(Tsol)
  z1R_bar = zero(Tsol); z2R_bar = zero(Tsol); z3R_bar = zero(Tsol); z4R_bar = zero(Tsol)
  pL_bar = zero(Tsol); pR_bar = zero(Tsol); z1avg_bar = zero(Tsol); z4avg_bar = zero(Tsol)

  mv_n = rho_hat*(dir[1]*u_hat + dir[2]*v_hat)  # normal momentum
  mv_n_bar = zero(Tres)

  # F[1]
  mv_n_bar += F_bar[1]

  # F[2]
  mv_n_bar  += u_hat*F_bar[2]
  u_hat_bar += mv_n*F_bar[2]
  p1_hat_bar += dir[1]*F_bar[2]

  # F[3]
  mv_n_bar += v_hat*F_bar[3]
  v_hat_bar += mv_n*F_bar[3]
  p1_hat_bar += dir[2]*F_bar[3]

  # F[4]
  mv_n_bar += h_hat*F_bar[4]
  h_hat_bar += mv_n*F_bar[4]

  # mv_n
  rho_hat_bar += (dir[1]*u_hat + dir[2]*v_hat)*mv_n_bar
  u_hat_bar += rho_hat*dir[1]*mv_n_bar
  v_hat_bar += rho_hat*dir[2]*mv_n_bar

  # h_hat
  p2_hat_bar += gamma*h_hat_bar/(rho_hat*gamma_1)
  rho_hat_bar += (-gamma*p2_hat/(rho_hat*rho_hat*gamma_1))*h_hat_bar
  u_hat_bar += u_hat*h_hat_bar
  v_hat_bar += v_hat*h_hat_bar

  # p2_hat
  p2tmp = (gamma_1/(2*gamma))*(z4L + z4R)/(z1L + z1R)
  z4avg_bar +=  ((gamma + 1)/(2*gamma))*p2_hat_bar/z1avg
  z1avg_bar += -(((gamma + 1)/(2*gamma))*z4avg/(z1avg*z1avg))*p2_hat_bar
  z4L_bar += (gamma_1/(2*gamma))*p2_hat_bar/(z1L + z1R)
  z4R_bar += (gamma_1/(2*gamma))*p2_hat_bar/(z1L + z1R)
  z1L_bar += -p2tmp/(z1L + z1R)*p2_hat_bar
  z1R_bar += -p2tmp/(z1L + z1R)*p2_hat_bar

  # p1_hat 
  z4L_bar += p1_hat_bar/(z1L + z1R)
  z4R_bar += p1_hat_bar/(z1L + z1R)
  z1L_bar += -p1_hat/(z1L + z1R)*p1_hat_bar
  z1R_bar += -p1_hat/(z1L + z1R)*p1_hat_bar

  # v_hat
  z3L_bar += v_hat_bar/(z1L + z1R)
  z3R_bar += v_hat_bar/(z1L + z1R)
  z1L_bar += -v_hat/(z1L + z1R)*v_hat_bar
  z1R_bar += -v_hat/(z1L + z1R)*v_hat_bar

  # u_hat
  z2L_bar += u_hat_bar/(z1L + z1R)
  z2R_bar += u_hat_bar/(z1L + z1R)
  z1L_bar += -u_hat/(z1L + z1R)*u_hat_bar
  z1R_bar += -u_hat/(z1L + z1R)*u_hat_bar

  # rho_hat
  z1L_bar += 0.5*z4avg*rho_hat_bar
  z1R_bar += 0.5*z4avg*rho_hat_bar
  z4avg_bar += 0.5*(z1L + z1R)*rho_hat_bar

  # log averages
  z4L_bar_tmp, z4R_bar_tmp = logavg_rev(z4L, z4R, z4avg_bar)
  z4L_bar += z4L_bar_tmp
  z4R_bar += z4R_bar_tmp

  z1L_bar_tmp, z1R_bar_tmp = logavg_rev(z1L, z1R, z1avg_bar)
  z1L_bar += z1L_bar_tmp
  z1R_bar += z1R_bar_tmp


  # z4L/R
  qL_bar[1] += 0.5*pL*z4L_bar/z4L; qR_bar[1] += 0.5*pR*z4R_bar/z4R
  pL_bar += 0.5*qL[1]*z4L_bar/z4L; pR_bar += 0.5*qR[1]*z4R_bar/z4R

  # z3L/R
  qL_bar[3] +=  z1L*z3L_bar/qL[1]; qR_bar[3] +=  z1R*z3R_bar/qR[1]
  qL_bar[1] += -z3L*z3L_bar/qL[1]; qR_bar[1] += -z3R*z3R_bar/qR[1]
  z1L_bar   +=  qL[3]*z3L_bar/qL[1]; z1R_bar +=  qR[3]*z3R_bar/qR[1]

  # z2L/R
  # z2L = z1L*qL[2]/qL[1]; z2R = z1R*qR[2]/qR[1]
  qL_bar[2] +=  z1L*z2L_bar/qL[1]; qR_bar[2] +=  z1R*z2R_bar/qR[1]
  qL_bar[1] += -z2L*z2L_bar/qL[1]; qR_bar[1] += -z2R*z2R_bar/qR[1]
  z1L_bar   +=  qL[2]*z2L_bar/qL[1]; z1R_bar +=  qR[2]*z2R_bar/qR[1]

  # z1L/R
  qL_bar[1] += (0.5/z1L)*z1L_bar/pL;  qR_bar[1] += (0.5/z1R)*z1R_bar/pR
  pL_bar += -(0.5/z1L)*(qL[1]/(pL*pL))*z1L_bar
  pR_bar += -(0.5/z1R)*(qR[1]/(pR*pR))*z1R_bar

  calcPressure_revq(params, qL, qL_bar, pL_bar)
  calcPressure_revq(params, qR, qR_bar, pR_bar)

  return nothing
end


function calcEulerFlux_IR_revm(params::ParamType{2, :conservative},
                  qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                  aux_vars::AbstractArray{Tres},
                  nrm::AbstractArray{Tmsh, 1}, nrm_bar::AbstractArray{Tmsh, 1},
                  F_bar::AbstractArray{Tres,1}) where {Tmsh, Tsol, Tres}

  @debug1 begin
    @assert size(nrm, 1) == 2
    @assert length(qL) == length(qR)
    @assert length(qL) == size(F_bar, 1)
    @assert size(nrm_bar, 1) == size(nrm, 1)
    @assert size(F_bar, 1) == length(qL)
  end


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
#=
  for i=1:2
    mv_n = rho_hat*(nrm[1]*u_hat + nrm[2]*v_hat)  # normal momentum
    F[1] = mv_n
    F[2] = mv_n*u_hat + nrm[1]*p1_hat
    F[3] = mv_n*v_hat + nrm[2]*p1_hat
    F[4] = mv_n*h_hat
  end
=#
  # reverse sweep
    mv_n = rho_hat*(nrm[1]*u_hat + nrm[2]*v_hat)  # normal momentum
    mv_n_bar = zero(mv_n)

    # F[4]
    mv_n_bar += h_hat*F_bar[4]

    # F[3]
    mv_n_bar += v_hat*F_bar[3]
    nrm_bar[2] += p1_hat*F_bar[3]

    # F[2]
    mv_n_bar += u_hat*F_bar[2]
    nrm_bar[1] += p1_hat*F_bar[2]

    # F[1]
    mv_n_bar += F_bar[1]

    # mv_n
    nrm_bar[1] += rho_hat*u_hat*mv_n_bar
    nrm_bar[2] += rho_hat*v_hat*mv_n_bar

  return nothing
end







#------------------------------------------------------------------------------
# 2D, multi-direction

"""
  Differentiated version of the multi-dimensional version of the IR flux,
  vector forward mode.

  The user is responsible for zeroing out the output arrays if needed!

  **Inputs**

   * Params: ParamType
   * qL: solution at the left node (numDofPerNode)
   * qg: solution at the right node (numDofPerNode)
   * aux_vars
   * nrm: mesh.dim x mesh.dim matrix of normal vectors, one per column
   * FL_dot: numDofPerNode x numDofPerNode x dim, jacobian of the flux
                with respect to qL, in each direction (overwritten)
   * FR_dot  similar to FR_dot, but jacobian with respect to qR (overwritten)

"""
function calcEulerFlux_IR_diff(params::ParamType{2, :conservative},
                   qL::AbstractArray{Tsol,1},
                   qR::AbstractArray{Tsol, 1},
                   aux_vars::AbstractArray{Tres, 1},
                   nrm::AbstractArray{Tmsh, 2},
                   FL_dot::AbstractArray{Tres, 3},
                   FR_dot::AbstractArray{Tres, 3}) where {Tmsh, Tsol, Tres}

  @debug1 begin
    @assert size(nrm, 1) == 2
    @assert size(nrm, 2) == 2
    @assert length(qL) == length(qR)
    @assert length(qL) == size(FL_dot, 1)
    @assert length(qL) == size(FL_dot, 2)
    @assert size(nrm, 2) == size(FL_dot, 3)
    @assert size(FL_dot, 1) == size(FR_dot, 1)
    @assert size(FL_dot, 2) == size(FR_dot, 2)
    @assert size(FR_dot, 3) == size(FR_dot, 3)
  end


  data = params.irfluxdata
  @unpack data pL_dot pR_dot
  @unpack data z1L_dot z2L_dot z3L_dot z4L_dot z1R_dot z2R_dot z3R_dot z4R_dot

  gamma = params.gamma
  gamma_1 = params.gamma_1

  pL = calcPressure_diff(params, qL, pL_dot)
  pR = calcPressure_diff(params, qR, pR_dot)
  uL = qL[2]/qL[1]; uR = qR[2]/qR[1]
  vL = qL[3]/qL[1]; vR = qR[3]/qR[1]
  z1L = sqrt(qL[1]/pL); z1R = sqrt(qR[1]/pR)
  z2L = z1L*uL; z2R = z1R*uR
  z3L = z1L*vL; z3R = z1R*vR
  z4L = sqrt(qL[1]*pL); z4R = sqrt(qR[1]*pR)

  fastzero!(z1L_dot); fastzero!(z1R_dot)
  fastzero!(z2L_dot); fastzero!(z2R_dot)
  fastzero!(z3L_dot); fastzero!(z3R_dot)
  fastzero!(z4L_dot); fastzero!(z4R_dot)

  # differentiate with respect to q (not including chain rule terms for p)
  z1L_dot[1] = (0.5/z1L)*1/pL; z1R_dot[1] = (0.5/z1R)*1/pR

  z2L_dot[1] = -z2L/qL[1]; z2R_dot[1] = -z2R/qR[1]
  z2L_dot[2] =  z1L/qL[1]; z2R_dot[2] =  z1R/qR[1]

  z3L_dot[1] = -z3L/qL[1]; z3R_dot[1] = -z3R/qR[1]
  z3L_dot[3] =  z1L/qL[1]; z3R_dot[3] =  z1R/qR[1]

  z4L_dot[1] =  (0.5/z4L)*pL; z4R_dot[1] = (0.5/z4R)*pR

  # do the pressure/z1L related terms
  @simd for i=1:4
    z1L_dot[i] += (0.5/z1L)*(-qL[1]/(pL*pL))*pL_dot[i]
    z1R_dot[i] += (0.5/z1R)*(-qR[1]/(pR*pR))*pR_dot[i]

    z2L_dot[i] += uL*z1L_dot[i]
    z2R_dot[i] += uR*z1R_dot[i]

    z3L_dot[i] += vL*z1L_dot[i]
    z3R_dot[i] += vR*z1R_dot[i]

    z4L_dot[i] += (0.5/z4L)*qL[1]*pL_dot[i]
    z4R_dot[i] += (0.5/z4R)*qR[1]*pR_dot[i]
  end

  @unpack data avgdata z4avg_dotL z4avg_dotR z1avg_dotL z1avg_dotR
  @unpack data rho_hat_dotL rho_hat_dotR u_hat_dotL u_hat_dotR
  @unpack data v_hat_dotL v_hat_dotR p1_hat_dotL p1_hat_dotR
  @unpack data h_hat_dotL h_hat_dotR

  # z4avg_dotL/r, z1avg_dotL/r, rho_hat, u_hat, v_hat, p1_hat, p2_hat
  z4avg = logavg_diff(avgdata, z4L, z4L_dot, z4R, z4R_dot, z4avg_dotL, z4avg_dotR)
  z1avg = logavg_diff(avgdata, z1L, z1L_dot, z1R, z1R_dot, z1avg_dotL, z1avg_dotR)

  z1 = z1L + z1R
  z1inv = 1/z1
  rho_hat = 0.5*(z1L + z1R)*z4avg
  u_hat = (z2L + z2R)*z1inv
  v_hat = (z3L + z3R)*z1inv
  p1_hat = (z4L + z4R)*z1inv
  p2_hat = ((gamma + 1)/(2*gamma) )*z4avg/z1avg + ( gamma_1/(2*gamma) )*p1_hat
  h_hat = gamma*p2_hat/(rho_hat*gamma_1) + 0.5*(u_hat*u_hat + v_hat*v_hat)

  @simd for i=1:4
    rho_hat_dotL[i] = 0.5*(z4avg*z1L_dot[i] + z1*z4avg_dotL[i])
    rho_hat_dotR[i] = 0.5*(z4avg*z1R_dot[i] + z1*z4avg_dotR[i])

    u_hat_dotL[i] = z2L_dot[i]*z1inv - u_hat*z1inv*z1L_dot[i]
    u_hat_dotR[i] = z2R_dot[i]*z1inv - u_hat*z1inv*z1R_dot[i]

    v_hat_dotL[i] = z3L_dot[i]*z1inv - v_hat*z1inv*z1L_dot[i]
    v_hat_dotR[i] = z3R_dot[i]*z1inv - v_hat*z1inv*z1R_dot[i]

    p1_hat_dotL[i] = z4L_dot[i]*z1inv - p1_hat*z1inv*z1L_dot[i]
    p1_hat_dotR[i] = z4R_dot[i]*z1inv - p1_hat*z1inv*z1R_dot[i]

    # p2_hat is an intermediate variable for h_hat below
    p2_hat_dotL = ((gamma + 1)/(2*gamma))*(z4avg_dotL[i]/z1avg +
                      -z4avg/(z1avg*z1avg)*z1avg_dotL[i]) + 
                      ( gamma_1/(2*gamma))*p1_hat_dotL[i]
    p2_hat_dotR = ((gamma + 1)/(2*gamma))*(z4avg_dotR[i]/z1avg +
                      -z4avg/(z1avg*z1avg)*z1avg_dotR[i]) +
                      ( gamma_1/(2*gamma))*p1_hat_dotR[i]

    h_hat_dotL[i] = (gamma/gamma_1)*(p2_hat_dotL/rho_hat +
                     -p2_hat/(rho_hat*rho_hat)*rho_hat_dotL[i]) +
                     u_hat*u_hat_dotL[i] + v_hat*v_hat_dotL[i]
    h_hat_dotR[i] = (gamma/gamma_1)*(p2_hat_dotR/rho_hat +
                     -p2_hat/(rho_hat*rho_hat)*rho_hat_dotR[i]) +
                      u_hat*u_hat_dotR[i] + v_hat*v_hat_dotR[i]
  end


  for j=1:2
    mv_n = rho_hat*(nrm[1, j]*u_hat + nrm[2, j]*v_hat)  # normal momentum
    #F[1, j] = mv_n
    #F[2, j] = mv_n*u_hat + nrm[1, j]*p1_hat
    #F[3, j] = mv_n*v_hat + nrm[2, j]*p1_hat
    #F[4, j] = mv_n*h_hat

    @simd for i=1:4
      mv_n_dotL = (nrm[1, j]*u_hat + nrm[2, j]*v_hat)*rho_hat_dotL[i] + 
                  rho_hat*(nrm[1, j]*u_hat_dotL[i] + nrm[2, j]*v_hat_dotL[i])
      mv_n_dotR = (nrm[1, j]*u_hat + nrm[2, j]*v_hat)*rho_hat_dotR[i] +
                  rho_hat*(nrm[1, j]*u_hat_dotR[i] + nrm[2, j]*v_hat_dotR[i])

      FL_dot[1, i, j] += mv_n_dotL
      FL_dot[2, i, j] += u_hat*mv_n_dotL + mv_n*u_hat_dotL[i] + 
                        nrm[1, j]*p1_hat_dotL[i]
      FL_dot[3, i, j] += v_hat*mv_n_dotL + mv_n*v_hat_dotL[i] +
                        nrm[2, j]*p1_hat_dotL[i]
      FL_dot[4, i, j] += h_hat*mv_n_dotL + mv_n*h_hat_dotL[i]
      
      FR_dot[1, i, j] += mv_n_dotR
      FR_dot[2, i, j] += u_hat*mv_n_dotR + mv_n*u_hat_dotR[i] +
                        nrm[1, j]*p1_hat_dotR[i]
      FR_dot[3, i, j] += v_hat*mv_n_dotR + mv_n*v_hat_dotR[i] +
                        nrm[2, j]*p1_hat_dotR[i]
      FR_dot[4, i, j] += h_hat*mv_n_dotR + mv_n*h_hat_dotR[i]
      
    end
  end

  return nothing
end

"""
  Reverse mode of IR flux with respect to q

  **Inputs**

   * params: ParamType
   * qL: solution at the left state
   * qR: solution at the right state
   * aux_vars
   * dir: `dim` x `dim` array of normal vectors, one per column
   * F_bar: `numDofPerNode` x `dim` seed vector for the flux

  **Inputs/Outputs**

   * qL_bar: vector to sum the result for the left state into (not overwritten)
   * qR_bar: vector to sum the result for the right state into (not overwritten)
"""
function calcEulerFlux_IR_revq(params::ParamType{2, :conservative},
                      qL::AbstractArray{Tsol,1}, qL_bar::AbstractArray{Tsol, 1},
                      qR::AbstractArray{Tsol, 1}, qR_bar::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tres}, dir::AbstractArray{Tmsh, 2},  
                      F_bar::AbstractArray{Tres, 2}) where {Tmsh, Tsol, Tres}


  @debug1 begin
    @assert length(qL) == length(qL_bar)
    @assert length(qR) == length(qR_bar)
    @assert size(F_bar, 1) == length(qL_bar)
    @assert size(dir, 1) == 2
    @assert size(F_bar, 2) == size(dir, 2)
  end

  gamma = params.gamma
  gamma_1 = params.gamma_1
  pL = calcPressure(params, qL); pR = calcPressure(params, qR)
  z1L = sqrt(qL[1]/pL); z1R = sqrt(qR[1]/pR)
  z2L = z1L*qL[2]/qL[1]; z2R = z1R*qR[2]/qR[1]
  z3L = z1L*qL[3]/qL[1]; z3R = z1R*qR[3]/qR[1]
  z4L = sqrt(qL[1]*pL); z4R = sqrt(qR[1]*pR)

  z4avg = logavg(z4L, z4R)
  z1avg = logavg(z1L, z1R)
  rho_hat = 0.5*(z1L + z1R)*z4avg
  u_hat = (z2L + z2R)/(z1L + z1R)
  v_hat = (z3L + z3R)/(z1L + z1R)
  p1_hat = (z4L + z4R)/(z1L + z1R)
  p2_hat = ((gamma + 1)/(2*gamma) )*z4avg/z1avg + ( gamma_1/(2*gamma) )*(z4L + z4R)/(z1L + z1R)
  h_hat = gamma*p2_hat/(rho_hat*gamma_1) + 0.5*(u_hat*u_hat + v_hat*v_hat)

#=
  for i=1:2
    mv_n = rho_hat*(dir[1, i]*u_hat + dir[2, i]*v_hat)  # normal momentum
    F[1, i] = mv_n
    F[2, i] = mv_n*u_hat + dir[1, i]*p1_hat
    F[3, i] = mv_n*v_hat + dir[2, i]*p1_hat
    F[4, i] = mv_n*h_hat
  end
=#
  # reverse sweep
  rho_hat_bar = zero(Tsol)
  u_hat_bar = zero(Tsol); v_hat_bar = zero(Tsol); p1_hat_bar = zero(Tsol)
  p2_hat_bar = zero(Tsol); h_hat_bar = zero(Tsol);
  z4avg_bar = zero(Tsol); z1avg_bar = zero(Tsol)

  z1L_bar = zero(Tsol); z2L_bar = zero(Tsol); z3L_bar = zero(Tsol); z4L_bar = zero(Tsol)
  z1R_bar = zero(Tsol); z2R_bar = zero(Tsol); z3R_bar = zero(Tsol); z4R_bar = zero(Tsol)
  pL_bar = zero(Tsol); pR_bar = zero(Tsol); z1avg_bar = zero(Tsol); z4avg_bar = zero(Tsol)

  for i=1:2
    mv_n = rho_hat*(dir[1, i]*u_hat + dir[2, i]*v_hat)  # normal momentum
    mv_n_bar = zero(Tres)

    # F[1, i]
    mv_n_bar += F_bar[1, i]

    # F[2, i]
    mv_n_bar  += u_hat*F_bar[2, i]
    u_hat_bar += mv_n*F_bar[2, i]
    p1_hat_bar += dir[1, i]*F_bar[2, i]

    # F[3, i]
    mv_n_bar += v_hat*F_bar[3, i]
    v_hat_bar += mv_n*F_bar[3, i]
    p1_hat_bar += dir[2, i]*F_bar[3, i]

    # F[4, i]
    mv_n_bar += h_hat*F_bar[4, i]
    h_hat_bar += mv_n*F_bar[4, i]

    # mv_n
    rho_hat_bar += (dir[1, i]*u_hat + dir[2, i]*v_hat)*mv_n_bar
    u_hat_bar += rho_hat*dir[1, i]*mv_n_bar
    v_hat_bar += rho_hat*dir[2, i]*mv_n_bar
  end

  # h_hat
  p2_hat_bar += gamma*h_hat_bar/(rho_hat*gamma_1)
  rho_hat_bar += (-gamma*p2_hat/(rho_hat*rho_hat*gamma_1))*h_hat_bar
  u_hat_bar += u_hat*h_hat_bar
  v_hat_bar += v_hat*h_hat_bar

  # p2_hat
  p2tmp = (gamma_1/(2*gamma))*(z4L + z4R)/(z1L + z1R)
  z4avg_bar +=  ((gamma + 1)/(2*gamma))*p2_hat_bar/z1avg
  z1avg_bar += -(((gamma + 1)/(2*gamma))*z4avg/(z1avg*z1avg))*p2_hat_bar
  z4L_bar += (gamma_1/(2*gamma))*p2_hat_bar/(z1L + z1R)
  z4R_bar += (gamma_1/(2*gamma))*p2_hat_bar/(z1L + z1R)
  z1L_bar += -p2tmp/(z1L + z1R)*p2_hat_bar
  z1R_bar += -p2tmp/(z1L + z1R)*p2_hat_bar

  # p1_hat 
  z4L_bar += p1_hat_bar/(z1L + z1R)
  z4R_bar += p1_hat_bar/(z1L + z1R)
  z1L_bar += -p1_hat/(z1L + z1R)*p1_hat_bar
  z1R_bar += -p1_hat/(z1L + z1R)*p1_hat_bar

  # v_hat
  z3L_bar += v_hat_bar/(z1L + z1R)
  z3R_bar += v_hat_bar/(z1L + z1R)
  z1L_bar += -v_hat/(z1L + z1R)*v_hat_bar
  z1R_bar += -v_hat/(z1L + z1R)*v_hat_bar

  # u_hat
  z2L_bar += u_hat_bar/(z1L + z1R)
  z2R_bar += u_hat_bar/(z1L + z1R)
  z1L_bar += -u_hat/(z1L + z1R)*u_hat_bar
  z1R_bar += -u_hat/(z1L + z1R)*u_hat_bar

  # rho_hat
  z1L_bar += 0.5*z4avg*rho_hat_bar
  z1R_bar += 0.5*z4avg*rho_hat_bar
  z4avg_bar += 0.5*(z1L + z1R)*rho_hat_bar

  # log averages
  z4L_bar_tmp, z4R_bar_tmp = logavg_rev(z4L, z4R, z4avg_bar)
  z4L_bar += z4L_bar_tmp
  z4R_bar += z4R_bar_tmp

  z1L_bar_tmp, z1R_bar_tmp = logavg_rev(z1L, z1R, z1avg_bar)
  z1L_bar += z1L_bar_tmp
  z1R_bar += z1R_bar_tmp


  # z4L/R
  qL_bar[1] += 0.5*pL*z4L_bar/z4L; qR_bar[1] += 0.5*pR*z4R_bar/z4R
  pL_bar += 0.5*qL[1]*z4L_bar/z4L; pR_bar += 0.5*qR[1]*z4R_bar/z4R

  # z3L/R
  qL_bar[3] +=  z1L*z3L_bar/qL[1]; qR_bar[3] +=  z1R*z3R_bar/qR[1]
  qL_bar[1] += -z3L*z3L_bar/qL[1]; qR_bar[1] += -z3R*z3R_bar/qR[1]
  z1L_bar   +=  qL[3]*z3L_bar/qL[1]; z1R_bar +=  qR[3]*z3R_bar/qR[1]

  # z2L/R
  # z2L = z1L*qL[2]/qL[1]; z2R = z1R*qR[2]/qR[1]
  qL_bar[2] +=  z1L*z2L_bar/qL[1]; qR_bar[2] +=  z1R*z2R_bar/qR[1]
  qL_bar[1] += -z2L*z2L_bar/qL[1]; qR_bar[1] += -z2R*z2R_bar/qR[1]
  z1L_bar   +=  qL[2]*z2L_bar/qL[1]; z1R_bar +=  qR[2]*z2R_bar/qR[1]

  # z1L/R
  qL_bar[1] += (0.5/z1L)*z1L_bar/pL;  qR_bar[1] += (0.5/z1R)*z1R_bar/pR
  pL_bar += -(0.5/z1L)*(qL[1]/(pL*pL))*z1L_bar
  pR_bar += -(0.5/z1R)*(qR[1]/(pR*pR))*z1R_bar

  calcPressure_revq(params, qL, qL_bar, pL_bar)
  calcPressure_revq(params, qR, qR_bar, pR_bar)

  return nothing
end


"""
  Reverse mode of [`calcIRFLux`](@ref) with respect to the normal vector.

  **Inputs**

   * params: ParamType
   * qL: solution at the left state
   * qR: solution at the right state
   * aux_vars
   * nrm: `dim` x `dim array of normal vectors, one per column
   * F_bar: seed vector for the flux, `numDofPerNode` x `dim`

  **Inputs**

   * nrm_bar: array, same size as `nrm` to sum the result into (not overwritten)

"""
function calcEulerFlux_IR_revm(params::ParamType{2, :conservative},
                  qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                  aux_vars::AbstractArray{Tres},
                  nrm::AbstractArray{Tmsh, 2}, nrm_bar::AbstractArray{Tmsh, 2},
                  F_bar::AbstractArray{Tres,2}) where {Tmsh, Tsol, Tres}

  @debug1 begin
    @assert size(F_bar, 1) == length(qL)
    @assert size(nrm, 1) == 2
    @assert size(nrm, 2) == 2
    @assert size(nrm_bar, 1) == size(nrm_bar, 2)
    @assert size(nrm_bar, 2) == size(nrm_bar, 2)
    @assert size(F_bar, 2) == size(nrm, 2)
  end


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
#=
  for i=1:2
    mv_n = rho_hat*(nrm[1, i]*u_hat + nrm[2, i]*v_hat)  # normal momentum
    F[1, i] = mv_n
    F[2, i] = mv_n*u_hat + nrm[1, i]*p1_hat
    F[3, i] = mv_n*v_hat + nrm[2, i]*p1_hat
    F[4, i] = mv_n*h_hat
  end
=#
  # reverse sweep
  for i=1:2
    mv_n = rho_hat*(nrm[1, i]*u_hat + nrm[2, i]*v_hat)  # normal momentum
    mv_n_bar = zero(mv_n)

    # F[4, i]
    mv_n_bar += h_hat*F_bar[4, i]

    # F[3, i]
    mv_n_bar += v_hat*F_bar[3, i]
    nrm_bar[2, i] += p1_hat*F_bar[3, i]

    # F[2, i]
    mv_n_bar += u_hat*F_bar[2, i]
    nrm_bar[1, i] += p1_hat*F_bar[2, i]

    # F[1, i]
    mv_n_bar += F_bar[1, i]

    # mv_n
    nrm_bar[1, i] += rho_hat*u_hat*mv_n_bar
    nrm_bar[2, i] += rho_hat*v_hat*mv_n_bar
  end

  return nothing
end



#------------------------------------------------------------------------------
# 3D, single direction

"""
  3D, single direction version
"""
function calcEulerFlux_IR_diff(params::ParamType{3, :conservative},
                   qL::AbstractArray{Tsol,1},
                   qR::AbstractArray{Tsol, 1},
                   aux_vars::AbstractArray{Tres, 1},
                   nrm::AbstractArray{Tmsh, 1},
                   FL_dot::AbstractArray{Tres, 2},
                   FR_dot::AbstractArray{Tres, 2}) where {Tmsh, Tsol, Tres}

  @debug1 begin
    @assert size(nrm, 1) == 3
    @assert length(qL) == length(qR)
    @assert length(qL) == size(FL_dot, 1)
    @assert length(qL) == size(FL_dot, 2)
    @assert size(nrm, 1) == 3
    @assert size(FL_dot, 1) == size(FR_dot, 1)
    @assert size(FL_dot, 2) == size(FR_dot, 2)
  end

  data = params.irfluxdata
  @unpack data pL_dot pR_dot
  @unpack data z1L_dot z2L_dot z3L_dot z4L_dot z5L_dot
  @unpack data z1R_dot z2R_dot z3R_dot z4R_dot z5R_dot

  gamma = params.gamma
  gamma_1 = params.gamma_1

  pL = calcPressure_diff(params, qL, pL_dot)
  pR = calcPressure_diff(params, qR, pR_dot)
  uL = qL[2]/qL[1]; uR = qR[2]/qR[1]
  vL = qL[3]/qL[1]; vR = qR[3]/qR[1]
  wL = qL[4]/qL[1]; wR = qR[4]/qR[1]
  z1L = sqrt(qL[1]/pL); z1R = sqrt(qR[1]/pR)
  z2L = z1L*uL; z2R = z1R*uR
  z3L = z1L*vL; z3R = z1R*vR
  z4L = z1L*wL; z4R = z1R*wR
  z5L = sqrt(qL[1]*pL); z5R = sqrt(qR[1]*pR)

  fastzero!(z1L_dot); fastzero!(z1R_dot)
  fastzero!(z2L_dot); fastzero!(z2R_dot)
  fastzero!(z3L_dot); fastzero!(z3R_dot)
  fastzero!(z4L_dot); fastzero!(z4R_dot)
  fastzero!(z5L_dot); fastzero!(z5R_dot)

  # differentiate with respect to q (not including chain rule terms for p)
  z1L_dot[1] = (0.5/z1L)*1/pL; z1R_dot[1] = (0.5/z1R)*1/pR

  z2L_dot[1] = -z2L/qL[1]; z2R_dot[1] = -z2R/qR[1]
  z2L_dot[2] =  z1L/qL[1]; z2R_dot[2] =  z1R/qR[1]

  z3L_dot[1] = -z3L/qL[1]; z3R_dot[1] = -z3R/qR[1]
  z3L_dot[3] =  z1L/qL[1]; z3R_dot[3] =  z1R/qR[1]

  z4L_dot[1] = -z4L/qL[1]; z4R_dot[1] = -z4R/qR[1]
  z4L_dot[4] =  z1L/qL[1]; z4R_dot[4] =  z1R/qR[1]

  z5L_dot[1] =  (0.5/z5L)*pL; z5R_dot[1] = (0.5/z5R)*pR

  # do the pressure/z1L related terms
  @simd for i=1:5
    z1L_dot[i] += (0.5/z1L)*(-qL[1]/(pL*pL))*pL_dot[i]
    z1R_dot[i] += (0.5/z1R)*(-qR[1]/(pR*pR))*pR_dot[i]

    z2L_dot[i] += uL*z1L_dot[i]
    z2R_dot[i] += uR*z1R_dot[i]

    z3L_dot[i] += vL*z1L_dot[i]
    z3R_dot[i] += vR*z1R_dot[i]

    z4L_dot[i] += wL*z1L_dot[i]
    z4R_dot[i] += wR*z1R_dot[i]

    z5L_dot[i] += (0.5/z5L)*qL[1]*pL_dot[i]
    z5R_dot[i] += (0.5/z5R)*qR[1]*pR_dot[i]
  end

  @unpack data avgdata z5avg_dotL z5avg_dotR z1avg_dotL z1avg_dotR
  @unpack data rho_hat_dotL rho_hat_dotR u_hat_dotL u_hat_dotR
  @unpack data v_hat_dotL v_hat_dotR w_hat_dotL w_hat_dotR 
  @unpack data p1_hat_dotL p1_hat_dotR h_hat_dotL h_hat_dotR

  # z4avg_dotL/r, z1avg_dotL/r, rho_hat, u_hat, v_hat, p1_hat, p2_hat
  z5avg = logavg_diff(avgdata, z5L, z5L_dot, z5R, z5R_dot, z5avg_dotL, z5avg_dotR)
  z1avg = logavg_diff(avgdata, z1L, z1L_dot, z1R, z1R_dot, z1avg_dotL, z1avg_dotR)

  z1 = z1L + z1R
  z1inv = 1/z1
  rho_hat = 0.5*(z1L + z1R)*z5avg
  u_hat = (z2L + z2R)*z1inv
  v_hat = (z3L + z3R)*z1inv
  w_hat = (z4L + z4R)*z1inv
  p1_hat = (z5L + z5R)*z1inv
  p2_hat = ((gamma + 1)/(2*gamma) )*z5avg/z1avg + ( gamma_1/(2*gamma) )*p1_hat
  h_hat = gamma*p2_hat/(rho_hat*gamma_1) + 0.5*(u_hat*u_hat + v_hat*v_hat + w_hat*w_hat)

  @simd for i=1:5
    rho_hat_dotL[i] = 0.5*(z5avg*z1L_dot[i] + z1*z5avg_dotL[i])
    rho_hat_dotR[i] = 0.5*(z5avg*z1R_dot[i] + z1*z5avg_dotR[i])

    u_hat_dotL[i] = z2L_dot[i]*z1inv - u_hat*z1inv*z1L_dot[i]
    u_hat_dotR[i] = z2R_dot[i]*z1inv - u_hat*z1inv*z1R_dot[i]

    v_hat_dotL[i] = z3L_dot[i]*z1inv - v_hat*z1inv*z1L_dot[i]
    v_hat_dotR[i] = z3R_dot[i]*z1inv - v_hat*z1inv*z1R_dot[i]

    w_hat_dotL[i] = z4L_dot[i]*z1inv - w_hat*z1inv*z1L_dot[i]
    w_hat_dotR[i] = z4R_dot[i]*z1inv - w_hat*z1inv*z1R_dot[i]

    p1_hat_dotL[i] = z5L_dot[i]*z1inv - p1_hat*z1inv*z1L_dot[i]
    p1_hat_dotR[i] = z5R_dot[i]*z1inv - p1_hat*z1inv*z1R_dot[i]

    # p2_hat is an intermediate variable for h_hat below
    p2_hat_dotL = ((gamma + 1)/(2*gamma))*(z5avg_dotL[i]/z1avg +
                      -z5avg/(z1avg*z1avg)*z1avg_dotL[i]) + 
                      ( gamma_1/(2*gamma))*p1_hat_dotL[i]
    p2_hat_dotR = ((gamma + 1)/(2*gamma))*(z5avg_dotR[i]/z1avg +
                      -z5avg/(z1avg*z1avg)*z1avg_dotR[i]) +
                      ( gamma_1/(2*gamma))*p1_hat_dotR[i]

    h_hat_dotL[i] = (gamma/gamma_1)*(p2_hat_dotL/rho_hat +
                     -p2_hat/(rho_hat*rho_hat)*rho_hat_dotL[i]) +
                    u_hat*u_hat_dotL[i] + v_hat*v_hat_dotL[i] +
                    w_hat*w_hat_dotL[i]
    h_hat_dotR[i] = (gamma/gamma_1)*(p2_hat_dotR/rho_hat +
                     -p2_hat/(rho_hat*rho_hat)*rho_hat_dotR[i]) +
                      u_hat*u_hat_dotR[i] + v_hat*v_hat_dotR[i] + 
                      w_hat*w_hat_dotR[i]
  end

  nx = nrm[1]; ny = nrm[2]; nz = nrm[3]
  mv_n = rho_hat*(nx*u_hat + ny*v_hat + nz*w_hat)  # normal momentum
  #F[1] = mv_n
  #F[2] = mv_n*u_hat + nx*p1_hat
  #F[3] = mv_n*v_hat + ny*p1_hat
  #F[4] = mv_n*h_hat

  @simd for i=1:5
    mv_n_dotL = (nx*u_hat + ny*v_hat + nz*w_hat)*rho_hat_dotL[i] + 
                rho_hat*(nx*u_hat_dotL[i] + ny*v_hat_dotL[i] +
                         nz*w_hat_dotL[i])
    mv_n_dotR = (nx*u_hat + ny*v_hat + nz*w_hat)*rho_hat_dotR[i] +
                rho_hat*(nx*u_hat_dotR[i] + ny*v_hat_dotR[i] +
                         nz*w_hat_dotR[i])

    FL_dot[1, i] += mv_n_dotL
    FL_dot[2, i] += u_hat*mv_n_dotL + mv_n*u_hat_dotL[i] + 
                       nx*p1_hat_dotL[i]
    FL_dot[3, i] += v_hat*mv_n_dotL + mv_n*v_hat_dotL[i] +
                       ny*p1_hat_dotL[i]
    FL_dot[4, i] += w_hat*mv_n_dotL + mv_n*w_hat_dotL[i] +
                       nz*p1_hat_dotL[i]
    FL_dot[5, i] += h_hat*mv_n_dotL + mv_n*h_hat_dotL[i]
    
    FR_dot[1, i] += mv_n_dotR
    FR_dot[2, i] += u_hat*mv_n_dotR + mv_n*u_hat_dotR[i] +
                       nx*p1_hat_dotR[i]
    FR_dot[3, i] += v_hat*mv_n_dotR + mv_n*v_hat_dotR[i] +
                       ny*p1_hat_dotR[i]
    FR_dot[4, i] += w_hat*mv_n_dotR + mv_n*w_hat_dotR[i] +
                       nz*p1_hat_dotR[i]
    FR_dot[5, i] += h_hat*mv_n_dotR + mv_n*h_hat_dotR[i]    
  end

  return nothing
end


"""
  3D, single direction version
"""
function calcEulerFlux_IR_revq(params::ParamType{3, :conservative},
                      qL::AbstractArray{Tsol,1}, qL_bar::AbstractArray{Tsol, 1},
                      qR::AbstractArray{Tsol, 1}, qR_bar::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tres}, dir::AbstractArray{Tmsh, 1},  
                      F_bar::AbstractArray{Tres, 1}) where {Tmsh, Tsol, Tres}


  @debug1 begin
    @assert length(qL) == length(qL_bar)
    @assert length(qR) == length(qR_bar)
    @assert size(F_bar, 1) == length(qL_bar)
    @assert size(dir, 1) == 3
  end

  gamma = params.gamma
  gamma_1 = params.gamma_1
  pL = calcPressure(params, qL); pR = calcPressure(params, qR)
  z1L = sqrt(qL[1]/pL); z1R = sqrt(qR[1]/pR)
  z2L = z1L*qL[2]/qL[1]; z2R = z1R*qR[2]/qR[1]
  z3L = z1L*qL[3]/qL[1]; z3R = z1R*qR[3]/qR[1]
  z4L = z1L*qL[4]/qL[1]; z4R = z1R*qR[4]/qR[1]
  z5L = sqrt(qL[1]*pL); z5R = sqrt(qR[1]*pR)

  z5avg = logavg(z5L, z5R)
  z1avg = logavg(z1L, z1R)
  rho_hat = 0.5*(z1L + z1R)*z5avg
  u_hat = (z2L + z2R)/(z1L + z1R)
  v_hat = (z3L + z3R)/(z1L + z1R)
  w_hat = (z4L + z4R)/(z1L + z1R)
  p1_hat = (z5L + z5R)/(z1L + z1R)
  p2_hat = ((gamma + 1)/(2*gamma) )*z5avg/z1avg + ( gamma_1/(2*gamma) )*(z5L + z5R)/(z1L + z1R)
  h_hat = gamma*p2_hat/(rho_hat*gamma_1) + 0.5*(u_hat*u_hat + v_hat*v_hat + w_hat*w_hat)

#=
    mv_n = rho_hat*(dir[1]*u_hat + dir[2]*v_hat)  # normal momentum
    F[1] = mv_n
    F[2] = mv_n*u_hat + dir[1]*p1_hat
    F[3] = mv_n*v_hat + dir[2]*p1_hat
    F[4] = mv_n*w_hat + dir[3]*p1_hat
    F[4] = mv_n*h_hat
=#
  # reverse sweep
  rho_hat_bar = zero(Tsol)
  u_hat_bar = zero(Tsol); v_hat_bar = zero(Tsol); w_hat_bar = zero(Tsol)
  p1_hat_bar = zero(Tsol); p2_hat_bar = zero(Tsol); h_hat_bar = zero(Tsol);
  z5avg_bar = zero(Tsol); z1avg_bar = zero(Tsol)

  z1L_bar = zero(Tsol); z2L_bar = zero(Tsol); z3L_bar = zero(Tsol);
  z4L_bar = zero(Tsol); z5L_bar = zero(Tsol)
  z1R_bar = zero(Tsol); z2R_bar = zero(Tsol); z3R_bar = zero(Tsol);
  z4R_bar = zero(Tsol); z5R_bar = zero(Tsol)
  pL_bar = zero(Tsol); pR_bar = zero(Tsol); z1avg_bar = zero(Tsol); z5avg_bar = zero(Tsol)

  # normal momentum
  mv_n = rho_hat*(dir[1]*u_hat + dir[2]*v_hat + dir[3]*w_hat)
  mv_n_bar = zero(Tres)

  # F[1]
  mv_n_bar += F_bar[1]

  # F[2]
  mv_n_bar  += u_hat*F_bar[2]
  u_hat_bar += mv_n*F_bar[2]
  p1_hat_bar += dir[1]*F_bar[2]

  # F[3]
  mv_n_bar += v_hat*F_bar[3]
  v_hat_bar += mv_n*F_bar[3]
  p1_hat_bar += dir[2]*F_bar[3]

  # F[4]
  mv_n_bar += w_hat*F_bar[4]
  w_hat_bar += mv_n*F_bar[4]
  p1_hat_bar += dir[3]*F_bar[4]

  # F[5]
  mv_n_bar += h_hat*F_bar[5]
  h_hat_bar += mv_n*F_bar[5]

  # mv_n
  rho_hat_bar += (dir[1]*u_hat + dir[2]*v_hat + dir[3]*w_hat)*mv_n_bar
  u_hat_bar += rho_hat*dir[1]*mv_n_bar
  v_hat_bar += rho_hat*dir[2]*mv_n_bar
  w_hat_bar += rho_hat*dir[3]*mv_n_bar

  # h_hat
  p2_hat_bar += gamma*h_hat_bar/(rho_hat*gamma_1)
  rho_hat_bar += (-gamma*p2_hat/(rho_hat*rho_hat*gamma_1))*h_hat_bar
  u_hat_bar += u_hat*h_hat_bar
  v_hat_bar += v_hat*h_hat_bar
  w_hat_bar += w_hat*h_hat_bar

  # p2_hat
  p2tmp = (gamma_1/(2*gamma))*(z5L + z5R)/(z1L + z1R)
  z5avg_bar +=  ((gamma + 1)/(2*gamma))*p2_hat_bar/z1avg
  z1avg_bar += -(((gamma + 1)/(2*gamma))*z5avg/(z1avg*z1avg))*p2_hat_bar
  z5L_bar += (gamma_1/(2*gamma))*p2_hat_bar/(z1L + z1R)
  z5R_bar += (gamma_1/(2*gamma))*p2_hat_bar/(z1L + z1R)
  z1L_bar += -p2tmp/(z1L + z1R)*p2_hat_bar
  z1R_bar += -p2tmp/(z1L + z1R)*p2_hat_bar

  # p1_hat 
  z5L_bar += p1_hat_bar/(z1L + z1R)
  z5R_bar += p1_hat_bar/(z1L + z1R)
  z1L_bar += -p1_hat/(z1L + z1R)*p1_hat_bar
  z1R_bar += -p1_hat/(z1L + z1R)*p1_hat_bar

  # w_hat
  z4L_bar += w_hat_bar/(z1L + z1R)
  z4R_bar += w_hat_bar/(z1L + z1R)
  z1L_bar += -w_hat/(z1L + z1R)*w_hat_bar
  z1R_bar += -w_hat/(z1L + z1R)*w_hat_bar

  # v_hat
  z3L_bar += v_hat_bar/(z1L + z1R)
  z3R_bar += v_hat_bar/(z1L + z1R)
  z1L_bar += -v_hat/(z1L + z1R)*v_hat_bar
  z1R_bar += -v_hat/(z1L + z1R)*v_hat_bar

  # u_hat
  z2L_bar += u_hat_bar/(z1L + z1R)
  z2R_bar += u_hat_bar/(z1L + z1R)
  z1L_bar += -u_hat/(z1L + z1R)*u_hat_bar
  z1R_bar += -u_hat/(z1L + z1R)*u_hat_bar

  # rho_hat
  z1L_bar += 0.5*z5avg*rho_hat_bar
  z1R_bar += 0.5*z5avg*rho_hat_bar
  z5avg_bar += 0.5*(z1L + z1R)*rho_hat_bar

  # log averages
  z5L_bar_tmp, z5R_bar_tmp = logavg_rev(z5L, z5R, z5avg_bar)
  z5L_bar += z5L_bar_tmp
  z5R_bar += z5R_bar_tmp

  z1L_bar_tmp, z1R_bar_tmp = logavg_rev(z1L, z1R, z1avg_bar)
  z1L_bar += z1L_bar_tmp
  z1R_bar += z1R_bar_tmp


  # z5L/R
  qL_bar[1] += 0.5*pL*z5L_bar/z5L; qR_bar[1] += 0.5*pR*z5R_bar/z5R
  pL_bar += 0.5*qL[1]*z5L_bar/z5L; pR_bar += 0.5*qR[1]*z5R_bar/z5R

  # z4L/R
  qL_bar[4] +=  z1L*z4L_bar/qL[1]; qR_bar[4] +=  z1R*z4R_bar/qR[1]
  qL_bar[1] += -z4L*z4L_bar/qL[1]; qR_bar[1] += -z4R*z4R_bar/qR[1]
  z1L_bar   +=  qL[4]*z4L_bar/qL[1]; z1R_bar +=  qR[4]*z4R_bar/qR[1]


  # z3L/R
  qL_bar[3] +=  z1L*z3L_bar/qL[1]; qR_bar[3] +=  z1R*z3R_bar/qR[1]
  qL_bar[1] += -z3L*z3L_bar/qL[1]; qR_bar[1] += -z3R*z3R_bar/qR[1]
  z1L_bar   +=  qL[3]*z3L_bar/qL[1]; z1R_bar +=  qR[3]*z3R_bar/qR[1]

  # z2L/R
  qL_bar[2] +=  z1L*z2L_bar/qL[1]; qR_bar[2] +=  z1R*z2R_bar/qR[1]
  qL_bar[1] += -z2L*z2L_bar/qL[1]; qR_bar[1] += -z2R*z2R_bar/qR[1]
  z1L_bar   +=  qL[2]*z2L_bar/qL[1]; z1R_bar +=  qR[2]*z2R_bar/qR[1]

  # z1L/R
  qL_bar[1] += (0.5/z1L)*z1L_bar/pL;  qR_bar[1] += (0.5/z1R)*z1R_bar/pR
  pL_bar += -(0.5/z1L)*(qL[1]/(pL*pL))*z1L_bar
  pR_bar += -(0.5/z1R)*(qR[1]/(pR*pR))*z1R_bar

  calcPressure_revq(params, qL, qL_bar, pL_bar)
  calcPressure_revq(params, qR, qR_bar, pR_bar)

  return nothing
end


"""
  3D, single direction version
"""
function calcEulerFlux_IR_revm(params::ParamType{3, :conservative},
                  qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                  aux_vars::AbstractArray{Tres},
                  nrm::AbstractArray{Tmsh, 1}, nrm_bar::AbstractArray{Tmsh, 1},
                  F_bar::AbstractArray{Tres, 1}) where {Tmsh, Tsol, Tres}
  @debug1 begin
    @assert length(qL) == length(qR)
    @assert size(F_bar, 1) == length(qL)
    @assert size(nrm, 1) == 3
    @assert size(nrm_bar, 1) == size(nrm, 1)
  end


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

#=
  mv_n = rho_hat*(nrm[1]*u_hat + nrm[2]*v_hat + nrm[3]*w_hat)  # normal momentum
  F[1] = mv_n
  F[2] = mv_n*u_hat + nrm[1]*p1_hat
  F[3] = mv_n*v_hat + nrm[2]*p1_hat
  F[4] = mv_n*w_hat + nrm[3]*p1_hat
  F[5] = mv_n*h_hat
=#
  # reverse sweep
    mv_n = rho_hat*(nrm[1]*u_hat + nrm[2]*v_hat + nrm[3]*w_hat)  # normal momentum
    mv_n_bar = zero(mv_n)

    # F[5]
    mv_n_bar += h_hat*F_bar[5]

    # F[4]
    mv_n_bar += w_hat*F_bar[4]
    nrm_bar[3] += p1_hat*F_bar[4]

    # F[3]
    mv_n_bar += v_hat*F_bar[3]
    nrm_bar[2] += p1_hat*F_bar[3]

    # F[2]
    mv_n_bar += u_hat*F_bar[2]
    nrm_bar[1] += p1_hat*F_bar[2]

    # F[1]
    mv_n_bar += F_bar[1]

    # mv_n
    nrm_bar[1] += rho_hat*u_hat*mv_n_bar
    nrm_bar[2] += rho_hat*v_hat*mv_n_bar
    nrm_bar[3] += rho_hat*w_hat*mv_n_bar

  return nothing
end



#------------------------------------------------------------------------------
# 3D, multi direction

"""
  3D, multi direction version
"""
function calcEulerFlux_IR_diff(params::ParamType{3, :conservative},
                   qL::AbstractArray{Tsol,1},
                   qR::AbstractArray{Tsol, 1},
                   aux_vars::AbstractArray{Tres, 1},
                   nrm::AbstractArray{Tmsh, 2},
                   FL_dot::AbstractArray{Tres, 3},
                   FR_dot::AbstractArray{Tres, 3}) where {Tmsh, Tsol, Tres}

  @debug1 begin
    @assert size(nrm, 1) == 3
    @assert size(nrm, 2) == 3
    @assert length(qL) == length(qR)
    @assert length(qL) == size(FL_dot, 1)
    @assert length(qL) == size(FL_dot, 2)
    @assert size(nrm, 2) == size(FL_dot, 3)
    @assert size(FL_dot, 1) == size(FR_dot, 1)
    @assert size(FL_dot, 2) == size(FR_dot, 2)
    @assert size(FR_dot, 3) == size(FR_dot, 3)
  end

  data = params.irfluxdata
  @unpack data pL_dot pR_dot
  @unpack data z1L_dot z2L_dot z3L_dot z4L_dot z5L_dot
  @unpack data z1R_dot z2R_dot z3R_dot z4R_dot z5R_dot

  gamma = params.gamma
  gamma_1 = params.gamma_1

  pL = calcPressure_diff(params, qL, pL_dot)
  pR = calcPressure_diff(params, qR, pR_dot)
  uL = qL[2]/qL[1]; uR = qR[2]/qR[1]
  vL = qL[3]/qL[1]; vR = qR[3]/qR[1]
  wL = qL[4]/qL[1]; wR = qR[4]/qR[1]
  z1L = sqrt(qL[1]/pL); z1R = sqrt(qR[1]/pR)
  z2L = z1L*uL; z2R = z1R*uR
  z3L = z1L*vL; z3R = z1R*vR
  z4L = z1L*wL; z4R = z1R*wR
  z5L = sqrt(qL[1]*pL); z5R = sqrt(qR[1]*pR)

  fastzero!(z1L_dot); fastzero!(z1R_dot)
  fastzero!(z2L_dot); fastzero!(z2R_dot)
  fastzero!(z3L_dot); fastzero!(z3R_dot)
  fastzero!(z4L_dot); fastzero!(z4R_dot)
  fastzero!(z5L_dot); fastzero!(z5R_dot)

  # differentiate with respect to q (not including chain rule terms for p)
  z1L_dot[1] = (0.5/z1L)*1/pL; z1R_dot[1] = (0.5/z1R)*1/pR

  z2L_dot[1] = -z2L/qL[1]; z2R_dot[1] = -z2R/qR[1]
  z2L_dot[2] =  z1L/qL[1]; z2R_dot[2] =  z1R/qR[1]

  z3L_dot[1] = -z3L/qL[1]; z3R_dot[1] = -z3R/qR[1]
  z3L_dot[3] =  z1L/qL[1]; z3R_dot[3] =  z1R/qR[1]

  z4L_dot[1] = -z4L/qL[1]; z4R_dot[1] = -z4R/qR[1]
  z4L_dot[4] =  z1L/qL[1]; z4R_dot[4] =  z1R/qR[1]

  z5L_dot[1] =  (0.5/z5L)*pL; z5R_dot[1] = (0.5/z5R)*pR

  # do the pressure/z1L related terms
  @simd for i=1:5
    z1L_dot[i] += (0.5/z1L)*(-qL[1]/(pL*pL))*pL_dot[i]
    z1R_dot[i] += (0.5/z1R)*(-qR[1]/(pR*pR))*pR_dot[i]

    z2L_dot[i] += uL*z1L_dot[i]
    z2R_dot[i] += uR*z1R_dot[i]

    z3L_dot[i] += vL*z1L_dot[i]
    z3R_dot[i] += vR*z1R_dot[i]

    z4L_dot[i] += wL*z1L_dot[i]
    z4R_dot[i] += wR*z1R_dot[i]

    z5L_dot[i] += (0.5/z5L)*qL[1]*pL_dot[i]
    z5R_dot[i] += (0.5/z5R)*qR[1]*pR_dot[i]
  end

  @unpack data avgdata z5avg_dotL z5avg_dotR z1avg_dotL z1avg_dotR
  @unpack data rho_hat_dotL rho_hat_dotR u_hat_dotL u_hat_dotR
  @unpack data v_hat_dotL v_hat_dotR w_hat_dotL w_hat_dotR 
  @unpack data p1_hat_dotL p1_hat_dotR h_hat_dotL h_hat_dotR

  # z4avg_dotL/r, z1avg_dotL/r, rho_hat, u_hat, v_hat, p1_hat, p2_hat
  z5avg = logavg_diff(avgdata, z5L, z5L_dot, z5R, z5R_dot, z5avg_dotL, z5avg_dotR)
  z1avg = logavg_diff(avgdata, z1L, z1L_dot, z1R, z1R_dot, z1avg_dotL, z1avg_dotR)

  z1 = z1L + z1R
  z1inv = 1/z1
  rho_hat = 0.5*(z1L + z1R)*z5avg
  u_hat = (z2L + z2R)*z1inv
  v_hat = (z3L + z3R)*z1inv
  w_hat = (z4L + z4R)*z1inv
  p1_hat = (z5L + z5R)*z1inv
  p2_hat = ((gamma + 1)/(2*gamma) )*z5avg/z1avg + ( gamma_1/(2*gamma) )*p1_hat
  h_hat = gamma*p2_hat/(rho_hat*gamma_1) + 0.5*(u_hat*u_hat + v_hat*v_hat + w_hat*w_hat)

  @simd for i=1:5
    rho_hat_dotL[i] = 0.5*(z5avg*z1L_dot[i] + z1*z5avg_dotL[i])
    rho_hat_dotR[i] = 0.5*(z5avg*z1R_dot[i] + z1*z5avg_dotR[i])

    u_hat_dotL[i] = z2L_dot[i]*z1inv - u_hat*z1inv*z1L_dot[i]
    u_hat_dotR[i] = z2R_dot[i]*z1inv - u_hat*z1inv*z1R_dot[i]

    v_hat_dotL[i] = z3L_dot[i]*z1inv - v_hat*z1inv*z1L_dot[i]
    v_hat_dotR[i] = z3R_dot[i]*z1inv - v_hat*z1inv*z1R_dot[i]

    w_hat_dotL[i] = z4L_dot[i]*z1inv - w_hat*z1inv*z1L_dot[i]
    w_hat_dotR[i] = z4R_dot[i]*z1inv - w_hat*z1inv*z1R_dot[i]

    p1_hat_dotL[i] = z5L_dot[i]*z1inv - p1_hat*z1inv*z1L_dot[i]
    p1_hat_dotR[i] = z5R_dot[i]*z1inv - p1_hat*z1inv*z1R_dot[i]

    # p2_hat is an intermediate variable for h_hat below
    p2_hat_dotL = ((gamma + 1)/(2*gamma))*(z5avg_dotL[i]/z1avg +
                      -z5avg/(z1avg*z1avg)*z1avg_dotL[i]) + 
                      ( gamma_1/(2*gamma))*p1_hat_dotL[i]
    p2_hat_dotR = ((gamma + 1)/(2*gamma))*(z5avg_dotR[i]/z1avg +
                      -z5avg/(z1avg*z1avg)*z1avg_dotR[i]) +
                      ( gamma_1/(2*gamma))*p1_hat_dotR[i]

    h_hat_dotL[i] = (gamma/gamma_1)*(p2_hat_dotL/rho_hat +
                     -p2_hat/(rho_hat*rho_hat)*rho_hat_dotL[i]) +
                    u_hat*u_hat_dotL[i] + v_hat*v_hat_dotL[i] +
                    w_hat*w_hat_dotL[i]
    h_hat_dotR[i] = (gamma/gamma_1)*(p2_hat_dotR/rho_hat +
                     -p2_hat/(rho_hat*rho_hat)*rho_hat_dotR[i]) +
                      u_hat*u_hat_dotR[i] + v_hat*v_hat_dotR[i] + 
                      w_hat*w_hat_dotR[i]
  end

  for j=1:3
    nx = nrm[1, j]; ny = nrm[2, j]; nz = nrm[3, j]
    mv_n = rho_hat*(nx*u_hat + ny*v_hat + nz*w_hat)  # normal momentum
    #F[1, j] = mv_n
    #F[2, j] = mv_n*u_hat + nx*p1_hat
    #F[3, j] = mv_n*v_hat + ny*p1_hat
    #F[4, j] = mv_n*h_hat

    @simd for i=1:5
      mv_n_dotL = (nx*u_hat + ny*v_hat + nz*w_hat)*rho_hat_dotL[i] + 
                  rho_hat*(nx*u_hat_dotL[i] + ny*v_hat_dotL[i] +
                           nz*w_hat_dotL[i])
      mv_n_dotR = (nx*u_hat + ny*v_hat + nz*w_hat)*rho_hat_dotR[i] +
                  rho_hat*(nx*u_hat_dotR[i] + ny*v_hat_dotR[i] +
                           nz*w_hat_dotR[i])

      FL_dot[1, i, j] += mv_n_dotL
      FL_dot[2, i, j] += u_hat*mv_n_dotL + mv_n*u_hat_dotL[i] + 
                         nx*p1_hat_dotL[i]
      FL_dot[3, i, j] += v_hat*mv_n_dotL + mv_n*v_hat_dotL[i] +
                         ny*p1_hat_dotL[i]
      FL_dot[4, i, j] += w_hat*mv_n_dotL + mv_n*w_hat_dotL[i] +
                         nz*p1_hat_dotL[i]
      FL_dot[5, i, j] += h_hat*mv_n_dotL + mv_n*h_hat_dotL[i]
      
      FR_dot[1, i, j] += mv_n_dotR
      FR_dot[2, i, j] += u_hat*mv_n_dotR + mv_n*u_hat_dotR[i] +
                         nx*p1_hat_dotR[i]
      FR_dot[3, i, j] += v_hat*mv_n_dotR + mv_n*v_hat_dotR[i] +
                         ny*p1_hat_dotR[i]
      FR_dot[4, i, j] += w_hat*mv_n_dotR + mv_n*w_hat_dotR[i] +
                         nz*p1_hat_dotR[i]
      FR_dot[5, i, j] += h_hat*mv_n_dotR + mv_n*h_hat_dotR[i]
      
    end
  end

  return nothing
end

"""
  3D, multi direction version
"""
function calcEulerFlux_IR_revq(params::ParamType{3, :conservative},
                      qL::AbstractArray{Tsol,1}, qL_bar::AbstractArray{Tsol, 1},
                      qR::AbstractArray{Tsol, 1}, qR_bar::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tres}, dir::AbstractArray{Tmsh, 2},  
                      F_bar::AbstractArray{Tres, 2}) where {Tmsh, Tsol, Tres}


  @debug1 begin
    @assert length(qL) == length(qL_bar)
    @assert length(qR) == length(qR_bar)
    @assert size(F_bar, 1) == length(qL_bar)
    @assert size(dir, 1) == 3
    @assert size(F_bar, 2) == size(dir, 2)
  end

  gamma = params.gamma
  gamma_1 = params.gamma_1
  pL = calcPressure(params, qL); pR = calcPressure(params, qR)
  z1L = sqrt(qL[1]/pL); z1R = sqrt(qR[1]/pR)
  z2L = z1L*qL[2]/qL[1]; z2R = z1R*qR[2]/qR[1]
  z3L = z1L*qL[3]/qL[1]; z3R = z1R*qR[3]/qR[1]
  z4L = z1L*qL[4]/qL[1]; z4R = z1R*qR[4]/qR[1]
  z5L = sqrt(qL[1]*pL); z5R = sqrt(qR[1]*pR)

  z5avg = logavg(z5L, z5R)
  z1avg = logavg(z1L, z1R)
  rho_hat = 0.5*(z1L + z1R)*z5avg
  u_hat = (z2L + z2R)/(z1L + z1R)
  v_hat = (z3L + z3R)/(z1L + z1R)
  w_hat = (z4L + z4R)/(z1L + z1R)
  p1_hat = (z5L + z5R)/(z1L + z1R)
  p2_hat = ((gamma + 1)/(2*gamma) )*z5avg/z1avg + ( gamma_1/(2*gamma) )*(z5L + z5R)/(z1L + z1R)
  h_hat = gamma*p2_hat/(rho_hat*gamma_1) + 0.5*(u_hat*u_hat + v_hat*v_hat + w_hat*w_hat)

#=
  for i=1:2
    mv_n = rho_hat*(dir[1, i]*u_hat + dir[2, i]*v_hat)  # normal momentum
    F[1, i] = mv_n
    F[2, i] = mv_n*u_hat + dir[1, i]*p1_hat
    F[3, i] = mv_n*v_hat + dir[2, i]*p1_hat
    F[4, i] = mv_n*h_hat
  end
=#
  # reverse sweep
  rho_hat_bar = zero(Tsol)
  u_hat_bar = zero(Tsol); v_hat_bar = zero(Tsol); w_hat_bar = zero(Tsol)
  p1_hat_bar = zero(Tsol); p2_hat_bar = zero(Tsol); h_hat_bar = zero(Tsol);
  z5avg_bar = zero(Tsol); z1avg_bar = zero(Tsol)

  z1L_bar = zero(Tsol); z2L_bar = zero(Tsol); z3L_bar = zero(Tsol);
  z4L_bar = zero(Tsol); z5L_bar = zero(Tsol)
  z1R_bar = zero(Tsol); z2R_bar = zero(Tsol); z3R_bar = zero(Tsol);
  z4R_bar = zero(Tsol); z5R_bar = zero(Tsol)
  pL_bar = zero(Tsol); pR_bar = zero(Tsol); z1avg_bar = zero(Tsol); z5avg_bar = zero(Tsol)

  for i=1:3
    mv_n = rho_hat*(dir[1, i]*u_hat + dir[2, i]*v_hat + dir[3, i]*w_hat)  # normal momentum
    mv_n_bar = zero(Tres)

    # F[1, i]
    mv_n_bar += F_bar[1, i]

    # F[2, i]
    mv_n_bar  += u_hat*F_bar[2, i]
    u_hat_bar += mv_n*F_bar[2, i]
    p1_hat_bar += dir[1, i]*F_bar[2, i]

    # F[3, i]
    mv_n_bar += v_hat*F_bar[3, i]
    v_hat_bar += mv_n*F_bar[3, i]
    p1_hat_bar += dir[2, i]*F_bar[3, i]

    # F[4, i]
    mv_n_bar += w_hat*F_bar[4, i]
    w_hat_bar += mv_n*F_bar[4, i]
    p1_hat_bar += dir[3, i]*F_bar[4, i]

    # F[5, i]
    mv_n_bar += h_hat*F_bar[5, i]
    h_hat_bar += mv_n*F_bar[5, i]

    # mv_n
    rho_hat_bar += (dir[1, i]*u_hat + dir[2, i]*v_hat + dir[3, i]*w_hat)*mv_n_bar
    u_hat_bar += rho_hat*dir[1, i]*mv_n_bar
    v_hat_bar += rho_hat*dir[2, i]*mv_n_bar
    w_hat_bar += rho_hat*dir[3, i]*mv_n_bar
  end

  # h_hat
  p2_hat_bar += gamma*h_hat_bar/(rho_hat*gamma_1)
  rho_hat_bar += (-gamma*p2_hat/(rho_hat*rho_hat*gamma_1))*h_hat_bar
  u_hat_bar += u_hat*h_hat_bar
  v_hat_bar += v_hat*h_hat_bar
  w_hat_bar += w_hat*h_hat_bar

  # p2_hat
  p2tmp = (gamma_1/(2*gamma))*(z5L + z5R)/(z1L + z1R)
  z5avg_bar +=  ((gamma + 1)/(2*gamma))*p2_hat_bar/z1avg
  z1avg_bar += -(((gamma + 1)/(2*gamma))*z5avg/(z1avg*z1avg))*p2_hat_bar
  z5L_bar += (gamma_1/(2*gamma))*p2_hat_bar/(z1L + z1R)
  z5R_bar += (gamma_1/(2*gamma))*p2_hat_bar/(z1L + z1R)
  z1L_bar += -p2tmp/(z1L + z1R)*p2_hat_bar
  z1R_bar += -p2tmp/(z1L + z1R)*p2_hat_bar

  # p1_hat 
  z5L_bar += p1_hat_bar/(z1L + z1R)
  z5R_bar += p1_hat_bar/(z1L + z1R)
  z1L_bar += -p1_hat/(z1L + z1R)*p1_hat_bar
  z1R_bar += -p1_hat/(z1L + z1R)*p1_hat_bar

  # w_hat
  z4L_bar += w_hat_bar/(z1L + z1R)
  z4R_bar += w_hat_bar/(z1L + z1R)
  z1L_bar += -w_hat/(z1L + z1R)*w_hat_bar
  z1R_bar += -w_hat/(z1L + z1R)*w_hat_bar

  # v_hat
  z3L_bar += v_hat_bar/(z1L + z1R)
  z3R_bar += v_hat_bar/(z1L + z1R)
  z1L_bar += -v_hat/(z1L + z1R)*v_hat_bar
  z1R_bar += -v_hat/(z1L + z1R)*v_hat_bar

  # u_hat
  z2L_bar += u_hat_bar/(z1L + z1R)
  z2R_bar += u_hat_bar/(z1L + z1R)
  z1L_bar += -u_hat/(z1L + z1R)*u_hat_bar
  z1R_bar += -u_hat/(z1L + z1R)*u_hat_bar

  # rho_hat
  z1L_bar += 0.5*z5avg*rho_hat_bar
  z1R_bar += 0.5*z5avg*rho_hat_bar
  z5avg_bar += 0.5*(z1L + z1R)*rho_hat_bar

  # log averages
  z5L_bar_tmp, z5R_bar_tmp = logavg_rev(z5L, z5R, z5avg_bar)
  z5L_bar += z5L_bar_tmp
  z5R_bar += z5R_bar_tmp

  z1L_bar_tmp, z1R_bar_tmp = logavg_rev(z1L, z1R, z1avg_bar)
  z1L_bar += z1L_bar_tmp
  z1R_bar += z1R_bar_tmp


  # z5L/R
  qL_bar[1] += 0.5*pL*z5L_bar/z5L; qR_bar[1] += 0.5*pR*z5R_bar/z5R
  pL_bar += 0.5*qL[1]*z5L_bar/z5L; pR_bar += 0.5*qR[1]*z5R_bar/z5R

  # z4L/R
  qL_bar[4] +=  z1L*z4L_bar/qL[1]; qR_bar[4] +=  z1R*z4R_bar/qR[1]
  qL_bar[1] += -z4L*z4L_bar/qL[1]; qR_bar[1] += -z4R*z4R_bar/qR[1]
  z1L_bar   +=  qL[4]*z4L_bar/qL[1]; z1R_bar +=  qR[4]*z4R_bar/qR[1]


  # z3L/R
  qL_bar[3] +=  z1L*z3L_bar/qL[1]; qR_bar[3] +=  z1R*z3R_bar/qR[1]
  qL_bar[1] += -z3L*z3L_bar/qL[1]; qR_bar[1] += -z3R*z3R_bar/qR[1]
  z1L_bar   +=  qL[3]*z3L_bar/qL[1]; z1R_bar +=  qR[3]*z3R_bar/qR[1]

  # z2L/R
  qL_bar[2] +=  z1L*z2L_bar/qL[1]; qR_bar[2] +=  z1R*z2R_bar/qR[1]
  qL_bar[1] += -z2L*z2L_bar/qL[1]; qR_bar[1] += -z2R*z2R_bar/qR[1]
  z1L_bar   +=  qL[2]*z2L_bar/qL[1]; z1R_bar +=  qR[2]*z2R_bar/qR[1]

  # z1L/R
  qL_bar[1] += (0.5/z1L)*z1L_bar/pL;  qR_bar[1] += (0.5/z1R)*z1R_bar/pR
  pL_bar += -(0.5/z1L)*(qL[1]/(pL*pL))*z1L_bar
  pR_bar += -(0.5/z1R)*(qR[1]/(pR*pR))*z1R_bar

  calcPressure_revq(params, qL, qL_bar, pL_bar)
  calcPressure_revq(params, qR, qR_bar, pR_bar)

  return nothing
end


"""
  3D, multi direction version
"""
function calcEulerFlux_IR_revm(params::ParamType{3, :conservative},
                  qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                  aux_vars::AbstractArray{Tres},
                  nrm::AbstractArray{Tmsh, 2}, nrm_bar::AbstractArray{Tmsh, 2},
                  F_bar::AbstractArray{Tres, 2}) where {Tmsh, Tsol, Tres}

  @debug1 begin
    @assert size(F_bar, 1) == length(qL)
    @assert size(nrm, 1) == 3
    @assert size(nrm, 2) == 3
    @assert size(nrm_bar, 1) == size(nrm_bar, 2)
    @assert size(nrm_bar, 2) == size(nrm_bar, 2)
    @assert size(F_bar, 2) == size(nrm, 2)
  end


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

#=
  for i=1:3
    mv_n = rho_hat*(nrm[1, i]*u_hat + nrm[2, i]*v_hat + nrm[3, i]*w_hat)  # normal momentum
    F[1, i] = mv_n
    F[2, i] = mv_n*u_hat + nrm[1, i]*p1_hat
    F[3, i] = mv_n*v_hat + nrm[2, i]*p1_hat
    F[4, i] = mv_n*w_hat + nrm[3, i]*p1_hat
    F[5, i] = mv_n*h_hat
  end
=#
  # reverse sweep
  for i=1:3
    mv_n = rho_hat*(nrm[1, i]*u_hat + nrm[2, i]*v_hat + nrm[3, i]*w_hat)  # normal momentum
    mv_n_bar = zero(mv_n)

    # F[5, i]
    mv_n_bar += h_hat*F_bar[5, i]

    # F[4, i]
    mv_n_bar += w_hat*F_bar[4, i]
    nrm_bar[3, i] += p1_hat*F_bar[4, i]

    # F[3, i]
    mv_n_bar += v_hat*F_bar[3, i]
    nrm_bar[2, i] += p1_hat*F_bar[3, i]

    # F[2, i]
    mv_n_bar += u_hat*F_bar[2, i]
    nrm_bar[1, i] += p1_hat*F_bar[2, i]

    # F[1, i]
    mv_n_bar += F_bar[1, i]

    # mv_n
    nrm_bar[1, i] += rho_hat*u_hat*mv_n_bar
    nrm_bar[2, i] += rho_hat*v_hat*mv_n_bar
    nrm_bar[3, i] += rho_hat*w_hat*mv_n_bar
  end

  return nothing
end




#------------------------------------------------------------------------------
# log average (needed by IR flux)

"""
  Differentiated version of logarithmic average (forward vector mode)

  **Inputs**

   * aL: left state (scalar)
   * aL_dot: vector of derivatives of aL (length arbitrary)
   * aR: right state (scalar)
   * aR_dot: vector of derivatives of aR (length same as aL_dot)

  **Inputs/Outputs

   * a_avg_dotL: d logavg/daL * aL_dot (forward mode propagation of aL_dot)
   * a_avg_dotR: d logavg/daR * aR_dot (forward mode propagation of aR_dot)

"""
function logavg_diff(data::LogAvgData, aL, aL_dot, aR, aR_dot, a_avg_dotL, a_avg_dotR)
# calculate the logarithmic average needed by the IR flux
  @assert length(aL_dot) == length(aR_dot)
  nd = length(aL_dot)

  # unpack args
  xi_dotL = data.xi_dotL
  xi_dotR = data.xi_dotR
  f_dotL = data.f_dotL
  f_dotR = data.f_dotR
  u_dotL = data.u_dotL
  u_dotR = data.u_dotR
  F_dotL = data.F_dotL
  F_dotR = data.F_dotR

  xi = aL/aR
  for i=1:nd
    xi_dotL[i] = aL_dot[i]/aR
    xi_dotR[i] = -aL/(aR*aR)*aR_dot[i]
  end

  f = (xi - 1)/(xi + 1)
  f_dotxi = 2/((xi + 1)*(xi + 1))
  for i=1:nd
    f_dotL[i] = f_dotxi*xi_dotL[i]
    f_dotR[i] = f_dotxi*xi_dotR[i]
  end
  u = f*f
  for i=1:nd
    u_dotL[i] = 2*f*f_dotL[i]
    u_dotR[i] = 2*f*f_dotR[i]
  end

  eps = 1e-3
  if u < eps
    F = @evalpoly( u, 1, 1/3, 1/5, 1/7, 1/9)
    F_dotu = @evalpoly(u, 1/3, 2/5, 3/7, 4/9)
    for i=1:nd
      F_dotL[i] = F_dotu*u_dotL[i]
      F_dotR[i] = F_dotu*u_dotR[i]
    end
#    F = 1.0 + u/3.0 + u*u/5.0 + u*u*u/7.0 + u*u*u*u/9.0
  else
    F = (log(xi)/2.0)/f
    for i=1:nd
      F_dotL[i] = 1/(xi*2*f)*xi_dotL[i] + (-F/f)*f_dotL[i]
      F_dotR[i] = 1/(xi*2*f)*xi_dotR[i] + (-F/f)*f_dotR[i]
    end
  end

  a_avg = (aL + aR)/(2*F)
  for i=1:nd
    a_avg_dotL[i] = aL_dot[i]/(2*F) + -(aL + aR)/(2*F*F)*F_dotL[i]
    a_avg_dotR[i] = aR_dot[i]/(2*F) + -(aL + aR)/(2*F*F)*F_dotR[i]
  end

  return a_avg
end


"""
  Reverse mode of [`logavg`](@ref)

  **Inputs**

   * aL: left variable (scalar)
   * aR: right variable (scalar)
   * a_avg_bar: seed value for the log average

  **Outputs**

   * aL_bar: adjoint part of aL
   * aR_bar: adjoint part of aR
"""
function logavg_rev(aL, aR, a_avg_bar)
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

  a_avg = (aL + aR)/(2*F)

  # reverse sweep
  aL_bar = zero(aL); aR_bar = zero(aR)
  F_bar = zero(F); u_bar = zero(u); f_bar = zero(F); xi_bar = zero(xi)

  # a_avg
  aL_bar += a_avg_bar/(2*F)
  aR_bar += a_avg_bar/(2*F)
  F_bar += -(a_avg/F)*a_avg_bar

  # F
  if u < eps
    F_dotu = @evalpoly(u, 1/3, 2/5, 3/7, 4/9)
    u_bar += F_dotu*F_bar
  else
    xi_bar += F_bar/(2*xi*f)
    f_bar += -(F/f)*F_bar
  end

  f_bar += 2*f*u_bar
  xi_bar += 2/( (xi + 1)*(xi + 1) )*f_bar
  aL_bar += xi_bar/aR
  aR_bar += -xi*xi_bar/aR

  return aL_bar, aR_bar
end



#------------------------------------------------------------------------------
# calcEulerFlux_* methods

function calcEulerFlux_IRSLF_diff(params::ParamType{Tdim, :conservative},
                   qL::AbstractArray{Tsol,1},
                   qR::AbstractArray{Tsol, 1},
                   aux_vars::AbstractArray{Tres, 1},
                   nrm::AbstractArray{Tmsh, 1},
                   FL_dot::AbstractArray{Tres, 2},
                   FR_dot::AbstractArray{Tres, 2}) where {Tmsh, Tsol, Tres, Tdim}


  kernel = params.entropy_lf_kernel

  calcEulerFlux_IR_diff(params, qL, qR, aux_vars, nrm, FL_dot, FR_dot)
  applyEntropyKernel_diagE_diff(params, kernel, qL, qR, aux_vars, nrm, FL_dot, FR_dot)

  return nothing
end

