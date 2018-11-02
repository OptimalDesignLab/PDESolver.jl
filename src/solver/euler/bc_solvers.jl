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

"""->
function RoeSolver(params::ParamType{2},
                   q::AbstractArray{Tsol,1},
                   qg::AbstractArray{Tsol, 1},
                   aux_vars::AbstractArray{Tres, 1},
                   nrm::AbstractArray{Tmsh,1},
                   flux::AbstractArray{Tres, 1}) where {Tmsh, Tsol, Tres}

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
  calcSAT(params, roe_vars, dq, nrm, sat)

  # calculate Euler flux in wall normal directiona
  # because edge numbering is rather arbitary, any memory access is likely to
  # be a cache miss, so we recalculate the Euler flux
  nrm2[1] = nx   # why are we assigning to nrm2?
  nrm2[2] = ny

  convertFromNaturalToWorkingVars(params, q, v_vals)
  calcEulerFlux(params, v_vals, aux_vars, nrm2, euler_flux)

  for i=1:4  # ArrayViews does not support flux[:] = .
    flux[i] = (sat_fac*sat[i] + euler_flux[i])
  end

  return nothing

end # ends the function eulerRoeSAT




"""
  The main Roe solver.  Populates `flux` with the computed flux.
"""
function RoeSolver(params::ParamType{3},
                   q::AbstractArray{Tsol,1},
                   qg::AbstractArray{Tsol, 1},
                   aux_vars::AbstractArray{Tres, 1},
                   nrm::AbstractArray{Tmsh,1},
                   flux::AbstractArray{Tres, 1}) where {Tmsh, Tsol, Tres}

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

  for i=1:length(dq)
    dq[i] = q[i] - qg[i]
  end

  roe_vars[1] = u
  roe_vars[2] = v
  roe_vars[3] = w
  roe_vars[4] = H

  calcSAT(params, roe_vars, dq, nrm, sat)

  # calculate Euler flux in wall normal directiona
  # because edge numbering is rather arbitary, any memory access is likely to
  # be a cache miss, so we recalculate the Euler flux
  convertFromNaturalToWorkingVars(params, q, v_vals)
  calcEulerFlux(params, v_vals, aux_vars, nrm, euler_flux)

  for i=1:5  # ArrayViews does not support flux[:] = .
    flux[i] = (sat_fac*sat[i] + euler_flux[i])
  end

  return nothing

end # ends the function eulerRoeSAT


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
function calcSAT(params::ParamType{2},
                 roe_vars::AbstractArray{Tsol, 1},
                 dq::AbstractArray{Tsol,1},
                 nrm::AbstractArray{Tmsh,1},
                 sat::AbstractArray{Tsol,1}) where {Tmsh, Tsol}
# roe_vars = [u, v, H] at Roe average 

  data = params.calcsatdata
  @unpack data E1dq E2dq

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


function calcSAT(params::ParamType{3},
                 roe_vars::AbstractArray{Tsol, 1},
                 dq::AbstractArray{Tsol,1},
                 nrm::AbstractArray{Tmsh,1},
                 sat::AbstractArray{Tsol,1}) where {Tsol, Tmsh}
  # roe_vars = [u, v, w, H] at Roe average 

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




"""
  Calculates the Lax-Friedrich flux function on the conservative variables
"""
function calcLFFlux(params::ParamType{Tdim, :conservative},
                    qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                    aux_vars::AbstractArray{Tsol, 1},
                    dir::AbstractArray{Tmsh, 1},  F::AbstractArray{Tres,1}) where {Tmsh, Tsol, Tres, Tdim}

  # compute Euler flux of left and right states
  data = params.lffluxdata
  @unpack data fluxL fluxR

  calcEulerFlux(params, qL, aux_vars, dir, fluxL)
  calcEulerFlux(params, qR, aux_vars, dir, fluxR)
  lambda_max = getLambdaMaxSimple(params, qL, qR, dir)

  for i=1:length(F)
    F[i] = 0.5*(fluxL[i] + fluxR[i] - lambda_max*(qR[i] - qL[i]))
  end

  return nothing
end


function calcEulerFlux_standard(
                      params::ParamType{2, :conservative},
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tsol, 1},
                      dir::AbstractArray{Tmsh, 1},  F::AbstractArray{Tres,1}) where {Tmsh, Tsol, Tres}
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
function calcEulerFlux_standard(
                      params::ParamType{2, :conservative},
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tsol, 1},
                      dir::AbstractArray{Tmsh, 2},  F::AbstractArray{Tres,2}) where {Tmsh, Tsol, Tres}
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


function calcEulerFlux_standard(
                      params::ParamType{3, :conservative},
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tres, 1},
                      dir::AbstractArray{Tmsh, 1},  F::AbstractArray{Tres,1}) where {Tmsh, Tsol, Tres}
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

function calcEulerFlux_standard(
                      params::ParamType{3, :conservative},
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tres, 1},
                      dir::AbstractArray{Tmsh, 2},  F::AbstractArray{Tres,2}) where {Tmsh, Tsol, Tres}
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
function calcEulerFlux_Ducros(params::ParamType{2, :conservative},
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tres},
                      dir::AbstractArray{Tmsh},  F::AbstractArray{Tres,1}) where {Tmsh, Tsol, Tres}
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

function calcEulerFlux_Ducros(params::ParamType{3, :conservative},
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tres},
                      dir::AbstractArray{Tmsh},  F::AbstractArray{Tres,1}) where {Tmsh, Tsol, Tres}
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
function calcEulerFlux_IR(params::ParamType{2, :conservative},
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tres},
                      dir::AbstractArray{Tmsh, 1},  F::AbstractArray{Tres,1}) where {Tmsh, Tsol, Tres}

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


  mn = rho_hat*(dir[1]*u_hat + dir[2]*v_hat)
  F[1] = mn
  F[2] = u_hat*mn + dir[1]*p1_hat
  F[3] = v_hat*mn + dir[2]*p1_hat
  F[4] = h_hat*mn
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
function calcEulerFlux_IR(params::ParamType{2, :conservative},
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tres},
                      dir::AbstractArray{Tmsh, 2},  F::AbstractArray{Tres,2}) where {Tmsh, Tsol, Tres}

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

  return nothing
end


function calcEulerFlux_IR(params::ParamType{3, :conservative},
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tres},
                      dir::AbstractArray{Tmsh, 1},  F::AbstractArray{Tres,1}) where {Tmsh, Tsol, Tres}

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
  
  mv_n = rho_hat*(dir[1]*u_hat + dir[2]*v_hat + dir[3]*w_hat)
  F[1] = mv_n
  F[2] = mv_n*u_hat + dir[1]*p1_hat
  F[3] = mv_n*v_hat + dir[2]*p1_hat
  F[4] = mv_n*w_hat + dir[3]*p1_hat
  F[5] = mv_n*h_hat

  return nothing
end


# multi-dimension version
function calcEulerFlux_IR(params::ParamType{3, :conservative},
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tres},
                      dir::AbstractArray{Tmsh, 2},  F::AbstractArray{Tres,2}) where {Tmsh, Tsol, Tres}

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


"""

  This function calculates the flux across an interface using the IR
  numerical flux function and a Lax-Friedrich type of entropy dissipation.

  Currently this is implemented for conservative variables only.

  This is the second method that takes in a normal vector directly.

  Inputs
    qL, qR: vectors conservative variables at left and right states
    aux_vars: aux_vars for qL
    nrm: a normal vector in xy space

  Inputs/Outputs
    F: vector to be updated with the result

  Alising restrictions:
    See getEntropyLFStab

"""
function calcEulerFlux_IRSLF(
                      params::ParamType{Tdim, :conservative},
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractVector{Tres},
                      dir::AbstractVector{Tmsh},  F::AbstractArray{Tres,1}) where {Tmsh, Tsol, Tres, Tdim}

  kernel = params.entropy_lf_kernel
  calcEulerFlux_IR(params, qL, qR, aux_vars, dir, F)
  applyEntropyKernel_diagE(params, kernel, qL, qR, aux_vars, dir, F)

  return nothing
end


function calcEulerFlux_IRSWF(
                      params::ParamType{Tdim, :conservative},
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractVector{Tres},
                      dir::AbstractVector{Tmsh},  F::AbstractArray{Tres,1}) where {Tmsh, Tsol, Tres, Tdim}

  kernel = params.entropy_lw2_kernel
  calcEulerFlux_IR(params, qL, qR, aux_vars, dir, F)
  applyEntropyKernel_diagE(params, kernel, qL, qR, aux_vars, dir, F)

  return nothing
end



"""
  Computes the logarithmic average required by the IR flux, in a numerically
  stable manner.  Uses one extra term in the polynomial as suggested by
  "On Discretely Entropy Conservative and Entropy Stable DG methods", Chan
  2018.

  **Inputs**

   * aL: left state
   * aR: right state

  **Outputs**

   * a_avg: logarithmic average of aL and aR
"""
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


"""
  Entropy conservative flux computed by solving an optimization problem
"""
function calcECOptFlux(params::ParamType{Tdim, :conservative},
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractVector{Tres},
                      nrm::AbstractVector{Tmsh},  F::AbstractArray{Tres,1}) where {Tmsh, Tsol, Tres, Tdim}

  numDofPerNode = length(qL)
  data = params.ecoptfluxdata
  @unpack data q_avg wL wR

  for i=1:numDofPerNode
    q_avg[i] = 0.5*(qL[i] + qR[i])
  end

  calcEulerFlux(params, q_avg, aux_vars, nrm, F)

  convertToIR(params, qL, wL)
  convertToIR(params, qR, wR)

  delta_psi = calcPotentialFluxIR(params, qL, nrm) -
              calcPotentialFluxIR(params, qR, nrm)

  num = zero(Tres)
  den = zero(Tres)

  for i=1:numDofPerNode
    delta_w = wL[i] - wR[i]
    num += delta_w*F[i]
    den += delta_w*delta_w
  end
  num -= delta_psi

  if den > 1e-13
    for i=1:numDofPerNode
      F[i] -= num*(wL[i] - wR[i])/den
    end
  end

  return nothing
end


function calcECOptFlux(params::ParamType{Tdim, :conservative},
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractVector{Tres},
                      nrm::AbstractMatrix{Tmsh},  F::AbstractArray{Tres,2}) where {Tmsh, Tsol, Tres, Tdim}

  #TODO: implement a specialized method that doesn't require converting to
  #      entropy variables repeatedly

  for i=1:size(nrm, 2)
    nrm_i = sview(nrm, :, i)
    F_i = sview(F, :, i)
    calcECOptFlux(params, qL, qR, aux_vars, nrm_i, F_i)
  end

  return nothing
end
