# differentiated version of functions in bc_solvers.jl


"""
  Differentiated version of the [`RoeSolver`](@ref).  Computes the jacobian
  of the flux with respect to `q` and `qg`.  Methods are available for 2D
  and 3D

  **Inputs**

   * params: ParamType, conservative variables only
   * q: vector of conservative variables for the left state
   * qg: vector of conservative variables for the right state
   * aux_vars: auxiliary variables for the left state
   * nrm: scaled normal vector in x-y space (outward wrt the element q lives on

  **Inputs/Outputs**

   * fluxL_dot: flux jacobian wrt `q`, numDofPerNode x numDofPerNode
   * fluxR_dot: flux jacobian wrt `qg`, numDofPerNode x numDofPerNode

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


  dq = params.v_vals2 # zeros(Tsol, 4)
  for i=1:length(dq)
    dq[i] = q[i] - qg[i]
  end

  # dq_dotL* = 1, dq_dotR* = -1, so omit them

  roe_vars = params.roe_vars
  roe_vars[1] = u
  roe_vars[2] = v
  roe_vars[3] = H

  roe_vars_dot = params.roe_vars_dot
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
  

  sat = params.sat_vals
#  sat_jacL = params.sat_jacL
#  sat_jacR = params.sat_jacR
  fill!(fluxL_dot, 0.0)
  fill!(fluxR_dot, 0.0)
  # pass in fluxL_dot and fluxR_dot here, then add the Euler flux
  # contribution below
  calcSAT_diff(params, roe_vars, roe_vars_dot,  dq, nrm, fluxL_dot, fluxR_dot)

#  calcSAT(params, nrm, dq, sat, u, v, H, use_efix)
  
  #euler_flux = params.flux_vals1
  euler_fluxjac = params.euler_fluxjac

  v_vals = params.q_vals
  nrm2 = params.nrm
  nrm2[1] = nx   # why are we assigning to nrm2?
  nrm2[2] = ny

#  convertFromNaturalToWorkingVars(params, q, v_vals)
  calcEulerFlux_diff(params, q, aux_vars, nrm2, euler_fluxjac)

  @simd for i=1:4
    @simd for j=1:4
      fluxL_dot[j, i] += euler_fluxjac[j, i]
    end
  end

  return nothing

end # ends the function RoeSolver

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


  dq = params.v_vals2 # zeros(Tsol, 4)
  for i=1:length(dq)
    dq[i] = q[i] - qg[i]
  end

  # dq_dotL* = 1, dq_dotR* = -1, so omit them

  roe_vars = params.roe_vars
  roe_vars[1] = u
  roe_vars[2] = v
  roe_vars[3] = w
  roe_vars[4] = H

  roe_vars_dot = params.roe_vars_dot
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
  

  sat = params.sat_vals
#  sat_jacL = params.sat_jacL
#  sat_jacR = params.sat_jacR
  fill!(fluxL_dot, 0.0)
  fill!(fluxR_dot, 0.0)
  # pass in fluxL_dot and fluxR_dot here, then add the Euler flux
  # contribution below
  calcSAT_diff(params, roe_vars, roe_vars_dot,  dq, nrm, fluxL_dot, fluxR_dot)

  euler_fluxjac = params.euler_fluxjac

  calcEulerFlux_diff(params, q, aux_vars, nrm, euler_fluxjac)

  @simd for i=1:5
    @simd for j=1:5
      fluxL_dot[j, i] += euler_fluxjac[j, i]
    end
  end

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

#  @printit u v H

#  @printit u_dotL1 u_dotR1 u_dotL2 u_dotR2 v_dotL1 v_dotR1 v_dotL3 v_dotR3 H_dotL1 H_dotR1 H_dotL2 H_dotR2 H_dotL3 H_dotR3 H_dotL4 H_dotR4

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


#  @printit Un_dotL1 phi_dotL1

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

#  @printit a_dotL1

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

#  @printit lambda1 lambda2 lambda3 lambda1_dotR1 lambda2_dotR1 lambda3_dotR1

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

#  @printit rhoA rhoA_dotR1 sat_Vn

  # Compute Eigen Values of the Flux Jacobian
  # The eigen values calculated above cannot be used directly. Near stagnation
  # points lambda3 approaches zero while near sonic lines lambda1 and lambda2
  # approach zero. This has a possibility of creating numerical difficulties.
  # As a result, the eigen values are limited by the following expressions.


  # see lambda1 expression below
  if absvalue(lambda1) > sat_Vn*rhoA
#    println("lambda1 is used")
    if lambda1 > 0
      fac = 1
    else
      fac = -1
    end
#    @printit fac tau lambda1 lambda2 lambda3
#    println("fac = ", fac)

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
#    println("not using lambda1")
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

#  println("after entropy fix")
#  @printit lambda1 lambda2 lambda3 lambda1_dotR1 lambda2_dotR1 lambda3_dotR1

                    
  dq1 = dq[1]
  dq2 = dq[2]
  dq3 = dq[3]
  dq4 = dq[4]

  #TODO: see if the sat values are needed
  # sat[1] = lambda3*dq1
  sat_jacL[1, 1] = lambda3 + dq1*lambda3_dotL1
  sat_jacL[1, 2] =         + dq1*lambda3_dotL2
  sat_jacL[1, 3] =         + dq1*lambda3_dotL3
  sat_jacL[1, 4] =         + dq1*lambda3_dotL4
  #TODO: zero out unused entries?

  sat_jacR[1, 1] = -lambda3 + dq1*lambda3_dotR1
  sat_jacR[1, 2] =          + dq1*lambda3_dotR2
  sat_jacR[1, 3] =          + dq1*lambda3_dotR3
  sat_jacR[1, 4] =          + dq1*lambda3_dotR4

  # sat[2] = lambda3*dq2
  sat_jacL[2, 1] =         + dq2*lambda3_dotL1
  sat_jacL[2, 2] = lambda3 + dq2*lambda3_dotL2
  sat_jacL[2, 3] =         + dq2*lambda3_dotL3
  sat_jacL[2, 4] =         + dq2*lambda3_dotL4

  sat_jacR[2, 1] =          + dq2*lambda3_dotR1
  sat_jacR[2, 2] = -lambda3 + dq2*lambda3_dotR2
  sat_jacR[2, 3] =          + dq2*lambda3_dotR3
  sat_jacR[2, 4] =          + dq2*lambda3_dotR4

  # sat[3] = lambda3*dq3
  sat_jacL[3, 1] =         + dq3*lambda3_dotL1
  sat_jacL[3, 2] =         + dq3*lambda3_dotL2
  sat_jacL[3, 3] = lambda3 + dq3*lambda3_dotL3
  sat_jacL[3, 4] =         + dq3*lambda3_dotL4

  sat_jacR[3, 1] =          + dq3*lambda3_dotR1
  sat_jacR[3, 2] =          + dq3*lambda3_dotR2
  sat_jacR[3, 3] = -lambda3 + dq3*lambda3_dotR3
  sat_jacR[3, 4] =          + dq3*lambda3_dotR4

  # sat[4] = lambda3*dq4
  sat_jacL[4, 1] =           dq4*lambda3_dotL1
  sat_jacL[4, 2] =           dq4*lambda3_dotL2
  sat_jacL[4, 3] =           dq4*lambda3_dotL3
  sat_jacL[4, 4] = lambda3 + dq4*lambda3_dotL4

  sat_jacR[4, 1] =            dq4*lambda3_dotR1
  sat_jacR[4, 2] =            dq4*lambda3_dotR2
  sat_jacR[4, 3] =            dq4*lambda3_dotR3
  sat_jacR[4, 4] = -lambda3 + dq4*lambda3_dotR4

#  @printit sat_jacR[2, 1]

  E1dq = params.res_vals1
  E2dq = params.res_vals2

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

#  @printit E1dq1_dotR1 E1dq2_dotR1 E1dq3_dotR1 E1dq4_dotR1 E2dq2_dotR1

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

#  @printit tmp1 tmp2 E1dq1_dotR1 E1dq[1] tmp1_dotR1 tmp2_dotR1 tmp3 E2dq[1]

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

#  println("after first summation")
#  @printit sat_jacR[2, 1]

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

#  @printit E1dq1_dotR1 E1dq2_dotR1 E1dq3_dotR1 E1dq4_dotR1 E2dq2_dotR1 E2dq3_dotR1 E2dq4_dotR1

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


#  println("after second summation")
#  @printit sat_jacR[2, 1]

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
  sat_jacL[1, 1] = lambda3 + dq1*lambda3_dotL1
  sat_jacL[1, 2] =         + dq1*lambda3_dotL2
  sat_jacL[1, 3] =         + dq1*lambda3_dotL3
  sat_jacL[1, 4] =         + dq1*lambda3_dotL4
  sat_jacL[1, 5] =         + dq1*lambda3_dotL5

  sat_jacR[1, 1] = -lambda3 + dq1*lambda3_dotR1
  sat_jacR[1, 2] =          + dq1*lambda3_dotR2
  sat_jacR[1, 3] =          + dq1*lambda3_dotR3
  sat_jacR[1, 4] =          + dq1*lambda3_dotR4
  sat_jacR[1, 5] =          + dq1*lambda3_dotR5

  # sat[2] = lambda3*dq2
  sat_jacL[2, 1] =         + dq2*lambda3_dotL1
  sat_jacL[2, 2] = lambda3 + dq2*lambda3_dotL2
  sat_jacL[2, 3] =         + dq2*lambda3_dotL3
  sat_jacL[2, 4] =         + dq2*lambda3_dotL4
  sat_jacL[2, 5] =         + dq2*lambda3_dotL5

  sat_jacR[2, 1] =          + dq2*lambda3_dotR1
  sat_jacR[2, 2] = -lambda3 + dq2*lambda3_dotR2
  sat_jacR[2, 3] =          + dq2*lambda3_dotR3
  sat_jacR[2, 4] =          + dq2*lambda3_dotR4
  sat_jacR[2, 5] =          + dq2*lambda3_dotR5

  # sat[3] = lambda3*dq3
  sat_jacL[3, 1] =         + dq3*lambda3_dotL1
  sat_jacL[3, 2] =         + dq3*lambda3_dotL2
  sat_jacL[3, 3] = lambda3 + dq3*lambda3_dotL3
  sat_jacL[3, 4] =         + dq3*lambda3_dotL4
  sat_jacL[3, 5] =         + dq3*lambda3_dotL5

  sat_jacR[3, 1] =          + dq3*lambda3_dotR1
  sat_jacR[3, 2] =          + dq3*lambda3_dotR2
  sat_jacR[3, 3] = -lambda3 + dq3*lambda3_dotR3
  sat_jacR[3, 4] =          + dq3*lambda3_dotR4
  sat_jacR[3, 5] =          + dq3*lambda3_dotR5

  sat_jacL[4, 1] =         + dq4*lambda3_dotL1
  sat_jacL[4, 2] =         + dq4*lambda3_dotL2
  sat_jacL[4, 3] =         + dq4*lambda3_dotL3
  sat_jacL[4, 4] = lambda3 + dq4*lambda3_dotL4
  sat_jacL[4, 5] =         + dq4*lambda3_dotL5

  sat_jacR[4, 1] =          + dq4*lambda3_dotR1
  sat_jacR[4, 2] =          + dq4*lambda3_dotR2
  sat_jacR[4, 3] =          + dq4*lambda3_dotR3
  sat_jacR[4, 4] = -lambda3 + dq4*lambda3_dotR4
  sat_jacR[4, 5] =          + dq4*lambda3_dotR5

  # sat[5] = lambda3*dq5
  sat_jacL[5, 1] =           dq5*lambda3_dotL1
  sat_jacL[5, 2] =           dq5*lambda3_dotL2
  sat_jacL[5, 3] =           dq5*lambda3_dotL3
  sat_jacL[5, 4] =           dq5*lambda3_dotL4
  sat_jacL[5, 5] = lambda3 + dq5*lambda3_dotL5

  sat_jacR[5, 1] =            dq5*lambda3_dotR1
  sat_jacR[5, 2] =            dq5*lambda3_dotR2
  sat_jacR[5, 3] =            dq5*lambda3_dotR3
  sat_jacR[5, 4] =            dq5*lambda3_dotR4
  sat_jacR[5, 5] = -lambda3 + dq5*lambda3_dotR5

  E1dq = params.res_vals1
  E2dq = params.res_vals2

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




function calcLFFlux_diff(
                      params::ParamType{Tdim, :conservative},
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tsol, 1},
                      dir::AbstractArray{Tmsh, 1},
                      F_dotL::AbstractMatrix{Tres}, F_dotR::AbstractMatrix{Tres}) where {Tmsh, Tsol, Tres, Tdim}

  numDofPerNode = length(qL)
#  lambda_dotL = zeros(Tres, numDofPerNode)
#  lambda_dotR = zeros(Tres, numDofPerNode)

  lambda_dotL = params.lambda_dotL
  lambda_dotR = params.lambda_dotR

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
    end
  end

  return nothing
end

"""
  This struct holds all the temporary arrays needed to calculate the IR flux
"""
struct IRFluxData{Tsol}
  pL_dot::Vector{Tsol}
  pR_dot::Vector{Tsol}
  z1L_dot::Vector{Tsol}
  z2L_dot::Vector{Tsol}
  z3L_dot::Vector{Tsol}
  z4L_dot::Vector{Tsol}

  z1R_dot::Vector{Tsol}
  z2R_dot::Vector{Tsol}
  z3R_dot::Vector{Tsol}
  z4R_dot::Vector{Tsol}

  z4avg_dotL::Vector{Tsol}
  z4avg_dotR::Vector{Tsol}
  z1avg_dotL::Vector{Tsol}
  z1avg_dotR::Vector{Tsol}

  rho_hat_dotL::Vector{Tsol}
  rho_hat_dotR::Vector{Tsol}
  u_hat_dotL::Vector{Tsol}
  u_hat_dotR::Vector{Tsol}
  v_hat_dotL::Vector{Tsol}
  v_hat_dotR::Vector{Tsol}
  p1_hat_dotL::Vector{Tsol}
  p1_hat_dotR::Vector{Tsol}
  h_hat_dotL::Vector{Tsol}
  h_hat_dotR::Vector{Tsol}
  logdata::LogAvgData{Tsol, Tsol}

  function IRFluxData{Tsol}(nd::Integer) where {Tsol}

    #TODO: consider making these views of an array to get spatial locality
    pL_dot = zeros(Tsol, nd)
    pR_dot = zeros(Tsol, nd)

    z1L_dot = zeros(Tsol, nd)
    z2L_dot = zeros(Tsol, nd)
    z3L_dot = zeros(Tsol, nd)
    z4L_dot = zeros(Tsol, nd)

    z1R_dot = zeros(Tsol, nd)
    z2R_dot = zeros(Tsol, nd)
    z3R_dot = zeros(Tsol, nd)
    z4R_dot = zeros(Tsol, nd)

    z4avg_dotL = zeros(Tsol, nd)
    z4avg_dotR = zeros(Tsol, nd)
    z1avg_dotL = zeros(Tsol, nd)
    z1avg_dotR = zeros(Tsol, nd)

    rho_hat_dotL = zeros(Tsol, nd)
    rho_hat_dotR = zeros(Tsol, nd)
    u_hat_dotL = zeros(Tsol, nd)
    u_hat_dotR = zeros(Tsol, nd)
    v_hat_dotL = zeros(Tsol, nd)
    v_hat_dotR = zeros(Tsol, nd)
    p1_hat_dotL = zeros(Tsol, nd)
    p1_hat_dotR = zeros(Tsol, nd)
    h_hat_dotL = zeros(Tsol, nd)
    h_hat_dotR = zeros(Tsol, nd)

    logdata = LogAvgData{Tsol, Tsol}(nd)

    return new(pL_dot, pR_dot, z1L_dot, z2L_dot, z2L_dot, z4L_dot,
               z1R_dot, z2R_dot, z3R_dot, z4R_dot,
               z4avg_dotL, z4avg_dotR, z1avg_dotL, z1avg_dotR,
               rho_hat_dotL, rho_hat_dotR, u_hat_dotL, u_hat_dotR,
               v_hat_dotL, v_hat_dotR, p1_hat_dotL, p1_hat_dotR,
               h_hat_dotL,h_hat_dotR,
               logdata)
  end
end


"""
  Differentiated version of the multi-dimensional version of the IR flux

  **Inputs**

   * Params: ParamType
   * qL: solution at the left node (numDofPerNode)
   * qg: solution at the right node (numDofPerNode)
   * aux_vars
   * nrm: mesh.dim x mesh.dim matrix of normal vectors, one per column
   * fluxL_dot: numDofPerNode x numDofPerNode x dim, jacobian of the flux
                with respect to q, in each direction
   * fluxR_dot  similar to fluxR_dot, but jacobian with respect to qg

"""
function calcEulerFlux_IR_diff(params::ParamType{2, :conservative},
                   qL::AbstractArray{Tsol,1},
                   qR::AbstractArray{Tsol, 1},
                   aux_vars::AbstractArray{Tres, 1},
                   nrm::AbstractArray{Tmsh, 2},
                   fluxL_dot::AbstractArray{Tres, 3},
                   fluxR_dot::AbstractArray{Tres, 3}) where {Tmsh, Tsol, Tres}


  # pL_dot, pR_dot, z1L_dot - z4L_dot, same for zR

  @unpack data z1L_dot z2L_dot z3L_dot z4L_dot z1R_dot z2R_dot z3R_dot z4R_dot

  gamma = params.gamma
  gamma_1 = params.gamma_1
  z1L = sqrt(qL[1]/pL); z1R = sqrt(qR[1]/pR)
  z2L = z1L*qL[2]/qL[1]; z2R = z1R*qR[2]/qR[1]
  z3L = z1L*qL[3]/qL[1]; z3R = z1R*qR[3]/qR[1]
  z4L = sqrt(qL[1]*pL); z4R = sqrt(qR[1]*pR)

  fastzero(z1L_dot); fastzero(z1R_dot)
  fastzero(z2L_dot); fastzero(z2R_dot)
  fastzero(z3L_dot); fastzero(z3R_dot)
  fastzero(z4L_dot); fastzero(z4L_dot)

  # differentiate with respect to q (not including chain rule terms for p)
  z1L_dot[1] = (-0.5/z1L)*1/pL; z1R_dot[1] = (-0.5/z1R)*1/pR

  z2L_dot[1] = -z2L/qL[1]; z2R_dot[1] = -z2R/qR[1]
  z2L_dot[2] =  z1L/qL[1]; z2R_dot[2] =  z1R/qR[1]

  z3L_dot[1] = -z3L/qL[1]; z3L_dot[1] = -z3R/qR[1]
  z3L_dot[3] =  z1L/qL[1]; z3R_dot[3] =  z1R/qR[1]

  z4L_dot[1] =  (-0.5/z4L)*pL; z4R_dot[1] = -(0.5/z4R)*pR

  # do the pressure/z1L related terms
  for i=1:4
    z1L_dot[i] += (-0.5/z1L)*(-qL[1]/(pL*pL))*pL_dot[i]
    z1R_dot[i] += (-0.5/z1R)*(-qR[1]/(pR*pR))*pR_dot[i]

    z2L_dot[i] += (qL[2]/qL[1])*z1L_dot[i]
    z2R_dot[i] += (qR[3]/qR[1])*z1R_dot[i]

    z3L_dot[i] += (qL[3]/qL[1])*z1L_dot[i]
    z3R_dot[i] += (qR[3]/qR[1])*z1R_dot[i]

    z4L_dot[i] += (-0.5/z4L)*qL[1]*pL_dot[i]
    z4R_dot[i] += (-0.5/z4R)*aR[1]*pR_dot[i]
  end

  @unpack data avgdata z4avg_dotL z4avg_dotR z1avg_dotL z1avg_dot
  @unpack data rho_hat_dotL rho_hat_dotR u_hat_dotL u_hat_dotR
  @unpack data v_hat_dotL v_hat_dotR p1_hat_dotL p1_hat_dotR
  @unpack data h_hat_dotL h_hat_dotR

  # z4avg_dotL/r, z1avg_dotL/r, rho_hat, u_hat, v_hat, p1_hat, p2_hat
  z4avg = logavg(avgdata, z4L, z4L_dot, z4R, z4R_dot, z4avg_dotL, z4avg_dotR)
  z1avg = logavg(avgdata, z1L, z1L_dot, z1R, z1R_dot, z1avg_dotL, z1avg_dotR)

  rho_hat = 0.5*(z1L + z1R)*z4avg
  u_hat = (z2L + z2R)/(z1L + z1R)
  v_hat = (z3L + z3R)/(z1L + z1R)
  p1_hat = (z4L + z4R)/(z1L + z1R)
  p2_hat = ((gamma + 1)/(2*gamma) )*z4avg/z1avg + ( gamma_1/(2*gamma) )*p1_hat
  h_hat = gamma*p2_hat/(rho_hat*gamma_1) + 0.5*(u_hat*u_hat + v_hat*v_hat)

  for i=1:4
    rho_hat_dotL[i] = 0.5*(z4avg*z1L_dot[i] + (z1L + z1R)*z4avg_dotL[i])
    rho_hat_dotR[i] = 0.5*(z4avg*z1R_dot[i] + (z1L + z1R)*z4avg_dotR[i])

    u_hat_dotL[i] = z2L_dot[i]/(z1L + z1R) - u_hat/(z1L + z1R)*z1L_dot[i]
    u_hat_dotR[i] = z2R_dot[i]/(z1L + z1R) - u_hat/(z1L + z1R)*z1R_dot[i]

    v_hat_dotL[i] = z3L_dot[i]/(z1L + z1R) - v_hat/(z1L + z1R)*z1L_dot[i]
    v_hat_dotR[i] = z3R_dot[i]/(z1L + z1R) - v_hat/(z1L + z1R)*z1R_dot[i]

    p1_hat_dotL[i] = z4L_dot[i]/(z1L + z1R) - p1_hat/(z1L + z1R)*z1L_dot[i]
    p1_hat_dotR[i] = z4L_dot[i]/(z1L + z1R) - p1_hat/(z1L + z1R)*z1R_dot[i]

    # p2_hat is an intermediate variable for h_hat below
    p2_hat_dotL = ((gamma + 1)/(2*gamma))*(z4avg_dotL[i]/z1avg +
                      -z4avg/(z1avg*z1avg)*z1avg_dotL[i]) + 
                      ( gamma_1/(2*gamma))*p1_hat_dotL[i]
    p2_hat_dotR = ((gamma + 1)/(2*gamma))*(z4avg_dotR[i]/z1avg +
                      -z4avg/(z1avg*z1avg)*z1avg_dotR[i]) +
                      ( gamma_1/(2*gamma))*p1_hat_dotR[i]

    h_hat_dotL[i] = (gamma/gamma_1)*(p2_hat_dotL[i]/rho_hat +
                     -p2_hat/(rho_hat*rho_hat)*rho_hat_dotL[i]) +
                     u_hat*u_hat_dotL[i] + v_hat*v_hat_dotL[i]
    h_hat_dotR[i] = (gamma/gamma_1)*(p2_hat_dotR[i]/rho_hat +
                     -p2_hat/(rho_hat*rho_hat)*rho_hat_dotR[i]) +
                      u_hat*u_hat_dotR[i] + v_hat*v_hat_dotR[i]
  end
                      

  for j=1:2
    mv_n = rho_hat*(dir[1, j]*u_hat + dir[2, j]*v_hat)  # normal momentum
    F[1, j] = mv_n
    F[2, j] = mv_n*u_hat + dir[1, j]*p1_hat
    F[3, j] = mv_n*v_hat + dir[2, j]*p1_hat
    F[4, j] = mv_n*h_hat

    for i=1:4
      mv_n_dotL = (dir[1, j]*u_hat + dir[2, j]*v_hat)*rho_hat_dotL[i] + 
                  rho_hat*(dir[1, j]*u_hat_dotL[i] + dir[2, j]*v_hat_dotL[j])
      mv_n_dotR = (dir[1, j]*u_hat + dir[2, j]*v_hat)*rho_hat_dor[i] +
                  rho_hat*(dir[1, j]*u_hat_dotR[i] + v_hat_doR[i])

      F_dotL[1, i, j] = mv_n_dotL
      F_dotL[2, i, j] = u_hat*mv_n_dotL + mv_n*u_hat_dotL[i] + 
                        dir[1, j]*p1_hat_dotL[i]
      F_dotL[3, i, j] = v_hat*mv_n_dotL + mv_n*v_hat_dotL[i] +
                        dir[2, j]*p1_hat_dotL[i]
      F_dotR[4, i, j] = h_hat*mv_n_dotL[i] + mv_n*h_hat_dotL[i]
      
      F_dotR[1, i, j] = mv_n_dotR
      F_dotR[2, i, j] = u_hat*mv_n_dotR + mv_n*u_hat_dotR[i] +
                        dir[1, j]*p1_hat_dotR[i]
      F_dotR[3, i, j] = v_hat*mv_n_dotR + mv_n*v_hat_dotR[i] +
                        dir[2, j]*p1_hat_dotR[i]
      F_dotR[4, i, j] = h_hat*mv_n_dotR + mv_n*h_hat_dotR[i]
      
    end
  end

  return nothing
end

"""
  Data needed by [`logavg_diff`](@ref)

  **Static Parameters**

   * Tl: datatype of left state
   * Tr: datatype of right state
"""
struct LogAvgData{Tl, Tr}
  xi_dotL::Vector{Tl}
  xi_dotR::Vector{Tr}
  f_dotL::Vector{Tl}
  f_dotR::Vector{Tr}
  u_dotL::Vector{Tl}
  u_dotR::Vector{Tr}
  F_dotL::Vector{Tl}
  F_dotR::Vector{Tr}

  function LogAvgData{Tl, Tr}(nd::Integer) where {Tl, Tr}
    xi_dotL = zeros(Tl, nd)
    xi_dotR = zeros(Tr, nd)
    f_dotL = zeros(Tl, nd)
    f_dotR = zeros(Tr, nd)
    u_dotL = zeros(Tl, nd)
    u_dotR = zeros(Tr, nd)
    F_dotL = zeros(Tl, nd)
    F_dotR = zeros(Tr, nd)

    return new(xi_dotL, xi_dotR, f_dotL, f_dotR, u_dotL, u_dotR, F_dotL, F_dotR)
  end
end




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
  @unpack data xi_dotL xi_dotR f_dotL f_dotR u_dotL u_dotR F_dotL F_dotR

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
