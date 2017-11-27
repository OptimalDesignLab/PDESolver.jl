type nonslipBC <: BCType
end
# low level function
function call{Tmsh, Tsol, Tres}(obj::nonslipBC, 
                                params::ParamType,
                                q::AbstractArray{Tsol,1},  
                                aux_vars::AbstractArray{Tres, 1},  
                                x::AbstractArray{Tmsh,1}, 
                                nrm_xy::AbstractArray{Tmsh,1}, 
                                bndryflux::AbstractArray{Tres, 1})

  dim = length(nrm_xy)
	qg = params.qg
  # adiabatic wall
	qg[1] = q[1]
	qg[2:dim+1] = 0.0
	qg[dim+2] = q[dim+2]
  # isothermal wall
	# qg[1] = q[1]
	# rhoV2 = (q[2]*q[2] + q[3]*q[3])/q[1]
	# qg[2:dim+1] = 0.0
	# qg[dim+2] = q[4] - 0.5*rhoV2

	v_vals = params.v_vals
	convertFromNaturalToWorkingVars(params, qg, v_vals)
	# this is a problem: q is in conservative variables even if
	# params says we are using entropy variables
	calcEulerFlux(params, v_vals, aux_vars, nrm_xy, bndryflux)

	return nothing
end


type ExactChannelBC <: BCType
end
# low level function
function call{Tmsh, Tsol, Tres}(obj::ExactChannelBC, 
                                params::ParamType{3},
                                q::AbstractArray{Tsol,1},  
                                aux_vars::AbstractArray{Tres, 1},  
                                xyz::AbstractArray{Tmsh,1}, 
                                nrm_xy::AbstractArray{Tmsh,1}, 
                                bndryflux::AbstractArray{Tres, 1})

  sigma = 0.01
  gamma = params.gamma
  gamma_1 = params.gamma - 1
  aoa = params.aoa
  beta = params.sideslip_angle
  rhoInf = 1.0
  uInf = params.Ma * cos(beta) * cos(aoa)
  vInf = params.Ma * sin(beta) * -1
  wInf = params.Ma * cos(beta) * sin(aoa)
  TInf = 1.0
  x = xyz[1]
  y = xyz[2]
  z = xyz[3]

  rho = rhoInf * (1 + sigma*x*y*z)
  ux = sin(pi*x) + 1
  uy = sin(pi*y) + 1
  uz = sin(pi*z) + 1
  u  = (1 + sigma*ux * uy * uz )* uInf
  vx = sin(pi*x) + 1
  vy = sin(pi*y) + 1
  vz = sin(pi*z) + 1
  v  = (1 + sigma*vx * vy * vz )* vInf
  wx = sin(pi*x) + 1
  wy = sin(pi*y) + 1
  wz = sin(pi*z) + 1
  w  = (1 + sigma*wx * wy * wz) * wInf
  T  = TInf 

  if !params.isViscous
    u += 0.2 * uInf
    v += 0.2 * vInf
    w += 0.2 * wInf
  end

  qg = Array(Tsol, 5)
	qg[1] = rho
	qg[2] = rho*u
	qg[3] = rho*v
	qg[4] = rho*w
  qg[5] = T/(gamma * gamma_1) + 0.5 * (u*u + v*v + w*w)
  qg[5] *= rho

	v_vals = params.v_vals
	convertFromNaturalToWorkingVars(params, qg, v_vals)
	# this is a problem: q is in conservative variables even if
	# params says we are using entropy variables
	# calcEulerFlux(params, v_vals, aux_vars, [nx2, ny2], bndryflux)
  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)

	return nothing
end

type ExactChannelBC <: BCType
end
# low level function
function call{Tmsh, Tsol, Tres}(obj::ExactChannelBC, 
                                params::ParamType{2},
                                q::AbstractArray{Tsol,1},  
                                aux_vars::AbstractArray{Tres, 1},  
                                x::AbstractArray{Tmsh,1}, 
                                nrm_xy::AbstractArray{Tmsh,1}, 
                                bndryflux::AbstractArray{Tres, 1})

  pi = 3.14159265358979323846264338
  gamma = params.gamma
  gamma_1 = params.gamma - 1
  # qInf = zeros(Tsol, 4)
  # calcFreeStream(params, x, qInf)
  aoa = params.aoa
  rhoInf = 1.0
  uInf = params.Ma*cos(aoa)
  vInf = params.Ma*sin(aoa)
  TInf = 1.0
  rho = rhoInf
  # rho = rhoInf * (0.1*sin(2*pi*x[1]) + 0.1*x[2] +  1.0)
  # u   = uInf * (-4.0 * y * (y-1.0)) + 0.1*uInf
  # u   = uInf * (-4.0 * y * (y-1.0)) 
  ux = (0.1*sin(2*pi*x[1]) + 0.2) * uInf
  # uy = -4.0 * x[2] * (x[2]-1.0)
  uy = sin(pi*x[2]) 
  u  = ux * uy
  v  = vInf 
  T  = TInf 
  if !params.isViscous
    u += 0.2 * uInf
  end

  qg = Array(Tsol, 4)
	qg[1] = rho
	qg[2] = rho*u
	qg[3] = rho*v
  qg[4] = T/(gamma * gamma_1) + 0.5 * (u*u + v*v)
  qg[4] *= rho

	v_vals = params.v_vals
	convertFromNaturalToWorkingVars(params, qg, v_vals)
	# this is a problem: q is in conservative variables even if
	# params says we are using entropy variables
  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)

	return nothing
end

type zeroPressGradientBC <: BCType
end

# low level function
function call{Tmsh, Tsol, Tres}(obj::zeroPressGradientBC,
                                params::ParamType,
                                q::AbstractArray{Tsol,1},
                                aux_vars::AbstractArray{Tres, 1},
                                x::AbstractArray{Tmsh,1},
                                nrm_xy::AbstractArray{Tmsh,1},
                                bndryflux::AbstractArray{Tres, 1})


  dim = length(nrm_xy)

	gamma = params.gamma
	gamma_1 = params.gamma_1
	qg = params.qg
	dim = 2
  rhoV2 = (norm(view(q, 2:dim+1))) / q[1]
	# rhoV2 = (q[2]*q[2] + q[3]*q[3]) / q[1]
	pinf = 1./gamma
	qg[1:dim+1] = q[1:dim+1]
	qg[dim+2] = pinf/gamma_1 + 0.5*rhoV2

	v_vals = params.v_vals
	convertFromNaturalToWorkingVars(params, qg, v_vals)
	# this is a problem: q is in conservative variables even if
	# params says we are using entropy variables
	calcEulerFlux(params, v_vals, aux_vars, nrm_xy, bndryflux)

	return nothing
end
