mutable struct nonslipBC <: BCType
end
# low level function
function (obj::nonslipBC)(
              params::ParamType,
              q::AbstractArray{Tsol,1},  
              aux_vars::AbstractArray{Tres, 1},  
              x::AbstractArray{Tmsh,1}, 
              nrm_xy::AbstractArray{Tmsh,1}, 
              bndryflux::AbstractArray{Tres, 1}) where {Tmsh, Tsol, Tres}

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


mutable struct ExactChannelBC <: BCType
end
# low level function
function (obj::ExactChannelBC)(
              params::ParamType{3},
              q::AbstractArray{Tsol,1},  
              aux_vars::AbstractArray{Tres, 1},  
              xyz::AbstractArray{Tmsh,1}, 
              nrm_xy::AbstractArray{Tmsh,1}, 
              bndryflux::AbstractArray{Tres, 1}) where {Tmsh, Tsol, Tres}

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

  qg = Array{Tsol}(5)
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

# low level function
function (obj::ExactChannelBC)(
              params::ParamType{2, :conservative},
              q::AbstractArray{Tsol,1},  
              aux_vars::AbstractArray{Tres, 1},  
              x::AbstractArray{Tmsh,1}, 
              nrm_xy::AbstractArray{Tmsh,1}, 
              bndryflux::AbstractArray{Tres, 1}) where {Tmsh, Tsol, Tres}

  # functor ExactChannel takes varibales on multiple nodes, so we need to reshape some variables
  xy = reshape(x, length(x), 1)
  norm = reshape(nrm_xy, length(nrm_xy), 1)
  q_in = reshape(q, length(q), 1)
  q_bnd = zeros(Tsol, 4, 1)
  bnd_functor = ExactChannel()
  bnd_functor(q_in, xy, norm, params, q_bnd)
  qg = reshape(q_bnd, length(q_bnd))

	v_vals = params.v_vals
	convertFromNaturalToWorkingVars(params, qg, v_vals)
  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)

	return nothing
end

mutable struct zeroPressGradientBC <: BCType
end

# low level function
function (obj::zeroPressGradientBC)(
                                params::ParamType,
                                q::AbstractArray{Tsol,1},
                                aux_vars::AbstractArray{Tres, 1},
                                x::AbstractArray{Tmsh,1},
                                nrm_xy::AbstractArray{Tmsh,1},
                                bndryflux::AbstractArray{Tres, 1}) where {Tmsh, Tsol, Tres}


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
