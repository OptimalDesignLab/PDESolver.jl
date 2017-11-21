#
#TODO: As far as we know `q`, `q_x`, `q_xx`, we can move
# the rest into a function, which is less error prone, and 
# will save a lot of space
#
@doc """
### EulerEquationMod.applySourceTerm

  This function updates eqn.res with the source term.  

  Inputs: 
    mesh
    sbp
    eqn
    opts
    src_func:  the functor that returns the value of the source term at a node
               This functor must have the signature:
               src_func(q, coords, params, t)
               where coords is a vector of length 2 or 3 containing the x and y 
               coordinates of the node, params is the ParamType, t is the 
               current time, and q is the vector to be populated with the 
               source term values.

  Outputs: none

  Aliasing restrictions: params.q_vals cannot be in use

"""->
function applySourceTerm(mesh,sbp, eqn, opts, src_func::SRCType)
#TODO: check that the k loop vectorizes
  weights = sbp.w
  q_vals = eqn.params.q_vals
  t = eqn.params.t

  for i=1:mesh.numEl
    jac_i = ro_sview(mesh.jac, :, i)
    res_i = sview(eqn.res, :, :, i)
    for j=1:mesh.numNodesPerElement
      coords_j = ro_sview(mesh.coords, :, j, i)
      src_func(q_vals, coords_j, eqn.params, t)
      fac = weights[j]/jac_i[j]
      for k=1:mesh.numDofPerNode
        res_i[k, j] += fac*q_vals[k]
      end
    end
  end

  return nothing
end


type SRCPolynomial <: SRCType
end
function call(obj::SRCPolynomial, 
              src::AbstractVector,
              xyz::AbstractVector, 
              params::ParamType{3}, 
              t)
  sigma = 0.01
	gamma = params.gamma
	gamma_1 = params.gamma_1
	aoa = params.aoa
  beta = eqn.params.sideslip_angle
	rhoInf = 1.0
  velInt = zeros(Tsol, 3)
  uInf = eqn.params.Ma * cos(beta) * cos(aoa)
  vInf = eqn.params.Ma * sin(beta) * -1
  wInf = eqn.params.Ma * cos(beta) * sin(aoa)
	TInf = 1.0

  x = xyz[1]
  y = xyz[2]
  z = xyz[3]

  rho = (x-x*x)*(y-y*y)* (z - z*z)
  u   = (x-x*x)*(y-y*y)* (z - z*z)
	v   = (x-x*x)*(y-y*y)* (z - z*z)
	w   = (x-x*x)*(y-y*y)* (z - z*z)
	T   = (x-x*x)*(y-y*y)* (z - z*z)

	rho_x = (1.0 - 2.0*x) * (y - y*y) * (z - z*z)
	u_x   = (1.0 - 2.0*x) * (y - y*y) * (z - z*z)
	v_x   = (1.0 - 2.0*x) * (y - y*y) * (z - z*z)
	w_x   = (1.0 - 2.0*x) * (y - y*y) * (z - z*z)
	T_x   = (1.0 - 2.0*x) * (y - y*y) * (z - z*z)

	rho_y = (x - x*x) * (1.0 - 2.0*y) * (z - z*z)
	u_y   = (x - x*x) * (1.0 - 2.0*y) * (z - z*z)
	v_y   = (x - x*x) * (1.0 - 2.0*y) * (z - z*z)
	w_y   = (x - x*x) * (1.0 - 2.0*y) * (z - z*z)
	T_y   = (x - x*x) * (1.0 - 2.0*y) * (z - z*z)

	rho_z = (x - x*x) * (y - y*y) * (1.0 - 2.0*z) 
	u_z   = (x - x*x) * (y - y*y) * (1.0 - 2.0*z) 
	v_z   = (x - x*x) * (y - y*y) * (1.0 - 2.0*z) 
	w_z   = (x - x*x) * (y - y*y) * (1.0 - 2.0*z) 
	T_z   = (x - x*x) * (y - y*y) * (1.0 - 2.0*z) 

	rho = (sigma*rho + 1.0)*rhoInf 
	u   = (sigma*u + 1.0)*uInf
	v   = (sigma*v + 1.0)*vInf
	w   = (sigma*w + 1.0)*wInf
	T   = (sigma*T + 1.0)*TInf

	rho_x = rho_x*rhoInf*sigma
	rho_y = rho_y*rhoInf*sigma
	rho_z = rho_z*rhoInf*sigma

	u_x = u_x*uInf*sigma
	u_y = u_y*uInf*sigma
	u_z = u_z*uInf*sigma

	v_x = v_x*vInf*sigma
	v_y = v_y*vInf*sigma
	v_z = v_z*vInf*sigma

	w_x = w_x*wInf*sigma
	w_y = w_y*wInf*sigma
	w_z = w_z*wInf*sigma

	T_x = T_x*TInf*sigma
	T_y = T_y*TInf*sigma
	T_z = T_z*TInf*sigma

	p   = rho*T/gamma
	E   = T/(gamma*gamma_1) + 0.5*(u*u + v*v + w*w)

	p_x   = 1.0/gamma*(rho_x*T + rho*T_x)
	p_y   = 1.0/gamma*(rho_y*T + rho*T_y)
	p_z   = 1.0/gamma*(rho_z*T + rho*T_z)
	E_x   = T_x/(gamma*gamma_1) + (u*u_x + v*v_x + w*w_x)	
	E_y   = T_y/(gamma*gamma_1) + (u*u_y + v*v_y + w*w_y)
	E_z   = T_z/(gamma*gamma_1) + (u*u_z + v*v_z + w*w_z)
	src[:] = 0.0
  #
  # contribution from inviscid terms
  #
  src[1]  = rho_x*u + rho*u_x + rho_y*v + rho*v_y + rho_z*w + rho*w_z
  
  src[2]  = rho_x*u*u + 2*rho*u*u_x + p_x
  src[2] += rho_y*u*v + rho*u_y*v + rho*u*v_y
  src[2] += rho_z*u*w + rho*u_z*w + rho*u*w_z

  src[3]  = rho_x*u*v + rho*u_x*v + rho*u*v_x
  src[3] += rho_y*v*v + 2*rho*v*v_y + p_y
  src[3] += rho_z*v*w + rho*v_z*w + rho*v*w_z

  src[4]  = rho_x*u*w + rho*u_x*w + rho*u*w_x
  src[4] += rho_y*v*w + rho*v_y*w + rho*v*w_y
  src[4] += rho_z*w*w + 2*rho*w*w_z + p_z

  src[5]  = rho_x*E*u + rho*E_x*u + rho*E*u_x + p_x*u + p*u_x
  src[5] += rho_y*E*v + rho*E_y*v + rho*E*v_y + p_y*v + p*v_y
  src[5] += rho_z*E*w + rho*E_z*w + rho*E*w_z + p_z*w + p*w_z

	if !params.isViscous 
		return nothing
	end
	#
	# contribution from viscous terms
	#
  u_xx = -2.0 * (y - y*y) * (z - z*z) 
	v_xx = -2.0 * (y - y*y) * (z - z*z) 
	w_xx = -2.0 * (y - y*y) * (z - z*z) 
	T_xx = -2.0 * (y - y*y) * (z - z*z) 

	u_yy = -2.0 * (x - x*x) * (z - z*z)
	v_yy = -2.0 * (x - x*x) * (z - z*z)
	w_yy = -2.0 * (x - x*x) * (z - z*z)
	T_yy = -2.0 * (x - x*x) * (z - z*z)

	u_zz = -2.0 * (x - x*x) * (y - y*y)
	v_zz = -2.0 * (x - x*x) * (y - y*y)
	w_zz = -2.0 * (x - x*x) * (y - y*y)
	T_zz = -2.0 * (x - x*x) * (y - y*y)

  u_xy = (1.0 - 2.0*x) * (1.0 - 2.0*y) * (z - z*z)
  v_xy = (1.0 - 2.0*x) * (1.0 - 2.0*y) * (z - z*z)
  w_xy = (1.0 - 2.0*x) * (1.0 - 2.0*y) * (z - z*z)

  u_xz = (1.0 - 2.0*x) * (1.0 - 2.0*z) * (y - y*y)
  v_xz = (1.0 - 2.0*x) * (1.0 - 2.0*z) * (y - y*y)
  w_xz = (1.0 - 2.0*x) * (1.0 - 2.0*z) * (y - y*y)

  u_yz = (1.0 - 2.0*y) * (1.0 - 2.0*z) * (x - x*x)
  v_yz = (1.0 - 2.0*y) * (1.0 - 2.0*z) * (x - x*x)
  w_yz = (1.0 - 2.0*y) * (1.0 - 2.0*z) * (x - x*x)

	u_xx *= sigma*uInf
	u_yy *= sigma*uInf
  u_zz *= sigma*uInf
	u_xy *= sigma*uInf
	u_xz *= sigma*uInf
	u_yz *= sigma*uInf
  u_yx = u_xy
  u_zx = u_xz
  u_zy = u_yz

	v_xx *= sigma*vInf
	v_yy *= sigma*vInf
	v_zz *= sigma*vInf
	v_xy *= sigma*vInf
	v_xz *= sigma*vInf
	v_yz *= sigma*vInf
  v_yx = v_xy
  v_zx = v_xz
  v_zy = v_yz

	w_xx *= sigma*wInf
	w_yy *= sigma*wInf
	w_zz *= sigma*wInf
	w_xy *= sigma*wInf
	w_xz *= sigma*wInf
	w_yz *= sigma*wInf
  w_yx = w_xy
  w_zx = w_xz
  w_zy = w_yz

	T_xx *= sigma*TInf
	T_yy *= sigma*TInf
	T_zz *= sigma*TInf

	muK = Array(typeof(coords[1]), 2)
	getMuK(T, muK)
	rmu = muK[1]
	rK = muK[2]

	txx = rmu * (4./3.*u_x - 2.0/3.0*v_y - 2.0/3.0*w_z)
	tyy = rmu * (4./3.*v_y - 2.0/3.0*u_x - 2.0/3.0*w_z)
	tzz = rmu * (4./3.*w_z - 2.0/3.0*u_x - 2.0/3.0*v_y)
	txy = rmu * (u_y + v_x) 
	txz = rmu * (u_z + w_x) 
	tyz = rmu * (w_y + v_z) 
  tyx = txy
  tzx = txz
  tzy = tyz


	txx_x = rmu*(4./3.*u_xx - 2.0/3.0*v_yx - 2.0/3.0*w_zx)
	txx_y = rmu*(4./3.*u_xy - 2.0/3.0*v_yy - 2.0/3.0*w_zy)
	txx_z = rmu*(4./3.*u_xz - 2.0/3.0*v_yz - 2.0/3.0*w_zz)

	tyy_x = rmu*(4./3.*v_yx - 2.0/3.0*u_xx - 2.0/3.0*w_zx)
	tyy_y = rmu*(4./3.*v_yy - 2.0/3.0*u_xy - 2.0/3.0*w_zy)
	tyy_z = rmu*(4./3.*v_yz - 2.0/3.0*u_xz - 2.0/3.0*w_zz)

	tzz_x = rmu*(4./3.*w_zx - 2.0/3.0*u_xx - 2.0/3.0*v_yx)
	tzz_y = rmu*(4./3.*w_zy - 2.0/3.0*u_xy - 2.0/3.0*v_yy)
	tzz_z = rmu*(4./3.*w_zz - 2.0/3.0*u_xz - 2.0/3.0*v_yz)

	txy_x = rmu*(u_yx + v_xx)
	txy_y = rmu*(u_yy + v_xy)
	txy_z = rmu*(u_yz + v_xz)
  txz_x = rmu*(u_zx + w_xx)
  txz_y = rmu*(u_zy + w_xy)
  txz_z = rmu*(u_zz + w_xz)
  tyz_x = rmu*(w_yx + v_zx)
  tyz_y = rmu*(w_yy + v_zy)
  tyz_z = rmu*(w_yz + v_zz)
	tyx_x = txy_x
	tyx_y = txy_y
	tyx_z = txy_z
	tzx_x = txz_x
	tzx_y = txz_y
	tzx_z = txz_z
	tzy_x = tyz_x
	tzy_y = tyz_y
	tzy_z = tyz_z
	
	Pr = 0.72
	c1 = params.Ma/params.Re
	c2 = c1/(Pr*gamma_1)
	src[2] -= c1*(txx_x + txy_y + txz_z)
	src[3] -= c1*(tyx_x + tyy_y + tyz_z)
  src[4] -= c1*(tzx_x + tzy_y + tzz_z)
	src[5] -= c1*(txx_x*u + txx*u_x + txy_x*v + txy*v_x + txz_x*w + txz*w_x) 
	src[5] -= c1*(tyx_y*u + tyx*u_y + tyy_y*v + tyy*v_y + tyz_y*w + tyz*w_y) 
	src[5] -= c1*(tzx_z*u + tzx*u_z + tzy_z*v + tzy*v_z + tzz_z*w + tzz*w_z) 
	src[5] -= c2*rK*(T_xx + T_yy + T_zz)

	return nothing
end

type SRCPolynomial <: SRCType
end
function call(obj::SRCPolynomial, 
              src::AbstractVector,
              coords::AbstractVector, 
              params::ParamType{2}, 
              t)
  sigma = 0.01
	gamma = params.gamma
	gamma_1 = params.gamma_1
	aoa = params.aoa
	rhoInf = 1.0
	uInf = params.Ma*cos(aoa)
	vInf = params.Ma*sin(aoa)
	TInf = 1.0
	x = coords[1]
	y = coords[2]

	rho = (x-x*x)*(y-y*y) 
	u   = (x-x*x)*(y-y*y)
	v   = (x-x*x)*(y-y*y)
	T   = (x-x*x)*(y-y*y)
	rho_x = (1.0 - 2.0*x) * (y - y*y) 
	u_x   = (1.0 - 2.0*x) * (y - y*y)
	v_x   = (1.0 - 2.0*x) * (y - y*y)
	T_x   = (1.0 - 2.0*x) * (y - y*y)
	rho_y = (x - x*x) * (1.0 - 2.0*y)
	u_y   = (x - x*x) * (1.0 - 2.0*y)
	v_y   = (x - x*x) * (1.0 - 2.0*y)
	T_y   = (x - x*x) * (1.0 - 2.0*y)

	rho = (sigma*rho + 1.0)*rhoInf 
	u   = (sigma*u + 1.0)*uInf
	v   = (sigma*v + 1.0)*vInf
	T   = (sigma*T + 1.0)*TInf
	rho_x = rho_x*rhoInf*sigma
	rho_y = rho_y*rhoInf*sigma
	u_x = u_x*uInf*sigma
	u_y = u_y*uInf*sigma
	v_x = v_x*vInf*sigma
	v_y = v_y*vInf*sigma
	T_x = T_x*TInf*sigma
	T_y = T_y*TInf*sigma

	p   = rho*T/gamma
	E   = T/(gamma*gamma_1) + 0.5*(u*u + v*v)

	p_x   = 1.0/gamma*(rho_x*T + rho*T_x)
	p_y   = 1.0/gamma*(rho_y*T + rho*T_y)
	E_x   = T_x/(gamma*gamma_1) + (u*u_x + v*v_x)	
	E_y   = T_y/(gamma*gamma_1) + (u*u_y + v*v_y)
	src[:] = 0.0
  #
  # contribution from inviscid terms
  #
  src[1]  = rho_x*u + rho*u_x + rho_y*v + rho*v_y 
  src[2]  = rho_x*u*u + 2*rho*u*u_x + p_x
  src[2] += rho_y*u*v + rho*u_y*v + rho*u*v_y
  src[3]  = rho_x*u*v + rho*u_x*v + rho*u*v_x
  src[3] += rho_y*v*v + 2*rho*v*v_y + p_y
  src[4]  = rho_x*E*u + rho*E_x*u + rho*E*u_x + p_x*u + p*u_x
  src[4] += rho_y*E*v + rho*E_y*v + rho*E*v_y + p_y*v + p*v_y

	if !params.isViscous 
		return nothing
	end
	#
	# contribution from viscous terms
	#
	u_xx = -2.0 * (y - y*y)
	v_xx = -2.0 * (y - y*y) 
	T_xx = -2.0 * (y - y*y) 
	u_yy = -2.0 * (x - x*x)
	v_yy = -2.0 * (x - x*x)
	T_yy = -2.0 * (x - x*x)
	u_xy = (1.0 - 2.0*x) * (1.0 - 2.0*y) 
	v_xy = (1.0 - 2.0*x) * (1.0 - 2.0*y)

	u_xx *= sigma*uInf
	u_xy *= sigma*uInf
	u_yy *= sigma*uInf

	v_xx *= sigma*vInf
	v_xy *= sigma*vInf
	v_yy *= sigma*vInf

	T_xx *= sigma*TInf
	T_yy *= sigma*TInf

	muK = Array(typeof(coords[1]), 2)
	getMuK(T, muK)
	rmu = muK[1]
	rK = muK[2]
	txx = rmu * (4./3.*u_x - 2.0/3.0*v_y)
	txy = rmu * (u_y + v_x) 
	tyx = rmu * (u_y + v_x) 
	tyy = rmu * (4./3.*v_y - 2.0/3.0*u_x)


	txx_x = rmu*(4./3.*u_xx - 2.0/3.0*v_xy)
	txx_y = rmu*(4./3.*u_xy - 2.0/3.0*v_yy)
	txy_x = rmu*(u_xy + v_xx)
	txy_y = rmu*(u_yy + v_xy)
	tyx_x = txy_x
	tyx_y = txy_y
	tyy_x = rmu*(4./3.*v_xy - 2.0/3.0*u_xx)
	tyy_y = rmu*(4./3.*v_yy - 2.0/3.0*u_xy)
	
	Pr = 0.72
	c1 = params.Ma/params.Re
	c2 = c1/(Pr*gamma_1)
	src[2] -= c1*(txx_x + txy_y)
	src[3] -= c1*(tyx_x + tyy_y)
	src[4] -= c1*(txx_x*u + txx*u_x + txy_x*v + txy*v_x) 
	src[4] -= c1*(txy_y*u + txy*u_y + tyy_y*v + tyy*v_y) 
	src[4] -= c2*rK*(T_xx + T_yy)

	return nothing
end


type SRCChannel <: SRCType
end
function call(obj::SRCChannel, 
              src::AbstractVector,
              coords::AbstractVector, 
              params::ParamType{2}, 
              t)
  pi = 3.14159265358979323846264338
	gamma = 1.4
	gamma_1 = gamma - 1.0
	aoa = params.aoa
	rhoInf = 1.0
	uInf = params.Ma*cos(aoa)
	vInf = params.Ma*sin(aoa)
	TInf = 1.0
	x = coords[1]
	y = coords[2]
	#
	# Exact solution in form of primitive variables
  rho = rhoInf
  # rho = rhoInf * (0.1*sin(2*pi*x) + 0.1*y +  1.0)
  # u   = uInf * (-4.0 * y * (y-1.0)) + 0.1*uInf
  ux = (0.1*sin(2*pi*x) + 0.2) * uInf
  # uy = -4.0 * y * (y-1.0)
  uy = sin(pi*y) 
  u  = ux * uy

  if !params.isViscous
    u += 0.2 * uInf
  end
  v  = vInf 
  T  = TInf 

  #
  # contribution from inviscid terms
  #
  rho_x = 0.0
  rho_y = 0.0
  # rho_x = rhoInf * 0.2*pi*cos(2*pi*x)
  # rho_y = rhoInf * 0.1
  u_x = uInf * 0.2*pi*cos(2*pi*x) * uy
  u_y = ux * pi * cos(pi*y) 
  v_x = 0.0 
  v_y = 0.0 
  T_x = 0.0 
	T_y = 0.0 

  p   = rho*T/gamma
	E   = T/(gamma*gamma_1) + 0.5*(u*u + v*v)

	p_x   = 1.0/gamma*(rho_x*T + rho*T_x)
	p_y   = 1.0/gamma*(rho_y*T + rho*T_y)
	E_x   = T_x/(gamma*gamma_1) + (u*u_x + v*v_x)	
	E_y   = T_y/(gamma*gamma_1) + (u*u_y + v*v_y)

	src[:] = 0.0
  src[1]  = rho_x*u + rho*u_x + rho_y*v + rho*v_y
  src[2]  = rho_x*u*u + 2*rho*u*u_x +  p_x
  src[2] += rho_y*u*v + rho*u_y*v + rho*u*v_y
  src[3]  = rho_x*u*v + rho*u_x*v + rho*u*v_x
  src[3] += rho_y*v*v + 2*rho*v*v_y + p_y
  src[4]  = rho_x*E*u + rho*E_x*u + rho*E*u_x + p_x*u + p*u_x
  src[4] += rho_y*E*v + rho*E_y*v + rho*E*v_y + p_y*v + p*v_y

  # println(real(src))

	if !params.isViscous 
    return nothing
	end
	
	#
	# contribution from viscous terms
	#
	muK = Array(typeof(coords[1]), 2)
	getMuK(T, muK)
	rmu = muK[1]
	rK  = muK[2]
  u_xx = uInf * -0.4 * pi * pi * sin(2*pi*x) * uy
  # u_xy = uInf * 0.2*pi* cos(2*pi*x) * (-8*y + 4)  
  u_xy = uInf * 0.2*pi* cos(2*pi*x) * pi* cos(pi*y)
  u_yy = ux * -pi*pi*sin(pi*y)
  # u_xx = uInf * -0.4 * pi * pi * sin(2*pi*x) * uy
  # u_xy = uInf * 0.2*pi* cos(2*pi*x) * (-8*y + 4)  
  # u_yy = ux * -8.0 
  u_yx = u_xy
  v_xx = 0.0 
  v_xy = 0.0
  v_yy = 0.0
  v_yx = v_xy
  T_xx = 0.0
  T_yy = 0.0
  
	txx = rmu * (4./3.*u_x - 2.0/3.0*v_y)
	txy = rmu * (u_y + v_x) 
	tyx = rmu * (u_y + v_x) 
	tyy = rmu * (4./3.*v_y - 2.0/3.0*u_x)

	txx_x = rmu*(4./3.*u_xx - 2.0/3.0*v_xy)
	txx_y = rmu*(4./3.*u_xy - 2.0/3.0*v_yy)
	txy_x = rmu*(u_xy + v_xx)
	txy_y = rmu*(u_yy + v_xy)
	tyx_x = txy_x
	tyx_y = txy_y
	tyy_x = rmu*(4./3.*v_xy - 2.0/3.0*u_xx)
	tyy_y = rmu*(4./3.*v_yy - 2.0/3.0*u_xy)

	Pr = 0.72
	c1 = params.Ma/params.Re
	c2 = c1/(Pr*gamma_1)
	src[2] -= c1*(txx_x + txy_y)
	src[3] -= c1*(tyx_x + tyy_y)
	src[4] -= c1*(txx_x*u + txx*u_x + txy_x*v + txy*v_x) 
	src[4] -= c1*(tyx_y*u + tyx*u_y + tyy_y*v + tyy*v_y) 
	src[4] -= c2*rK*(T_xx + T_yy)
	
	return nothing
end

#######################
#  ___________
# |           |
# |    ---    |
# |   |   |   |
# |    ---    |
# |           |
#  -----------
#
#######################
# inner block = nonslip wall
# outer block = freestream
type SRCDoubleSquare <: SRCType
end
function call(obj::SRCDoubleSquare, 
              src::AbstractVector,
              coords::AbstractVector, 
              params::ParamType{2}, 
              t)
	gamma = 1.4
	gamma_1 = gamma - 1.0
	aoa = params.aoa
	rhoInf = 1.0
	uInf = params.Ma*cos(aoa)
	vInf = params.Ma*sin(aoa)
	TInf = 1.0
	x = coords[1]
	y = coords[2]
	#
	# Exact solution in form of primitive variables
  #
  si = [0.5, 1.5]
  a = [-2.375, 16.875, -45.0, 55.0, -30.0, 6.0]
  gx    = 0.0
  gy    = 0.0
  gx_x  = 0.0
  gy_y  = 0.0
  gx_xx = 0.0
  gy_yy = 0.0

  if x >= si[1] && x < si[2]
    gx = a[1] + a[2]*x + a[3]*x*x + a[4]*x^3 + a[5]*x^4 + a[6]*x^5
    gx_x = a[2] + 2*a[3]*x + 3*a[4]*x*x + 4*a[5]*x*x*x + 5*a[6]*x^4
    gx_xx = 2*a[3] + 6*a[4]*x + 12*a[5]*x*x + 20*a[6]*x*x*x 
  elseif x >= si[2] 
    gx = 1.0
  elseif x <= -si[1] && x > -si[2]
    gx = a[1] - a[2]*x + a[3]*x*x - a[4]*x^3 + a[5]*x^4 - a[6]*x^5
    gx_x = -a[2] + 2*a[3]*x - 3*a[4]*x*x + 4*a[5]*x*x*x - 5*a[6]*x^4
    gx_xx = 2*a[3] - 6*a[4]*x + 12*a[5]*x*x - 20*a[6]*x*x*x 
  elseif x <= -si[2]
    gx = 1.0
  end  
  if y >= si[1] && y < si[2]
    gy = a[1] + a[2]*y + a[3]*y*y + a[4]*y^3 + a[5]*y^4 + a[6]*y^5
    gy_y = a[2] + 2*a[3]*y + 3*a[4]*y*y + 4*a[5]*y*y*y + 5*a[6]*y^4
    gy_yy = 2*a[3] + 6*a[4]*y + 12*a[5]*y*y + 20*a[6]*y*y*y 
  elseif y >= si[2] 
    gy = 1.0
  elseif y <= -si[1] && y > -si[2]
    gy = a[1] - a[2]*y + a[3]*y*y - a[4]*y^3 + a[5]*y^4 - a[6]*y^5
    gy_y = -a[2] + 2*a[3]*y - 3*a[4]*y*y + 4*a[5]*y*y*y - 5*a[6]*y^4
    gy_yy = 2*a[3] - 6*a[4]*y + 12*a[5]*y*y - 20*a[6]*y*y*y 
  elseif y <= -si[2]
    gy = 1.0
  end  

  rho = rhoInf
  u   = uInf * (gx + gy - gx*gy) 
  v   = vInf * (gx + gy - gx*gy) 
  T   = TInf

  #
  # contribution from inviscid terms
  #
  rho_x = 0.0
  rho_y = 0.0
  u_x = uInf * (gx_x - gx_x * gy) 
  u_y = uInf * (gy_y - gx * gy_y) 
  v_x = vInf * (gx_x - gx_x * gy) 
  v_y = vInf * (gy_y - gx * gy_y) 
  T_x = 0.0 
	T_y = 0.0 

  p   = rho*T/gamma
	E   = T/(gamma*gamma_1) + 0.5*(u*u + v*v)

	p_x   = 1.0/gamma*(rho_x*T + rho*T_x)
	p_y   = 1.0/gamma*(rho_y*T + rho*T_y)
	E_x   = T_x/(gamma*gamma_1) + (u*u_x + v*v_x)	
	E_y   = T_y/(gamma*gamma_1) + (u*u_y + v*v_y)

	src[:] = 0.0
  src[1]  = rho_x*u + rho*u_x + rho_y*v + rho*v_y
  src[2]  = rho_x*u*u + 2*rho*u*u_x +  p_x
  src[2] += rho_y*u*v + rho*u_y*v + rho*u*v_y
  src[3]  = rho_x*u*v + rho*u_x*v + rho*u*v_x
  src[3] += rho_y*v*v + 2*rho*v*v_y + p_y
  src[4]  = rho_x*E*u + rho*E_x*u + rho*E*u_x + p_x*u + p*u_x
  src[4] += rho_y*E*v + rho*E_y*v + rho*E*v_y + p_y*v + p*v_y

	if !params.isViscous 
    return nothing
	end
	
	#
	# contribution from viscous terms
	#
	muK = Array(typeof(coords[1]), 2)
	getMuK(T, muK)
	rmu = muK[1]
	rK  = muK[2]
  u_xx = uInf * (gx_xx - gx_xx * gy) 
  u_xy = uInf * (-gx_x * gy_y) 
  u_yy = uInf * (gy_yy - gx * gy_yy)
  u_yx = u_xy
  v_xx = vInf * (gx_xx - gx_xx * gy) 
  v_xy = vInf * (-gx_x * gy_y) 
  v_yy = vInf * (gy_yy - gx * gy_yy)
  v_yx = v_xy
  T_xx = 0.0
  T_yy = 0.0
  
	txx = rmu * (4./3.*u_x - 2.0/3.0*v_y)
	txy = rmu * (u_y + v_x) 
	tyx = rmu * (u_y + v_x) 
	tyy = rmu * (4./3.*v_y - 2.0/3.0*u_x)

	txx_x = rmu*(4./3.*u_xx - 2.0/3.0*v_xy)
	txx_y = rmu*(4./3.*u_xy - 2.0/3.0*v_yy)
	txy_x = rmu*(u_xy + v_xx)
	txy_y = rmu*(u_yy + v_xy)
	tyx_x = txy_x
	tyx_y = txy_y
	tyy_x = rmu*(4./3.*v_xy - 2.0/3.0*u_xx)
	tyy_y = rmu*(4./3.*v_yy - 2.0/3.0*u_xy)

	Pr = 0.72
	c1 = params.Ma/params.Re
	c2 = c1/(Pr*gamma_1)
	src[2] -= c1*(txx_x + txy_y)
	src[3] -= c1*(tyx_x + tyy_y)
	src[4] -= c1*(txx_x*u + txx*u_x + txy_x*v + txy*v_x) 
	src[4] -= c1*(tyx_y*u + tyx*u_y + tyy_y*v + tyy*v_y) 
	src[4] -= c2*rK*(T_xx + T_yy)
	
	return nothing
end

type SRCTrigWall <: SRCType
end
function call(obj::SRCTrigWall, 
>>>>>>> 3d_wing2
              src::AbstractVector,
              coords::AbstractVector, 
              params::ParamType{2}, 
              t)
  sigma = 0.01
  pi = 3.14159265358979323846264338
	gamma = 1.4
	gamma_1 = gamma - 1.0
	aoa = params.aoa
	rhoInf = 1.0
	uInf = params.Ma*cos(aoa)
	vInf = params.Ma*sin(aoa)
	TInf = 1.0
	x = coords[1]*pi
	y = coords[2]*pi
	x2 = 2*x
	y2 = 2*y 
	x4 = 4*x
	y4 = 4*y
	sx  = sin(x)
	sy  = sin(y)
	cx  = cos(x)
	cy  = cos(y)
	sx2 = sin(x2)
	sy2 = sin(y2)
	cx2 = cos(x2)
	cy2 = cos(y2)
	sx4 = sin(x4)
	sy4 = sin(y4)
	cx4 = cos(x4)
	cy4 = cos(y4)
	#
	# Exact solution in form of primitive variables
	#
  rho = rhoInf * ( 1.0 + sigma*exp(x)*exp(y) )
  u   = uInf * sx2 * sy2 
  v   = vInf * sx2 * sy2 
  T   = TInf 
  rho_x = rhoInf * sigma * exp(x) * exp(y)
  rho_y = rhoInf * sigma * exp(x) * exp(y)
  u_x = uInf * 2*pi * cx2 * sy2
  u_y = uInf * 2*pi * sx2 * cy2 
  v_x = vInf * 2*pi * cx2 * sy2
  v_y = vInf * 2*pi * sx2 * cy2 
	T_x = 0  
	T_y = 0 


	p   = rho*T/gamma
	E   = T/(gamma*gamma_1) + 0.5*(u*u + v*v)

	p_x   = 1.0/gamma*(rho_x*T + rho*T_x)
	p_y   = 1.0/gamma*(rho_y*T + rho*T_y)
	E_x   = T_x/(gamma*gamma_1) + (u*u_x + v*v_x)	
	E_y   = T_y/(gamma*gamma_1) + (u*u_y + v*v_y)

	src[:] = 0.0
  #
  # contribution from inviscid terms
  #
  src[1]  = rho_x*u + rho*u_x + rho_y*v + rho*v_y
  src[2]  = rho_x*u*u + 2*rho*u*u_x +  p_x
  src[2] += rho_y*u*v + rho*u_y*v + rho*u*v_y
  src[3]  = rho_x*u*v + rho*u_x*v + rho*u*v_x
  src[3] += rho_y*v*v + 2*rho*v*v_y + p_y
  src[4]  = rho_x*E*u + rho*E_x*u + rho*E*u_x + p_x*u + p*u_x
  src[4] += rho_y*E*v + rho*E_y*v + rho*E*v_y + p_y*v + p*v_y

	if !params.isViscous 
		return nothing
	end
	
	#
	# contribution from viscous terms
	#
	muK = Array(typeof(coords[1]), 2)
	getMuK(T, muK)
	rmu = muK[1]
	rK  = muK[2]

	u_xx = -4*pi*pi * sx2 * sy2 * uInf
	u_xy =  4*pi*pi * cx2 * cy2 * uInf
	u_yy = -4*pi*pi * sx2 * sy2 * uInf
	v_xx = -4*pi*pi * sx2 * sy2 * vInf
	v_xy =  4*pi*pi * cx2 * cy2 * vInf
	v_yy = -4*pi*pi * sx2 * sy2 * vInf
	T_xx = 0.0 
	T_yy = 0.0 

	txx = rmu * (4./3.*u_x - 2.0/3.0*v_y)
	txy = rmu * (u_y + v_x) 
	tyx = rmu * (u_y + v_x) 
	tyy = rmu * (4./3.*v_y - 2.0/3.0*u_x)

	txx_x = rmu*(4./3.*u_xx - 2.0/3.0*v_xy)
	txx_y = rmu*(4./3.*u_xy - 2.0/3.0*v_yy)
	txy_x = rmu*(u_xy + v_xx)
	txy_y = rmu*(u_yy + v_xy)
	tyx_x = txy_x
	tyx_y = txy_y
	tyy_x = rmu*(4./3.*v_xy - 2.0/3.0*u_xx)
	tyy_y = rmu*(4./3.*v_yy - 2.0/3.0*u_xy)

	Pr = 0.72
	c1 = params.Ma/params.Re
	c2 = c1/(Pr*gamma_1)
	src[2] -= c1*(txx_x + txy_y)
	src[3] -= c1*(tyx_x + tyy_y)
	src[4] -= c1*(txx_x*u + txx*u_x + txy_x*v + txy*v_x) 
	src[4] -= c1*(tyx_y*u + tyx*u_y + tyy_y*v + tyy*v_y) 
	src[4] -= c2*rK*(T_xx + T_yy)
	
	return nothing
end

type SRCTrigonometric <: SRCType
end
function call(obj::SRCTrigonometric, 
              src::AbstractVector,
              xyz::AbstractVector, 
              params::ParamType{3}, 
              t)
  sigma = 0.01
	gamma = params.gamma
	gamma_1 = params.gamma_1
	aoa = params.aoa
  beta = eqn.params.sideslip_angle
	rhoInf = 1.0
  velInt = zeros(Tsol, 3)
  uInf = eqn.params.Ma * cos(beta) * cos(aoa)
  vInf = eqn.params.Ma * sin(beta) * -1
  wInf = eqn.params.Ma * cos(beta) * sin(aoa)
	TInf = 1.0

  xyz2 = 2 * pi * xyz
  xyz4 = 4 * pi * xyz
  sin_val_1 = sin.(xyz)
  cos_val_1 = cos.(xyz)
  sin_val_2 = sin.(xyz2)
  cos_val_2 = cos.(xyz2)
  sin_val_4 = sin.(xyz4)
  cos_val_4 = cos.(xyz4)

  rho = 0.125 * sin_val_2[1] * sin_val_2[2] * sin_val_2[3] 
  u   = 0.125 * sin_val_4[1] * sin_val_4[2] * sin_val_4[3]
  v   = 0.125 * sin_val_2[1] * sin_val_2[2] * sin_val_2[3]
  w   = 0.125 * sin_val_1[1] * sin_val_1[2] * sin_val_1[3] 
  T   = 0.125 * (1.0 - cos_val_4[1]) * (1.0 - cos_val_4[2]) * (1.0 - cos_val_4[3])

  rho_x = 0.25*pi * cos_val_2[1] * sin_val_2[2] * sin_val_2[3] 
  u_x   = 0.50*pi * cos_val_4[1] * sin_val_4[2] * sin_val_4[3] 
  v_x   = 0.25*pi * cos_val_2[1] * sin_val_2[2] * sin_val_2[3] 
  w_x   = 0.125*pi* cos_val_1[1] * sin_val_1[2] * sin_val_1[3] 
  T_x   = 0.50*pi * sin_val_4[1] * (1.0 - cos_val_4[2]) * (1.0 - cos_val_4[3])  

  rho_y = 0.25*pi * cos_val_2[2] * sin_val_2[1] * sin_val_2[3] 
  u_y   = 0.50*pi * cos_val_4[2] * sin_val_4[1] * sin_val_4[3] 
  v_y   = 0.25*pi * cos_val_2[2] * sin_val_2[1] * sin_val_2[3] 
  w_y   = 0.125*pi* cos_val_1[2] * sin_val_1[1] * sin_val_1[3] 
  T_y   = 0.50*pi * sin_val_4[2] * (1.0 - cos_val_4[1]) * (1.0 - cos_val_4[3])

  rho_z = 0.25*pi * cos_val_2[3] * sin_val_2[1] * sin_val_2[2] 
  u_z   = 0.50*pi * cos_val_4[3] * sin_val_4[1] * sin_val_4[2] 
  v_z   = 0.25*pi * cos_val_2[3] * sin_val_2[1] * sin_val_2[2] 
  w_z   = 0.125*pi* cos_val_1[3] * sin_val_1[1] * sin_val_1[2] 
  T_z   = 0.50*pi * sin_val_4[3] * (1.0 - cos_val_4[1]) * (1.0 - cos_val_4[2])


	rho = (sigma*rho + 1.0)*rhoInf 
	u   = (sigma*u + 1.0)*uInf
	v   = (sigma*v + 1.0)*vInf
	w   = (sigma*w + 1.0)*wInf
	T   = (sigma*T + 1.0)*TInf

	rho_x = rho_x*rhoInf*sigma
	rho_y = rho_y*rhoInf*sigma
	rho_z = rho_z*rhoInf*sigma

	u_x = u_x*uInf*sigma
	u_y = u_y*uInf*sigma
	u_z = u_z*uInf*sigma

	v_x = v_x*vInf*sigma
	v_y = v_y*vInf*sigma
	v_z = v_z*vInf*sigma

	w_x = w_x*wInf*sigma
	w_y = w_y*wInf*sigma
	w_z = w_z*wInf*sigma

	T_x = T_x*TInf*sigma
	T_y = T_y*TInf*sigma
	T_z = T_z*TInf*sigma

	p   = rho*T/gamma
	E   = T/(gamma*gamma_1) + 0.5*(u*u + v*v + w*w)

	p_x   = 1.0/gamma*(rho_x*T + rho*T_x)
	p_y   = 1.0/gamma*(rho_y*T + rho*T_y)
	p_z   = 1.0/gamma*(rho_z*T + rho*T_z)
	E_x   = T_x/(gamma*gamma_1) + (u*u_x + v*v_x + w*w_x)	
	E_y   = T_y/(gamma*gamma_1) + (u*u_y + v*v_y + w*w_y)
	E_z   = T_z/(gamma*gamma_1) + (u*u_z + v*v_z + w*w_z)
	src[:] = 0.0
  #
  # contribution from inviscid terms
  #
  src[1]  = rho_x*u + rho*u_x + rho_y*v + rho*v_y + rho_z*w + rho*w_z
  
  src[2]  = rho_x*u*u + 2*rho*u*u_x + p_x
  src[2] += rho_y*u*v + rho*u_y*v + rho*u*v_y
  src[2] += rho_z*u*w + rho*u_z*w + rho*u*w_z

  src[3]  = rho_x*u*v + rho*u_x*v + rho*u*v_x
  src[3] += rho_y*v*v + 2*rho*v*v_y + p_y
  src[3] += rho_z*v*w + rho*v_z*w + rho*v*w_z

  src[4]  = rho_x*u*w + rho*u_x*w + rho*u*w_x
  src[4] += rho_y*v*w + rho*v_y*w + rho*v*w_y
  src[4] += rho_z*w*w + 2*rho*w*w_z + p_z

  src[5]  = rho_x*E*u + rho*E_x*u + rho*E*u_x + p_x*u + p*u_x
  src[5] += rho_y*E*v + rho*E_y*v + rho*E*v_y + p_y*v + p*v_y
  src[5] += rho_z*E*w + rho*E_z*w + rho*E*w_z + p_z*w + p*w_z

	if !params.isViscous 
		return nothing
	end
	#
	# contribution from viscous terms
	#
  u_xx = -2.0*pi*pi  * sin_val_4[1] * sin_val_4[2] * sin_val_4[3]
  v_xx = -0.5*pi*pi  * sin_val_2[1] * sin_val_2[2] * sin_val_2[3] 
  w_xx = -0.125*pi*pi* sin_val_1[1] * sin_val_1[2] * sin_val_1[3] 
  T_xx =  2.0*pi*pi  * cos_val_4[1] * (1.0 - cos_val_4[2]) * (1.0 - cos_val_4[3])

  u_yy = -2.0*pi*pi  * sin_val_4[2] * sin_val_4[1] * sin_val_4[3]
  v_yy = -0.5*pi*pi  * sin_val_2[2] * sin_val_2[1] * sin_val_2[3] 
  w_yy = -0.125*pi*pi* sin_val_1[2] * sin_val_1[1] * sin_val_1[3] 
  T_yy =  2.0*pi*pi  * cos_val_4[2] * (1.0 - cos_val_4[1]) * (1.0 - cos_val_4[3])

  u_zz = -2.0*pi*pi  * sin_val_4[3] * sin_val_4[1] * sin_val_4[2]
  v_zz = -0.5*pi*pi  * sin_val_2[3] * sin_val_2[1] * sin_val_2[2] 
  w_zz = -0.125*pi*pi* sin_val_1[3] * sin_val_1[1] * sin_val_1[2] 
  T_zz =  2.0*pi*pi  * cos_val_4[3] * (1.0 - cos_val_4[1]) * (1.0 - cos_val_4[2])

  u_xy = 2.0*pi*pi  * cos_val_4[1] * cos_val_4[2] * sin_val_4[3] 
  v_xy = 0.5*pi*pi  * cos_val_2[1] * cos_val_2[2] * sin_val_2[3]
  w_xy = 0.125*pi*pi* cos_val_1[1] * cos_val_1[2] * sin_val_1[3]

  u_xz = 2.0*pi*pi  * cos_val_4[1] * cos_val_4[3] * sin_val_4[2] 
  v_xz = 0.5*pi*pi  * cos_val_2[1] * cos_val_2[3] * sin_val_2[2]
  w_xz = 0.125*pi*pi* cos_val_1[1] * cos_val_1[3] * sin_val_1[2]

  u_yz = 2.0*pi*pi  * cos_val_4[2] * cos_val_4[3] * sin_val_4[1] 
  v_yz = 0.5*pi*pi  * cos_val_2[2] * cos_val_2[3] * sin_val_2[1]
  w_yz = 0.125*pi*pi* cos_val_1[2] * cos_val_1[3] * sin_val_1[1]

	u_xx *= sigma*uInf
	u_yy *= sigma*uInf
  u_zz *= sigma*uInf
	u_xy *= sigma*uInf
	u_xz *= sigma*uInf
	u_yz *= sigma*uInf
  u_yx = u_xy
  u_zx = u_xz
  u_zy = u_yz

	v_xx *= sigma*vInf
	v_yy *= sigma*vInf
	v_zz *= sigma*vInf
	v_xy *= sigma*vInf
	v_xz *= sigma*vInf
	v_yz *= sigma*vInf
  v_yx = v_xy
  v_zx = v_xz
  v_zy = v_yz

	w_xx *= sigma*wInf
	w_yy *= sigma*wInf
	w_zz *= sigma*wInf
	w_xy *= sigma*wInf
	w_xz *= sigma*wInf
	w_yz *= sigma*wInf
  w_yx = w_xy
  w_zx = w_xz
  w_zy = w_yz

	T_xx *= sigma*TInf
	T_yy *= sigma*TInf
	T_zz *= sigma*TInf

	muK = Array(typeof(coords[1]), 2)
	getMuK(T, muK)
	rmu = muK[1]
	rK = muK[2]

	txx = rmu * (4./3.*u_x - 2.0/3.0*v_y - 2.0/3.0*w_z)
	tyy = rmu * (4./3.*v_y - 2.0/3.0*u_x - 2.0/3.0*w_z)
	tzz = rmu * (4./3.*w_z - 2.0/3.0*u_x - 2.0/3.0*v_y)
	txy = rmu * (u_y + v_x) 
	txz = rmu * (u_z + w_x) 
	tyz = rmu * (w_y + v_z) 
  tyx = txy
  tzx = txz
  tzy = tyz


	txx_x = rmu*(4./3.*u_xx - 2.0/3.0*v_yx - 2.0/3.0*w_zx)
	txx_y = rmu*(4./3.*u_xy - 2.0/3.0*v_yy - 2.0/3.0*w_zy)
	txx_z = rmu*(4./3.*u_xz - 2.0/3.0*v_yz - 2.0/3.0*w_zz)

	tyy_x = rmu*(4./3.*v_yx - 2.0/3.0*u_xx - 2.0/3.0*w_zx)
	tyy_y = rmu*(4./3.*v_yy - 2.0/3.0*u_xy - 2.0/3.0*w_zy)
	tyy_z = rmu*(4./3.*v_yz - 2.0/3.0*u_xz - 2.0/3.0*w_zz)

	tzz_x = rmu*(4./3.*w_zx - 2.0/3.0*u_xx - 2.0/3.0*v_yx)
	tzz_y = rmu*(4./3.*w_zy - 2.0/3.0*u_xy - 2.0/3.0*v_yy)
	tzz_z = rmu*(4./3.*w_zz - 2.0/3.0*u_xz - 2.0/3.0*v_yz)

	txy_x = rmu*(u_yx + v_xx)
	txy_y = rmu*(u_yy + v_xy)
	txy_z = rmu*(u_yz + v_xz)
  txz_x = rmu*(u_zx + w_xx)
  txz_y = rmu*(u_zy + w_xy)
  txz_z = rmu*(u_zz + w_xz)
  tyz_x = rmu*(w_yx + v_zx)
  tyz_y = rmu*(w_yy + v_zy)
  tyz_z = rmu*(w_yz + v_zz)
	tyx_x = txy_x
	tyx_y = txy_y
	tyx_z = txy_z
	tzx_x = txz_x
	tzx_y = txz_y
	tzx_z = txz_z
	tzy_x = tyz_x
	tzy_y = tyz_y
	tzy_z = tyz_z
	
	Pr = 0.72
	c1 = params.Ma/params.Re
	c2 = c1/(Pr*gamma_1)
	src[2] -= c1*(txx_x + txy_y + txz_z)
	src[3] -= c1*(tyx_x + tyy_y + tyz_z)
  src[4] -= c1*(tzx_x + tzy_y + tzz_z)
	src[5] -= c1*(txx_x*u + txx*u_x + txy_x*v + txy*v_x + txz_x*w + txz*w_x) 
	src[5] -= c1*(tyx_y*u + tyx*u_y + tyy_y*v + tyy*v_y + tyz_y*w + tyz*w_y) 
	src[5] -= c1*(tzx_z*u + tzx*u_z + tzy_z*v + tzy*v_z + tzz_z*w + tzz*w_z) 
	src[5] -= c2*rK*(T_xx + T_yy + T_zz)

	return nothing
end

type SRCTrigonometric <: SRCType
end
function call(obj::SRCTrigonometric, 
              src::AbstractVector,
              coords::AbstractVector, 
              params::ParamType{2}, 
              t)
  sigma = 0.01
  pi = 3.14159265358979323846264338
	gamma = 1.4
	gamma_1 = gamma - 1.0
	aoa = params.aoa
	rhoInf = 1.0
	uInf = params.Ma*cos(aoa)
	vInf = params.Ma*sin(aoa)
	TInf = 1.0
	x = coords[1]*pi
	y = coords[2]*pi
	x2 = 2*x
	y2 = 2*y 
	x4 = 4*x
	y4 = 4*y
	sx  = sin(x)
	sy  = sin(y)
	cx  = cos(x)
	cy  = cos(y)
	sx2 = sin(x2)
	sy2 = sin(y2)
	cx2 = cos(x2)
	cy2 = cos(y2)
	sx4 = sin(x4)
	sy4 = sin(y4)
	cx4 = cos(x4)
	cy4 = cos(y4)
	#
	# Exact solution in form of primitive variables
	#
	rho = 0.25 * sx2 * sy2
	u   = 0.25 * sx4 * sy4
	v   = 0.25 * (cx4  + 1.0) * (cy4 + 1.0)
	T   = 0.25 * (1.0 - cx4) * (1.0 - cy4)
	rho_x = 0.5*pi * cx2 * sy2
	rho_y = 0.5*pi * sx2 * cy2
	u_x = pi* cx4 * sy4
	u_y = pi* sx4 * cy4
	v_x = -pi* sx4 * (cy4 + 1.0)
	v_y = -pi* (cx4  + 1) * sy4
	T_x = pi* sx4 * (1.0 - cy4)
	T_y = pi* (1.0 - cx4) * sy4

	rho = (sigma*rho + 1.0)*rhoInf 
	u   = (sigma*u + 1.0)*uInf
	v   = (sigma*v + 1.0)*vInf
	T   = (sigma*T + 1.0)*TInf
	rho_x = rho_x*rhoInf*sigma
	rho_y = rho_y*rhoInf*sigma
	u_x = u_x*uInf*sigma
	u_y = u_y*uInf*sigma
	v_x = v_x*vInf*sigma
	v_y = v_y*vInf*sigma
	T_x = T_x*TInf*sigma
	T_y = T_y*TInf*sigma

	p   = rho*T/gamma
	E   = T/(gamma*gamma_1) + 0.5*(u*u + v*v)

	p_x   = 1.0/gamma*(rho_x*T + rho*T_x)
	p_y   = 1.0/gamma*(rho_y*T + rho*T_y)
	E_x   = T_x/(gamma*gamma_1) + (u*u_x + v*v_x)	
	E_y   = T_y/(gamma*gamma_1) + (u*u_y + v*v_y)

	src[:] = 0.0
  #
  # contribution from inviscid terms
  #
  src[1]  = rho_x*u + rho*u_x + rho_y*v + rho*v_y
  src[2]  = rho_x*u*u + 2*rho*u*u_x +  p_x
  src[2] += rho_y*u*v + rho*u_y*v + rho*u*v_y
  src[3]  = rho_x*u*v + rho*u_x*v + rho*u*v_x
  src[3] += rho_y*v*v + 2*rho*v*v_y + p_y
  src[4]  = rho_x*E*u + rho*E_x*u + rho*E*u_x + p_x*u + p*u_x
  src[4] += rho_y*E*v + rho*E_y*v + rho*E*v_y + p_y*v + p*v_y

	if !params.isViscous 
		return nothing
	end
	
	#
	# contribution from viscous terms
	#
	muK = Array(typeof(coords[1]), 2)
	getMuK(T, muK)
	rmu = muK[1]
	rK  = muK[2]
	u_xx = -4*pi*pi * sx4 * sy4
	u_xy =  4*pi*pi * cx4 * cy4
	u_yy = -4*pi*pi * sx4 * sy4
	v_xx = -4*pi*pi * cx4 * (cy4 + 1.0)
	v_xy =  4*pi*pi * sx4 * sy4
	v_yy = -4*pi*pi * (cx4 + 1) * cy4
	T_xx =  4*pi*pi * cx4 * (1.0 - cy4)
	T_yy =  4*pi*pi * (1.0 - cx4) * cy4
	u_xx = u_xx * uInf * sigma
	u_xy = u_xy * uInf * sigma
	u_yy = u_yy * uInf * sigma
	v_xx = v_xx * vInf * sigma
	v_xy = v_xy * vInf * sigma
	v_yy = v_yy * vInf * sigma
	T_xx = T_xx * TInf * sigma
	T_yy = T_yy * TInf * sigma

	txx = rmu * (4./3.*u_x - 2.0/3.0*v_y)
	txy = rmu * (u_y + v_x) 
	tyx = rmu * (u_y + v_x) 
	tyy = rmu * (4./3.*v_y - 2.0/3.0*u_x)

	txx_x = rmu*(4./3.*u_xx - 2.0/3.0*v_xy)
	txx_y = rmu*(4./3.*u_xy - 2.0/3.0*v_yy)
	txy_x = rmu*(u_xy + v_xx)
	txy_y = rmu*(u_yy + v_xy)
	tyx_x = txy_x
	tyx_y = txy_y
	tyy_x = rmu*(4./3.*v_xy - 2.0/3.0*u_xx)
	tyy_y = rmu*(4./3.*v_yy - 2.0/3.0*u_xy)

	Pr = 0.72
	c1 = params.Ma/params.Re
	c2 = c1/(Pr*gamma_1)
	src[2] -= c1*(txx_x + txy_y)
	src[3] -= c1*(tyx_x + tyy_y)
	src[4] -= c1*(txx_x*u + txx*u_x + txy_x*v + txy*v_x) 
	src[4] -= c1*(tyx_y*u + tyx*u_y + tyy_y*v + tyy*v_y) 
	src[4] -= c2*rK*(T_xx + T_yy)
	
	return nothing
end


type SRC0 <: SRCType  # dummy source functor, it should nevery actually be called
end

@doc """
### AdvectionEquationMod.SRC0

  This is the zero source term.  This is the default of source term
  is specified
"""->
type SRCExp <: SRCType
end

function call(obj::SRCExp, q::AbstractVector, coords::AbstractVector, params::ParamType{2}, t)
  x = coords[1]
  y = coords[2]
  gamma_1 = params.gamma_1
  # a and b are parameters determine the scale and offset of the solution
  # the must be the same here and in the boundary condition
  af = 1/5  # a = 1/af
  b = 0.01

  q[1] = 2*y*af*exp(2*x*y*af + b) + 3*x*af*exp(3*x*y*af + b)
  q[2] = 3*y*af*exp(3*x*y*af + b) + 5*y*af*exp(5*x*y*af + b) + 4*x*af*exp(4*x*y*af+b)
  q[3] = 4*y*af*exp(4*x*y*af + b) + 10*x*af*exp(5*x*y*af + b)
  q[4] = 6*y*af*(1/gamma_1 + 1.5)*exp(6*x*y*af + b) + 2*y*af*exp(4*x*y*af + b) + 7*x*af*(1/gamma_1 + 1.5)*exp(7*x*y*af + b) + 2.5*x*af*exp(5*x*y*af + b)

  # make this work for other variables
  convertFromNaturalToWorkingVars(params, q, q)
end

# declare some global variables
# this is somewhat namespace polluting
global const MMSExp_a = 1/500
global const MMSExp_b = 0.01
global const MMSExp_c1 = 1
global const MMSExp_c2 = 2
global const MMSExp_c3 = 3
global const MMSExp_c4 = 4
global const MMSExp_c5 = 20
global const MMSExp_d1 = 1
global const MMSExp_d2 = 0.05
global const MMSExp_d3 = 0.15
global const MMSExp_d4 = 0.25
global const MMSExp_d5 = 1

function call(obj::SRCExp, q::AbstractVector, coords::AbstractVector, params::ParamType{3}, t)

  x = coords[1]
  y = coords[2]
  z = coords[3]

  # constant parameters
  gamma_1 = params.gamma_1
  a = MMSExp_a
  b = MMSExp_b
  c1 = MMSExp_c1 
  c2 = MMSExp_c2
  c3 = MMSExp_c3
  c4 = MMSExp_c4
  c5 = MMSExp_c5
  d1 = MMSExp_d1
  d2 = MMSExp_d2
  d3 = MMSExp_d3
  d4 = MMSExp_d4
  d5 = MMSExp_d5

  t2 = exp(b);
  t3 = a*c2*x*y*z;
  t4 = exp(t3);
  t5 = a*c4*x*y*z;
  t6 = exp(t5);
  t7 = c4*d4*t6*x*y;
  t8 = a*c3*x*y*z;
  t9 = exp(t8);
  t10 = c3*d3*t9*x*z;
  t11 = a*c5*x*y*z;
  t12 = exp(t11);
  t13 = 1.0/d1;
  t16 = a*c1*x*y*z;
  t14 = exp(-t16);
  t15 = c2*d2*t4*y*z;
  t17 = b-t16;
  t18 = exp(t17);
  t19 = d2*d2;
  t20 = a*c2*x*y*z*2.0;
  t21 = exp(t20);
  t22 = d3*d3;
  t23 = a*c3*x*y*z*2.0;
  t24 = exp(t23);
  t25 = d4*d4;
  t26 = a*c4*x*y*z*2.0;
  t27 = exp(t26);
  t28 = b+t11;
  t29 = exp(t28);
  t30 = c2*t19*t21;
  t31 = c3*t22*t24;
  t32 = c4*t25*t27;
  t33 = t30+t31+t32;
  t34 = t19*t21;
  t35 = t22*t24;
  t36 = t25*t27;
  t37 = t34+t35+t36;
  t38 = 1.0/gamma_1;
  t39 = c1-c4;
  t40 = exp(-a*t39*x*y*z);
  t41 = c1-c3;
  t42 = exp(-a*t41*x*y*z);
  t43 = d5*t29;
  t44 = d5*t29*t38;
  t45 = t13*t18*t37*(1.0/2.0);
  t46 = t43+t44+t45;
  t47 = c1-c2;
  t48 = exp(-a*t47*x*y*z);
  q[1] = a*t2*(t7+t10+t15);

  q[2] = a*c5*d5*t2*t12*y*z + a*d2*t2*t4*t13*t14*(t7+t10 - c1*d4*t6*x*y + c2*d4*t6*x*y - c1*d3*t9*x*z + c2*d3*t9*x*z - c1*d2*t4*y*z + c2*d2*t4*y*z*2.0);

  q[3] = a*c5*d5*t2*t12*x*z + a*d3*t2*t9*t13*t14*( (t7+t15 - c1*d4*t6*x*y + c3*d4*t6*x*y - c1*d3*t9*x*z) + (c3*d3*t9*x*z*2.0 - c1*d2*t4*y*z + c3*d2*t4*y*z) );

  q[4] = a*c5*d5*t2*t12*x*y + a*d4*t2*t6*t13*t14*( (t10+t15 - c1*d4*t6*x*y + c4*d4*t6*x*y*2.0 - c1*d3*t9*x*z) + (c4*d3*t9*x*z - c1*d2*t4*y*z + c4*d2*t4*y*z) );

  q[5] = d4*t13*t40*(a*t13*t18*t33*x*y + a*c5*d5*t29*x*y + a*c5*d5*t29*t38*x*y - a*c1*t13*t18*t37*x*y*(1.0/2.0)) + d3*t13*t42*(a*t13*t18*t33*x*z + a*c5*d5*t29*x*z + a*c5*d5*t29*t38*x*z - a*c1*t13*t18*t37*x*z*(1.0/2.0)) + d2*t13*t48*(a*t13*t18*t33*y*z + a*c5*d5*t29*y*z + a*c5*d5*t29*t38*y*z - a*c1*t13*t18*t37*y*z*(1.0/2.0)) - a*d4*t13*t39*t40*t46*x*y - a*d3*t13*t41*t42*t46*x*z - a*d2*t13*t46*t47*t48*y*z;

  return nothing
end

"""
  Functor for source term corresponding to ICPeriodicMMS
"""
type SRCPeriodicMMS <: SRCType
end

function call(obj::SRCPeriodicMMS, q::AbstractVector, coords::AbstractVector, 
              params::ParamType{2}, t)

  x = coords[1]
  y = coords[2]
  gamma_1 = params.gamma_1

  t4 = t*2.0;
  t2 = -t4+x+y;
  t3 = 3.141592653589793*t2;
  t5 = cos(t3);
  t6 = sin(t3);
  t7 = t6+1.5E1;
  t8 = 3.141592653589793*gamma_1*t5*t7*(1.0/5.0E1);
  q[1] = 0
  q[2] = t8;
  q[3] = t8;
  q[4] = 3.141592653589793*gamma_1*t5*t7*(1.0/2.5E1);

  return nothing
end

function call(obj::SRCPeriodicMMS, q::AbstractVector, coords::AbstractVector, 
              params::ParamType{3}, t)

  x = coords[1]
  y = coords[2]
  z = coords[3]
  gamma_1 = params.gamma_1
  gamma = params.gamma
#=
  t2 = x+y+z;
  t3 = 3.141592653589793*t2;
  t4 = cos(t3);
  t5 = sin(t3);
  t6 = t5*(1.0/1.0E1);
  t7 = t6+2.0;
  t8 = 3.141592653589793*gamma_1*t4*t7*(1.0/5.0);
=#
#  tmp = pi*gamma_1*cos(pi*(x + y + z))*(sin(pi*(x + y + z))/10 + 2)/5
#=
  q[1] = 0;
  q[2] = 0
  q[3] = 0
  q[4] = 0
  q[5] = 0;
=#

  t5 = t*2.0;
  t2 = -t5+x+y+z;
  t3 = 3.141592653589793*t2;
  t4 = cos(t3);
  t6 = gamma_1*2.5E1;
  t7 = sin(t3);
  t8 = gamma_1*t7*2.0;
  t9 = t6+t8+1.0E1;
  t10 = 3.141592653589793*t4*t9*(1.0/1.0E2);
  q[1] = 3.141592653589793*t4*(1.0/1.0E1);
  q[2] = t10;
  q[3] = t10;
  q[4] = t10;
  q[5] = 3.141592653589793*t4*(gamma_1*7.5E1+t7*2.0+gamma_1*t7*6.0+4.0E1)*(1.0/1.0E2);

  return nothing
end

include("source_viscous.jl")

@doc """
### EulerEquationMod.SRCDict

  Stores the functors that evaluate source terms at a node.  Every new 
  functor should be added to this dictonary

  All functors must have the signature:

  src_func(q, coords, params::ParamType, t)

  where coords is the vector of length 2 containing the x and y coordinates
  of the node, t is the current time, and q is the vector to be populated with
  the source term values.
"""->
global const SRCDict = Dict{ASCIIString, SRCType}(
"SRCExp" => SRCExp(),
"SRCPeriodicMMS" => SRCPeriodicMMS(),
"SRC0" => SRC0(),
"SRCTrigonometric" => SRCTrigonometric(),
"SRCTrigWall" => SRCTrigWall(),
"SRCDoubleSquare" => SRCDoubleSquare(),
"SRCPolynomial" => SRCPolynomial(),
"SRCChannel" => SRCChannel(),
)


@doc """
### EulerEquationMod.getSRCFunctors

  This function gets the functor specified by opts["SRCname"] and stores
  it to the equation object.  Currently one 1 source functor is allowed.

"""->
function getSRCFunctors(mesh::AbstractMesh, sbp::AbstractSBP, 
                        eqn::EulerData, opts)

  # currently we only allow 1 source functor
  eqn.src_func = SRCDict[opts["SRCname"]]
  println("using source term functor ", eqn.src_func)
  return nothing
end


