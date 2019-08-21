@doc """
### EulerEquationMod.calcIsentropicVortex

  This function calculates the isentropic vortex solution to the Euler
  equations at a node.  This function uses an inner radius of 1, an inner
  radius density of 2, an inner radius Mach number of 0.95, and an inner radius
  pressure of 1/gamma.  The denstiy as a function of radius r can be found in,
  for example,
  "Output Error Estimates for Summation-by-parts Finite-difference Schemes",
  JE Hicken.

  This function uses conservative variables regardless of the static parameter
  of params.

  Inputs:
    coords: a vector of length 2 containing the x and y coordinates of the point
    params: the params object.

  Inputs/Outputs:
    sol: vector of length 4 to be populated with the solution

  Aliasing restrictions: none

"""->
function calcIsentropicVortex(params::ParamType2,
                  coords::AbstractArray{Tmsh},
                  sol::AbstractVector{Tsol}) where {Tmsh, Tsol}
  # calculates the solution at a point of the isentropic vortex
  # 2D only


  # unpack arguments
  x = coords[1]
  y = coords[2]

  # get some values out of eqn
  cv = params.cv
  R = params.R
  gamma = params.gamma

  # the (hard coded) parameters are
  r_in = 1  # inner radius of sector of circle
  rho_in = 2 # density at inner radius
  M_in = 0.95  # Mach number
  p_in =  1/gamma


  # calculate r, theta coordinates from x,y
  r = sqrt(x*x + y*y)
  theta = atan2(y, x)  # angle in radians

  # calculate values at r radius
  tmp1 = ((gamma-1)/2)*M_in*M_in
  rho_r = rho_in*(1 + tmp1*(1- (r_in*r_in)/(r*r)))^(1/(gamma-1))

  p_r = p_in*(rho_r/rho_in)^gamma

  a_r = sqrt( gamma*p_r/rho_r )

  M_r = sqrt( (2/(gamma-1))*((rho_in/rho_r)^(gamma-1))*(1 + tmp1) - 2/(gamma-1) )
  U_r = M_r*a_r  # velocity magnitude

  u_r = U_r*sin(theta)
  v_r = -U_r*cos(theta)
  e_r = cv*p_r/(rho_r*R)
  E_r = rho_r*e_r + 0.5*rho_r*U_r*U_r

  # save solution to sol
  sol[1] = rho_r
  sol[2] = rho_r*u_r
  sol[3] = rho_r*v_r
  sol[4] = E_r

  return nothing

end


"""
  Reverse mode
"""
function calcIsentropicVortex_rev(params::ParamType2,
                  coords::AbstractArray{Tmsh}, coords_bar::AbstractVector{Tres},
                  sol_bar::AbstractVector{Tsol}) where {Tmsh, Tsol, Tres}
  # calculates the solution at a point of the isentropic vortex
  # 2D only


  # unpack arguments
  x = coords[1]
  y = coords[2]

  # get some values out of eqn
  cv = params.cv
  R = params.R
  gamma = params.gamma
  gamma_1 = params.gamma_1

  # the (hard coded) parameters are
  r_in = 1  # inner radius of sector of circle
  rho_in = 2 # density at inner radius
  M_in = 0.95  # Mach number
  p_in =  1/gamma


  # calculate r, theta coordinates from x,y
  r = sqrt(x*x + y*y)
  theta = atan2(y, x)  # angle in radians

  # calculate values at r radius
  tmp1 = ((gamma-1)/2)*M_in*M_in
  rho_r = rho_in*(1 + tmp1*(1- (r_in*r_in)/(r*r)))^(1/(gamma-1))

  p_r = p_in*(rho_r/rho_in)^gamma

  a_r = sqrt( gamma*p_r/rho_r )

  M_r = sqrt( (2/(gamma-1))*((rho_in/rho_r)^(gamma-1))*(1 + tmp1) - 2/(gamma-1) )
  U_r = M_r*a_r  # velocity magnitude

  u_r = U_r*sin(theta)
  v_r = -U_r*cos(theta)
  e_r = cv*p_r/(rho_r*R)
  E_r = rho_r*e_r + 0.5*rho_r*U_r*U_r

  # save solution to sol
  #sol[1] = rho_r
  #sol[2] = rho_r*u_r
  #sol[3] = rho_r*v_r
  #sol[4] = E_r

  #----------------------------------------------------------------------------
  # reverse sweep

  E_r_bar = zero(Tres)
  rho_r_bar = zero(Tres)
  v_r_bar = zero(Tres)
  u_r_bar = zero(Tres)
  U_r_bar = zero(Tres)
  e_r_bar = zero(Tres)
  p_r_bar = zero(Tres)
  a_r_bar = zero(Tres)
  M_r_bar = zero(Tres)
  r_bar   = zero(Tres)
  x_bar   = zero(Tres)
  y_bar   = zero(Tres)
  theta_bar = zero(Tres)

  E_r_bar   += sol_bar[4]

  rho_r_bar += v_r*sol_bar[3]
  v_r_bar   += rho_r*sol_bar[3]

  rho_r_bar += u_r*sol_bar[2]
  u_r_bar   += rho_r*sol_bar[2]

  rho_r_bar += sol_bar[1]


  # E_r through u_r
  rho_r_bar += (e_r + 0.5*U_r*U_r)*E_r_bar
  e_r_bar   += rho_r*E_r_bar
  U_r_bar   += rho_r*U_r*E_r_bar

  p_r_bar   += cv/(rho_r*R)*e_r_bar
  rho_r_bar += -e_r/rho_r*e_r_bar

  U_r_bar +=  sin(theta)*u_r_bar
  U_r_bar +=  -cos(theta)*v_r_bar
  theta_bar += U_r*sin(theta)*v_r_bar + U_r*cos(theta)*u_r_bar


  # calculating values at r radius
  a_r_bar += M_r*U_r_bar
  M_r_bar += a_r*U_r_bar

  rho_r_bar += (M_r_bar/M_r)*((rho_in/rho_r)^(gamma - 2))*(-rho_in/(rho_r*rho_r))*(1 + tmp1)
#  tmp1_bar   = M_r_bar/(M_r*gamma_1)*(rho_in/rho_r)^(gamma_1)

  p_r_bar   += (0.5/a_r)*gamma*a_r_bar/rho_r
  rho_r_bar += (0.5/a_r)*(-gamma*p_r/(rho_r*rho_r))*a_r_bar

  rho_r_bar += p_in*gamma*((rho_r/rho_in)^gamma_1)*p_r_bar/rho_in

  t2 = (rho_in/gamma_1)*(1 + tmp1*(1 - (r_in*r_in)/(r*r)))^(1/gamma_1 - 1)
#  tmp1_bar  += t2*(1 - (r_in*r_in)/(r*r))
  r_bar += rho_r_bar*t2*2*tmp1*(r_in*r_in)/(r*r*r)


  # r, theta
  theta, y_bar, x_bar = atan2_rev(y, x, theta_bar)
  x_bar += (0.5/r)*2*x*r_bar
  y_bar += (0.5/r)*2*y*r_bar

  coords_bar[1] += x_bar
  coords_bar[2] += y_bar


  return nothing
end
 


function calcIsentropicVortex(params::ParamType3,
                  coords::AbstractArray{Tmsh},
                  sol::AbstractVector{Tsol}) where {Tmsh, Tsol}
# calculates the solution at a point of the isentropic vortex


  # unpack arguments
  x = coords[1]; x_orig = x
  y = coords[2]; y_orig = y
  z = coords[3]; z_orig = z;

  #println("x = ", x, ", y = ", y, ", z = ", z)

  phi_z = pi/4  # angle of axis relative to z axis
  sgn = 1.0

  # calculate new x coordinate in rotated coordinate system
  theta1 = atan2(z, x)
  #println("theta1 = ", theta1)
  phi2 = 0.5*pi - theta1
  #println("phi2 = ", phi2)
  r_xz = sqrt(x*x + z*z)  # distance from origin to point in xz plane
  #println("r_xz = ", r_xz)
  x = r_xz*sin(phi_z + phi2)  # x coordinate in rotated system
  #println("x = ", x)
  #println("diff = ", abs(x - z_orig))

  theta3 = theta1 + phi_z + phi2 - 0.5*pi

  # get some values out of params
  cv = params.cv
  R = params.R
  gamma = params.gamma

  # the (hard coded) parameters are
  r_in = 1  # inner radius of sector of circle
  rho_in = 2 # density at inner radius
  M_in = 0.95  # Mach number
  p_in =  1/gamma

  # calculate r, theta coordinates from x,y
  r = sqrt(x*x + y*y)
  if phi_z > 0
    theta = atan2(x, y)
  else
    theta = atan2(y, x)
  end

  # calculate values at r radius
  tmp1 = ((gamma-1)/2)*M_in*M_in
  rho_r = rho_in*(1 + tmp1*(1- (r_in*r_in)/(r*r)))^(1/(gamma-1))
  p_r = p_in*(rho_r/rho_in)^gamma
  a_r = sqrt( gamma*p_r/rho_r )
  M_r = sqrt( (2/(gamma-1))*((rho_in/rho_r)^(gamma-1))*(1 + tmp1) - 2/(gamma-1) )
  U_r = M_r*a_r  # velocity magnitude

  if phi_z > 0
    v_r = U_r*sin(theta)
    u_r = -U_r*cos(theta)
  else
    u_r = U_r*sin(theta)
    v_r = -U_r*cos(theta)
  end
  e_r = cv*p_r/(rho_r*R)
  E_r = rho_r*e_r + 0.5*rho_r*U_r*U_r

  # convert velocities back to base frame of reference
  w_r = u_r*sin(theta3)
  u_r = u_r*cos(theta3)

  # save solution to sol
  sol[1] = rho_r
  sol[2] = rho_r*u_r
  sol[3] = rho_r*v_r
  sol[4] = rho_r*w_r
  sol[5] = E_r

  return nothing

end


@doc """
### EulerEquationMod.calcFreeStream

  This function calculates the free stream solution for an airfoil problem
  based on the angle of attack and Mach number in nondimensionalized variables.

  Density and energy are set to params.rho_free (usually 1.0) and params.E_free,
  (usually 1/(gamma*gamma_1) + 0.5*Ma*Ma), and the x and y momenta as

  rho*Ma*cos(angle of attack)  and rho*Ma*sin(angle of attack).

  The angle of attack must be in radians.

  This function uses conservative variables regardless of the static parameter
  of params.

  Inputs:
    coords: a vector of length 2 containing the x and y coordinates of the point
    params: the params object.

  Inputs/Outputs:
    sol: vector of length 4 to be populated with the solution

  Aliasing restrictions: none

"""->
function calcFreeStream(params::ParamType2,
            coords::AbstractArray{Tmsh, 1},
            sol::AbstractArray{Tsol, 1}) where {Tmsh, Tsol}
# calculate the free stream conditions using the fields of params


  rho = sol[1] = params.rho_free
  E = sol[4] = params.E_free

  Ma = params.Ma

  sol[2] = rho*Ma*cos(params.aoa)
  sol[3] = rho*Ma*sin(params.aoa)

  return nothing
end

function calcFreeStream(params::ParamType3,
            coords::AbstractArray{Tmsh, 1},
            sol::AbstractArray{Tsol, 1}) where {Tmsh, Tsol}
# calculate the free stream conditions using the fields of params

  # calculate the free stream conditions using the fields of params


  rho = sol[1] = params.rho_free
  E = sol[5] = params.E_free

  Ma = params.Ma

  sol[2] = rho*Ma*cos(params.aoa)
  sol[3] = 0.0
  sol[4] = -rho*Ma*sin(params.aoa)

  return nothing
end

"""
  Like [`calcFreeStream`](@ref), but assumes zero angle of attack, regardless
  of what the options dictionary says.
"""
function calcFreeStream0(params::ParamType2,
            coords::AbstractArray{Tmsh, 1},
            sol::AbstractArray{Tsol, 1}) where {Tmsh, Tsol}
# calculate the free stream conditions using the fields of params


  rho = sol[1] = params.rho_free
  E = sol[4] = params.E_free

  Ma = params.Ma

  sol[2] = rho*Ma
  sol[3] = 0

  return nothing
end



@doc """
### EulerEquationMod.calcFreeStream_daoa

  This function calculates the free stream solution for an airfoil problem
  based on the angle of attack and Mach number in nondimensionalized variables.

  Density and energy are set to params.rho_free (usually 1.0) and params.E_free,
  (usually 1/(gamma*gamma_1) + 0.5*Ma*Ma), and the x and y momenta as

  rho*Ma*cos(angle of attack)  and rho*Ma*sin(angle of attack).

  The angle of attack must be in radians.

  This function uses conservative variables regardless of the static parameter
  of params.

  Inputs:
    coords: a vector of length 2 containing the x and y coordinates of the point
    params: the params object.

  Inputs/Outputs:
    sol: vector of length 4 to be populated with the solution

  Aliasing restrictions: none

"""->

function calcFreeStream_dAlpha(params::ParamType2,
                   coords::AbstractArray{Tmsh, 1},
                   sol::AbstractArray{Tsol, 1}) where {Tmsh, Tsol}
# calculate the free stream conditions using the fields of params



  rho = params.rho_free
  Ma = params.Ma

  sol[2] = -rho*Ma*sin(params.aoa)
  sol[3] = -rho*Ma*cos(params.aoa)

  return nothing
end

function calcFreeStream_daoa(params::ParamType3,
                 coords::AbstractArray{Tmsh, 1},
                 sol::AbstractArray{Tsol, 1}) where {Tmsh, Tsol}

  # calculate the free stream conditions using the fields of params


  rho = params.rho_free
  Ma = params.Ma

  sol[2] = -rho*Ma*sin(params.aoa)
  sol[3] = 0.0
  sol[4] = -rho*Ma*cos(params.aoa)

  return nothing
end



@doc """
### EulerEquationMod.calcUnsteadyVortex

  This function calculates the unsteady vortex solution to the Euler equations
  at time params.t, where the vortex was centered at x = params.vortex_x0 at
  time t=0.  The analytical solution can be found in, for example,
  K. Mattsson et al./ Computers & Fluxs 36 (2007) 636-649

  This function uses conservative variables regardless of the static parameter
  of params.

  Inputs:
    coords: a vector of length 2 containing the x and y coordinates of the point
    params: the params object.

  Inputs/Outputs:
    sol: vector of length 4 to be populated with the solution

  Aliasing restrictions: none

"""->
function calcUnsteadyVortex(params::ParamType{Tdim},
          coords::AbstractArray{Tmsh, 1},
          sol::AbstractArray{Tsol, 1}) where {Tmsh, Tsol, Tdim}

  fval = vortex_f(coords, params)
  t = params.t
  epsilon = params.vortex_strength
  Ma = params.Ma
  gamma = params.gamma
  gamma_1 = params.gamma_1
  x = coords[1]
  y = coords[2]
  x0 = params.vortex_x0

  coeff1 = epsilon*epsilon*gamma_1*Ma*Ma/(8*pi*pi)
  rho = (1 -coeff1* e^(fval))^(1/gamma_1)

  coeff2 = epsilon*y/(2*pi)
  u = 1 - coeff2*e^(fval/2)

  coeff3 = epsilon*(x - x0 - t)/(2*pi)
  v = coeff3*e^(fval/2)

  q2 = rho*u
  q3 = rho*v

  # calculate E
  # the paper gives pressure, so calculate E from that and the
  # pressure expression in terms of E
  term1 = (rho^gamma)/(gamma_1*gamma*Ma*Ma)
  term2 = 0.5*(q2*q2 + q3*q3)/rho
  E = term1 + term2

  sol[1] = rho
  sol[2] = q2
  sol[3] = q3

  if Tdim == 2
    sol[4] = E
  else
    sol[4] = 0
    sol[5] = E
  end

  return nothing
end

function vortex_f(coords, params)
  t = params.t
  x = coords[1]
  y = coords[2]
  x0 = params.vortex_x0
  return 1 - ( (  (x-x0) - t)^2 + y*y)
end

"""
  Unsteady vortex for Carpenters paper "Entropy Stable Staggered Grid Spectral
  Collocation for hte Burgers and Compressible Navier-Stokes Equations

  It is similar to the other unsteady vortex, except the vortex travels
  at an angle

  Does *not* use params.Ma or params.vortex_strength, or vortex.x0

"""
function calcUnsteadyVortex2(params::ParamType{Tdim},
          coords::AbstractArray{Tmsh, 1},
          sol::AbstractArray{Tsol, 1}) where {Tmsh, Tsol, Tdim}
  t = params.t
  gamma = params.gamma
  gamma_1 = params.gamma_1
  x = coords[1]
  y = coords[2]
  x0 = 0.0
  y0 = 0
  epsilon = 5.0
  Ma = 0.5
#  alpha = 45.0  # degrees
  sind_alpha = 0.7071067811865476
  cosd_alpha = 0.7071067811865476
#  cinf = 1.0  # assumption, I think this is a free parameter
#  Uinf = Ma*cinf

  Uinf = 1.0
  cinf = Ma/Uinf

  ycoeff = y - y0 - Uinf*sind_alpha*t
  xcoeff = x - x0 - Uinf*cosd_alpha*t
  f = 1 - ( xcoeff^2 + ycoeff^2 )

  coeff1 = epsilon*epsilon*Ma*Ma*gamma_1/(8*pi*pi)
  T = 1 - coeff1*exp(f)
#  T = 1 - coeff1*e^f
  rho = T^(1/gamma_1)
 
#  u1 = Uinf - (epsilon*y/(2*pi))*e^(f/2)
  u1 = Uinf*cosd_alpha - (epsilon*ycoeff/(2*pi))*exp(f/2)
#  u2 = (epsilon*xcoeff/(2*pi))*e^(f/2)
  u2 = Uinf*sind_alpha + ((epsilon*xcoeff)/(2*pi))*exp(f/2)
 
  q2 = rho*u1
  q3 = rho*u2

#  p = rho*params.R*T
#  p = (rho^gamma)/(gamma*Ma*ma)
#  term1 = (p)/(gamma_1*gamma*Ma*Ma)
  term1 = (rho^gamma)/(gamma_1*gamma*Ma*Ma)
  term2 = 0.5*(q2*q2 + q3*q3)/rho
  E = term1 + term2

  sol[1] = rho
  sol[2] = q2
  sol[3] = q3
  sol[4] = E

  return nothing
end







  


@doc """
### EulerEquationMod.calcRho1Energy2

  This function sets the density at a node to 1, energy to 2 and the momenta to
  zero.

  This function uses conservative variables regardless of the static parameter
  of params.

  Inputs:
    coords: a vector of length 2 containing the x and y coordinates of the point
    params: the params object.

  Inputs/Outputs:
    sol: vector of length 4 to be populated with the solution

  Aliasing restrictions: none

"""->
function calcRho1Energy2(params::ParamType2,
             coords::AbstractArray{Tmsh, 1},
             sol::AbstractArray{Tsol,1}) where {Tmsh, Tsol}
  # for square test case with rho = 1, everything else  = 0

  sol[1] = 1.0
  sol[2] = 0.0
  sol[3] = 0.0
  sol[4] = 2.0

  return nothing
end
function calcRho1Energy2(params::ParamType3,
             coords::AbstractArray{Tmsh, 1},
             sol::AbstractArray{Tsol,1}) where {Tmsh, Tsol}
  # for square test case with rho = 1, everything else  = 0

  sol[1] = 1.0
  sol[2] = 0.0
  sol[3] = 0.0
  sol[4] = 0.0
  sol[5] = 2.0

  return nothing
end



@doc """
### EulerEquationMod.calcOnes

  This function sets all the solution variables at a node to 1.0

  This function uses conservative variables regardless of the static parameter
  of params.

  Inputs:
    coords: a vector of length 2 containing the x and y coordinates of the point
    params: the params object.

  Inputs/Outputs:
    sol: vector of length 4 to be populated with the solution

  Aliasing restrictions: none

"""->
function calcOnes(params::ParamType2,
      coords::AbstractArray{Tmsh, 1},
      sol::AbstractArray{Tsol,1}) where {Tmsh, Tsol}

  fill!(sol, 1.0)

  return nothing
end  # end function calcOnes

@doc """
### EulerEquationMod.calcZeros

  This function sets all the solution variables at a node to 0.0

  This function uses conservative variables regardless of the static parameter
  of params.

  Inputs:
    coords: a vector of length 2 containing the x and y coordinates of the point
    params: the params object.

  Inputs/Outputs:
    sol: vector of length 4 to be populated with the solution

  Aliasing restrictions: none

"""->

function calcZeros(params::ParamType2,
       coords::AbstractArray{Tmsh, 1},
       sol::AbstractArray{Tsol,1}) where {Tmsh, Tsol}

  fill!(sol, 0.0)

  return nothing
end  # end function calcZeros


@doc """
### EulerEquationMod.calcRho1Energy2U1VW0

  Sets the density values 1.0, x momentum to 1.0, 
  v & w momenta to 0.0, and energy to 2.0 at a node.

  It should work for 2D and 3D meshes.

  This function uses conservative variables regardless of the static parameter
  of params.

  Inputs:
    coords: a vector of length 2 containing the x and y coordinates of the point
    params: the params object.

  Inputs/Outputs:
    sol: vector of length 4 to be populated with the solution

  Aliasing restrictions: none

"""->
function calcRho1Energy2U1VW0(params::ParamType2,
                           coords::AbstractArray{Tmsh},
                           sol::AbstractArray{Tsol, 1}) where {Tmsh, Tsol}
  # for square test case with rho = 1, digonal momentum, energy

  sol[1] = 1.0
  sol[2] = 1.0
  sol[3] = 0.0
  sol[4] = 2.0

  return nothing
end

function calcRho1Energy2U1VW0(params::ParamType3,
                           coords::AbstractArray{Tmsh},
                           sol::AbstractArray{Tsol, 1}) where {Tmsh, Tsol}
  # for square test case with rho = 1, digonal momentum, energy

  sol[1] = 1.0
  sol[2] = 1.0
  sol[3] = 0.0
  sol[4] = 0.0
  sol[5] = 2.0

  return nothing
end

@doc """
### EulerEquationMod.calcRho1Energy2U3

  Sets the density values 1.0, momenta to 0.35355, and
  energy to 2.0 at a node.



  This function uses conservative variables regardless of the static parameter
  of params.

  Inputs:
    coords: a vector of length 2 containing the x and y coordinates of the point
    params: the params object.

  Inputs/Outputs:
    sol: vector of length 4 to be populated with the solution

  Aliasing restrictions: none

"""->
function calcRho1Energy2U3(params::ParamType2,
               coords::AbstractArray{Tmsh},
               sol::AbstractArray{Tsol, 1}) where {Tmsh, Tsol}
  # for square test case with rho = 1, digonal momentum, energy

  sol[1] = 1.0
  sol[2] = 0.35355
  sol[3] = 0.35355
  sol[4] = 2.0

  return nothing
end

function calcRho1Energy2U3(params::ParamType3,
               coords::AbstractArray{Tmsh},
               sol::AbstractArray{Tsol, 1}) where {Tmsh, Tsol}
  # for square test case with rho = 1, digonal momentum, energy

  sol[1] = 1.0
  sol[2] = 0.35355
  sol[3] = 0.35355
  sol[4] = 0.35355
  sol[5] = 2.0

  return nothing
end


"""
  Alias for [`calcRho1Energy2U3`](@ref) used by initial conditions.
"""
const calcRho1E2U3 = calcRho1Energy2U3



@doc """
### EulerEquationMod.calcVortex

  Sets the density 1.0, energy to 2.0 at a node.  The momenta are calculated
  according to solid body rotation with an angular velocity of 0.5 centered
  at x = 0.


  This function uses conservative variables regardless of the static parameter
  of params.

  Inputs:
    coords: a vector of length 2 containing the x and y coordinates of the point
    params: the params object.

  Inputs/Outputs:
    sol: vector of length 4 to be populated with the solution

  Aliasing restrictions: none

"""->
function calcVortex(params::ParamType2,
        coords::AbstractArray{Tmsh,1},
        sol::AbstractArray{Tsol,1}) where {Tmsh, Tsol}
# solid body rotation
  x = coords[1]
  y = coords[2]

  r = sqrt(x*x + y*y)
  theta = atan2(y, x)

  omega = 0.5

  u_norm = omega*r

  u = -u_norm*sin(theta)
  v = u_norm*cos(theta)

  rho = 1.0
  E = 2.0

  sol[1] = rho
  sol[2] = rho*u
  sol[3] = rho*v
  sol[4] = E

  return nothing
end

"""
  Calculates a manufactured solution based on exponentials
"""
function calcExp(params::ParamType2, coords::AbstractArray{Tmsh,1},
                 q::AbstractArray{Tsol,1}) where {Tmsh, Tsol}

  x = coords[1]
  y = coords[2]
  # a and b are parameters determine the scale and offset of the solution
  # the must be the same here and in the source term
  af = 1/5  # a = 1/af
  b = 0.01
  gamma_1 = params.gamma_1
  q[1] = exp(af*x*y + b)
  q[2] = exp(af*2*x*y + b)
  q[3] = exp(af*3*x*y + b)
  q[4] = (1/gamma_1 + 0.5)*exp(af*5*x*y + b) + 0.5*exp(af*3*x*y + b)

  return nothing
end


function calcExp_rev(params::ParamType2, coords::AbstractArray{Tmsh, 1},
                 coords_bar::AbstractArray{Tmsh, 1},
                 q_bar::AbstractArray{Tsol, 1}) where {Tmsh, Tsol}

  af = 1/5  # a = 1/af
  b = 0.01
  gamma_1 = params.gamma_1

  x = coords[1]
  y = coords[2]

  x_bar = zero(Tmsh)
  y_bar = zero(Tmsh)

  # some reused quantities
  t1 = (1/gamma_1 + 0.5)
  e1 = exp(af*x*y + b)
  e2 = exp(af*2*x*y + b)
  e3 = exp(af*3*x*y + b)
  e5 = exp(af*5*x*y + b)

  x_bar += (t1*af*5*y*e5 + 0.5*af*3*y*e3)*q_bar[4]
  y_bar += (t1*af*5*x*e5 + 0.5*af*3*x*e3)*q_bar[4]

  x_bar += af*3*y*e3*q_bar[3]
  y_bar += af*3*x*e3*q_bar[3]

  x_bar += af*2*y*e2*q_bar[2]
  y_bar += af*2*x*e2*q_bar[2]

  x_bar += af*y*e1*q_bar[1]
  y_bar += af*x*e1*q_bar[1]

  coords_bar[1] += x_bar
  coords_bar[2] += y_bar

  return nothing
end

function calcExp(params::ParamType3, coords::AbstractArray{Tmsh,1},
                 q::AbstractArray{Tsol,1}) where {Tmsh, Tsol}
  x = coords[1]
  y = coords[2]
  z = coords[3]

  # constant parameters
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

  # af = 1/50, b = 0.01, c = 8, d = 0.5, f = 10 works
  gamma_1 = params.gamma_1

  t2 = exp(b);
  t3 = a*c1*x*y*z;
  q[1] = d1*t2*exp(t3);
  q[2] = d2*t2*exp(a*c2*x*y*z);
  q[3] = d3*t2*exp(a*c3*x*y*z);
  q[4] = d4*t2*exp(a*c4*x*y*z);
  q[5] = (t2*exp(-t3)*((d2*d2)*exp(a*c2*x*y*z*2.0)+(d3*d3)*exp(a*c3*x*y*z*2.0)+(d4*d4)*exp(a*c4*x*y*z*2.0))*(1.0/2.0))/d1+(d5*t2*exp(a*c5*x*y*z))/gamma_1;

  return nothing
end

"""
  Calculates a manufactured solution from Gassner, Winters, Kopriva: Split
  Form Nodal DG Schemes with SBP Propertiy for the Compressible Euler Equations.
  This is typically used with a mesh that spans [-1, 1] in all directions

"""
function calcPeriodicMMS(params::ParamType2,
             coords::AbstractArray{Tmsh,1},
             q::AbstractArray{Tsol,1}) where {Tmsh, Tsol}

  x = coords[1]
  y = coords[2]
  t = params.t

  t7 = t*2.0;
  t2 = -t7+x+y;
  t3 = 3.141592653589793*t2;
  t4 = sin(t3);
  t5 = t4*(1.0/1.0E1);
  t6 = t5+2.0;
  q[1] = t6;
  q[2] = t6;
  q[3] = t6;
  q[4] = t6*t6;

  return nothing
end

function calcPeriodicMMS(params::ParamType3,
             coords::AbstractArray{Tmsh,1},
             q::AbstractArray{Tsol,1}) where {Tmsh, Tsol}

  x = coords[1]
  y = coords[2]
  z = coords[3]
  t = params.t
#=
  rho = 2 + 0.1*sin(pi*(x + y + z))
#  rho = 2 + 0.5*(x*x + y*y + z*z)
  q[1] = rho
  q[2] = 0
  q[3] = 0
  q[4] = 0
  q[5] = 10
=#
  t7 = t*2.0;
  t2 = -t7+x+y+z;
  t3 = 3.141592653589793*t2;
  t4 = sin(t3);
  t5 = t4*(1.0/1.0E1);
  t6 = t5+2.0;

  q[1] = t6
  q[2] = t6
  q[3] = t6
  q[4] = t6
  q[5] = t6*t6



  return nothing
end

"""
  Manufactured condition for a channel with y in [0, 1].
  No source term.  2D only
"""
function calcChannelMMS(params::ParamType2,
            coords::AbstractArray{Tmsh,1},
            q::AbstractArray{Tsol,1}) where {Tmsh, Tsol}

  R = params.R
  gamma_1 = params.gamma_1

  rho_inf = 1
  u_inf = 1
  offset = 0.1
  T_inf = 1  # make this smaller to reduce Mach number

  y = coords[2]

  q[1] = rho_inf
  q[2] = rho_inf*u_inf*( y*(1-y) + offset)
  q[3] = 0
  q[4] = 0.5*rho_inf*u_inf*u_inf*( offset - y*(y-1))^2 + (R*T_inf*rho_inf)/gamma_1

  return nothing
end

"""
  One dimensional square wave, 2D only

  The base flow is uniform, and if  -0.05 < x < 0.05, then it is perturbed
  by a small amount.
"""
function calcSquare1D(params::ParamType2,
            coords::AbstractArray{Tmsh,1},
            q::AbstractArray{Tsol,1}) where {Tmsh, Tsol}

  # apply base state
  q[1] = 1.0
  q[2] = 0.3
  q[3] = 0.3
  q[4] = 5.0

  x = coords[1]
  if x > -0.05 && x < 0.05
    # apply pertubration
    q[1] += 0.1
    q[2] += 0.1
    q[3] += 0.1
    q[4] += 0.1
  end

  return nothing
end

"""
  Two dimensional square wave, 2D only

  The base flow is uniform, and if  -0.05 < x, y < 0.05, then it is perturbed
  by a small amount.
"""
function calcSquare2D(params::ParamType2,
            coords::AbstractArray{Tmsh,1},
            q::AbstractArray{Tsol,1}) where {Tmsh, Tsol}

  # apply base state
  q[1] = 1.0
  q[2] = 0.3
  q[3] = 0.3
  q[4] = 5.0

  x = coords[1]
  y = coords[2]
  if x > 7.5 && x < 12.5 && y > -1.25 && y < 1.25
    # apply pertubration
    q[1] += 0.1
    q[2] += 0.1
    q[3] += 0.1
    q[4] += 0.1
  end

  return nothing
end

"""
  Sedov explosion, 2D only

  Uniform fluid at rest with energetic region in the circle of radius 
  0.05 centered at the origin

"""
function calcSedovExplosion(params::ParamType2,
                        coords::AbstractArray{Tmsh,1},
                        q::AbstractArray{Tsol,1}) where {Tmsh, Tsol}


  # parameters
  E = 1.0 # total amount of energy contained in the blast (not the usual E)
  nu = 2  # cylindrical blast wave
  dr = 0.05  # radius of initial blast

  q[1] = 1.0
  q[2] = 0.0
  q[3] = 0.0
  q[4] = 100*(1e-5)/params.gamma_1

  x = coords[1]
  y = coords[2]
  if x*x + y*y < dr*dr
    q[4] = 100*3*params.gamma_1*E/( (nu+1)*Float64(pi)*(dr^nu) )
  end

  return nothing
end

"""
  Initial condition for SU2 inviscid channel
"""
function calcInvChannelIC(params::ParamType2,
             coords::AbstractArray{Tmsh,1},
             q::AbstractArray{Tsol,1}) where {Tmsh, Tsol}


  rho = 1.22531
  rhou = rho*170.104
  rhov = 0
  E = 101300/params.gamma_1 + 0.5*rho*(rhou*rhou + rhov*rhov)

  q[1] = rho
  q[2] = rhou
  q[3] = rhov
  q[4] = E

  return nothing
end

"""
  Free stream conditions for SU2 inviscid channel
"""
function calcInvChannelFreeStream(params::ParamType2,
                         coords::AbstractArray{Tmsh,1},
                         q::AbstractArray{Tsol,1}) where {Tmsh, Tsol}


  rho = 1.22531
  rhou = rho*170.104
  rhov = 0
  E = 101100/params.gamma_1 + 0.5*rho*(rhou*rhou + rhov*rhov)

  q[1] = rho
  q[2] = rhou
  q[3] = rhov
  q[4] = E

  return nothing
end

"""
  Used for testing shock capturing scheme
"""
function calcLaplaceSolution(params::ParamType,
                         coords::AbstractArray{Tmsh,1},
                         q::AbstractArray{Tsol,1}) where {Tmsh, Tsol}

  x = coords[1]; y = coords[2]
  val = x*x*x + y*y*y
  val = sin(x)*sin(2*y)
  for i=1:length(q)
    q[i] = val + i
  end

  return nothing
end


global const wedge20_tanbeta = tand(53.422940527228654)
global const tan20 = tand(20)
global const cos20 = cosd(20)
global const sin20 = sind(20)

"""
  Compute the exact solution for a supersonic wedge where the semi-angle is
  20 degrees.  The vertex of the wedge must be at x, y = (-1, 0), and the
  free stream Mach number must be 2.0
"""
function calcWedge20(params::ParamType2,
                     coords::AbstractArray{Tmsh,1},
                     q::AbstractArray{Tsol,1}) where {Tmsh, Tsol}

#  @assert params.Ma == 2.0
#  @assert params.aoa == 0.0

  # figure out if the given point is post-shock or pre-shock

  x = coords[1]; y = coords[2]; tanbeta = wedge20_tanbeta

  if (y > 0 && y > tanbeta*x + tanbeta) || (y <= 0 && y < -tanbeta*x - tanbeta)
    # pre-shock
    calcFreeStream(params, coords, q)
  else
    # computed from the oblique shock relations in Anderson's Aerodynamics
    M2 = 1.2102184008268027
    a2 = 1.1799115865336027

    q[1] = 2.0420057206059625
    q[2] = q[1]*M2*a2*cos20
    q[3] = q[1]*M2*a2*sin20
    q[4] = 7.1584095248438615

    if y < 0
      q[3] = -q[3]
    end
  end

  return nothing
end


"""
  Computes the potential flow solution for flow over a 20 degree wedge.
  Potential flow assumes incompressible flow, so the Mach number should be
  small when running this case
"""
function calcWedge20Potential(params::ParamType2,
                     coords::AbstractArray{Tmsh,1},
                     q::AbstractArray{Tsol,1}) where {Tmsh, Tsol}

  n = 9/8  # this corresponds to a 20 degree semi-angle wedge
  rprime = 2*sqrt(1 + 1.5*1.5)  # define pressure = free stream pressure at this
                                # point
  Vprime = n*params.Ma*rprime^(n-1)
  Vprime2 = Vprime*Vprime

  # Bernoulli's equation constant at rprime.  Use this for *all* streamlines
  c1 = params.p_free/params.rho_free + 0.5*Vprime2

  # the potential flow solution has the wedge geometry as:
  #  \                                    /
  #   \             rather than          /
  #    \                                /
  #     \________                ______/
  # so we transform the coordinates at the beginning of the computation
  # and then transform the velocity at the end.
  # Also, the potential flow solution as the origin at the corner of the wedge,
  # not the lower right corner

  x = -(coords[1] + 1.0); y = coords[2]

  # The potential flow solution si f(z) = z^n, where z = x + i*y, but computing
  # it this way would not be complex-stepable,  Instead apply Eulers identity

  # compute r-theta coordinates
  r = sqrt(x*x + y*y); nrn = n*params.Ma*r^(n-1)
  theta = atan2(y, x); nm1_theta = (n-1)*theta

  # compute velocity
  u = nrn*cos(nm1_theta); v = -nrn*sin(nm1_theta)

  # transform velocity from the potential flow geometry to the real geometry
  v = -v

  # compute pressure using Bernoulli's equation, then energy from it (using
  # ideal gas law)
  v_squared = u*u + v*v
  p = params.rho_free*(c1 - 0.5*v_squared)
  E = p/params.gamma_1 + 0.5*params.rho_free*v_squared

  q[1] = params.rho_free
  q[2] = params.rho_free*u
  q[3] = params.rho_free*v
  q[4] = E

  return nothing
end



