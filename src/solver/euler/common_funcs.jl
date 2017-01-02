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
function calcIsentropicVortex{Tmsh, Tsol}(coords::AbstractArray{Tmsh}, 
                              params::ParamType{2}, sol::AbstractVector{Tsol})
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
theta = atan2(y,x)  # angle in radians

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


function calcIsentropicVortex{Tmsh, Tsol}(coords::AbstractArray{Tmsh}, 
                              params::ParamType{3}, sol::AbstractVector{Tsol})
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
function calcFreeStream{Tmsh, Tsol}(coords::AbstractArray{Tmsh, 1}, 
                        params::ParamType{2}, sol::AbstractArray{Tsol, 1})
# calculate the free stream conditions using the fields of params

  
  rho = sol[1] = params.rho_free
  E = sol[4] = params.E_free

  Ma = params.Ma

  sol[2] = rho*Ma*cos(params.aoa)
  sol[3] = -rho*Ma*sin(params.aoa)

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
function calcUnsteadyVortex{Tmsh, Tsol}(coords::AbstractArray{Tmsh, 1}, 
                            params::ParamType{2}, sol::AbstractArray{Tsol, 1})

  function f(coords, params)
    t = params.t
    x = coords[1]
    y = coords[2]
    x0 = params.vortex_x0
    return 1 - ( (  (x-x0) - t)^2 + y*y)
  end
  
  fval = f(coords, params)
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
function calcRho1Energy2{Tmsh, Tsol}(coords::AbstractArray{Tmsh, 1}, 
                         params::ParamType{2}, sol::AbstractArray{Tsol,1})
  # for square test case with rho = 1, everything else  = 0

  sol[1] = 1.0
  sol[2] = 0.0
  sol[3] = 0.0
  sol[4] = 2.0

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
function calcOnes{Tmsh, Tsol}(coords::AbstractArray{Tmsh, 1}, 
                  params::ParamType{2}, sol::AbstractArray{Tsol,1})
  
  fill!(sol, 1.0)

  return nothing
end  # end function calcOnes


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
function calcRho1Energy2U3{Tmsh, Tsol}(coords::AbstractArray{Tmsh}, 
                           params::ParamType{2}, sol::AbstractArray{Tsol, 1})
  # for square test case with rho = 1, digonal momentum, energy

  sol[1] = 1.0
  sol[2] = 0.35355
  sol[3] = 0.35355
  sol[4] = 2.0

  return nothing
end

function calcRho1Energy2U3{Tmsh, Tsol}(coords::AbstractArray{Tmsh}, 
                           params::ParamType{3}, sol::AbstractArray{Tsol, 1})
  # for square test case with rho = 1, digonal momentum, energy

  sol[1] = 1.0
  sol[2] = 0.35355
  sol[3] = 0.35355
  sol[4] = 0.35355
  sol[5] = 2.0

  return nothing
end



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
function calcVortex{Tmsh, Tsol}(coords::AbstractArray{Tmsh,1}, 
                    params::ParamType{2}, sol::AbstractArray{Tsol,1})
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
function calcExp{Tmsh, Tsol}(coords::AbstractArray{Tmsh,1}, 
                    params::ParamType{2}, q::AbstractArray{Tsol,1})

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

function calcExp{Tmsh, Tsol}(coords::AbstractArray{Tmsh,1}, 
                    params::ParamType{3}, q::AbstractArray{Tsol,1})
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
function calcPeriodicMMS{Tmsh, Tsol}(coords::AbstractArray{Tmsh,1}, 
                    params::ParamType{2}, q::AbstractArray{Tsol,1})

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

function calcPeriodicMMS{Tmsh, Tsol}(coords::AbstractArray{Tmsh,1}, 
                    params::ParamType{3}, q::AbstractArray{Tsol,1})

  x = coords[1]
  y = coords[2]
  z = coords[3]
  t = params.t

  rho = 2 + 0.1*sin(pi*(x + y + z))
#  rho = 2 + 0.5*(x*x + y*y + z*z)
  q[1] = rho
  q[2] = 0
  q[3] = 0
  q[4] = 0
  q[5] = 10

  return nothing
end
