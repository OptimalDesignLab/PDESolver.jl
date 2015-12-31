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
  K. Mattsson et al./ Computers & Fluixs 36 (2007) 636-649

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

  Sets the density values 1.0, x and y momenta to 0.35355, and
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



