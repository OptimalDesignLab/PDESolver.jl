function calcIsentropicVortex{Tmsh, Tsol}(coords::AbstractArray{Tmsh}, params::ParamType{2}, sol::AbstractVector{Tsol})
# calculates the solution at a point of the isentripic vortex
# 2D only
# coords contains xy coordinates
# sol is a vector to be populated with the solution at the point

#sol = zeros(4)  # storage for solution

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
# println("r = ", r)
theta = atan2(y,x)  # angle in radians
# println("theta = ", theta)

#println("r = ", r)
#println("theta = ", theta)

tmp1 = ((gamma-1)/2)*M_in*M_in
# println("tmp1 = ", tmp1)

# calculate values at r radius
rho_r = rho_in*(1 + tmp1*(1- (r_in*r_in)/(r*r)))^(1/(gamma-1))
# println("rho_r = ", rho_r)

p_r = p_in*(rho_r/rho_in)^gamma
# println("p_r = ", p_r)

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



function calcFreeStream{Tmsh, Tsol}(coords::AbstractArray{Tmsh, 1}, params::ParamType{2}, sol::AbstractArray{Tsol, 1})
# calculate the free stream conditions using the fields of params

  
  rho = sol[1] = params.rho_free
  E = sol[4] = params.E_free

  cv = params.cv
  gamma = params.gamma
  R = params.R
  Ma = params.Ma

#  num = gamma*R*E/(rho*cv)
#  denom = 1/(Ma*Ma) + gamma*R/(2*cv)


#  u_norm = sqrt(num/denom)  # magnitude of free stream velocity

  sol[2] = rho*Ma*cos(params.aoa)
  sol[3] = -rho*Ma*sin(params.aoa)

  return nothing
end


function calcUnsteadyVortex{Tmsh, Tsol}(coords::AbstractArray{Tmsh, 1}, params::ParamType{2}, sol::AbstractArray{Tsol, 1})
# calculate the unsteady vortex at a given point
# coords specifies the xy coordinates
# the time is contained in params
# from K. Mattsson et al./ Computers & Fluixs 36 (2007) 636-649

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




function calcRho1Energy2{Tmsh, Tsol}(coords::AbstractArray{Tmsh, 1}, params::ParamType{2}, sol::AbstractArray{Tsol,1})
  # for square test case with rho = 1, everything else  = 0

  sol[1] = 1.0
  sol[4] = 2.0

  return nothing
end


function calcRho1Energy2U3{Tmsh, Tsol}(coords::AbstractArray{Tmsh}, params::ParamType{2}, sol::AbstractArray{Tsol, 1})
  # for square test case with rho = 1, digonal momentum, energy

  sol[1] = 1.0
  sol[2] = 0.35355
  sol[3] = 0.35355
  sol[4] = 2.0

#=
  sol[1] = 2.0
  sol[2] = 0.2
  sol[3] = 0.2
  sol[4] = 4.0
=#

  return nothing
end


function calcVortex{Tmsh, Tsol}(coords::AbstractArray{Tmsh,1}, params::ParamType{2}, sol::AbstractArray{Tsol,1})
# solid body rotation
  x = coords[1]
  y = coords[2]

#  r = sqrt(x*x + y*y)
#  theta = atan2(y, x)

  omega = 0.5

  u_norm = omega*r

#  u = -u_norm*sin(theta)
#  v = u_norm*cos(theta)
  u0 = 0.1


  rho = 1.0
  E = 2.0

  sol[1] = rho
  sol[2] = rho*u0*x
  sol[3] = 0.0
  sol[4] = E

  return nothing
end



