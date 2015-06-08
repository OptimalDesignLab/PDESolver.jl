
function calcIsentropicVortex(coords::AbstractVector, eqn::EulerEquation, sol::AbstractVector)
# calculates the solution at a point of the isentripic vortex
# 2D only
# coords contains xy coordinates
# sol is a vector to be populated with the solution at the point

#sol = zeros(4)  # storage for solution

# unpack arguments
x = coords[1]
y = coords[2]

# get some values out of eqn
cv = eqn.cv
R = eqn.R
gamma = eqn.gamma

# the (hard coded) parameters are
r_in = 1  # inner radius of sector of circle
rho_in = 2 # density at inner radius
M_in = 0.95  # Mach number
p_in =  1/gamma


# calculate r, theta coordinates from x,y
r = sqrt(x*x + y*y)
theta = atan2(y,x)  # angle in radians
# println("theta = ", theta)

# calculate values at inner radius
a_in = sqrt(gamma*p_in/rho_in)  # speed of sound
u_norm_in = M_in*a_in  # magnitude of velocity

# calculate values at r
rho_r = rho_in*(1 + (gamma-1)*M_in*M_in*(1 - (r_in*r_in)/(r*r))/2)^(1/(gamma-1))
# println("rho_r = ", rho_r)
p_r = p_in*( (rho_r/rho_in)^gamma)  # isentropic relation
a_r = sqrt(gamma*p_r/rho_r)
u_norm_r = M_in*a_r  # M_in is constant
# println("u_norm_r = ", u_norm_r)
u_r = u_norm_r*sin(theta)
v_r = -u_norm_r*cos(theta)
# println("============================================ U_R: ",u_r,"  V_R: ",v_r)

# println("u_r = ", u_r, " v_r = ", v_r)

T_r = p_r/(rho_r*R)  # ideal gas law
E_r = rho_r*cv*T_r + 0.5*M_in*M_in*gamma*p_r

# save solution to sol
sol[1] = rho_r
sol[2] = rho_r*u_r
sol[3] = rho_r*v_r
sol[4] = E_r

return nothing

end

function calcRho1Energy2(coords::AbstractVector, eqn::EulerEquation, sol::AbstractVector)
  # for square test case with rho = 1, everything else  = 0

  sol[1] = 1.0
  sol[4] = 2.0

  return nothing
end


function calcRho1Energy2U3(coords::AbstractVector, eqn::EulerEquation, sol::AbstractVector)
  # for square test case with rho = 1, everything else  = 0

  sol[1] = 1.0
  sol[2] = 0.49
  sol[3] = 0.49
  sol[4] = 2.0

  return nothing
end


function calcVortex(coords::AbstractVector, eqn::EulerEquation, sol::AbstractVector)
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



