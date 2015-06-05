# functions that populate the initial conditions
# List of functions:
#   ICZero (all zeros)
#   ICRho1E2 (all zeros, except rho = 1, E = 2)
#   ICLinear
#   ICsmoothHeavisideder
#   ICsmoothHeaviside
#   ICIsentropicVortex

export ICZero, ICRho1E2, ICLinear, ICsmoothHeavisideder, ICsmoothHeaviside, ICIsentropicVortex

function ICZero(mesh::AbstractMesh, operator::SBPOperator, eqn::EulerEquation, u0::AbstractVector)
# populate u0 with initial values
# this is a template for all other initial conditions

numEl = getNumEl(mesh)
nnodes = operator.numnodes
dofpernode = getNumDofPerNode(mesh)
for i=1:numEl
  dofnums_i = getGlobalNodeNumbers(mesh, i)  # get dof nums for this element
  coords = getElementVertCoords(mesh, [i])

  for j=1:nnodes
      # get dof numbers for each variable
      dofnum_rho = dofnums_i[1,j]
      dofnum_rhou = dofnums_i[2,j]
      dofnum_rhov = dofnums_i[3,j]
      dofnum_e = dofnums_i[4,j]

      # coordinates of this node (must be a vertex)
      x = coords[1,j]
      y = coords[2,j]
      z = coords[3,j]

      # apply initial conditions here
      u0[dofnum_rho] = 0.0
      u0[dofnum_rhou] = 0.0
      u0[dofnum_rhov] = 0.0
      u0[dofnum_e] = 0.0
  end
end

return nothing

end  # end function

function ICRho1E2(mesh::AbstractMesh, operator::SBPOperator, eqn::EulerEquation, u0::AbstractVector)
# populate u0 with initial values
# this is a template for all other initial conditions

numEl = getNumEl(mesh)
nnodes = operator.numnodes
dofpernode = getNumDofPerNode(mesh)
for i=1:numEl
  dofnums_i = getGlobalNodeNumbers(mesh, i)  # get dof nums for this element
  coords = getElementVertCoords(mesh, [i])

  for j=1:nnodes
      # get dof numbers for each variable
      dofnum_rho = dofnums_i[1,j]
      dofnum_rhou = dofnums_i[2,j]
      dofnum_rhov = dofnums_i[3,j]
      dofnum_e = dofnums_i[4,j]

      # coordinates of this node (must be a vertex)
      x = coords[1,j]
      y = coords[2,j]
      z = coords[3,j]

      # apply initial conditions here
      u0[dofnum_rho] = 1.0
      u0[dofnum_rhou] = 0.0
      u0[dofnum_rhov] = 0.0
      u0[dofnum_e] = 2.0
  end
end

return nothing

end  # end function


function ICRho1E2U3(mesh::AbstractMesh, operator::SBPOperator, eqn::EulerEquation, u0::AbstractVector)
# populate u0 with initial values
# this is a template for all other initial conditions

numEl = getNumEl(mesh)
nnodes = operator.numnodes
dofpernode = getNumDofPerNode(mesh)
sol = zeros(4)
for i=1:numEl
  dofnums_i = getGlobalNodeNumbers(mesh, i)  # get dof nums for this element
  coords = getElementVertCoords(mesh, [i])

  for j=1:nnodes
      # get dof numbers for each variable
      dofnum_rho = dofnums_i[1,j]
      dofnum_rhou = dofnums_i[2,j]
      dofnum_rhov = dofnums_i[3,j]
      dofnum_e = dofnums_i[4,j]

      # coordinates of this node (must be a vertex)
      x = coords[1,j]
      y = coords[2,j]
      z = coords[3,j]

      calcRho1Energy2U3(coords[:,j], eqn, sol)


      # apply initial conditions here
#      u0[dofnum_rho] = 1.0
#      u0[dofnum_rhou] = 3.0
#      u0[dofnum_rhov] = 0.0
#      u0[dofnum_e] = 2.0

      u0[dofnums_i[:,j]] = sol
  end
end

return nothing

end  # end function



function ICVortex(mesh::AbstractMesh, operator::SBPOperator, eqn::EulerEquation, u0::AbstractVector)
# populate u0 with initial values
# this is a template for all other initial conditions

numEl = getNumEl(mesh)
nnodes = operator.numnodes
dofpernode = getNumDofPerNode(mesh)
sol = zeros(4)
for i=1:numEl
  dofnums_i = getGlobalNodeNumbers(mesh, i)  # get dof nums for this element
  coords = getElementVertCoords(mesh, [i])

  for j=1:nnodes
      # get dof numbers for each variable
      dofnum_rho = dofnums_i[1,j]
      dofnum_rhou = dofnums_i[2,j]
      dofnum_rhov = dofnums_i[3,j]
      dofnum_e = dofnums_i[4,j]

      # coordinates of this node (must be a vertex)
      x = coords[1,j]
      y = coords[2,j]
      z = coords[3,j]

      calcVortex(coords[:,j], eqn, sol)


      # apply initial conditions here
#      u0[dofnum_rho] = 1.0
#      u0[dofnum_rhou] = 3.0
#      u0[dofnum_rhov] = 0.0
#      u0[dofnum_e] = 2.0

      u0[dofnums_i[:,j]] = sol
  end
end

return nothing

end  # end function




function ICLinear(mesh::AbstractMesh, operator::SBPOperator, eqn::EulerEquation, u0::AbstractVector)
# populate u0 with initial values
# this is a template for all other initial conditions

nnodes = operator.numnodes
dofpernode = getNumDofPerNode(mesh)
dofnums_i = zeros(dofpernode)

cntr = 1
for i=1:mesh.numVert
  for j=1:dofpernode
    dofnums_i[j] = getNumberJ(mesh.dofnums_Nptr, mesh.verts[i], 0, j-1)
  end

      dofnum_rho = dofnums_i[1]
      dofnum_rhou = dofnums_i[2]
      dofnum_rhov = dofnums_i[3]
      dofnum_e = dofnums_i[4]


      # apply initial conditions here
      u0[dofnum_rho] = cntr
      u0[dofnum_rhou] = cntr+1
      u0[dofnum_rhov] = cntr+2
      u0[dofnum_e] = cntr+3

      cntr += 4
end

return nothing

end  # end function


function ICsmoothHeavisideder(mesh::AbstractMesh, operator::SBPOperator, eqn::EulerEquation, u0::AbstractVector)
# calculate the value of the smooth heaviside function derivative at a location x
# x0 is specified within this function

# smooth heaviside  parameters
  x0 = 0
  L = 5
  k = 5



numEl = getNumEl(mesh)
nnodes = operator.numnodes
dofpernode = getNumDofPerNode(mesh)
for i=1:numEl
  dofnums_i = getGlobalNodeNumbers(mesh, i)  # get dof nums for this element
  coords = getElementVertCoords(mesh, [i])

  for j=1:nnodes
      # get dof numbers for each variable
      dofnum_rho = dofnums_i[1,j]
      dofnum_rhou = dofnums_i[2,j]
      dofnum_rhov = dofnums_i[3,j]
      dofnum_e = dofnums_i[4,j]

      # coordinates of this node (must be a vertex)
      x = coords[1,j]
      y = coords[2,j]
      z = coords[3,j]

      # apply initial conditions here
      u0[dofnum_rho] = L*(2*k*e^(-2*k*x))/(e^(-2*k*x) +1 )^2
      u0[dofnum_rhou] = 0.0
      u0[dofnum_rhov] = 0.0
      u0[dofnum_e] = 0.0
  end
end

return nothing



end

function ICsmoothHeaviside(mesh::AbstractMesh, operator::SBPOperator, eqn::EulerEquation, u0::AbstractVector)
# calculate the value of the smooth heaviside function at a location x
# x0 is specified within this function

# smooth heaviside  parameters
  x0 = 0
  L = 5
  k = 5



numEl = getNumEl(mesh)
nnodes = operator.numnodes
dofpernode = getNumDofPerNode(mesh)
for i=1:numEl
  dofnums_i = getGlobalNodeNumbers(mesh, i)  # get dof nums for this element
  coords = getElementVertCoords(mesh, [i])

  for j=1:nnodes
      # get dof numbers for each variable
      dofnum_rho = dofnums_i[1,j]
      dofnum_rhou = dofnums_i[2,j]
      dofnum_rhov = dofnums_i[3,j]
      dofnum_e = dofnums_i[4,j]

      # coordinates of this node (must be a vertex)
      x = coords[1,j]
      y = coords[2,j]
      z = coords[3,j]

      # apply initial conditions here
      u0[dofnum_rho] = L/(1 + e^(-k*(x-x0)))
      u0[dofnum_rhou] = 0.0
      u0[dofnum_rhov] = 0.0
      u0[dofnum_e] = 0.0
  end
end

return nothing



end

function ICIsentropicVortex(mesh::AbstractMesh, operator::SBPOperator, eqn::EulerEquation, u0::AbstractVector)
# populate u0 with initial values
# this is a template for all other initial conditions

numEl = getNumEl(mesh)
nnodes = operator.numnodes
dofpernode = getNumDofPerNode(mesh)
sol = zeros(4)
for i=1:numEl
  dofnums_i = getGlobalNodeNumbers(mesh, i)  # get dof nums for this element
  coords = getElementVertCoords(mesh, [i])

  for j=1:nnodes

      # coordinates of this node (must be a vertex)
      coords_j = coords[:,j]
      calcIsentropicVortex(coords_j, eqn, sol)

      # apply initial conditions here
      u0[dofnums_i[:,j]] = sol
  end
end

return nothing

end  # end function

function ICIsentropicVortexWithNoise(mesh::AbstractMesh, operator::SBPOperator, eqn::EulerEquation, u0::AbstractVector)
# populate u0 with initial values
# this is a template for all other initial conditions

numEl = getNumEl(mesh)
nnodes = operator.numnodes
dofpernode = getNumDofPerNode(mesh)
sol = zeros(4)
for i=1:numEl
  dofnums_i = getGlobalNodeNumbers(mesh, i)  # get dof nums for this element
  coords = getElementVertCoords(mesh, [i])

  for j=1:nnodes

      # coordinates of this node (must be a vertex)
      coords_j = coords[:,j]
      calcIsentropicVortex(coords_j, eqn, sol)

      # apply initial conditions here
#       u0[dofnums_i[:,j]] = sol
      u0[dofnums_i[:,j]] = sol+0.1*rand(4)
  end
end

return nothing

end  # end function


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



