# functions that populate the initial conditions
# List of functions:

@doc """
### EulerEquationMod.ICTrigonometric

Sets all components of the solution to the free stream condition according
to the angle of attack and and Mach number.

Inputs:
mesh
sbp
eqn
opts

Inputs/Outputs: 
u0: vector to populate with the solution

Aliasing restrictions: none.

"""->

function ICPolynomial{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                                        sbp::AbstractSBP{Tsbp}, 
                                        eqn::EulerData{Tsol, Tsol, 3}, 
                                        opts, 
                                        u0::AbstractVector{Tsol})
  sigma = 0.01
  params = eqn.params
  gamma = params.gamma
	gamma_1 = params.gamma_1
	aoa = params.aoa
  beta = eqn.params.sideslip_angle
	rhoInf = 1.0
  uInf = eqn.params.Ma * cos(beta) * cos(aoa)
  vInf = eqn.params.Ma * sin(beta) * -1
  wInf = eqn.params.Ma * cos(beta) * sin(aoa)
	TInf = 1.0

  numEl = mesh.numEl
  nnodes = mesh.numNodesPerElement
  dofpernode = mesh.numDofPerNode
  sol = zeros(Tsol, 5)

  for i=1:numEl
    for j=1:nnodes
      coords_j = sview(mesh.coords, :, j, i)
      dofnums_j = sview(mesh.dofs, :, j, i)
      x = coords_j[1]
      y = coords_j[2]
      z = coords_j[3]

      rho = (x-x*x)*(y-y*y)* (z - z*z)
      u   = (x-x*x)*(y-y*y)* (z - z*z)
      v   = (x-x*x)*(y-y*y)* (z - z*z)
      w   = (x-x*x)*(y-y*y)* (z - z*z)
      T   = (x-x*x)*(y-y*y)* (z - z*z)
      rho = (sigma*rho + 1.0)*rhoInf 
      u   = (sigma*u + 1.0)*uInf
      v   = (sigma*v + 1.0)*vInf
      w   = (sigma*w + 1.0)*wInf
      T   = (sigma*T + 1.0)*TInf
      p   = rho*T/gamma
      E   = T/(gamma*gamma_1) + 0.5*(u*u + v*v + w*w)

      u0[dofnums_j[1]] = rho
      u0[dofnums_j[2]] = rho*u
      u0[dofnums_j[3]] = rho*v
      u0[dofnums_j[4]] = rho*v
      u0[dofnums_j[5]] = T/(gamma*gamma_1) + 0.5*(u*u + v*v + w*w)
      u0[dofnums_j[5]] *= rho
    end
  end
end

function ICChannel{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                                     operator::AbstractSBP{Tsbp}, 
                                     eqn::EulerData{Tsol, Tsol, 3}, 
                                     opts, 
                                     u0::AbstractVector{Tsol})
  # populate u0 with initial values
  # this is a template for all other initial conditions
  Tdim = 3
  sigma = 0.01
  pi = 3.14159265358979323846264338
  params = eqn.params
  gamma = params.gamma
	gamma_1 = params.gamma_1
	aoa = params.aoa
  beta = params.sideslip_angle
	rhoInf = 1.0
  uInf = params.Ma * cos(beta) * cos(aoa)
  vInf = params.Ma * sin(beta) * -1
  wInf = params.Ma * cos(beta) * sin(aoa)
	TInf = 1.0

  numEl = mesh.numEl
  nnodes = mesh.numNodesPerElement

  for i=1:numEl
    for j=1:nnodes
      coords_j = sview(mesh.coords, :, j, i)
      dofnums_j = sview(mesh.dofs, :, j, i)

      x = coords_j[1]
      y = coords_j[2]
      z = coords_j[3]

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

      u0[dofnums_j[1]] = rho
      u0[dofnums_j[2]] = rho*u
      u0[dofnums_j[3]] = rho*v
      u0[dofnums_j[4]] = rho*w
      u0[dofnums_j[5]] = T/(gamma*gamma_1) + 0.5*(u*u + v*v + w*w)
      u0[dofnums_j[5]] *= rho
    end
  end

  return nothing

end  # end function

function ICChannel{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                                     operator::AbstractSBP{Tsbp}, 
                                     eqn::EulerData{Tsol, Tsol, 2}, 
                                     opts, 
                                     u0::AbstractVector{Tsol})
  # populate u0 with initial values
  # this is a template for all other initial conditions
  sigma = 0.01
  pi = 3.14159265358979323846264338
  gamma = 1.4
  gamma_1 = gamma - 1.0
  aoa = eqn.params.aoa
  rhoInf = 1.0
  uInf = eqn.params.Ma*cos(aoa)
  vInf = eqn.params.Ma*sin(aoa)
  TInf = 1.0

  numEl = mesh.numEl
  nnodes = mesh.numNodesPerElement
  dofpernode = mesh.numDofPerNode
  sol = zeros(Tsol, 4)

  for i=1:numEl
    for j=1:nnodes
      coords_j = sview(mesh.coords, :, j, i)
      dofnums_j = sview(mesh.dofs, :, j, i)

      x = coords_j[1]
      y = coords_j[2]

      calcFreeStream(eqn.params, coords_j, sol)

      rho = rhoInf
      # rho = rhoInf * (0.1*sin(2*pi*x) + 0.1*y +  1.0)
      # u   = uInf * (-4.0 * y * (y-1.0)) + 0.1*uInf
      # u   = uInf * (-4.0 * y * (y-1.0)) 
      ux = (0.1*sin(2*pi*x) + 0.2) * uInf
      uy = sin(pi*y) 
      # uy = -4.0 * y * (y-1.0)
      u  = ux * uy
      v  = vInf 
      T  = TInf 

      if !eqn.params.isViscous
        u += 0.2 * uInf
      end

      u0[dofnums_j[1]] = rho
      u0[dofnums_j[2]] = rho*u
      u0[dofnums_j[3]] = rho*v
      u0[dofnums_j[4]] = T/(gamma*gamma_1) + 0.5*(u*u + v*v)
      u0[dofnums_j[4]] *= rho
    end
  end

  return nothing

end  # end function

function ICDoubleSquare{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                                          operator::AbstractSBP{Tsbp}, 
                                          eqn::EulerData{Tsol, Tsol, 2}, 
                                          opts, 
                                          u0::AbstractVector{Tsol})
  # populate u0 with initial values
  # this is a template for all other initial conditions
  sigma = 0.01
  pi = 3.14159265358979323846264338
  gamma = 1.4
  gamma_1 = gamma - 1.0
  aoa = eqn.params.aoa
  rhoInf = 1.0
  uInf = eqn.params.Ma*cos(aoa)
  vInf = eqn.params.Ma*sin(aoa)
  TInf = 1.0

  numEl = mesh.numEl
  nnodes = mesh.numNodesPerElement
  dofpernode = mesh.numDofPerNode
  sol = zeros(Tsol, 4)

  si = [0.5, 1.5]
  a = [-2.375, 16.875, -45.0, 55.0, -30.0, 6.0]
  # si = [0.5, 2.5]
  # a = [-0.220703125,1.46484375,-3.515625,3.59375,-1.40625,0.1875]
  # si = [0.5, 5.0]
  # a = [-0.01610526850581724,0.10161052685058161,-0.2235431590712794,0.19102779047909357,-0.04470863181425595,0.003251536859218615]
  for i=1:numEl
    for j=1:nnodes
      coords_j = sview(mesh.coords, :, j, i)
      dofnums_j = sview(mesh.dofs, :, j, i)

      x = coords_j[1]
      y = coords_j[2]

      calcFreeStream(eqn.params, coords_j, sol)
      E = sol[4]/sol[1]
      V2 = (sol[2]*sol[2] + sol[3]*sol[3]) / (sol[1]*sol[1])

      gx = 0.0
      gy = 0.0
      if x >= si[1] && x < si[2]
        gx = a[1] + a[2]*x + a[3]*x*x + a[4]*x^3 + a[5]*x^4 + a[6]*x^5
      elseif x >= si[2] 
        gx = 1.0
      elseif x <= -si[1] && x > -si[2]
        gx = a[1] - a[2]*x + a[3]*x*x - a[4]*x^3 + a[5]*x^4 - a[6]*x^5
      elseif x <= -si[2]
        gx = 1.0
      end  

      if y >= si[1] && y < si[2]
        gy = a[1] + a[2]*y + a[3]*y*y + a[4]*y^3 + a[5]*y^4 + a[6]*y^5
      elseif y >= si[2] 
        gy = 1.0
      elseif y <= -si[1] && y > -si[2]
        gy = a[1] - a[2]*y + a[3]*y*y - a[4]*y^3 + a[5]*y^4 - a[6]*y^5
      elseif y <= -si[2] 
        gy = 1.0
      end  

      rho = rhoInf
      u   = uInf * (gx + gy - gx*gy) 
      v   = vInf * (gx + gy - gx*gy) 
      T   = TInf 

      u0[dofnums_j[1]] = rho
      u0[dofnums_j[2]] = rho*u
      u0[dofnums_j[3]] = rho*v
      u0[dofnums_j[4]] = T/(gamma*gamma_1) + 0.5*(u*u + v*v)
      u0[dofnums_j[4]] *= rho
    end
  end

  return nothing

end  # end function

function ICTrigonometric{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                                           operator::AbstractSBP{Tsbp}, 
                                           eqn::EulerData{Tsol, Tsol, 2}, 
                                           opts, 
                                           u0::AbstractVector{Tsol})
  # populate u0 with initial values
  # this is a template for all other initial conditions
  sigma = 0.01
  pi = 3.14159265358979323846264338
  gamma = 1.4
  gamma_1 = gamma - 1.0
  aoa = eqn.params.aoa
  rhoInf = 1.0
  uInf = eqn.params.Ma*cos(aoa)
  vInf = eqn.params.Ma*sin(aoa)
  TInf = 1.0

  numEl = mesh.numEl
  nnodes = mesh.numNodesPerElement
  dofpernode = mesh.numDofPerNode
  for i=1:numEl
    for j=1:nnodes
      coords_j = sview(mesh.coords, :, j, i)
      dofnums_j = sview(mesh.dofs, :, j, i)

      x = coords_j[1]
      y = coords_j[2]

      x2 = 2*x*pi
      y2 = 2*y*pi
      x4 = 4*x*pi
      y4 = 4*y*pi
      sx2 = sin(x2)
      sy2 = sin(y2)
      sx4 = sin(x4)
      sy4 = sin(y4)
      cx2 = cos(x2)
      cx4 = cos(x4)
      cy2 = cos(y2)
      cy4 = cos(y4)
      #
      # Exact solution in form of primitive variables
      #
      rho = 0.25 * sx2 * sy2
      u   = 0.25 * sx4 * sy4
      v   = 0.25 * (cx4  + 1.0) * (cy4 + 1.0)
      T   = 0.25 * (1.0 - cx4) * (1.0 - cy4)
      rho = (sigma*rho + 1.0)*rhoInf 
      u   = (sigma*u + 1.0)*uInf
      v   = (sigma*v + 1.0)*vInf
      T   = (sigma*T + 1.0)*TInf

      u0[dofnums_j[1]] = rho
      u0[dofnums_j[2]] = rho*u
      u0[dofnums_j[3]] = rho*v
      u0[dofnums_j[4]] = T/(gamma*gamma_1) + 0.5*(u*u + v*v)
      u0[dofnums_j[4]]  = u0[dofnums_j[4]] * rho
    end
  end

  return nothing

end  # end function

function ICTrigonometric{Tmsh, Tsbp, Tsol, Tres}(mesh::AbstractMesh{Tmsh}, 
                                                 operator::AbstractSBP{Tsbp}, 
                                                 eqn::EulerData{Tsol, Tres, 3}, 
                                                 opts, 
                                                 u0::AbstractVector{Tsol})
  # populate u0 with initial values
  # this is a template for all other initial conditions
  sigma = 0.01
  pi = 3.14159265358979323846264338
  gamma = 1.4
  gamma_1 = gamma - 1.0
  aoa = eqn.params.aoa
  beta = eqn.params.sideslip_angle
  rhoInf = 1.0
  uInf = eqn.params.Ma * cos(beta) * cos(aoa)
  vInf = eqn.params.Ma * sin(beta) * -1
  wInf = eqn.params.Ma * cos(beta) * sin(aoa)
  TInf = 1.0

  numEl = mesh.numEl
  nnodes = mesh.numNodesPerElement
  dofpernode = mesh.numDofPerNode
  sol = zeros(Tsol, 5)
  for i=1:numEl
    for j=1:nnodes
      xyz = sview(mesh.coords, :, j, i)
      dofnums_j = sview(mesh.dofs, :, j, i)

      calcFreeStream(eqn.params, xyz, sol)

      xyz2 = 2 * pi * xyz
      xyz4 = 4 * pi * xyz
      sin_val_1 = sin(xyz)
      cos_val_1 = cos(xyz)
      sin_val_2 = sin(xyz2)
      cos_val_2 = cos(xyz2)
      sin_val_4 = sin(xyz4)
      cos_val_4 = cos(xyz4)
      #
      # Exact solution in form of primitive variables
      #
      rho = 0.125 * sin_val_2[1] * sin_val_2[2] * sin_val_2[3] 
      u   = 0.125 * sin_val_4[1] * sin_val_4[2] * sin_val_4[3]
      v   = 0.125 * sin_val_2[1] * sin_val_2[2] * sin_val_2[3]
      w   = 0.125 * sin_val_1[1] * sin_val_1[2] * sin_val_1[3] 
      T   = 0.125 * (1.0 - cos_val_4[1]) * (1.0 - cos_val_4[2]) * (1.0 - cos_val_4[3])

      rho = (sigma*rho + 1.0)*rhoInf 
      u = (sigma*u + 1.0) * uInf
      v = (sigma*v + 1.0) * vInf
      w = (sigma*w + 1.0) * wInf
      T = (sigma*T + 1.0) * TInf
      vel2 = u*u + v*v + w*w

      u0[dofnums_j[1]] = rho
      u0[dofnums_j[2]] = rho*u
      u0[dofnums_j[3]] = rho*v
      u0[dofnums_j[4]] = rho*w
      u0[dofnums_j[5]] = rho*(T/(gamma*gamma_1) + 0.5*vel2)
    end
  end

  return nothing
end  # end function

@doc """
### EulerEquationMod.ICZero

  Sets all components of the solution to zero

  Inputs:
    mesh
    sbp
    eqn
    opts

  Inputs/Outputs: 
    u0: vector to populate with the solution

  Aliasing restrictions: none.

"""->
function ICZero{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                operator::AbstractSBP{Tsbp}, eqn::EulerData{Tsol}, opts, 
                u0::AbstractVector{Tsol})
# populate u0 with initial values
# this is a template for all other initial conditions

numEl = mesh.numEl
nnodes = operator.numnodes
dofpernode = mesh.numDofPerNode
for i=1:numEl

  for j=1:nnodes
      coords_j = sview(mesh.coords, :, j, i)
      dofnums_j = sview(mesh.dofs, :, j, i)
 
      # get dof numbers for each variable
      dofnum_rho = dofnums_j[1]
      dofnum_rhou = dofnums_j[2]
      dofnum_rhov = dofnums_j[3]
      dofnum_e = dofnums_j[4]

      x = coords_j[1]
      y = coords_j[2]

      # apply initial conditions here
      u0[dofnum_rho] = 0.0
      u0[dofnum_rhou] = 0.0
      u0[dofnum_rhov] = 0.0
      u0[dofnum_e] = 0.0
  end
end

return nothing

end  # end function


@doc """
### EulerEquationMod.ICOnes

  Sets all components of the solution to 1.0

  Inputs:
    mesh
    sbp
    eqn
    opts

  Inputs/Outputs: 
    u0: vector to populate with the solution

  Aliasing restrictions: none.

"""->

function ICOnes{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                operator::AbstractSBP{Tsbp}, eqn::EulerData{Tsol}, opts,
                u0::AbstractVector{Tsol})

  numEl = mesh.numEl
  nnodes = operator.numnodes
  dofpernode = mesh.numDofPerNode
  for i=1:numEl
    for j=1:nnodes
      coords_j = sview(mesh.coords, :, j, i)
      dofnums_j = sview(mesh.dofs, :, j, i)
 
      # get dof numbers for each variable
      dofnum_rho = dofnums_j[1]
      dofnum_rhou = dofnums_j[2]
      dofnum_rhov = dofnums_j[3]
      dofnum_e = dofnums_j[4]

      x = coords_j[1]
      y = coords_j[2]

      # apply initial conditions here
      u0[dofnum_rho] = 1.0
      u0[dofnum_rhou] = 1.0
      u0[dofnum_rhov] = 1.0
      u0[dofnum_e] = 1.0

      # u0 = 2*u0
    end
  end

  return nothing
end # end function ICOnes


@doc """
### EulerEquationMod.ICRho1E2

  Sets all density values to 1.0 and energy values to 2.0

  Inputs:
    mesh
    sbp
    eqn
    opts

  Inputs/Outputs: 
    u0: vector to populate with the solution

  Aliasing restrictions: none.

"""->


function ICRho1E2{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                  operator::AbstractSBP{Tsbp}, eqn::EulerData{Tsol}, opts, 
                  u0::AbstractVector{Tsol})
# populate u0 with initial values
# this is a template for all other initial conditions

numEl = mesh.numEl
nnodes = operator.numnodes
dofpernode = mesh.numDofPerNode
for i=1:numEl
  for j=1:nnodes
      coords_j = sview(mesh.coords, :, j, i)
      dofnums_j = sview(mesh.dofs, :, j, i)
 
      # get dof numbers for each variable
      dofnum_rho = dofnums_j[1]
      dofnum_rhou = dofnums_j[2]
      dofnum_rhov = dofnums_j[3]
      dofnum_e = dofnums_j[4]

      x = coords_j[1]
      y = coords_j[2]

      # apply initial conditions here
      u0[dofnum_rho] = 1.0
      u0[dofnum_rhou] = 0.0
      u0[dofnum_rhov] = 0.0
      u0[dofnum_e] = 2.0
  end
end

return nothing

end  # end function

@doc """
### EulerEquationMod.ICRho1E2U1VW0

  Sets the density values 1.0, x momentum to 1.0, 
  v & w momenta to 0.0, and energy to 2.0 at a node.

  It should work for 2D and 3D meshes.

  Inputs:
    mesh
    sbp
    eqn
    opts

  Inputs/Outputs: 
    u0: vector to populate with the solution

  Aliasing restrictions: none.

"""->
function ICRho1E2U1VW0{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                    operator::AbstractSBP{Tsbp}, eqn::EulerData{Tsol}, 
                    opts, u0::AbstractVector{Tsol})
# populate u0 with initial values
# this is a template for all other initial conditions


numEl = mesh.numEl
nnodes = mesh.numNodesPerElement
dofpernode = mesh.numDofPerNode
sol = zeros(Tsol, mesh.numDofPerNode)
for i=1:numEl
  for j=1:nnodes

      coords_j = sview(mesh.coords, :, j, i)
      dofnums_j = sview(mesh.dofs, :, j, i)
      # get dof numbers for each variable

      calcRho1Energy2U1VW0(eqn.params, coords_j, sol)

      for k=1:dofpernode
        u0[dofnums_j[k]] = sol[k]
      end
  end
end

return nothing

end  # end function


@doc """
### EulerEquationMod.ICRho1E2U3

  Sets all components density values to 1.0, x and y momenta to 0.35355, and
  energy to 2.0

  Inputs:
    mesh
    sbp
    eqn
    opts

  Inputs/Outputs: 
    u0: vector to populate with the solution

  Aliasing restrictions: none.

"""->
function ICRho1E2U3{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                    operator::AbstractSBP{Tsbp}, eqn::EulerData{Tsol}, 
                    opts, u0::AbstractVector{Tsol})
# populate u0 with initial values
# this is a template for all other initial conditions


numEl = mesh.numEl
nnodes = mesh.numNodesPerElement
dofpernode = mesh.numDofPerNode
sol = zeros(Tsol, mesh.numDofPerNode)
for i=1:numEl
  for j=1:nnodes

      coords_j = sview(mesh.coords, :, j, i)
      dofnums_j = sview(mesh.dofs, :, j, i)
      # get dof numbers for each variable

      calcRho1Energy2U3(eqn.params, coords_j, sol)

      for k=1:dofpernode
        u0[dofnums_j[k]] = sol[k]
      end
  end
end

return nothing

end  # end function

@doc """
### EulerEquationMod.ICFreeStream

  Sets all components of the solution to the free stream condition according
  to the angle of attack and and Mach number.

  Inputs:
    mesh
    sbp
    eqn
    opts

  Inputs/Outputs: 
    u0: vector to populate with the solution

  Aliasing restrictions: none.

"""->


function ICFreeStream{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                      operator::AbstractSBP{Tsbp}, eqn::EulerData{Tsol}, opts, 
                      u0::AbstractVector{Tsol})
# populate u0 with initial values
# this is a template for all other initial conditions


numEl = mesh.numEl
nnodes = mesh.numNodesPerElement
dofpernode = mesh.numDofPerNode
sol = zeros(Tsol, dofpernode)
for i=1:numEl
  for j=1:nnodes
      coords_j = sview(mesh.coords, :, j, i)
      dofnums_j = sview(mesh.dofs, :, j, i)
 
      calcFreeStream(eqn.params, coords_j, sol)

      for k=1:dofpernode
        u0[dofnums_j[k]] = sol[k]
      end

  end
end

return nothing

end  # end function




# what is this? how is it different than ICIsentropic Vortex?
function ICVortex{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                  operator::AbstractSBP{Tsbp}, eqn::EulerData{Tsol}, opts, 
                  u0::AbstractVector{Tsol})
# populate u0 with initial values
# this is a template for all other initial conditions

numEl = mesh.numEl
nnodes = operator.numnodes
dofpernode = mesh.numDofPerNode
sol = zeros(Tsol, 4)
for i=1:numEl
  for j=1:nnodes
      coords_j = sview(mesh.coords, :, j, i)
      dofnums_j = sview(mesh.dofs, :, j, i)
 
      # get dof numbers for each variable
      dofnum_rho = dofnums_j[1]
      dofnum_rhou = dofnums_j[2]
      dofnum_rhov = dofnums_j[3]
      dofnum_e = dofnums_j[4]

      x = coords_j[1]
      y = coords_j[2]

      calcVortex(eqn.params, coords_j, sol)


      for k=1:dofpernode
        u0[dofnums_j[k]] = sol[k]
      end

  end
end

return nothing

end  # end function

@doc """
### EulerEquationMod.ICsmoothHeavisideder

  Sets the density to the derivative of the smooth Heaviside function, all 
  other components to zero.

  Inputs:
    mesh
    sbp
    eqn
    opts

  Inputs/Outputs: 
    u0: vector to populate with the solution

  Aliasing restrictions: none.

"""->
function ICsmoothHeavisideder{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                              operator::AbstractSBP{Tsbp}, eqn::EulerData{Tsol},
                              opts, u0::AbstractVector{Tsol})
# calculate the value of the smooth heaviside function derivative at a location x
# x0 is specified within this function

# smooth heaviside  parameters
  x0 = 0
  L = 5
  k = 5



numEl = mesh.numEl
nnodes = operator.numnodes
dofpernode = mesh.numDofPerNode
for i=1:numEl
#  dofnums_i = sview(mesh, i)  # get dof nums for this element
#  coords = sview(mesh, [i])

  for j=1:nnodes
      coords_j = sview(mesh.coords, :, j, i)
      dofnums_j = sview(mesh.dofs, :, j, i)
 
      # get dof numbers for each variable
      dofnum_rho = dofnums_j[1]
      dofnum_rhou = dofnums_j[2]
      dofnum_rhov = dofnums_j[3]
      dofnum_e = dofnums_j[4]

      x = coords_j[1]
      y = coords_j[2]

      # apply initial conditions here
      u0[dofnum_rho] = L*(2*k*e^(-2*k*x))/(e^(-2*k*x) +1 )^2
      u0[dofnum_rhou] = 0.0
      u0[dofnum_rhov] = 0.0
      u0[dofnum_e] = 0.0
  end
end

return nothing



end


@doc """
### EulerEquationMod.ICZero

  Sets the density to the smooth Heaviside function, all other components to
  zero.

  Inputs:
    mesh
    sbp
    eqn
    opts

  Inputs/Outputs: 
    u0: vector to populate with the solution

  Aliasing restrictions: none.

"""->
function ICsmoothHeaviside{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                                             operator::AbstractSBP{Tsbp}, 
                                             eqn::EulerData{Tsol}, 
                                             opts, u0::AbstractArray{Tsol, 1})
  # calculate the value of the smooth heaviside function at a location x
  # x0 is specified within this function

  # smooth heaviside  parameters
  x0 = 0
  L = 5
  k = 5



  numEl = mesh.numEl
  nnodes = operator.numnodes
  dofpernode = mesh.numDofPerNode
  for i=1:numEl
    for j=1:nnodes
      coords = sview(mesh.coords, :, j, i)
      dofnums_j = sview(mesh.dofs, :, j, i)

      # get dof numbers for each variable
      dofnum_rho = dofnums_j[1]
      dofnum_rhou = dofnums_j[2]
      dofnum_rhov = dofnums_j[3]
      dofnum_e = dofnums_j[4]

      x = coords_j[1]
      y = coords_j[2]

      # apply initial conditions here
      u0[dofnum_rho] = L/(1 + e^(-k*(x-x0)))
      u0[dofnum_rhou] = 0.0
      u0[dofnum_rhov] = 0.0
      u0[dofnum_e] = 0.0
    end
  end

  return nothing
end

@doc """
### EulerEquationMod.ICIsentropicVortex

Sets the solution to the steady isentropic vortex solution.

Inputs:
mesh
    sbp
    eqn
    opts

  Inputs/Outputs: 
    u0: vector to populate with the solution

  Aliasing restrictions: none.

"""->
function ICIsentropicVortex{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                                              operator::AbstractSBP{Tsbp}, 
                                              eqn::EulerData{Tsol}, 
                                              opts, u0::AbstractArray{Tsol})
  # populate u0 with initial values
  # this is a template for all other initial conditions

  println("entered ICIsentropicVortex")

  numEl = mesh.numEl
  nnodes = operator.numnodes
  dofpernode = mesh.numDofPerNode
  sol = zeros(Tsol, dofpernode)
  for i=1:numEl
    #  println("i = ", i)
    #  coords = sview(mesh, [i])

    for j=1:nnodes
      dofnums_j = sview(mesh.dofs, :, j, i)
      coords_j = sview(mesh.coords, :, j, i)
      calcIsentropicVortex(eqn.params, coords_j, sol)

      # apply initial conditions here
      for k=1:dofpernode
        u0[dofnums_j[k]] = sol[k]
      end

    end
  end

  return nothing

end  # end function


@doc """
### EulerEquationMod.ICIsentropicVortexWithNoise

  Sets the solution to the steady isentropic vortex solution plus 
  a small random noise component.

  Inputs:
    mesh
    sbp
    eqn
    opts

  Inputs/Outputs: 
    u0: vector to populate with the solution

  Aliasing restrictions: none.

"""->
function ICIsentropicVortexWithNoise{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh},
                                     operator::AbstractSBP{Tsbp}, 
                                     eqn::EulerData{Tsol}, 
                                     opts, u0::AbstractVector{Tsol})
  # populate u0 with initial values
  # this is a template for all other initial conditions

  numEl = mesh.numEl
  nnodes = operator.numnodes
  dofpernode = mesh.numDofPerNode
  sol = zeros(Tsol, 4)
  for i=1:numEl
    for j=1:nnodes
      coords_j = sview(mesh.coords, :, j, i)
      dofnums_j = sview(mesh.dofs, :, j, i)
      calcIsentropicVortex(eqn.params, coords_j, sol)

      # apply initial conditions here
      for k=1:dofpernode
        u0[dofnums_j[k]] = sol[k] + 0.1*rand()
      end

    end
  end

  return nothing

end  # end function

@doc """
### EulerEquationMod.ICUnsteadyVortex

  Sets the solution to the unsteady vortex problem.  eqn.params.t is used to
  determine what time to use for the solution.

  Inputs:
    mesh
    sbp
    eqn
    opts

  Inputs/Outputs: 
    u0: vector to populate with the solution

  Aliasing restrictions: none.

"""->
function ICUnsteadyVortex{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                          operator::AbstractSBP{Tsbp}, eqn::EulerData{Tsol}, 
                          opts, u0::AbstractArray{Tsol})
  # populate u0 with initial values
  # this is a template for all other initial conditions

  numEl = mesh.numEl
  nnodes = operator.numnodes
  dofpernode = mesh.numDofPerNode
  sol = zeros(Tsol, mesh.numDofPerNode)


  for i=1:numEl
    for j=1:nnodes
      dofnums_j = sview(mesh.dofs, :, j, i)

      coords_j = sview(mesh.coords, :, j, i)
      calcUnsteadyVortex(eqn.params, coords_j, sol)

      for k=1:dofpernode
        u0[dofnums_j[k]] = sol[k]
      end

    end
  end

  return nothing

end  # end function

"""
  Vortex travelling at an angle
"""
function ICUnsteadyVortex2{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                          operator::AbstractSBP{Tsbp}, eqn::EulerData{Tsol}, 
                          opts, u0::AbstractArray{Tsol})
  # populate u0 with initial values
  # this is a template for all other initial conditions

  numEl = mesh.numEl
  nnodes = operator.numnodes
  dofpernode = mesh.numDofPerNode
  sol = zeros(Tsol, mesh.numDofPerNode)


  for i=1:numEl
    for j=1:nnodes
      dofnums_j = sview(mesh.dofs, :, j, i)

      coords_j = sview(mesh.coords, :, j, i)
      calcUnsteadyVortex2(eqn.params, coords_j, sol)

      for k=1:dofpernode
        u0[dofnums_j[k]] = sol[k]
      end

    end
  end

  return nothing

end  # end function



@doc """
### EulerEquationMod.ICFile

  This function reads a vector from a file on disk and set the solution to it.
  The vector must contain the same number of entries as there are degrees of 
  freedom in the mesh. 

  This function is useful for things like restarting from a checkpoint.
  In this case, the file should be the output of writedlm(eqn.q).  The degree 
  of freedom number must be the same for both simulation for this to work (the 
  file contains no degree of freedom number information).


  Inputs:
    mesh
    sbp
    eqn
    opts

  Inputs/Outputs: 
    u0: vector to populate with the solution

  Aliasing restrictions: none.

"""->
function ICFile{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                operator::AbstractSBP{Tsbp}, eqn::EulerData{Tsol}, opts, 
                u0::AbstractVector{Tsol})
  # populate u0 with initial values from a disk file
  # the file name comes from opts["ICfname"]

  fname = get_parallel_fname(opts["ICfname"], mesh.myrank)
  vals = readdlm(fname)

  @assert length(vals) == mesh.numDof

  for i=1:mesh.numDof
    u0[i] = vals[i]
  end

end

"""
  Assigns exp(k*x*y*z) as the initial condition, of each node, where k is 
  the index of the degree of freedom of the node
"""
function ICExp{Tmsh, Tsol,}(mesh::AbstractMesh{Tmsh}, sbp, eqn::EulerData{Tsol}, opts, u0::AbstractVector{Tsol})

  q = eqn.params.q_vals
  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      dofs = sview(mesh.dofs, :, j, i)
      coords = sview(mesh.coords, :, j, i)
      calcExp(eqn.params, coords, q)
      for k=1:mesh.numDofPerNode
        u0[dofs[k]] = q[k]
      end
    end
  end

  return nothing
end

"""
  Writes calcPeriodicMMS to the initial condition vector u0
"""
function ICPeriodicMMS{Tmsh, Tsol,}(mesh::AbstractMesh{Tmsh}, sbp, eqn::EulerData{Tsol}, opts, u0::AbstractVector{Tsol})

  q = eqn.params.q_vals
  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      dofs = sview(mesh.dofs, :, j, i)
      coords = sview(mesh.coords, :, j, i)
      calcPeriodicMMS(eqn.params, coords, q)
      for k=1:mesh.numDofPerNode
        u0[dofs[k]] = q[k]
      end
    end
  end

  return nothing
end

"""
  This function applies the initial condition for the Taylor Green vortex,
  using the constants in Gassner, Winters, and Kopriva's Split form Nodal
  DG paper
"""
function ICTaylorGreen{Tmsh, Tsol,}(mesh::AbstractMesh{Tmsh}, sbp, 
                       eqn::EulerData{Tsol}, opts, u0::AbstractVector{Tsol})

  # parameters
  M = 1.0  # Mach number
  gamma_1 = eqn.params.gamma_1
  gamma = eqn.params.gamma

  q = eqn.params.q_vals
  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      dofs = sview(mesh.dofs, :, j, i)
      coords = sview(mesh.coords, :, j, i)
      x = coords[1]
      y = coords[2]
      z = coords[3]

      p = 100/gamma + (1/16)*( cos(2*x)*cos(2*z) + 2*cos(2*y) + 2*cos(2*x) + cos(2*y)*cos(2*z) )
      q[1] = 1
      q[2] = q[1]*M*sin(x)*cos(y)*cos(z)
      q[3] = -q[1]*M*cos(x)*sin(y)*cos(z)
      q[4] = 0
      q[5] = p/gamma_1 + (q[2]*q[2] + q[3]*q[3] + q[4]*q[4])/(2*q[1])

      for k=1:mesh.numDofPerNode
        u0[dofs[k]] = q[k]
      end
    end
  end

  return nothing
end

"""
  Initial condition of channel MMS
"""
function ICChannelMMS{Tmsh, Tsol,}(mesh::AbstractMesh{Tmsh}, sbp, eqn::EulerData{Tsol}, opts, u0::AbstractVector{Tsol})

  q = eqn.params.q_vals
  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      dofs = sview(mesh.dofs, :, j, i)
      coords = sview(mesh.coords, :, j, i)
      calcChannelMMS(eqn.params, coords, q)
      for k=1:mesh.numDofPerNode
        u0[dofs[k]] = q[k]
      end
    end
  end

  return nothing
end

"""
  Initial for square wave in 1D
"""
function ICSquare1D{Tmsh, Tsol,}(mesh::AbstractMesh{Tmsh}, sbp, eqn::EulerData{Tsol}, opts, u0::AbstractVector{Tsol})

  q = eqn.params.q_vals
  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      dofs = sview(mesh.dofs, :, j, i)
      coords = sview(mesh.coords, :, j, i)
      calcSquare1D(eqn.params, coords, q)
      for k=1:mesh.numDofPerNode
        u0[dofs[k]] = q[k]
      end
    end
  end

  return nothing
end

 """
  Initial for square wave in 2D
"""
function ICSquare2D{Tmsh, Tsol,}(mesh::AbstractMesh{Tmsh}, sbp, eqn::EulerData{Tsol}, opts, u0::AbstractVector{Tsol})

  q = eqn.params.q_vals
  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      dofs = sview(mesh.dofs, :, j, i)
      coords = sview(mesh.coords, :, j, i)
      calcSquare2D(eqn.params, coords, q)
      for k=1:mesh.numDofPerNode
        u0[dofs[k]] = q[k]
      end
    end
  end

  return nothing
end

 """
  Initial for square wave in 2D
"""
function ICSedovExplosion{Tmsh, Tsol,}(mesh::AbstractMesh{Tmsh}, sbp,
                          eqn::EulerData{Tsol}, opts, u0::AbstractVector{Tsol})

  q = eqn.params.q_vals
  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      dofs = sview(mesh.dofs, :, j, i)
      coords = sview(mesh.coords, :, j, i)
      calcSedovExplosion(eqn.params, coords, q)
      for k=1:mesh.numDofPerNode
        u0[dofs[k]] = q[k]
      end
    end
  end

  return nothing
end

 



# declare a const dictionary here that maps strings to function (used for input arguments)
"""
  Map IC names to functions.  Generally the name is the same as the function
  name
"""
# global const ICDict = Dict{Any, Function}(
  # "ICZero" => ICZero,
  # "ICOnes" => ICOnes,
  # "ICRho1E2" => ICRho1E2,
  # "ICRho1E2U1VW0" => ICRho1E2U1VW0,
  # "ICRho1E2U3" => ICRho1E2U3,
  # "ICFreeStream" => ICFreeStream,
  # "ICVortex" => ICVortex,
  # #"ICLinear" => ICLinear,
  # "ICsmoothHeavisideder" => ICsmoothHeavisideder,
  # "ICsmoothHeaviside" => ICsmoothHeaviside,
  # "ICIsentropicVortex" => ICIsentropicVortex,
  # "ICUnsteadyVortex" => ICUnsteadyVortex,
  # "ICUnsteadyVortex2" => ICUnsteadyVortex2,
  # "ICIsentropicVortexWithNoise" => ICIsentropicVortexWithNoise,
  # "ICFile" => ICFile,
  # "ICExp" => ICExp,
  # "ICPeriodicMMS" => ICPeriodicMMS,
  # "ICTaylorGreen" => ICTaylorGreen,
  # "ICChannelMMS" => ICChannelMMS,
  # "ICSquare1D" => ICSquare1D,
  # "ICSquare2D" => ICSquare2D,
  # "ICSedovExplosion" => ICSedovExplosion,
  # "ICTrigonometric" => ICTrigonometric,
  # "ICPolynomial" => ICPolynomial,
  # "ICChannel" => ICChannel,
  # "ICDoubleSquare" => ICDoubleSquare,
# )


