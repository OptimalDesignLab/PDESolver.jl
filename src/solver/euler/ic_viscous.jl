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


