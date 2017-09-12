# functions that populate the initial conditions
# List of functions:

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
                           operator::AbstractSBP{Tsbp}, eqn::EulerData{Tsol}, 
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
                            operator::AbstractSBP{Tsbp}, eqn::EulerData{Tsol}, 
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

global const ICDict = Dict{Any, Function}(
"ICZero" => ICZero,
"ICOnes" => ICOnes,
"ICRho1E2" => ICRho1E2,
"ICRho1E2U3" => ICRho1E2U3,
"ICFreeStream" => ICFreeStream,
"ICVortex" => ICVortex,
#"ICLinear" => ICLinear,
"ICsmoothHeavisideder" => ICsmoothHeavisideder,
"ICsmoothHeaviside" => ICsmoothHeaviside,
"ICIsentropicVortex" => ICIsentropicVortex,
"ICUnsteadyVortex" => ICUnsteadyVortex,
"ICIsentropicVortexWithNoise" => ICIsentropicVortexWithNoise,
"ICFile" => ICFile,
"ICExp" => ICExp,
"ICPeriodicMMS" => ICPeriodicMMS,
"ICTaylorGreen" => ICTaylorGreen,
"ICChannelMMS" => ICChannelMMS,
"ICSquare1D" => ICSquare1D,
"ICSquare2D" => ICSquare2D,
"ICSedovExplosion" => ICSedovExplosion,
)


