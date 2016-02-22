# functions that populate the initial conditions
# List of functions:

export ICDict

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
                operator::SBPOperator{Tsbp}, eqn::EulerData{Tsol}, opts, 
                u0::AbstractVector{Tsol})
# populate u0 with initial values
# this is a template for all other initial conditions

numEl = getNumEl(mesh)
nnodes = operator.numnodes
dofpernode = getNumDofPerNode(mesh)
for i=1:numEl

  for j=1:nnodes
      coords_j = view(mesh.coords, :, j, i)
      dofnums_j = view(mesh.dofs, :, j, i)
 
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
                operator::SBPOperator{Tsbp}, eqn::EulerData{Tsol}, opts,
                u0::AbstractVector{Tsol})

  numEl = getNumEl(mesh)
  nnodes = operator.numnodes
  dofpernode = getNumDofPerNode(mesh)
  for i=1:numEl
    for j=1:nnodes
      coords_j = view(mesh.coords, :, j, i)
      dofnums_j = view(mesh.dofs, :, j, i)
 
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
                  operator::SBPOperator{Tsbp}, eqn::EulerData{Tsol}, opts, 
                  u0::AbstractVector{Tsol})
# populate u0 with initial values
# this is a template for all other initial conditions

numEl = getNumEl(mesh)
nnodes = operator.numnodes
dofpernode = getNumDofPerNode(mesh)
for i=1:numEl
  for j=1:nnodes
      coords_j = view(mesh.coords, :, j, i)
      dofnums_j = view(mesh.dofs, :, j, i)
 
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
                    operator::SBPOperator{Tsbp}, eqn::EulerData{Tsol}, 
                    opts, u0::AbstractVector{Tsol})
# populate u0 with initial values
# this is a template for all other initial conditions

numEl = mesh.numEl
nnodes = mesh.numNodesPerElement
dofpernode = mesh.numDofPerNode
sol = zeros(Tsol, 4)
for i=1:numEl
  for j=1:nnodes

      coords_j = view(mesh.coords, :, j, i)
      dofnums_j = view(mesh.dofs, :, j, i)
      # get dof numbers for each variable
      dofnum_rho = dofnums_j[1]
      dofnum_rhou = dofnums_j[2]
      dofnum_rhov = dofnums_j[3]
      dofnum_e = dofnums_j[4]

      x = coords_j[1]
      y = coords_j[2]

      calcRho1Energy2U3(coords_j, eqn.params, sol)

      sol[2] += 0*sin(x)  # add a perturbation

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
                      operator::SBPOperator{Tsbp}, eqn::EulerData{Tsol}, opts, 
                      u0::AbstractVector{Tsol})
# populate u0 with initial values
# this is a template for all other initial conditions

numEl = mesh.numEl
nnodes = mesh.numNodesPerElement
dofpernode = mesh.numDofPerNode
sol = zeros(Tsol, 4)
for i=1:numEl
  for j=1:nnodes
      coords_j = view(mesh.coords, :, j, i)
      dofnums_j = view(mesh.dofs, :, j, i)
 
      # get dof numbers for each variable
      dofnum_rho = dofnums_j[1]
      dofnum_rhou = dofnums_j[2]
      dofnum_rhov = dofnums_j[3]
      dofnum_e = dofnums_j[4]

      x = coords_j[1]
      y = coords_j[2]

      calcFreeStream(coords_j, eqn.params, sol)

      for k=1:dofpernode
        u0[dofnums_j[k]] = sol[k]
      end

  end
end

return nothing

end  # end function




# what is this? how is it different than ICIsentropic Vortex?
function ICVortex{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                  operator::SBPOperator{Tsbp}, eqn::EulerData{Tsol}, opts, 
                  u0::AbstractVector{Tsol})
# populate u0 with initial values
# this is a template for all other initial conditions

numEl = getNumEl(mesh)
nnodes = operator.numnodes
dofpernode = getNumDofPerNode(mesh)
sol = zeros(Tsol, 4)
for i=1:numEl
  for j=1:nnodes
      coords_j = view(mesh.coords, :, j, i)
      dofnums_j = view(mesh.dofs, :, j, i)
 
      # get dof numbers for each variable
      dofnum_rho = dofnums_j[1]
      dofnum_rhou = dofnums_j[2]
      dofnum_rhov = dofnums_j[3]
      dofnum_e = dofnums_j[4]

      x = coords_j[1]
      y = coords_j[2]

      calcVortex(coords_j, eqn.params, sol)


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
                              operator::SBPOperator{Tsbp}, eqn::EulerData{Tsol},
                              opts, u0::AbstractVector{Tsol})
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
#  dofnums_i = view(mesh, i)  # get dof nums for this element
#  coords = view(mesh, [i])

  for j=1:nnodes
      coords_j = view(mesh.coords, :, j, i)
      dofnums_j = view(mesh.dofs, :, j, i)
 
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
                           operator::SBPOperator{Tsbp}, eqn::EulerData{Tsol}, 
                           opts, u0::AbstractArray{Tsol, 1})
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
  for j=1:nnodes
      coords = view(mesh.coords, :, j, i)
      dofnums_j = view(mesh.dofs, :, j, i)
 
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
### EulerEquationMod.ICZero

  Sets the solution to the isentropic vortex solution.

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
                            operator::SBPOperator{Tsbp}, eqn::EulerData{Tsol}, 
                            opts, u0::AbstractArray{Tsol})
# populate u0 with initial values
# this is a template for all other initial conditions

println("entered ICIsentropicVortex")

numEl = getNumEl(mesh)
nnodes = operator.numnodes
dofpernode = getNumDofPerNode(mesh)
sol = zeros(Tsol, 4)
for i=1:numEl
#  println("i = ", i)
#  coords = view(mesh, [i])

  for j=1:nnodes
      dofnums_j = view(mesh.dofs, :, j, i)
      coords_j = view(mesh.coords, :, j, i)
      calcIsentropicVortex(coords_j, eqn.params, sol)

      # apply initial conditions here
      for k=1:dofpernode
        u0[dofnums_j[k]] = sol[k]
      end

  end
end

return nothing

end  # end function


@doc """
### EulerEquationMod.ICZero

  Sets the solutoin to the isentropic vortex solution plus a small random noise
  component.

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
                                     operator::SBPOperator{Tsbp}, 
                                     eqn::EulerData{Tsol}, 
                                     opts, u0::AbstractVector{Tsol})
# populate u0 with initial values
# this is a template for all other initial conditions

numEl = getNumEl(mesh)
nnodes = operator.numnodes
dofpernode = getNumDofPerNode(mesh)
sol = zeros(Tsol, 4)
for i=1:numEl
  for j=1:nnodes
      coords_j = view(mesh.coords, :, j, i)
      dofnums_j = view(mesh.dofs, :, j, i)
      calcIsentropicVortex(coords_j, eqn.params, sol)

      # apply initial conditions here
      for k=1:dofpernode
        u0[dofnums_j[k]] = sol[k] + 0.1*rand()
      end

  end
end

return nothing

end  # end function

@doc """
### EulerEquationMod.ICZero

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
                          operator::SBPOperator{Tsbp}, eqn::EulerData{Tsol}, 
                          opts, u0::AbstractArray{Tsol})
# populate u0 with initial values
# this is a template for all other initial conditions

println("entered ICIsentropicVortex")

numEl = getNumEl(mesh)
nnodes = operator.numnodes
dofpernode = getNumDofPerNode(mesh)
sol = zeros(Tsol, 4)
for i=1:numEl
  for j=1:nnodes
      dofnums_j = view(mesh.dofs, :, j, i)
 
      coords_j = view(mesh.coords, :, j, i)
      calcUnsteadyVortex(coords_j, eqn.params, sol)

      for k=1:dofpernode
        u0[dofnums_j[k]] = sol[k]
      end

  end
end

return nothing

end  # end function


@doc """
### EulerEquationMod.ICZero

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
                operator::SBPOperator{Tsbp}, eqn::EulerData{Tsol}, opts, 
                u0::AbstractVector{Tsol})
# populate u0 with initial values from a disk file
# the file name comes from opts["ICfname"]

fname = opts["ICfname"]
vals = readdlm(fname)

@assert length(vals) == mesh.numDof

for i=1:mesh.numDof
  u0[i] = vals[i]
end

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
"ICFile" => ICFile
)


