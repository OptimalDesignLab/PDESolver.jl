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
#  dofnums_i = view(mesh.dofs, :, :, i)  # get dof nums for this element
#  coords = view(mesh.coords, :, :, i)

  for j=1:nnodes
      coords = view(mesh.coords, :, j, i)
      dofnums_j = view(mesh.dofs, :, j, i)  # get dof nums for this element
 
      # get dof numbers for each variable
      dofnum_rho = dofnums_j[1]
      dofnum_rhou = dofnums_j[2]
      dofnum_rhov = dofnums_j[3]
      dofnum_e = dofnums_j[4]

      # coordinates of this node (must be a vertex)
      x = coords[1]
      y = coords[2]
#      z = coords[3,j]

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
### EulerEquationMod.ICZero

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
#  dofnums_i = view(mesh.dofs, :, :, i)  # get dof nums for this element
#  coords = view(mesh.coords, :, :, i)


#    dofnums_i = view(mesh, i)  # get dof nums for this element
#    coords = view(mesh, [i])

    for j=1:nnodes
      coords = view(mesh.coords, :, j, i)
      dofnums_j = view(mesh.dofs, :, j, i)  # get dof nums for this element
 
      # get dof numbers for each variable
      dofnum_rho = dofnums_j[1]
      dofnum_rhou = dofnums_j[2]
      dofnum_rhov = dofnums_j[3]
      dofnum_e = dofnums_j[4]

      # coordinates of this node (must be a vertex)
      x = coords[1]
      y = coords[2]
      z = coords[3,j]

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

function ICRho1E2{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, operator::SBPOperator{Tsbp}, eqn::EulerData{Tsol}, opts, u0::AbstractVector{Tsol})
# populate u0 with initial values
# this is a template for all other initial conditions

numEl = getNumEl(mesh)
nnodes = operator.numnodes
dofpernode = getNumDofPerNode(mesh)
for i=1:numEl
#  dofnums_i = view(mesh.dofs, :, :, i)  # get dof nums for this element
#  coords = view(mesh.coords, :, :, i)


#  dofnums_i = view(mesh, i)  # get dof nums for this element
#  coords = view(mesh, [i])

  for j=1:nnodes
      coords = view(mesh.coords, :, j, i)
      dofnums_j = view(mesh.dofs, :, j, i)  # get dof nums for this element
 
      # get dof numbers for each variable
      dofnum_rho = dofnums_j[1]
      dofnum_rhou = dofnums_j[2]
      dofnum_rhov = dofnums_j[3]
      dofnum_e = dofnums_j[4]

      # coordinates of this node (must be a vertex)
      x = coords[1]
      y = coords[2]

      # apply initial conditions here
      u0[dofnum_rho] = 1.0
      u0[dofnum_rhou] = 0.0
      u0[dofnum_rhov] = 0.0
      u0[dofnum_e] = 2.0
  end
end

return nothing

end  # end function


function ICRho1E2U3{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, operator::SBPOperator{Tsbp}, eqn::EulerData{Tsol}, opts, u0::AbstractVector{Tsol})
# populate u0 with initial values
# this is a template for all other initial conditions

numEl = mesh.numEl
nnodes = mesh.numNodesPerElement
dofpernode = mesh.numDofPerNode
sol = zeros(Tsol, 4)
for i=1:numEl


#  dofnums_i = mesh.dofs[:, :, i]
#  coords = mesh.coords[:, :, i]

  for j=1:nnodes

      coords = view(mesh.coords, :, j, i)
      dofnums_j = view(mesh.dofs, :, j, i)  # get dof nums for this element
      # get dof numbers for each variable
      dofnum_rho = dofnums_j[1]
      dofnum_rhou = dofnums_j[2]
      dofnum_rhov = dofnums_j[3]
      dofnum_e = dofnums_j[4]

      # coordinates of this node (must be a vertex)
      x = coords[1]
      y = coords[2]
#      z = coords[3,j]

      calcRho1Energy2U3(coords, eqn.params, sol)

      sol[2] += 0*sin(x)  # add a perturbation


      # apply initial conditions here
#      u0[dofnum_rho] = 1.0
#      u0[dofnum_rhou] = 3.0
#      u0[dofnum_rhov] = 0.0
#      u0[dofnum_e] = 2.0

      u0[dofnums_j] = sol
  end
end

return nothing

end  # end function


function ICFreeStream{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, operator::SBPOperator{Tsbp}, eqn::EulerData{Tsol}, opts, u0::AbstractVector{Tsol})
# populate u0 with initial values
# this is a template for all other initial conditions

numEl = mesh.numEl
nnodes = mesh.numNodesPerElement
dofpernode = mesh.numDofPerNode
sol = zeros(Tsol, 4)
for i=1:numEl
#  dofnums_i = view(mesh.dofs, :, :, i)  # get dof nums for this element
#  coords = view(mesh.coords, :, :, i)


#  dofnums_i = mesh.dofs[:, :, i]
#  coords = mesh.coords[:, :, i]

  for j=1:nnodes
      coords = view(mesh.coords, :, j, i)
      dofnums_j = view(mesh.dofs, :, j, i)  # get dof nums for this element
 
      # get dof numbers for each variable
      dofnum_rho = dofnums_j[1]
      dofnum_rhou = dofnums_j[2]
      dofnum_rhov = dofnums_j[3]
      dofnum_e = dofnums_j[4]

      # coordinates of this node (must be a vertex)
      x = coords[1]
      y = coords[2]
#      z = coords[3,j]

      calcFreeStream(coords, eqn.params, sol)

      # apply initial conditions here
#      u0[dofnum_rho] = 1.0
#      u0[dofnum_rhou] = 3.0
#      u0[dofnum_rhov] = 0.0
#      u0[dofnum_e] = 2.0

      u0[dofnums_j] = sol
  end
end

return nothing

end  # end function




# what is this? how is it different than ICIsentropic Vortex?
function ICVortex{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, operator::SBPOperator{Tsbp}, eqn::EulerData{Tsol}, opts, u0::AbstractVector{Tsol})
# populate u0 with initial values
# this is a template for all other initial conditions

numEl = getNumEl(mesh)
nnodes = operator.numnodes
dofpernode = getNumDofPerNode(mesh)
sol = zeros(Tsol, 4)
for i=1:numEl
#  dofnums_i = view(mesh.dofs, :, :, i)  # get dof nums for this element
#  coords = view(mesh.coords, :, :, i)


#  dofnums_i = view(mesh, i)  # get dof nums for this element
#  coords = view(mesh, [i])

  for j=1:nnodes
      coords = view(mesh.coords, :, j, i)
      dofnums_j = view(mesh.dofs, :, j, i)  # get dof nums for this element
 
      # get dof numbers for each variable
      dofnum_rho = dofnums_j[1]
      dofnum_rhou = dofnums_j[2]
      dofnum_rhov = dofnums_j[3]
      dofnum_e = dofnums_j[4]

      # coordinates of this node (must be a vertex)
      x = coords[1]
      y = coords[2]

      calcVortex(coords[:,j], eqn.params, sol)


      # apply initial conditions here
#      u0[dofnum_rho] = 1.0
#      u0[dofnum_rhou] = 3.0
#      u0[dofnum_rhov] = 0.0
#      u0[dofnum_e] = 2.0

      u0[dofnums_j] = sol
  end
end

return nothing

end  # end function



#=
function ICLinear{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, operator::SBPOperator{Tsbp}, eqn::EulerData{Tsol}, opts, u0::AbstractArray{Tsol,1})
# populate u0 with initial values
# this is a template for all other initial conditions

nnodes = operator.numnodes
dofpernode = getNumDofPerNode(mesh)
dofnums_i = zeros(Int, dofpernode)

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
=#

function ICsmoothHeavisideder{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, operator::SBPOperator{Tsbp}, eqn::EulerData{Tsol}, opts, u0::AbstractVector{Tsol})
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
      coords = view(mesh.coords, :, j, i)
      dofnums_j = view(mesh.dofs, :, j, i)  # get dof nums for this element
 
      # get dof numbers for each variable
      dofnum_rho = dofnums_j[1]
      dofnum_rhou = dofnums_j[2]
      dofnum_rhov = dofnums_j[3]
      dofnum_e = dofnums_j[4]

      # coordinates of this node (must be a vertex)
      x = coords[1]
      y = coords[2]

      # apply initial conditions here
      u0[dofnum_rho] = L*(2*k*e^(-2*k*x))/(e^(-2*k*x) +1 )^2
      u0[dofnum_rhou] = 0.0
      u0[dofnum_rhov] = 0.0
      u0[dofnum_e] = 0.0
  end
end

return nothing



end

function ICsmoothHeaviside{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, operator::SBPOperator{Tsbp}, eqn::EulerData{Tsol}, opts, u0::AbstractArray{Tsol, 1})
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
#  dofnums_i = view(mesh, i)  # get dof nums for this element
#  coords = view(mesh, [i])

  for j=1:nnodes
      coords = view(mesh.coords, :, j, i)
      dofnums_j = view(mesh.dofs, :, j, i)  # get dof nums for this element
 
      # get dof numbers for each variable
      dofnum_rho = dofnums_j[1]
      dofnum_rhou = dofnums_j[2]
      dofnum_rhov = dofnums_j[3]
      dofnum_e = dofnums_j[4]

      # coordinates of this node (must be a vertex)
      x = coords[1]
      y = coords[2]

      # apply initial conditions here
      u0[dofnum_rho] = L/(1 + e^(-k*(x-x0)))
      u0[dofnum_rhou] = 0.0
      u0[dofnum_rhov] = 0.0
      u0[dofnum_e] = 0.0
  end
end

return nothing



end

function ICIsentropicVortex{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, operator::SBPOperator{Tsbp}, eqn::EulerData{Tsol}, opts, u0::AbstractArray{Tsol})
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
#      coords = view(mesh.coords, :, j, i)
      dofnums_j = view(mesh.dofs, :, j, i)  # get dof nums for this element
 

      # coordinates of this node (must be a vertex)
#      coords_j = coords[:,j]
      coords_j = view(mesh.coords, :, j, i)
      calcIsentropicVortex(coords_j, eqn.params, sol)

#      println( "  j = ", j, " sol = ", sol, " coords_j = ", coords_j)

      # apply initial conditions here
      u0[dofnums_j] = sol
  end
end

return nothing

end  # end function

function ICIsentropicVortexWithNoise{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, operator::SBPOperator{Tsbp}, eqn::EulerData{Tsol}, opts, u0::AbstractVector{Tsol})
# populate u0 with initial values
# this is a template for all other initial conditions

numEl = getNumEl(mesh)
nnodes = operator.numnodes
dofpernode = getNumDofPerNode(mesh)
sol = zeros(Tsol, 4)
for i=1:numEl
#  dofnums_i = view(mesh, i)  # get dof nums for this element
#  coords = view(mesh, [i])

  for j=1:nnodes
      coords = view(mesh.coords, :, j, i)
      dofnums_j = view(mesh.dofs, :, j, i)  # get dof nums for this element
 
      # coordinates of this node (must be a vertex)
      calcIsentropicVortex(coords, eqn.params, sol)

      # apply initial conditions here
#       u0[dofnums_j] = sol
      u0[dofnums_j] = sol+0.1*rand(4)
  end
end

return nothing

end  # end function


function ICUnsteadyVortex{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, operator::SBPOperator{Tsbp}, eqn::EulerData{Tsol}, opts, u0::AbstractArray{Tsol})
# populate u0 with initial values
# this is a template for all other initial conditions

println("entered ICIsentropicVortex")

numEl = getNumEl(mesh)
nnodes = operator.numnodes
dofpernode = getNumDofPerNode(mesh)
sol = zeros(Tsol, 4)
for i=1:numEl
#  println("i = ", i)
#  dofnums_i = view(mesh, i)  # get dof nums for this element
#  coords = view(mesh, [i])

  for j=1:nnodes
#      coords = view(mesh.coords, :, j, i)
      dofnums_j = view(mesh.dofs, :, j, i)  # get dof nums for this element
 
      # coordinates of this node (must be a vertex)
#      coords_j = coords[:,j]
      coords = view(mesh.coords, :,j, i)
      calcUnsteadyVortex(coords, eqn.params, sol)

#      println( "  j = ", j, " sol = ", sol, " coords_j = ", coords_j)

      # apply initial conditions here
      u0[dofnums_j] = sol
  end
end

return nothing

end  # end function




function ICFile{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, operator::SBPOperator{Tsbp}, eqn::EulerData{Tsol}, opts, u0::AbstractVector{Tsol})
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

global const ICDict = Dict{Any, Function} (
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


