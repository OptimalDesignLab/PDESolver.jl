# functions that populate the initial conditions
# List of functions:

@doc """
### EllipticEquationMod.ICZero

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
function ICZero(mesh::AbstractMesh{Tmsh}, 
                operator::AbstractSBP{Tsbp}, 
                eqn::EllipticData{Tsol}, 
                opts, 
                u0::AbstractVector{Tsol}) where {Tmsh, Tsbp, Tsol}
  # populate u0 with initial values
  # this is a template for all other initial conditions

  numEl = mesh.numEl
  nnodes = operator.numnodes
  dofpernode = mesh.numDofPerNode
  for i=1:numEl

    for j=1:nnodes
      coords_j = sview(mesh.coords, :, j, i)
      dofnums_j = sview(mesh.dofs, :, j, i)
      for dof = 1: dofpernode
        dofnum = dofnums_j[dof]
        u0[dofnum] = 0.0
      end
    end
  end

  return nothing

end  # end function


@doc """
### EllipticEquationMod.ICOnes

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

function ICOnes(mesh::AbstractMesh{Tmsh}, 
                operator::AbstractSBP{Tsbp}, 
                eqn::EllipticData{Tsol}, 
                opts,
                u0::AbstractVector{Tsol}) where {Tmsh, Tsbp, Tsol}

  numEl = mesh.numEl
  nnodes = operator.numnodes
  dofpernode = mesh.numDofPerNode
  for i=1:numEl
    for j=1:nnodes
      coords_j = sview(mesh.coords, :, j, i)
      dofnums_j = sview(mesh.dofs, :, j, i)

      x = coords_j[1]
      y = coords_j[2]

      # apply initial conditions here
      u0[dofnums_j[:]] = 1.0
    end
  end

  return nothing
end # end function ICOnes

function ICRandom(mesh::AbstractMesh{Tmsh}, 
                  operator::AbstractSBP{Tsbp}, 
                  eqn::EllipticData{Tsol}, 
                  opts,
                  u0::AbstractVector{Tsol}) where {Tmsh, Tsbp, Tsol}

  numEl = mesh.numEl
  nnodes = operator.numnodes
  dofpernode = mesh.numDofPerNode
  for i=1:numEl
    for j=1:nnodes
      coords_j = sview(mesh.coords, :, j, i)
      dofnums_j = sview(mesh.dofs, :, j, i)

      x = coords_j[1]
      y = coords_j[2]

      r = rand(1:100)
      # apply initial conditions here
      u0[dofnums_j[:]] = r * 1.0e-1
    end
  end

  return nothing
end # end function ICOnes



# declare a const dictionary here that maps strings to function (used for input arguments)

global const ICDict = Dict{Any, Function}(
  "ICZero" => ICZero,
  "ICOnes" => ICOnes,
  "ICRandom" => ICRandom,
)


