# ic.jl
# Needed to initialize a problem.

@doc """
### AdvectionEquationMod.ICx5plusy5

Computes the initial conditions for the state variable
    u = x^5 + y^5

**Inputs**

*  `mesh` : AbstractMesh type
*  `sbp`  : Summation-by-parts operator
*  `eqn`  : Advection equation object
*  `opts` : Options dictionary
*  `u0`   : Array that stores inital state variables

**Outputs**

*  None

"""->

function ICx5plusy5{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                    sbp::SBPOperator{Tsbp}, eqn::AdvectionData{Tsol}, 
                    opts, u0::AbstractArray{Tsol})

  for i = 1:mesh.numEl
  	for j = 1:mesh.numNodesPerElement
  	  dofnums_j = view(mesh.dofs, :, j, i)
  	  x = mesh.coords[1,j,i]
  	  y = mesh.coords[2,j,i]
  	  u0[dofnums_j] = x^5 + y^5
  	end
  end

  return nothing
end # end function ICx5plusy5

@doc """
### AdvectionEquationMod.ICexp_xplusy

Computes the initial conditions for the state variable
    u = exp(x + y)

**Inputs**

*  `mesh` : AbstractMesh type
*  `sbp`  : Summation-by-parts operator
*  `eqn`  : Advection equation object
*  `opts` : Options dictionary
*  `u0`   : Array that stores inital state variables

**Outputs**

*  None

"""->

function ICexp_xplusy{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                    sbp::SBPOperator{Tsbp}, eqn::AdvectionData{Tsol}, 
                    opts, u0::AbstractArray{Tsol})

  for i = 1:mesh.numEl
  	for j = 1:mesh.numNodesPerElement
  	  dofnums_j = view(mesh.dofs, :, j, i)
  	  x = mesh.coords[1,j,i]
  	  y = mesh.coords[2,j,i]
  	  u0[dofnums_j] = exp(x + y)
  	end
  end

  return nothing
end # end function exp_xplusy


function ICsinwave{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                    sbp::SBPOperator{Tsbp}, eqn::AdvectionData{Tsol}, 
                    opts, u0::AbstractArray{Tsol})

  for i = 1:mesh.numEl
  	for j = 1:mesh.numNodesPerElement
  	  dofnums_j = view(mesh.dofs, :, j, i)
  	  x = mesh.coords[1,j,i]
  	  y = mesh.coords[2,j,i]
  	  u0[dofnums_j] = calc_sinwave(mesh.coords[:, j, i], eqn.t)
  	end
  end

  return nothing
end # end function exp_xplusy




@doc """
### AdvectionEquationMod.ICFile

This function reads a vector from a file on disk and set the solution to it.
The vector must contain the same number of entries as there are degrees of 
freedom in the mesh. 

This function is useful for things like restarting from a checkpoint.
In this case, the file should be the output of writedlm(eqn.q).  The degree 
of freedom number must be the same for both simulation for this to work (the 
file contains no degree of freedom number information).


**Inputs**

*  `mesh`
*  `sbp`
*  `eqn`
*  `opts`

**Inputs/Outputs**

*  `u0`: vector to populate with the solution

  Aliasing restrictions: none.

"""->
function ICFile{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                operator::SBPOperator{Tsbp}, eqn::AdvectionData{Tsol}, opts, 
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

global const ICDict = Dict{Any, Function} (
"ICx5plusy5" => ICx5plusy5,
"ICexp_xplusy" => ICexp_xplusy,
"ICsinwave" => ICsinwave,
"ICFile" => ICFile
)
