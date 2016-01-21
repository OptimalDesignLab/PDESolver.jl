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
                    sbp::SBPOperator{Tsbp}, eqn::EulerData{Tsol}, 
                    opts, u0::AbstractArray{Tsol})

  for i = 1:mesh.numEl
  	for j = 1:mesh.numNodenumNodesPerElement
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

function exp_ICxplusy{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, 
                    sbp::SBPOperator{Tsbp}, eqn::EulerData{Tsol}, 
                    opts, u0::AbstractArray{Tsol})

  for i = 1:mesh.numEl
  	for j = 1:mesh.numNodenumNodesPerElement
  	  dofnums_j = view(mesh.dofs, :, j, i)
  	  x = mesh.coords[1,j,i]
  	  y = mesh.coords[2,j,i]
  	  u0[dofnums_j] = exp(x + y)
  	end
  end

  return nothing
end # end function exp_xplusy