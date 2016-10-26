# ic.jl
# Needed to initialize a problem.

function ICConstant(mesh::AbstractMesh, sbp::AbstractSBP, eqn::SimpleODEData, opts,
                    u0::AbstractArray)

  for i=1:length(u0)
    u0[i] = 2
  end

  return nothing
end

@doc """
### SimpleODEMod.ICx2plust2

Computes the initial conditions for the state variable
    u = x^2 + t^2

**Inputs**

*  `mesh` : AbstractMesh type
*  `sbp`  : Summation-by-parts operator
*  `eqn`  : Simple ODE equation object
*  `opts` : Options dictionary
*  `u0`   : Array that stores inital state variables

**Outputs**

*  None

"""->

function ICx2plust2{Tmsh, Tsbp, Tsol, Tres}(mesh::AbstractMesh{Tmsh}, 
                    sbp::AbstractSBP{Tsbp}, eqn::SimpleODEData{Tsol, Tres, 2}, 
                    opts, u0::AbstractArray{Tsol})

  t = 0.0

  for i = 1:mesh.numEl
  	for j = 1:mesh.numNodesPerElement
  	  dofnums_j = sview(mesh.dofs, :, j, i)
#   	  x = mesh.coords[1, j, i]
#   	  y = mesh.coords[2, j, i]
#   	  u0[dofnums_j] = x^2 + t^2
  	  u0[dofnums_j] = calc_x2_t2(mesh.coords[:, j, i], eqn.params, t)
  	end
  end

  return nothing
end # end function ICx5plusy5

global const ICDict = Dict{Any, Function}(
"ICConstant" => ICConstant,
"ICx2plust2" => ICx2plust2,
)
