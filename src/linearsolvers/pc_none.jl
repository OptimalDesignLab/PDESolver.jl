# define the special PCNone type for use with direct solvers

"""
  Preconditioner type for direct solve.
  Do not use with PetscLinearOperator.

  **Public Fields**

   * none
"""
mutable struct PCNone <: AbstractPC
  is_shared::Bool
end

"""
  Outer constructor

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
"""
function PCNone(mesh::AbstractMesh, sbp::AbstractOperator,
                eqn::AbstractSolutionData, opts::Dict)

  is_shared = false
  return PCNone(is_shared)
end

function free(pc::PCNone)

  return nothing
end

function calcPC(pc::PCNone, mesh::AbstractMesh, sbp::AbstractOperator,
                eqn::AbstractSolutionData, opts::Dict, ctx_residual, t)

  return nothing
end


function applyPC(pc::PCNone, mesh::AbstractMesh, sbp::AbstractOperator,
                 eqn::AbstractSolutionData, opts::Dict, t, b::AbstractVector, 
                 x::AbstractVector)


  return nothing
end


function applyPCTranspose(pc::PCNone, mesh::AbstractMesh, sbp::AbstractOperator,
                 eqn::AbstractSolutionData, opts::Dict, t, b::AbstractVector, 
                 x::AbstractVector)


  return nothing
end

function getBasePC(pc::PCNone)

  # this is the bottom of the recursion tree
  return pc
end
