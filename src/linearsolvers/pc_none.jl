# define the special PCNONE type for use with direct solvers

"""
  Preconditioner type for direct solve.
  Do not use with PetscLinearOperator.
"""
type PCNone <: AbstractPC
end

function PCNone(pc::PCNone, mesh::AbstractMesh, sbp::AbstractSBP,
                eqn::AbstractSolutionData, opts::Dict)

  return PCNone()
end

function calcPC(pc::PCNone, mesh::AbstractMesh, sbp::AbstractSBP,
                eqn::AbstractSolutionData, opts::Dict, ctx_residual, t)

  return nothing
end


function applyPC(pc::PCNone, mesh::AbstractMesh, sbp::AbstractSBP,
                 eqn::AbstractSolutionData, opts::Dict, t, b::AbstractVector, 
                 x::AbstractVector)


  return nothing
end


function applyPCTranspose(pc::PCNone, mesh::AbstractMesh, sbp::AbstractSBP,
                 eqn::AbstractSolutionData, opts::Dict, t, b::AbstractVector, 
                 x::AbstractVector)


  return nothing
end

function getBasePC(pc::PCNone)

  # this is the bottom of the recursion tree
  return pc
end
