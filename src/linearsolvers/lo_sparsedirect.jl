# Linear operator for SparseMatrixCSC

"""
  LinearOperator type to be used with sparse direct solve.  Note that the
  sparsity pattern of the matrix A must be constant.

  **Fields**

   * A: the matrix
"""
type SparseDirectLO <: AbstractSparseDirectLO
  A::SparseMatrixCSC{Float64, Int64}
  fac::UmfpackLU{Float64, Int64}
  is_setup::Bool
end

function calcLinearOperator(lo::SparseDirectLO, mesh::AbstractMesh,
                            sbp::AbstractSBP, eqn::AbstractSolutionData,
                            opts::Dict, ctx_residual, t)

  physicsJac(lo, mesh, sbp, eqn, opts, lo.A, ctx_residual, t)
  lo.is_setup = false

  return nothing
end

function applyLinearOperator(lo::SparseDirectLO, mesh::AbstractMesh,
                             sbp::AbstractSBP, eqn::AbstractSolutionData,
                             opts::Dict, ctx_residual, t, x::AbstractVector, 
                             b::AbstractVector)

  A_mul_B!(1, lo.A, x, 0, b)

  return nothing
end


function applyLinearOperatorTranspose(lo::AbstractLinearOperator, 
                             mesh::AbstractMesh, sbp::AbstractSBP,
                             eqn::AbstractSolutionData, opts::Dict, 
                             ctx_residual, t, x::AbstractVector, 
                             b::AbstractVector)

  A_mul_Bt!(1, lo.A, x, 0, b)

  return nothing
end

function getBaseLinearOperator(lo::SparseDirectLO)

  # this is the bottom of the recursion tree
  return lo
end
