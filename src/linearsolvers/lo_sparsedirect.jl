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

#TODO: in the constructor, make colptr and rowval alias A
#      also, compute the symbolic factorization
#      also, create finalizer for symbolic factorizatrion

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


function applyLinearOperatorTranspose(lo::SparseDirectLO, 
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

#------------------------------------------------------------------------------
# helper functions

"""
  Take a 1-based SparseMatrixCSC and make it 0-based.  Throws an exception
  if the matrix is not 1-based.
"""
function make_zerobased(A::SparseMatrixCSC)

  @assert A.colptr[1] == 1

  rowval = A.rowval
  colptr = A.colptr
  @simd @inbound for i=1:length(rowval)
    rowval[i] -= 1
  end

  @simd @inbounds for i=1:length(colptr)
    colptr[i] -= 1
  end

  return nothing
end


"""
  Take a 0-based SparseMatrixCSC and make it 1-based.  Throws an exception if
  the matrix is not 1-based
"""
function make_onebased(A::SparseMatrixCSC)

  @assert A.colptr[1] == 0

  rowval = A.rowval
  colptr = A.colptr
  @simd @inbound for i=1:length(rowval)
    rowval[i] += 1
  end

  @simd @inbounds for i=1:length(colptr)
    colptr[i] += 1
  end

  return nothing
end


