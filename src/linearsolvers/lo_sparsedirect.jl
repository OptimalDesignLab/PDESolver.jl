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
  is_finalized::Bool
  nfactorizations::Int
  nsolves::Int
  ntsolves::Int
end

#TODO: in the constructor, make colptr and rowval alias A
#      also, compute the symbolic factorization
#      also, create finalizer for symbolic factorizatrion

function SparseDirectLO(pc::PCNone, mesh::AbstractMesh, sbp::AbstractSBP,
                        eqn::AbstractSolutionData, opts::Dict)

  if typeof(mesh) <: AbstractCGMesh
    jac = SparseMatrixCSC(mesh.sparsity_bnds, Float64)
  else
    jac = SparseMatrixCSC(mesh, Float64)
  end

  fac = UmfpackLU{Float64, Int}(C_NULL, C_NULL, mesh.numDof, mesh.numDof,
                                jac.colptr, jac.rowval, jac.nzval)

  make_zerobased(jac)
  umfpack_symbolic!(fac)
  make_onebased(jac)

  is_setup = false
  nfactorizations = 0
  nsolves = 0
  ntsolves = 0

  return SparseDirectLO(jac, fac, is_setup, false, nfactorizations, nsolves, ntsolves)
end


function free(lo::SparseDirectLO)

  if !lo.is_finalized
    umfpack_free_symbolic(lo.fac)
    umfpack_free_numeric(lo.fac)
  end

  lo.is_finalized = true

  return nothing
end


function calcLinearOperator(lo::SparseDirectLO, mesh::AbstractMesh,
                            sbp::AbstractSBP, eqn::AbstractSolutionData,
                            opts::Dict, ctx_residual, t)

#  physicsJac(lo, mesh, sbp, eqn, opts, lo.A, ctx_residual, t)
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

  At_mul_B!(1, lo.A, x, 0, b)

  return nothing
end

function getBaseLO(lo::SparseDirectLO)

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
  @inbounds @simd for i=1:length(rowval)
    rowval[i] -= 1
  end

  @inbounds @simd for i=1:length(colptr)
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
  @inbounds @simd for i=1:length(rowval)
    rowval[i] += 1
  end

  @inbounds @simd for i=1:length(colptr)
    colptr[i] += 1
  end

  return nothing
end


