# Linear operator for SparseMatrixCSC

"""
  LinearOperator type to be used with sparse direct solve.  Note that the
  sparsity pattern of the matrix A must be constant.

  **Fields**

   * A: the matrix
"""
mutable struct SparseDirectLO <: AbstractSparseDirectLO
  A::SparseMatrixCSC{Float64, Int64}
  fac::UmfpackLU{Float64, Int64}
  is_setup::Array{Bool, 1}
  is_shared::Bool
  is_finalized::Bool
  nfactorizations::Int
  nsolves::Int
  ntsolves::Int

  # MPI stuff
  comm::MPI.Comm
  myrank::Int
  commsize::Int
end

"""
  Outer constructor for SparseDirectLO

  **Inputs**

   * pc: a PCNone
   * mesh
   * sbp
   * eqn
   * opts
"""
function SparseDirectLO(pc::AbstractPCNone, mesh::AbstractMesh, sbp::AbstractOperator,
                        eqn::AbstractSolutionData, opts::Dict)

  if typeof(mesh) <: AbstractCGMesh
    jac = SparseMatrixCSC(mesh.sparsity_bnds, Float64)
  else
    face_type = getFaceType(mesh.sbpface)
    disc_type = INVISCID
    if opts["preallocate_jacobian_coloring"]
      disc_type = COLORING
    end
    jac = SparseMatrixCSC(mesh, Float64, disc_type, face_type)
  end

  # Note: colptr and rowval alias A
  #      also, compute the symbolic factorization

  fac = UmfpackLU{Float64, Int}(C_NULL, C_NULL, mesh.numDof, mesh.numDof,
                                jac.colptr, jac.rowval, jac.nzval)

  make_zerobased(jac)
  umfpack_symbolic!(fac)
  make_onebased(jac)

  is_setup = Bool[false]
  is_shared = false
  nfactorizations = 0
  nsolves = 0
  ntsolves = 0

  comm = eqn.comm
  myrank = eqn.myrank
  commsize = eqn.commsize

  return SparseDirectLO(jac, fac, is_setup, is_shared, false, nfactorizations,
                        nsolves, ntsolves, comm, myrank, commsize)
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
                            sbp::AbstractOperator, eqn::AbstractSolutionData,
                            opts::Dict, ctx_residual, t)

#  physicsJac(lo, mesh, sbp, eqn, opts, lo.A, ctx_residual, t)
  setIsSetup(lo, false)

  return nothing
end

function applyLinearOperator(lo::AbstractSparseDirectLO, mesh::AbstractMesh,
                             sbp::AbstractOperator, eqn::AbstractSolutionData,
                             opts::Dict, ctx_residual, t, x::AbstractVector, 
                             b::AbstractVector)

  lo2 = getBaseLO(lo)
  A_mul_B!(1, lo2.A, x, 0, b)

  return nothing
end


function applyLinearOperatorTranspose(lo::AbstractSparseDirectLO, 
                             mesh::AbstractMesh, sbp::AbstractOperator,
                             eqn::AbstractSolutionData, opts::Dict, 
                             ctx_residual, t, x::AbstractVector, 
                             b::AbstractVector)
  lo2 = getBaseLO(lo)
  At_mul_B!(1, lo2.A, x, 0, b)

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


