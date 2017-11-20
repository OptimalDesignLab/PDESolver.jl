# Dense matrix linear operator

"""
  Dense array linear operator.  Serial only.
"""
type DenseLO <: AbstractDenseLO
  A::Array{Float64, 2}
  ipiv::Array{BlasInt, 1}
  is_setup::Array{Bool, 1} # true if LO is already setup (ie. factored), false otherwise
  nfactorizations::Int
  nsolves::Int  # regular solves
  ntsolves::Int  # transpose solves

  # MPI stuff
  comm::MPI.Comm
  myrank::Int
  commsize::Int
end

"""
  Outer constructor for DenseLO

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
"""
function DenseLO(pc::PCNone, mesh::AbstractMesh, sbp::AbstractSBP,
                 eqn::AbstractSolutionData, opts::Dict)

  @assert mesh.commsize == 1  # serial only
  A = zeros(mesh.numDof, mesh.numDof)
  ipiv = zeros(BlasInt, mesh.numDof)
  is_setup = Bool[false]
  nfactorizations = 0
  nsolves = 0
  ntsolves = 0
  comm = eqn.comm
  myrank = eqn.myrank
  commsize = eqn.commsize

  return DenseLO(A, ipiv, is_setup, nfactorizations, nsolves, ntsolves,
                 comm, myrank, commsize)
end


function free(lo::DenseLO)
  # nothing to do here

  return nothing
end


function calcLinearOperator(lo::DenseLO, mesh::AbstractMesh,
                            sbp::AbstractSBP, eqn::AbstractSolutionData,
                            opts::Dict, ctx_residual, t)

#  physicsJac(lo, mesh, sbp, eqn, opts, lo.A, ctx_residual, t)
  setIsSetup(lo, false)

  return nothing
end


# note: this may be faster if lo is not yet set up
function applyLinearOperator(lo::AbstractDenseLO, mesh::AbstractMesh,
                             sbp::AbstractSBP, eqn::AbstractSolutionData,
                             opts::Dict, ctx_residual, t, x::AbstractVector, 
                             b::AbstractVector)

  lo2 = getBaseLO(lo)
  # BLAS operates in-place
  if getIsSetup(lo2)  # A is already factored
    throw(ErrorException("Multiplying factored matrix not supported"))
    # this doesn't work because the permutation is not correct for the vector
    # I think the problem is that ipiv is not a permuatation matrix, so really
    # what should happen is laswp(L)*U*b rather than laswp(L*U*b)
    #=
    println("multiplying in factored form")
    # A = P*L*U
    copy!(b, x)
    trmv!('U', 'N', 'N', lo2.A, b)
    trmv!('L', 'N', 'U', lo2.A, b)
    println("b = \n", b)
    println("ipiv = \n", lo2.ipiv)
    laswp!(b, 1, length(b), lo2.ipiv)
    println("b = \n", b)
    =#
  else
    println("standard multiply")
    gemv!('N', 1.0, lo2.A, x, 0.0, b)
  end

  return nothing
end


function applyLinearOperatorTranspose(lo::AbstractDenseLO,
                             mesh::AbstractMesh, sbp::AbstractSBP,
                             eqn::AbstractSolutionData, opts::Dict, 
                             ctx_residual, t, x::AbstractVector, 
                             b::AbstractVector)
  lo2 = getBaseLO(lo)
  if getIsSetup(lo2)  # A is already factored
    # A.'x = (P*L*U).'x = (U.'*L.'*P.')*x
    # there is no way to apply P.' in lapack (ipiv is not a true permutation)

    throw(ErrorException("applyLinearOperatorTranspose not supported for DenseLO in factored state: recompute the linear operator"))
  else
    gemv!('T', 1.0, lo2.A, x, 0.0, b)
  end

  return nothing
end


function getBaseLO(lo::DenseLO)

  # this is the bottom of the recursion tree
  return lo
end
