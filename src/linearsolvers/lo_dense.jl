# Dense matrix linear operator

"""
  Dense array linear operator.  Serial only.
"""
type DenseLO <: AbstractDenseLO
  A::Array{Float64, 2}
  ipiv::Array{BlasInt, 1}
  is_setup::Bool # true if LO is already setup (ie. factored), false otherwise
end

#TODO: outer constructor

function calcLinearOperator(lo::DenseLO, mesh::AbstractMesh,
                            sbp::AbstractSBP, eqn::AbstractSolutionData,
                            opts::Dict, ctx_residual, t)

  physicsJac(lo, mesh, sbp, eqn, opts, lo.A, ctx_residual, t)
  lo.is_setup = false

  return nothing
end


# note: this may be faster if lo is not yet set up
function applyLinearOperator(lo::DenseLO, mesh::AbstractMesh,
                             sbp::AbstractSBP, eqn::AbstractSolutionData,
                             opts::Dict, ctx_residual, t, x::AbstractVector, 
                             b::AbstractVector)

  # BLAS operates in-place
  if lo.is_setup  # A is already factored
    # A = P*L*U
    copy!(b, x)
    trmv!('U', 'N', 'N', lo.A, b)
    trmv!('L', 'N', 'U', lo.A, b)
    laswp!(b, 1, length(b), lo.ipiv)
  else
    gemv!('N', 1.0, lo.A, x, 0.0, b)
  end

  return nothing
end


function applyLinearOperatorTranspose(lo::DenseLO,
                             mesh::AbstractMesh, sbp::AbstractSBP,
                             eqn::AbstractSolutionData, opts::Dict, 
                             ctx_residual, t, x::AbstractVector, 
                             b::AbstractVector)

  if lo.is_setup  # A is already factored
    # A.'x = (P*L*U).'x = (U.'*L.'*P.')*x
    # there is no way to apply P.' in lapack (ipiv is not a true permutation)

    throw(ErrorException("applyLinearOperatorTranspose not supported for DenseLO in factored state: recompute the linear operator"))
  else
    gemv!('T', 1.0, lo.A, x, 0.0, b)
  end

  return nothing
end


function getBaseLinearOperator(lo::DenseLO)

  # this is the bottom of the recursion tree
  return lo
end
