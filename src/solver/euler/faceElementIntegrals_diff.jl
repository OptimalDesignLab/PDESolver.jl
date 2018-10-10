# differentiated version of faceElementIntegrals.jl



#------------------------------------------------------------------------------
# AbstractEntropyKernel 

"""
  Differentiated version of the LFKernel method.
"""
function applyEntropyKernel_diff(obj::LFKernel, params::ParamType, 
                          q_avg::AbstractVector{Tsol}, q_avg_dot::AbstractMatrix,
                          delta_w::AbstractVector, delta_w_dot::AbstractMatrix,
                          nrm::AbstractVector,
                          flux::AbstractVector{Tres}, flux_dot::AbstractMatrix)  where {Tsol, Tres}

  nd = size(q_avg_dot, 2)
  numDofPerNode = size(q_avg, 1)

  @assert nd == size(delta_w_dot, 2)
  @assert nd == size(flux_dot, 2)

  @unpack obj A0 A0_dot t1 t1_dot lambda_max_dot
  numDofPerNode = length(flux)

  #TODO: this is the dominant cost.  Perhaps compute 5 matrices: dA0/dq_i 
  # for i=1:5, then multiply each by q_avg_dot in each direction (as part of
  # triple for loop below)
  getIRA0_diff(params, q_avg, q_avg_dot, A0, A0_dot)

  # t1 = A0*delta_w
  smallmatvec!(A0, delta_w, t1)

  # t1_dot[:, i] = A0_dot[:, :, i]*delta_w + A0*delta_w_dot[:, i] for i=1:nd
  smallmatmat!(A0, delta_w_dot, sview(t1_dot, :, 1:nd))

  # this is basically repeated matrix-vector multiplication
  @simd for d=1:nd
    @simd for i=1:numDofPerNode
      @simd for j=1:numDofPerNode
        t1_dot[j, d] += A0_dot[j, i, d]*delta_w[i]
      end
    end
  end

  # now multiply with lambda_max
  lambda_max = getLambdaMax_diff(params, q_avg, nrm, lambda_max_dot)
  for i=1:numDofPerNode
    flux[i] = lambda_max*t1[i]
  end

  for d=1:nd
    # compute lambda_max_dot in direction d
    lambda_max_dotd = zero(Tres)
    for i=1:numDofPerNode
      lambda_max_dotd += lambda_max_dot[i]*q_avg_dot[i, d]
    end

    # compute flux_dot
    for i=1:numDofPerNode
      flux_dot[i, d] += lambda_max*t1_dot[i, d] + lambda_max_dotd*t1[i]
    end
  end

  return nothing
end
