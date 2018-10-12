# differentiated version of faceElementIntegrals.jl

"""
  Computes the Jacobian of [`calcECFaceIntegral`](@ref).  See that function
  for the meaning of the input arguments

  Note that the output arguments are updated, not overwritten.  The caller
  must zero them out when required!

  **Inputs**

   * params: ParamType
   * sbpface: SBPFace ojbect
   * iface: interface
   * qL: solution at the volume nodes of the left element
   * qR: solution at the volume nodes of the right element
   * aux_vars
   * nrm_xy: normal vector at each face node

  **Inputs/Outputs**

   * jacLL: jacobian of residual of left element wrt solution on left element
   * jacLR: jacobian of residual of left element wrt solution on right element
   * jacRL: jacobian of residual of right element wrt solution on left element
   * jacRR: jacobian of residual of right element wrt solution on left element

"""
function calcECFaceIntegral_diff(
     params::AbstractParamType{Tdim}, 
     sbpface::DenseFace, 
     iface::Interface,
     qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, 
     aux_vars::AbstractMatrix{Tres}, 
     nrm_xy::AbstractMatrix{Tmsh},
     functor::FluxType_diff, 
     jacLL::AbstractArray{Tres, 4}, jacLR::AbstractArray{Tres, 4},
     jacRL::AbstractArray{Tres, 4}, jacRR::AbstractArray{Tres, 4}
     ) where {Tdim, Tsol, Tres, Tmsh}


  data = params.calc_ec_face_integral_data
  @unpack data fluxD nrmD
  numDofPerNode = size(qL, 1)

  fL_dot = zeros(Tres, numDofPerNode, numDofPerNode, Tdim)
  fR_dot = zeros(Tres, numDofPerNode, numDofPerNode, Tdim)

  fill!(nrmD, 0.0)
  for d=1:Tdim
    nrmD[d, d] = 1
  end

  # loop over the nodes of "left" element that are in the stencil of interp
  for i = 1:sbpface.stencilsize
    p_i = sbpface.perm[i, iface.faceL]
    qi = ro_sview(qL, :, p_i)
    aux_vars_i = ro_sview(aux_vars, :, p_i)  # !!!! why no aux_vars_j???

    # loop over the nodes of "right" element that are in the stencil of interp
    for j = 1:sbpface.stencilsize
      p_j = sbpface.perm[j, iface.faceR]
      qj = ro_sview(qR, :, p_j)

      # compute flux and add contribution to left and right elements
      fill!(fL_dot, 0.0)
      fill!(fR_dot, 0.0)
      functor(params, qi, qj, aux_vars_i, nrmD, fL_dot, fR_dot)

      @simd for dim = 1:Tdim
        # accumulate entry p_i, p_j of E
        Eij = zero(Tres)
        @simd for k = 1:sbpface.numnodes
          # the computation of nrm_k could be moved outside i,j loops and saved
          # in an array of size [3, sbp.numnodes]
          nrm_k = nrm_xy[dim, k]
          kR = sbpface.nbrperm[k, iface.orient]
          Eij += sbpface.interp[i,k]*sbpface.interp[j,kR]*sbpface.wface[k]*nrm_k
        end  # end loop k
 

        # update jac
        @simd for q=1:numDofPerNode
          @simd for p=1:numDofPerNode
            valL = Eij*fL_dot[p, q, dim]
            valR = Eij*fR_dot[p, q, dim]

            jacLL[p, q, p_i, p_i] -= valL
            jacLR[p, q, p_i, p_j] -= valR
            jacRL[p, q, p_j, p_i] += valL
            jacRR[p, q, p_j, p_j] += valR
          end
        end

      end  # end loop dim
    end  # end loop j
  end  # end loop i


  return nothing
end



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
