# differentiated version of faceElementIntegrals.jl

#------------------------------------------------------------------------------
# calcECFaceIntegral

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
  @unpack data fluxD nrmD fL_dotD fR_dotD
  numDofPerNode = size(qL, 1)

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
      fill!(fL_dotD, 0.0)
      fill!(fR_dotD, 0.0)
      functor(params, qi, qj, aux_vars_i, nrmD, fL_dotD, fR_dotD)

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
            valL = Eij*fL_dotD[p, q, dim]
            valR = Eij*fR_dotD[p, q, dim]

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

# SparseFace method
function calcECFaceIntegral_diff(
     params::AbstractParamType{Tdim}, 
     sbpface::SparseFace, 
     iface::Interface,
     qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, 
     aux_vars::AbstractMatrix{Tres}, 
     nrm_xy::AbstractMatrix{Tmsh},
     functor::FluxType_diff, 
     jacLL::AbstractArray{Tres, 4}, jacLR::AbstractArray{Tres, 4},
     jacRL::AbstractArray{Tres, 4}, jacRR::AbstractArray{Tres, 4}
     ) where {Tdim, Tsol, Tres, Tmsh}

  data = params.calc_ec_face_integral_data
  @unpack data fL_dot fR_dot

  numDofPerNode = size(qL, 1)

  for i=1:sbpface.numnodes
    p_i = sbpface.perm[i, iface.faceL]
    q_i = ro_sview(qL, :, p_i)
    aux_vars_i = ro_sview(aux_vars, :, p_i)

    # get the corresponding node on faceR
    pnbr = sbpface.nbrperm[i, iface.orient]
    p_j = sbpface.perm[pnbr, iface.faceR]
#    p_j = sbpface.nbrperm[sbpface.perm[i, iface.faceR], iface.orient]
    q_j = ro_sview(qR, :, p_j)

    # compute flux in face normal direction
    nrm_i = ro_sview(nrm_xy, :, i)
    fill!(fL_dot, 0.0); fill!(fR_dot, 0.0)
    functor(params, q_i, q_j, aux_vars_i, nrm_i, fL_dot, fR_dot)

    w_i = sbpface.wface[i]
    @simd for q=1:numDofPerNode
      @simd for p=1:numDofPerNode
        valL = w_i*fL_dot[p, q]
        valR = w_i*fR_dot[p, q]

        jacLL[p, q, p_i, p_i] -= valL
        jacLR[p, q, p_i, p_j] -= valR
        jacRL[p, q, p_j, p_i] += valL
        jacRR[p, q, p_j, p_j] += valR
      end
    end


  end  # end loop i

  return nothing
end

#------------------------------------------------------------------------------
# calcEntropyPenaltyIntegral

function calcEntropyPenaltyIntegral_diff(
             params::ParamType{Tdim, :conservative, Tsol, Tres, Tmsh},
             sbpface::DenseFace, iface::Interface,
             kernel::AbstractEntropyKernel,
             qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, 
             aux_vars::AbstractMatrix{Tres}, nrm_face::AbstractArray{Tmsh, 2},
             jacLL::AbstractArray{Tres, 4}, jacLR::AbstractArray{Tres, 4},
             jacRL::AbstractArray{Tres, 4}, jacRR::AbstractArray{Tres, 4}
             ) where {Tdim, Tsol, Tres, Tmsh}

  numDofPerNode = size(qL, 1)
  nd = 2*numDofPerNode

  # convert qL and qR to entropy variables (only the nodes that will be used)
  data = params.calc_entropy_penalty_integral_data
  @unpack data wL wR wL_i wR_i qL_i qR_i delta_w q_avg flux

  # convert to IR entropy variables
  for i=1:sbpface.stencilsize
    # apply sbpface.perm here
    p_iL = sbpface.perm[i, iface.faceL]
    p_iR = sbpface.perm[i, iface.faceR]
    # these need to have different names from qL_i etc. below to avoid type
    # instability
    qL_itmp = ro_sview(qL, :, p_iL)
    qR_itmp = ro_sview(qR, :, p_iR)
    wL_itmp = sview(wL, :, i)
    wR_itmp = sview(wR, :, i)
    convertToIR(params, qL_itmp, wL_itmp)
    convertToIR(params, qR_itmp, wR_itmp)
  end


  @unpack data q_avg_dot delta_w_dot flux_dot_i flux_dotL flux_dotR
  @unpack data jacLL_tmp jacLR_tmp jacRL_tmp jacRR_tmp A0invL A0invR
  # flux_doti is passed into applyEntropyKernel_diff, flux_dotL, flux_dotR
  # are then used to reformat the data the way SBP wants it
#  q_avg_dot = zeros(Tsol, numDofPerNode, nd)
#  delta_w_dot = zeros(Tsol, numDofPerNode, nd)
#  flux_dot_i = zeros(Tres, numDofPerNode, nd)
  flux_dot_iL = sview(flux_dot_i, :, 1:numDofPerNode)  # wL part
  flux_dot_iR = sview(flux_dot_i, :, (numDofPerNode + 1):nd)  # wR part

  fill!(delta_w_dot, 0.0)
  fill!(jacLL_tmp, 0.0); fill!(jacLR_tmp, 0.0);
  fill!(jacRL_tmp, 0.0); fill!(jacRR_tmp, 0.0);

  # arrays needed by SBP function
#  flux_dotL = zeros(Tres, numDofPerNode, numDofPerNode, sbpface.numnodes)
#  flux_dotR = zeros(Tres, numDofPerNode, numDofPerNode, sbpface.numnodes)

 # jacLL_tmp = zeros(jacLL)
 # jacLR_tmp = zeros(jacLR)
 # jacRL_tmp = zeros(jacRL)
 # jacRR_tmp = zeros(jacRR)

 # A0invL = zeros(Tsol, numDofPerNode, numDofPerNode)
 # A0invR = zeros(Tsol, numDofPerNode, numDofPerNode)

 # flux_dot_iL = sview(flux_dot_i, :, 1:numDofPerNode)  # wL part
 # flux_dot_iR = sview(flux_dot_i, :, (numDofPerNode + 1):nd)  # wR part


  # accumulate wL at the node
  @simd for i=1:sbpface.numnodes  # loop over face nodes
    ni = sbpface.nbrperm[i, iface.orient]
    dir = ro_sview(nrm_face, :, i)
    fastzero!(wL_i)
    fastzero!(wR_i)

    # interpolate wL and wR to this node
    @simd for j=1:sbpface.stencilsize
      interpL = sbpface.interp[j, i]
      interpR = sbpface.interp[j, ni]

      @simd for k=1:numDofPerNode
        wL_i[k] += interpL*wL[k, j]
        wR_i[k] += interpR*wR[k, j]
      end
    end

    convertToConservativeFromIR_(params, wL_i, qL_i)
    convertToConservativeFromIR_(params, wR_i, qR_i)
    
    # compute average qL
    # also delta w (used later)
    # Approach for differentiating: compute d flux/d wL_i, d flux_d wR_i first
    # then compute flux jacobian at volume nodes
    @simd for j=1:numDofPerNode
      q_avg[j] = 0.5*(qL_i[j] + qR_i[j])
      delta_w[j] = wL_i[j] - wR_i[j]

      # set derivative part of wL_i, wR_i = I
      delta_w_dot[j, j] = 1; delta_w_dot[j, j + numDofPerNode] = -1
    end

    # get the derivative part of qL_i, qR_i, and q_avg all at the same time
    q_avg_dotL = sview(q_avg_dot, :, 1:numDofPerNode)
    q_avg_dotR = sview(q_avg_dot, :, (numDofPerNode + 1):nd)
    getIRA0(params, qL_i, q_avg_dotL)
    getIRA0(params, qR_i, q_avg_dotR)
    scale!(q_avg_dot, 0.5)



    # call kernel (apply symmetric semi-definite matrix)
    fill!(flux_dot_i, 0.0)
    applyEntropyKernel_diff(kernel, params, q_avg, q_avg_dot, delta_w, delta_w_dot, dir, flux, flux_dot_i)

    # put the data in the format required by SBP later
    # the copy could be avoided if we had a decent strided view implementation
    flux_dotL_i = sview(flux_dotL, :, :, i)
    flux_dotR_i = sview(flux_dotR, :, :, i)
    copy!(flux_dotL_i, flux_dot_iL); copy!(flux_dotR_i, flux_dot_iR)
  end  # end loop i

  # given flux jacobians wrt entropy variables at the face, compute jacobian
  # of residual contribution wrt entropy variables at solution nodes
  interiorFaceIntegrate_jac!(sbpface, iface, flux_dotL, flux_dotR,
                             jacLL_tmp, jacLR_tmp, jacRL_tmp, jacRR_tmp,
                             SummationByParts.Subtract())
                             


  # now jacLL_tmp and others contain dR/dw, where R and w live on the volume
  # nodes.  Now multiply by A0inv = dw/dq to transform them to dR/dq
  for i=1:sbpface.stencilsize
    i_pL = sbpface.perm[i, iface.faceL]
    i_pR = sbpface.perm[i, iface.faceR]

    # get A0inv for qL and qR
    qL_i = sview(qL, :, i_pL)
    qR_i = sview(qR, :, i_pR)

    getIRA0inv(params, qL_i, A0invL)
    getIRA0inv(params, qR_i, A0invR)

    # multiply all node level jacobians that depend on qL_i and qR_i by
    # A0inv computed using qL_i and qR_i
    # dRL/qL = A0InvL * dRL/wL  for i=1: number of residual nodes affected
    for j=1:sbpface.stencilsize
      j_pL = sbpface.perm[j, iface.faceL]
      j_pR = sbpface.perm[j, iface.faceR]

      # matrices to multiply against
      nodejacLL = sview(jacLL_tmp, :, :, j_pL, i_pL)
      nodejacLR = sview(jacLR_tmp, :, :, j_pL, i_pR)
      nodejacRL = sview(jacRL_tmp, :, :, j_pR, i_pL)
      nodejacRR = sview(jacRR_tmp, :, :, j_pR, i_pR)

      # output matrices (accumulated into)
      nodejac2LL = sview(jacLL, :, :, j_pL, i_pL)
      nodejac2LR = sview(jacLR, :, :, j_pL, i_pR)
      nodejac2RL = sview(jacRL, :, :, j_pR, i_pL)
      nodejac2RR = sview(jacRR, :, :, j_pR, i_pR)

      # transposed mat-mat is faster than regular mat-mat for column major
      # arrays, and A0inv is symmetric.
      # If  jacLL and others were in a compressed representation based on
      # the sparsity of sbpface.perm, the loop j could be replaced with a
      # single large mat-mat
      smallmatmat_kernel!(nodejacLL, A0invL, nodejac2LL, 1, 1)
      smallmatmat_kernel!(nodejacRL, A0invL, nodejac2RL, 1, 1)
      smallmatmat_kernel!(nodejacLR, A0invR, nodejac2LR, 1, 1)
      smallmatmat_kernel!(nodejacRR, A0invR, nodejac2RR, 1, 1)
    end
  end

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
