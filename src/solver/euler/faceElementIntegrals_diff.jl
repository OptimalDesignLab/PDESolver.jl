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

"""
  Reverse mode of [`calcECFaceIntegral`](@ref) wrt the metrics
"""
function calcECFaceIntegral_revm(
     params::AbstractParamType{Tdim}, 
     sbpface::DenseFace, 
     iface::Interface,
     qL::AbstractMatrix{Tsol}, 
     qR::AbstractMatrix{Tsol}, 
     aux_vars::AbstractMatrix{Tres}, 
     nrm_xy::AbstractMatrix{Tmsh},
     nrm_bar::AbstractMatrix{Tmsh},
     functor::FluxType, 
     resL_bar::AbstractMatrix{Tres}, 
     resR_bar::AbstractMatrix{Tres}) where {Tdim, Tsol, Tres, Tmsh}

  data = params.calc_ec_face_integral_data
  @unpack data fluxD nrmD
  numDofPerNode = size(fluxD, 1)

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
      functor(params, qi, qj, aux_vars_i, nrmD, fluxD)

      @simd for dim = 1:Tdim  # move this inside the j loop, at least
        # accumulate entry p_i, p_j of E
        #Eij = zero(Tres)
        Eij_bar = zero(Tres)
        #@simd for k = 1:sbpface.numnodes
        #  nrm_k = nrm_xy[dim, k]
        #  kR = sbpface.nbrperm[k, iface.orient]
        #  Eij += sbpface.interp[i,k]*sbpface.interp[j,kR]*sbpface.wface[k]*nrm_k
        #end  # end loop k
 
       
        @simd for p=1:numDofPerNode
          #resL[p, p_i] -= Eij*fluxD[p, dim]
          #resR[p, p_j] += Eij*fluxD[p, dim]

          #----------------------------------
          # reverse sweep
          Eij_bar -= fluxD[p, dim]*resL_bar[p, p_i]
          Eij_bar += fluxD[p, dim]*resR_bar[p, p_j]
        end

        @simd for k=1:sbpface.numnodes
          kR = sbpface.nbrperm[k, iface.orient]
          nrm_bar[dim, k] += sbpface.interp[i,k]*sbpface.interp[j,kR]*sbpface.wface[k]*Eij_bar
        end

        # none of the remaining computation depends on the metrics

      end  # end loop dim
    end  # end loop j
  end  # end loop i


  return nothing
end


"""
  Reverse mode of [`calcECFaceIntegral`](@ref) wrt the q
"""
function calcECFaceIntegral_revq(
     params::AbstractParamType{Tdim}, 
     sbpface::DenseFace, 
     iface::Interface,
     qL::AbstractMatrix{Tsol}, qL_bar::AbstractMatrix{Tres},
     qR::AbstractMatrix{Tsol}, qR_bar::AbstractMatrix{Tres},
     aux_vars::AbstractMatrix{Tres}, 
     nrm_xy::AbstractMatrix{Tmsh},
     functor_revq::FluxType_revq, 
     resL_bar::AbstractMatrix{Tres}, 
     resR_bar::AbstractMatrix{Tres}) where {Tdim, Tsol, Tres, Tmsh}

  data = params.calc_ec_face_integral_data
  @unpack data nrmD fluxD fluxD_bar
  numDofPerNode = size(fluxD, 1)


  fill!(nrmD, 0.0)
  for d=1:Tdim
    nrmD[d, d] = 1
  end

  # loop over the nodes of "left" element that are in the stencil of interp
  for i = 1:sbpface.stencilsize
    p_i = sbpface.perm[i, iface.faceL]
    qi = ro_sview(qL, :, p_i)
    qi_bar = sview(qL_bar, :, p_i)
    aux_vars_i = ro_sview(aux_vars, :, p_i)  # !!!! why no aux_vars_j???

    # loop over the nodes of "right" element that are in the stencil of interp
    for j = 1:sbpface.stencilsize
      p_j = sbpface.perm[j, iface.faceR]
      qj = ro_sview(qR, :, p_j)
      qj_bar = sview(qR_bar, :, p_j)

      # compute flux and add contribution to left and right elements
      #functor(params, qi, qj, aux_vars_i, nrmD, fluxD)

      fill!(fluxD_bar, 0)
      @simd for dim = 1:Tdim  # move this inside the j loop, at least
        # accumulate entry p_i, p_j of E
        Eij = zero(Tres)
        @simd for k = 1:sbpface.numnodes
          nrm_k = nrm_xy[dim, k]
          kR = sbpface.nbrperm[k, iface.orient]
          Eij += sbpface.interp[i,k]*sbpface.interp[j,kR]*sbpface.wface[k]*nrm_k
        end  # end loop k
 
       
        @simd for p=1:numDofPerNode
          #resL[p, p_i] -= Eij*fluxD[p, dim]
          #resR[p, p_j] += Eij*fluxD[p, dim]

          #----------------------------------
          # reverse sweep
          fluxD_bar[p, dim] -= Eij*resL_bar[p, p_i]
          fluxD_bar[p, dim] += Eij*resR_bar[p, p_j]
        end
      end  # end loop dim

      functor_revq(params, qi, qi_bar, qj, qj_bar, aux_vars, nrmD, fluxD_bar)
    end  # end loop j
  end  # end loop i


  return nothing
end




#------------------------------------------------------------------------------
# calcESFaceIntegral

"""
  Differentiated version of [`calcESFaceIntegral`](@ref)

  **Inputs**

   * `params`: `AbstractParamType`
   * `sbpface`: an `AbstractFace`.  Methods are available for both sparse
              and dense faces
   * `iface`: the [`Interface`](@ref) object for the given face
   * `kernel`: an [`AbstractEntropyKernel`](@ref) specifying what kind of
               dissipation to apply.
   * `qL`: the solution at the volume nodes of the left element (`numDofPerNode`
     x `numNodesPerElement)
   * `qR`: the solution at the volume nodes of the right element
   * `aux_vars`: the auxiliary variables for `qL`
   * `nrm_xy`: the normal vector at each face node, `dim` x `numNodesPerFace`
   * `functor: the differentiated flux function, of type [`FluxType_diff`](@ref)

   **Inputs/Outputs**

    * jacLL: jacobian of `resL` wrt `qL`
    * jacLR: jacobian of `resL` wrt `qR`
    * jacRL: jacobian of `resR` wrt `qL`
    * jacRR: jacobian of `resR` wrt `qR`
"""
function calcESFaceIntegral_diff(
     params::AbstractParamType{Tdim}, 
     sbpface::AbstractFace, 
     iface::Interface,
     kernel::AbstractEntropyKernel,
     qL::AbstractMatrix{Tsol}, 
     qR::AbstractMatrix{Tsol}, 
     aux_vars::AbstractMatrix{Tres}, 
     nrm_face::AbstractMatrix{Tmsh},
     functor::FluxType_diff,
     jacLL::AbstractArray{Tres, 4}, jacLR::AbstractArray{Tres, 4},
     jacRL::AbstractArray{Tres, 4}, jacRR::AbstractArray{Tres, 4}) where {Tdim, Tsol, Tres, Tmsh}

  calcECFaceIntegral_diff(params, sbpface, iface, qL, qR, aux_vars, nrm_face, 
                     functor, jacLL, jacLR, jacRL, jacRR)
  calcEntropyPenaltyIntegral_diff(params, sbpface, iface, kernel, qL, qR,
                                  aux_vars, nrm_face, jacLL, jacLR, jacRL, jacRR)

  return nothing
end

"""
  Reverse mode of [`calcESFaceIntegrals`](@ref)

  **Inputs**

   * params
   * sbpface
   * iface
   * kernel
   * qL
   * qR
   * aux_vars
   * nrm_face
   * functor: the same functor passed into the primal method (the reverse
              mode functor is not required)
   * resL_bar
   * resR_bar

  **Inputs/Outputs**

   * nrm_face_bar: updated with back-propigation of `resL_bar`, `resR_bar`,
                   same shape as `nrm_face`
"""
function calcESFaceIntegral_revm(
     params::AbstractParamType{Tdim}, 
     sbpface::AbstractFace, 
     iface::Interface,
     kernel::AbstractEntropyKernel,
     qL::AbstractMatrix{Tsol}, 
     qR::AbstractMatrix{Tsol}, 
     aux_vars::AbstractMatrix{Tres}, 
     nrm_face::AbstractMatrix{Tmsh},
     nrm_face_bar::AbstractMatrix{Tmsh},
     functor::FluxType, # interestingly, functor_revm is never required
     resL_bar::AbstractMatrix{Tres}, 
     resR_bar::AbstractMatrix{Tres}) where {Tdim, Tsol, Tres, Tmsh}

  calcECFaceIntegral_revm(params, sbpface, iface, qL, qR, aux_vars, nrm_face, 
                     nrm_face_bar, functor, resL_bar, resR_bar)
  calcEntropyPenaltyIntegral_revm(params, sbpface, iface, kernel, qL, qR, aux_vars, 
                               nrm_face, nrm_face_bar, resL_bar, resR_bar)

  return nothing
end


function calcESFaceIntegral_revq(
     params::AbstractParamType{Tdim}, 
     sbpface::AbstractFace, 
     iface::Interface,
     kernel::AbstractEntropyKernel,
     qL::AbstractMatrix{Tsol}, qL_bar::AbstractMatrix{Tres},
     qR::AbstractMatrix{Tsol}, qR_bar::AbstractMatrix{Tres},
     aux_vars::AbstractMatrix{Tres}, 
     nrm_face::AbstractMatrix{Tmsh},
     functor_revq::FluxType_revq,
     resL_bar::AbstractMatrix{Tres}, 
     resR_bar::AbstractMatrix{Tres}) where {Tdim, Tsol, Tres, Tmsh}

  calcECFaceIntegral_revq(params, sbpface, iface, qL, qL_bar, qR, qR_bar,
                     aux_vars, nrm_face, functor_revq, resL_bar, resR_bar)
  calcEntropyPenaltyIntegral_revq(params, sbpface, iface, kernel, qL, qL_bar,
                               qR, qR_bar, aux_vars, nrm_face,
                               resL_bar, resR_bar)

  return nothing
end



#------------------------------------------------------------------------------
# calcEntropyPenaltyIntegral

function calcEntropyPenaltyIntegral_diff(
             params::ParamType{Tdim, :conservative},
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

  flux_dot_iL = sview(flux_dot_i, :, 1:numDofPerNode)  # wL part
  flux_dot_iR = sview(flux_dot_i, :, (numDofPerNode + 1):nd)  # wR part

  fill!(delta_w_dot, 0.0)
  fill!(jacLL_tmp, 0.0); fill!(jacLR_tmp, 0.0);
  fill!(jacRL_tmp, 0.0); fill!(jacRR_tmp, 0.0);


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
  interiorFaceCombined_jac!(sbpface, iface, flux_dotL, flux_dotR,
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


"""
  Reverse mode of [`calcEntropyPenaltyIntegral`](@ref) wrt the metrics

  **Inputs**

   * params
   * sbpface
   * kernel
   * qL
   * qR
   * aux_vars
   * nrm_face
   * resL_bar, resR_bar

  **Inputs/Outputs**

   * nrm_bar: updated with back propigation from `resL_bar`, `resR_bar`
"""
function calcEntropyPenaltyIntegral_revm(
           params::ParamType{Tdim, :conservative},
           sbpface::DenseFace, iface::Interface,
           kernel::AbstractEntropyKernel,
           qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, 
           aux_vars::AbstractMatrix{Tres},
           nrm_face::AbstractArray{Tmsh, 2}, nrm_face_bar::AbstractArray{Tmsh, 2},
           resL_bar::AbstractMatrix{Tres}, resR_bar::AbstractMatrix{Tres}) where {Tdim, Tsol, Tres, Tmsh}

  numDofPerNode = size(qL, 1)

  # convert qL and qR to entropy variables (only the nodes that will be used)
  data = params.calc_entropy_penalty_integral_data
  @unpack data wL wR wL_i wR_i qL_i qR_i delta_w q_avg flux flux_bar

  # convert to IR entropy variables
  for i=1:sbpface.stencilsize
    # apply sbpface.perm here
    p_iL = sbpface.perm[i, iface.faceL]
    p_iR = sbpface.perm[i, iface.faceR]
    # these need to have different names from qL_i etc. below to avoid type
    # instability
    qL_itmp = ro_sview(qL, :, p_iL)
    qR_itmp = ro_sview(qR, :, p_iR)
    wL_itmp = sview(wL, :, i)  # wL is a compacted representation: it only
                               # stores the required values, in the order they
                               # will be used later
    wR_itmp = sview(wR, :, i)
    convertToIR(params, qL_itmp, wL_itmp)
    convertToIR(params, qR_itmp, wR_itmp)
  end



  # accumulate wL at the node
  @simd for i=1:sbpface.numnodes  # loop over face nodes
    ni = sbpface.nbrperm[i, iface.orient]
    dir = ro_sview(nrm_face, :, i)
    dir_bar = sview(nrm_face_bar, :, i)
    fastzero!(wL_i); fastzero!(wR_i)
    fastzero!(flux_bar)

    # interpolate wL and wR to this node
    @simd for j=1:sbpface.stencilsize
      interpL = sbpface.interp[j, i]
      interpR = sbpface.interp[j, ni]

      @simd for k=1:numDofPerNode
        wL_i[k] += interpL*wL[k, j]
        wR_i[k] += interpR*wR[k, j]
      end
    end

    #TODO: write getLambdaMaxSimple and getIRA0 in terms of the entropy
    #      variables to avoid the conversion
    convertToConservativeFromIR_(params, wL_i, qL_i)
    convertToConservativeFromIR_(params, wR_i, qR_i)
    
    # compute average qL
    # also delta w (used later)
    @simd for j=1:numDofPerNode
      q_avg[j] = 0.5*(qL_i[j] + qR_i[j])
      delta_w[j] = wL_i[j] - wR_i[j]
    end

    # call kernel (apply symmetric semi-definite matrix)
    #applyEntropyKernel(kernel, params, q_avg, delta_w, dir, flux)
    #for j=1:numDofPerNode
    #  flux[j] *= sbpface.wface[i]
    #end

    # interpolate back to volume nodes
    @simd for j=1:sbpface.stencilsize
      j_pL = sbpface.perm[j, iface.faceL]
      j_pR = sbpface.perm[j, iface.faceR]

      @simd for p=1:numDofPerNode
#        resL[p, j_pL] -= sbpface.interp[j, i]*flux[p]
#        resR[p, j_pR] += sbpface.interp[j, ni]*flux[p]

        #--------------------------
        # reverse sweep
        flux_bar[p] -= sbpface.interp[j, i ]*resL_bar[p, j_pL]
        flux_bar[p] += sbpface.interp[j, ni]*resR_bar[p, j_pR]
      end
    end

    for j=1:numDofPerNode
      flux_bar[j] *= sbpface.wface[i]
    end
    applyEntropyKernel_revm(kernel, params, q_avg, delta_w, dir, dir_bar, flux,
                            flux_bar)

    # nothing further depends on the metrics
  end  # end loop i



  return nothing
end

"""
  Reverse mode of [`calcEntropyPenaltyIntegral`](@ref) wrt the metrics

  **Inputs**

   * params
   * sbpface
   * kernel
   * qL
   * qR
   * aux_vars
   * nrm_face
   * resL_bar, resR_bar

  **Inputs/Outputs**

   * qL_bar, qR_bar: updated with back propigation from `resL_bar`, `resR_bar`
"""
function calcEntropyPenaltyIntegral_revq(
           params::ParamType{Tdim, :conservative},
           sbpface::DenseFace, iface::Interface,
           kernel::AbstractEntropyKernel,
           qL::AbstractMatrix{Tsol}, qL_bar::AbstractMatrix{Tres},
           qR::AbstractMatrix{Tsol}, qR_bar::AbstractMatrix{Tres},
           aux_vars::AbstractMatrix{Tres},
           nrm_face::AbstractArray{Tmsh, 2},
           resL_bar::AbstractMatrix{Tres}, resR_bar::AbstractMatrix{Tres}) where {Tdim, Tsol, Tres, Tmsh}

  numDofPerNode = size(qL, 1)

  # convert qL and qR to entropy variables (only the nodes that will be used)
  data = params.calc_entropy_penalty_integral_data
  @unpack data wL wR wL_i wR_i qL_i qR_i delta_w q_avg flux
  @unpack data A0 qL_bar_i qR_bar_i wL_bar_i wR_bar_i delta_w_bar q_avg_bar
  @unpack data wL_bar wR_bar

  fill!(wL_bar, 0)
  fill!(wR_bar, 0)


  flux_bar = zeros(flux)

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



  # accumulate wL at the node
  @simd for i=1:sbpface.numnodes  # loop over face nodes
    ni = sbpface.nbrperm[i, iface.orient]
    dir = ro_sview(nrm_face, :, i)
    fastzero!(wL_i); fastzero!(qL_bar_i); fastzero!(wL_bar_i)
    fastzero!(wR_i); fastzero!(qR_bar_i); fastzero!(wR_bar_i)
    fastzero!(delta_w_bar); fastzero!(q_avg_bar); fastzero!(flux_bar)


    # interpolate wL and wR to this node
    @simd for j=1:sbpface.stencilsize
      interpL = sbpface.interp[j, i]
      interpR = sbpface.interp[j, ni]

      @simd for k=1:numDofPerNode
        wL_i[k] += interpL*wL[k, j]
        wR_i[k] += interpR*wR[k, j]
      end
    end

    #TODO: write getLambdaMaxSimple and getIRA0 in terms of the entropy
    #      variables to avoid the conversion
    convertToConservativeFromIR_(params, wL_i, qL_i)
    convertToConservativeFromIR_(params, wR_i, qR_i)
    
    # compute average qL
    # also delta w (used later)
    @simd for j=1:numDofPerNode
      q_avg[j] = 0.5*(qL_i[j] + qR_i[j])
      delta_w[j] = wL_i[j] - wR_i[j]
    end

    # call kernel (apply symmetric semi-definite matrix)
    applyEntropyKernel(kernel, params, q_avg, delta_w, dir, flux)
    for j=1:numDofPerNode
      flux[j] *= sbpface.wface[i]
    end

    # interpolate back to volume nodes
    @simd for j=1:sbpface.stencilsize
      j_pL = sbpface.perm[j, iface.faceL]
      j_pR = sbpface.perm[j, iface.faceR]

      @simd for p=1:numDofPerNode
#        resL[p, j_pL] -= sbpface.interp[j, i]*flux[p]
#        resR[p, j_pR] += sbpface.interp[j, ni]*flux[p]

        #--------------------------
        # reverse sweep
        flux_bar[p] -= sbpface.interp[j, i ]*resL_bar[p, j_pL]
        flux_bar[p] += sbpface.interp[j, ni]*resR_bar[p, j_pR]
      end
    end

    for j=1:numDofPerNode
      flux_bar[j] *= sbpface.wface[i]
    end
    applyEntropyKernel_revq(kernel, params, q_avg, q_avg_bar, delta_w,
                            delta_w_bar, dir, flux, flux_bar)

    # none of the remaining calculation depends on the metricsa
    # compute average state
    @simd for j=1:numDofPerNode
      qL_bar_i[j] += 0.5*q_avg_bar[j]
      qR_bar_i[j] += 0.5*q_avg_bar[j]
      wL_bar_i[j] += delta_w_bar[j]
      wR_bar_i[j] -= delta_w_bar[j]
    end

    # reverse of convertToConservatvieFromIR
    getIRA0(params, qL_i, A0)
    smallmatTvec_kernel!(A0, qL_bar_i, wL_bar_i, 1, 1)
    getIRA0(params, qR_i, A0)
    smallmatTvec_kernel!(A0, qR_bar_i, wR_bar_i, 1, 1)

    # interpolate wL and wR to this node
    @simd for j=1:sbpface.stencilsize
      interpL = sbpface.interp[j, i]
      interpR = sbpface.interp[j, ni]

      @simd for k=1:numDofPerNode
        #wL_i[k] += interpL*wL[k, j]
        #wR_i[k] += interpR*wR[k, j]
        wL_bar[k, j] += interpL*wL_bar_i[k]
        wR_bar[k, j] += interpR*wR_bar_i[k]
      end
    end
    
  end  # end loop i

  # convert to IR entropy variables
  for i=1:sbpface.stencilsize
    # apply sbpface.perm here
    p_iL = sbpface.perm[i, iface.faceL]
    p_iR = sbpface.perm[i, iface.faceR]
    # these need to have different names from qL_i etc. below to avoid type
    # instability
    qL_itmp = ro_sview(qL, :, p_iL)
    qR_itmp = ro_sview(qR, :, p_iR)

    qL_bar_itmp = sview(qL_bar, :, p_iL)
    qR_bar_itmp = sview(qR_bar, :, p_iR)
    wL_bar_itmp = sview(wL_bar, :, i)
    wR_bar_itmp = sview(wR_bar, :, i)
    getIRA0inv(params, qL_itmp, A0)
    smallmatTvec_kernel!(A0, wL_bar_itmp, qL_bar_itmp, 1, 1)
    getIRA0inv(params, qR_itmp, A0)
    smallmatTvec_kernel!(A0, wR_bar_itmp, qR_bar_itmp, 1, 1)
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


function applyEntropyKernel_revq(obj::LFKernel, params::ParamType, 
                        q_avg::AbstractVector, q_bar::AbstractVector,
                        delta_w::AbstractVector, delta_w_bar::AbstractVector,
                        nrm::AbstractVector, flux::AbstractVector,
                        flux_bar::AbstractVector)

  @unpack obj t1 A0 t1_bar A0_bar

  numDofPerNode = length(q_avg)
  
  getIRA0(params, q_avg, A0)
  #for j=1:size(A0, 1)
  #  A0[j, j] = 1
  #end

  lambda_max = getLambdaMax(params, q_avg, nrm)
  # lambda_max * A0 * delta w

  smallmatvec!(A0, delta_w, t1)
  for i=1:length(flux)
    flux[i] = lambda_max*t1[i]
  end

  # reverse sweep

  # flux = t1*lambda_max
  lambda_max_bar = zero(lambda_max)
  for i=1:numDofPerNode
    lambda_max_bar += flux_bar[i]*t1[i]
    t1_bar[i] = flux_bar[i]*lambda_max
  end

  # t1 = A0*delta_w
  for i=1:numDofPerNode
    for j=1:numDofPerNode
      A0_bar[j, i] = t1_bar[j]*delta_w[i]
      delta_w_bar[i] += A0[j, i]*t1_bar[j]
    end
  end

  getLambdaMax_revq(params, q_avg, q_bar, nrm, lambda_max_bar)

  getIRA0_revq(params, q_avg, q_bar, A0, A0_bar)

  return nothing
end


function applyEntropyKernel_revm(obj::LFKernel, params::ParamType, 
                        q_avg::AbstractVector, delta_w::AbstractVector,
                        nrm::AbstractVector, nrm_bar::AbstractVector,
                        flux::AbstractVector,
                        flux_bar::AbstractVector)


  @unpack obj t1 A0 A0_bar
  
  numDofPerNode = length(q_avg)
  
  getIRA0(params, q_avg, A0)
  #for j=1:size(A0, 1)
  #  A0[j, j] = 1
  #end

  lambda_max = getLambdaMax(params, q_avg, nrm)
  # lambda_max * A0 * delta w

  smallmatvec!(A0, delta_w, t1)
  for i=1:length(flux)
    flux[i] = lambda_max*t1[i]
  end

  # reverse sweep

  # flux = t1*lambda_max
  lambda_max_bar = zero(lambda_max)
  for i=1:numDofPerNode
    lambda_max_bar += flux_bar[i]*t1[i]
  end

  getLambdaMax_revm(params, q_avg, nrm, nrm_bar, lambda_max_bar)

  return nothing
end


function applyEntropyKernel_revm(obj::IdentityKernel, params::ParamType, 
                        q_avg::AbstractVector, delta_w::AbstractVector,
                        nrm::AbstractVector, nrm_bar::AbstractVector,
                        flux::AbstractVector,
                        flux_bar::AbstractVector)

  # nothing to do

  return nothing
end

function applyEntropyKernel_revq(obj::IdentityKernel, params::ParamType, 
                        q_avg::AbstractVector, q_bar::AbstractVector,
                        delta_w::AbstractVector, delta_w_bar::AbstractVector,
                        nrm::AbstractVector, flux::AbstractVector,
                        flux_bar::AbstractVector)

  for i=1:length(q_avg)
    delta_w_bar[i] += flux_bar[i]
  end

  return nothing
end



#------------------------------------------------------------------------------
# apply entropy kernels to diagonal E operators

"""
  Differentiated version of [`applyEntropyKernel_diagE`](@ref)
"""
function applyEntropyKernel_diagE_diff(
                      params::ParamType{Tdim, :conservative},
                      kernel::AbstractEntropyKernel,
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tres},
                      dir::AbstractArray{Tmsh},
                      FL_dot::AbstractMatrix{Tres}, FR_dot::AbstractMatrix{Tres},
                      ) where {Tmsh, Tsol, Tres, Tdim}

#  nd = 2*params.numDofPerNode
  data = params.apply_entropy_kernel_diagE_data
  @unpack data q_avg q_avg_dot F F_dot

  fill!(q_avg_dot, 0.0); fill!(F_dot, 0.0)
  @debug1 begin
    @assert size(FL_dot, 2) <= params.numDofPerNode
    @assert size(FR_dot, 2) <= params.numDofPerNode
  end

  for i=1:length(q_avg)
    q_avg[i] = 0.5*(qL[i] + qR[i])
    q_avg_dot[i, i] = 0.5
    q_avg_dot[i, i + params.numDofPerNode] = 0.5
  end

  applyEntropyKernel_diagE_inner_diff(params, kernel, qL, qR, q_avg, q_avg_dot, aux_vars, dir, F, F_dot)

  # copy derivative into arrays in the format the caller provided
  @simd for i=1:params.numDofPerNode
    @simd for j=1:params.numDofPerNode
      FL_dot[j, i] += F_dot[j, i                       ]
      FR_dot[j, i] += F_dot[j, i + params.numDofPerNode]
    end
  end

  return nothing
end

"""
  Reverse mode wrt metrics of [`applyEntropyKernel_diagE`](@ref)

  **Inputs**
  
   * params
   * kernel
   * qL
   * qR
   * aux_var
   * dir
   * F_bar

  **Inputs/Outputs**

   * dir_bar: updated with back-propigation of `F_bar`

"""
function applyEntropyKernel_diagE_revm(
                      params::ParamType{Tdim, :conservative},
                      kernel::AbstractEntropyKernel,
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      aux_vars::AbstractArray{Tres},
                      dir::AbstractArray{Tmsh}, dir_bar::AbstractArray{Tmsh},
                      F_bar::AbstractArray{Tres,1}) where {Tmsh, Tsol, Tres, Tdim}

  q_avg = params.apply_entropy_kernel_diagE_data.q_avg

  for i=1:length(q_avg)
    q_avg[i] = 0.5*(qL[i] + qR[i])
  end

  applyEntropyKernel_diagE_inner_revm(params, kernel, qL, qR, q_avg, aux_vars,
                                      dir, dir_bar, F_bar)

  return nothing
end


"""
  Reverse mode wrt q of [`applyEntropyKernel_diagE`](@ref)

  **Inputs**
  
   * params
   * kernel
   * qL
   * qR
   * aux_var
   * dir
   * F_bar

  **Inputs/Outputs**

   * qL_bar, qR_bar: updated with back-propigation of `F_bar`

"""
function applyEntropyKernel_diagE_revq(
                      params::ParamType{Tdim, :conservative},
                      kernel::AbstractEntropyKernel,
                      qL::AbstractArray{Tsol,1}, qL_bar::AbstractArray{Tres, 1},
                      qR::AbstractArray{Tsol,1}, qR_bar::AbstractArray{Tres, 1},
                      aux_vars::AbstractArray{Tres},
                      dir::AbstractArray{Tmsh},
                      F_bar::AbstractArray{Tres,1}) where {Tmsh, Tsol, Tres, Tdim}

  data = params.apply_entropy_kernel_diagE_data
  @unpack data q_avg q_avg_bar

  fill!(q_avg_bar, 0)

  for i=1:length(q_avg)
    q_avg[i] = 0.5*(qL[i] + qR[i])
  end

  applyEntropyKernel_diagE_inner_revq(params, kernel, qL, qL_bar, qR, qR_bar,
                                      q_avg, q_avg_bar, aux_vars, dir, F_bar)

  for i=1:length(q_avg)
    qL_bar[i] += 0.5*q_avg_bar[i]
    qR_bar[i] += 0.5*q_avg_bar[i]
  end

  return nothing
end






function applyEntropyKernel_diagE_inner_diff(
                      params::ParamType{Tdim, :conservative}, 
                      kernel::AbstractEntropyKernel,
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      q_avg::AbstractArray{Tsol}, q_avg_dot::AbstractMatrix{Tsol},
                      aux_vars::AbstractArray{Tres},
                      dir::AbstractArray{Tmsh},
                      F::AbstractArray{Tres,1}, F_dot::AbstractMatrix{Tres}
                      ) where {Tmsh, Tsol, Tres, Tdim}

  @unpack params.apply_entropy_kernel_diagE_data vL vR F_tmp delta_w_dot
  gamma = params.gamma
  gamma_1inv = 1/params.gamma_1
#  p = calcPressure(params, q_avg)

  nd = 2*params.numDofPerNode
  delta_w_dotL = sview(delta_w_dot, :, 1:params.numDofPerNode)
  delta_w_dotR = sview(delta_w_dot, :, (params.numDofPerNode+1):nd)
   
  convertToIR(params, qL, vL)
  convertToIR(params, qR, vR)

  for i=1:length(vL)
    vL[i] = vL[i] - vR[i]
  end
  # A0inv = inv(dq/dw) = dw/dq = the derivative of convertToIR
  getIRA0inv(params, qL, delta_w_dotL)
  getIRA0inv(params, qR, delta_w_dotR)
  scale!(delta_w_dotR, -1)


  applyEntropyKernel_diff(kernel, params, q_avg, q_avg_dot, vL, delta_w_dot, dir, F_tmp, F_dot)

  for i=1:length(F_tmp)
    F[i] += F_tmp[i]
  end

  return nothing
end


function applyEntropyKernel_diagE_inner_revm(
                      params::ParamType{Tdim, :conservative}, 
                      kernel::AbstractEntropyKernel,
                      qL::AbstractArray{Tsol,1}, qR::AbstractArray{Tsol, 1},
                      q_avg::AbstractArray{Tsol}, aux_vars::AbstractArray{Tres},
                      dir::AbstractArray{Tmsh}, dir_bar::AbstractArray{Tmsh},
                      F_bar::AbstractArray{Tres,1}) where {Tmsh, Tsol, Tres, Tdim}
  @unpack params.apply_entropy_kernel_diagE_data vL vR F_tmp
  gamma = params.gamma
  gamma_1inv = 1/params.gamma_1
#  p = calcPressure(params, q_avg)

  convertToIR(params, qL, vL)
  convertToIR(params, qR, vR)

  for i=1:length(vL)
    vL[i] = vL[i] - vR[i]
  end

  applyEntropyKernel_revm(kernel, params, q_avg, vL, dir, dir_bar, F_tmp, F_bar)

  return nothing
end


function applyEntropyKernel_diagE_inner_revq(
                      params::ParamType{Tdim, :conservative}, 
                      kernel::AbstractEntropyKernel,
                      qL::AbstractArray{Tsol,1}, qL_bar::AbstractArray{Tres,1},
                      qR::AbstractArray{Tsol,1}, qR_bar::AbstractArray{Tres, 1},
                      q_avg::AbstractArray{Tsol},
                      q_avg_bar::AbstractArray{Tres, 1},
                      aux_vars::AbstractArray{Tres},
                      dir::AbstractArray{Tmsh},
                      F_bar::AbstractArray{Tres,1}) where {Tmsh, Tsol, Tres, Tdim}
  @unpack params.apply_entropy_kernel_diagE_data vL vR F_tmp delta_w delta_w_bar A0inv

  fill!(delta_w_bar, 0)

  gamma = params.gamma
  gamma_1inv = 1/params.gamma_1
#  p = calcPressure(params, q_avg)

  convertToIR(params, qL, vL)
  convertToIR(params, qR, vR)

  for i=1:length(vL)
    delta_w[i] = vL[i] - vR[i]
  end

  applyEntropyKernel_revq(kernel, params, q_avg, q_avg_bar, delta_w, delta_w_bar, dir, F_tmp, F_bar)

  # combine the delta_w and convertToIR steps
  getIRA0inv(params, qL, A0inv)
  smallmatTvec_kernel!(A0inv, delta_w_bar, qL_bar, 1, 1)
  getIRA0inv(params, qR, A0inv)
  smallmatTvec_kernel!(A0inv, delta_w_bar, qR_bar, -1, 1)

  return nothing
end




#------------------------------------------------------------------------------
# extend functors with differentiated method

"""
  Differentiated method for `ECFaceIntegral`
"""
function calcFaceElementIntegral_diff(obj::ECFaceIntegral,
              params::AbstractParamType{Tdim}, 
              sbpface::AbstractFace, iface::Interface,
              qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, 
              aux_vars::AbstractMatrix{Tres}, nrm_face::AbstractMatrix{Tmsh},
              functor::FluxType_diff, 
              jacLL::AbstractArray{Tres, 4}, jacLR::AbstractArray{Tres, 4},
              jacRL::AbstractArray{Tres, 4}, jacRR::AbstractArray{Tres, 4}) where {Tsol, Tres, Tmsh, Tdim}


  calcECFaceIntegral_diff(params, sbpface, iface, qL, qR, aux_vars, nrm_face, 
                      functor, jacLL, jacLR, jacRL, jacRR)

end

"""
  revm method for `ECFaceIntegral`
"""
function calcFaceElementIntegral_revm(obj::ECFaceIntegral,
              params::AbstractParamType{Tdim}, 
              sbpface::AbstractFace, iface::Interface,
              qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, 
              aux_vars::AbstractMatrix{Tres}, nrm_face::AbstractMatrix{Tmsh},
              nrm_face_bar::AbstractMatrix{Tmsh},
              functor::FluxType, 
              resL_bar::AbstractMatrix{Tres}, resR_bar::AbstractMatrix{Tres}) where {Tsol, Tres, Tmsh, Tdim}

  calcECFaceIntegral_revm(params, sbpface, iface, qL, qR, aux_vars, nrm_face,
                          nrm_face_bar, functor, resL_bar, resR_bar)
end

"""
  revq method for `ECFaceIntegral`
"""
function calcFaceElementIntegral_revq(obj::ECFaceIntegral,
              params::AbstractParamType{Tdim}, 
              sbpface::AbstractFace, iface::Interface,
              qL::AbstractMatrix{Tsol}, qL_bar::AbstractMatrix{Tres},
              qR::AbstractMatrix{Tsol}, qR_bar::AbstractMatrix{Tres},
              aux_vars::AbstractMatrix{Tres}, nrm_face::AbstractMatrix{Tmsh},
              functor_revq::FluxType_revq, 
              resL_bar::AbstractMatrix{Tres}, resR_bar::AbstractMatrix{Tres}) where {Tsol, Tres, Tmsh, Tdim}

  calcECFaceIntegral_revq(params, sbpface, iface, qL, qL_bar, qR, qR_bar,
                          aux_vars, nrm_face, functor_revq, resL_bar, resR_bar)
end

#--------------------------
# ESLFFaceIntegral

"""
  Differentiated method for `ESLFFaceIntegral`
"""
function calcFaceElementIntegral_diff(obj::ESLFFaceIntegral,
              params::AbstractParamType{Tdim}, 
              sbpface::AbstractFace, iface::Interface,
              qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, 
              aux_vars::AbstractMatrix{Tres}, nrm_face::AbstractMatrix{Tmsh},
              functor::FluxType_diff, 
              jacLL::AbstractArray{Tres, 4}, jacLR::AbstractArray{Tres, 4},
              jacRL::AbstractArray{Tres, 4}, jacRR::AbstractArray{Tres, 4}) where {Tsol, Tres, Tmsh, Tdim}


  calcESFaceIntegral_diff(params, sbpface, iface, obj.kernel, qL, qR, aux_vars,
                          nrm_face, functor, jacLL, jacLR, jacRL, jacRR)

end

"""
  revm method for `ESLFFaceIntegral`
"""
function calcFaceElementIntegral_revm(obj::ESLFFaceIntegral,
              params::AbstractParamType{Tdim}, 
              sbpface::AbstractFace, iface::Interface,
              qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, 
              aux_vars::AbstractMatrix{Tres}, nrm_face::AbstractMatrix{Tmsh},
              nrm_face_bar::AbstractMatrix{Tmsh},
              functor::FluxType, 
              resL_bar::AbstractMatrix{Tres}, resR_bar::AbstractMatrix{Tres}) where {Tsol, Tres, Tmsh, Tdim}

  calcESFaceIntegral_revm(params, sbpface, iface, obj.kernel, qL, qR, aux_vars,
                          nrm_face, nrm_face_bar, functor, resL_bar, resR_bar)
end

"""
  revq method for `ESLFFaceIntegral`
"""
function calcFaceElementIntegral_revq(obj::ESLFFaceIntegral,
              params::AbstractParamType{Tdim}, 
              sbpface::AbstractFace, iface::Interface,
              qL::AbstractMatrix{Tsol}, qL_bar::AbstractMatrix{Tres},
              qR::AbstractMatrix{Tsol}, qR_bar::AbstractMatrix{Tres},
              aux_vars::AbstractMatrix{Tres}, nrm_face::AbstractMatrix{Tmsh},
              functor_revq::FluxType_revq, 
              resL_bar::AbstractMatrix{Tres}, resR_bar::AbstractMatrix{Tres}) where {Tsol, Tres, Tmsh, Tdim}

  calcESFaceIntegral_revq(params, sbpface, iface, obj.kernel, qL, qL_bar,
                          qR, qR_bar, aux_vars, nrm_face, functor_revq,
                          resL_bar, resR_bar)
end


#------------------------
# ELFPenaltyIntegral

"""
  Differentiated method for `ELFPenaltyIntegral`
"""
function calcFaceElementIntegral_diff(obj::ELFPenaltyFaceIntegral,
              params::AbstractParamType{Tdim}, 
              sbpface::AbstractFace, iface::Interface,
              qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, 
              aux_vars::AbstractMatrix{Tres}, nrm_face::AbstractMatrix{Tmsh},
              functor::FluxType_diff, 
              jacLL::AbstractArray{Tres, 4}, jacLR::AbstractArray{Tres, 4},
              jacRL::AbstractArray{Tres, 4}, jacRR::AbstractArray{Tres, 4}) where {Tsol, Tres, Tmsh, Tdim}


  calcEntropyPenaltyIntegral_diff(params, sbpface, iface, obj.kernel, qL, qR,
                        aux_vars, nrm_face, functor, jacLL, jacLR, jacRL, jacRR)

end


function calcFaceElementIntegral_revm(obj::ELFPenaltyFaceIntegral,
              params::AbstractParamType{Tdim}, 
              sbpface::AbstractFace, iface::Interface,
              qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, 
              aux_vars::AbstractMatrix{Tres}, nrm_face::AbstractMatrix{Tmsh},
              nrm_face_bar::AbstractMatrix{Tmsh},
              functor::FluxType, 
              resL_bar::AbstractMatrix{Tres}, resR_bar::AbstractMatrix{Tres}) where {Tsol, Tres, Tmsh, Tdim}

  calcEntropyPenaltyIntegral_revm(params, sbpface, iface, obj.kernel, qL, qR,
                aux_vars, nrm_face, nrm_face_bar, resL_bar, resR_bar)
end


function calcFaceElementIntegral_revq(obj::ELFPenaltyFaceIntegral,
              params::AbstractParamType{Tdim}, 
              sbpface::AbstractFace, iface::Interface,
              qL::AbstractMatrix{Tsol}, qL_bar::AbstractMatrix{Tres},
              qR::AbstractMatrix{Tsol}, qR_bar::AbstractMatrix{Tres},
              aux_vars::AbstractMatrix{Tres}, nrm_face::AbstractMatrix{Tmsh},
              functor_revq::FluxType_revq, 
              resL_bar::AbstractMatrix{Tres}, resR_bar::AbstractMatrix{Tres}) where {Tsol, Tres, Tmsh, Tdim}

  calcEntropyPenaltyIntegral_revq(params, sbpface, iface, obj.kernel, qL,
                qL_bar, qR, qR_bar, aux_vars, nrm_face,
                resL_bar, resR_bar)
end


#---------------------
# entropy jump penalty

function calcFaceElementIntegral_revm(obj::EntropyJumpPenaltyFaceIntegral,
              params::AbstractParamType{Tdim}, 
              sbpface::AbstractFace, iface::Interface,
              qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, 
              aux_vars::AbstractMatrix{Tres}, nrm_face::AbstractMatrix{Tmsh},
              nrm_face_bar::AbstractMatrix{Tmsh},
              functor::FluxType, 
              resL_bar::AbstractMatrix{Tres}, resR_bar::AbstractMatrix{Tres}) where {Tsol, Tres, Tmsh, Tdim}

  calcEntropyPenaltyIntegral_revm(params, sbpface, iface, obj.kernel, qL, qR,
                aux_vars, nrm_face, nrm_face_bar, resL_bar, resR_bar)
end


function calcFaceElementIntegral_revq(obj::EntropyJumpPenaltyFaceIntegral,
              params::AbstractParamType{Tdim}, 
              sbpface::AbstractFace, iface::Interface,
              qL::AbstractMatrix{Tsol}, qL_bar::AbstractMatrix{Tres},
              qR::AbstractMatrix{Tsol}, qR_bar::AbstractMatrix{Tres},
              aux_vars::AbstractMatrix{Tres}, nrm_face::AbstractMatrix{Tmsh},
              functor_revq::FluxType_revq, 
              resL_bar::AbstractMatrix{Tres}, resR_bar::AbstractMatrix{Tres}) where {Tsol, Tres, Tmsh, Tdim}

  calcEntropyPenaltyIntegral_revq(params, sbpface, iface, obj.kernel, qL,
                qL_bar, qR, qR_bar, aux_vars, nrm_face,
                resL_bar, resR_bar)
end


