# functions that do face integral-like operations, but operate on data from
# the entire element
include("IR_stab.jl")  # stabilization for the IR flux

# naming convention
# EC -> entropy conservative
# ES -> entropy stable (ie. dissipative)
# LF -> Lax-Friedrich
# LW -> Lax-Wendroff
#
# so for example, ESLFFaceIntegral is an entropy stable face integral function
# that uses Lax-Friedrich type dissipation

#-----------------------------------------------------------------------------
# entry point functions
"""
  Calculate the face integrals in an entropy conservative manner for a given
  interface.  Unlike standard face integrals, this requires data from
  the entirety of both elements, not just data interpolated to the face

  resL and resR are updated with the results of the computation for the 
  left and right elements, respectively.

  Note that nrm_xy must contains the normal vector in x-y space at the
  face nodes.

  Aliasing restrictions: none, although its unclear what the meaning of this
                         function would be if resL and resR alias

  Performance note: the version in the tests is the same speed as this one
                    for p=1 Omega elements and about 10% faster for 
                    p=4 elements, but would not be able to take advantage of 
                    the sparsity of R for SBP Gamma elements
"""
function calcECFaceIntegral{Tdim, Tsol, Tres, Tmsh}(
                             params::AbstractParamType{Tdim}, 
                             sbpface::AbstractFace, 
                             iface::Interface,
                             qL::AbstractMatrix{Tsol}, 
                             qR::AbstractMatrix{Tsol}, 
                             aux_vars::AbstractMatrix{Tres}, 
                             nrm_xy::AbstractMatrix{Tmsh},
                             functor::FluxType, 
                             resL::AbstractMatrix{Tres}, 
                             resR::AbstractMatrix{Tres})


#  Flux_tmp = params.flux_vals1
  fluxD = params.flux_valsD
  numDofPerNode = size(fluxD, 1)
#  numDofPerNode = length(Flux_tmp)
#  nrm = params.nrm

  nrmD = params.nrmD
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
#        for d=1:Tdim
#          nrm[d] = 0
#        end
#        fill!(nrm, 0.0)
#        nrm[dim] = 1

        # accumulate entry p_i, p_j of E
        Eij = zero(Tres)  # should be Tres
        @simd for k = 1:sbpface.numnodes
          # the computation of nrm_k could be moved outside i,j loops and saved
          # in an array of size [3, sbp.numnodes]
          nrm_k = nrm_xy[dim, k]
          kR = sbpface.nbrperm[k, iface.orient]
          Eij += sbpface.interp[i,k]*sbpface.interp[j,kR]*sbpface.wface[k]*nrm_k
        end  # end loop k
 
       
        @simd for p=1:numDofPerNode
          resL[p, p_i] -= Eij*fluxD[p, dim]
          resR[p, p_j] += Eij*fluxD[p, dim]
        end

      end  # end loop dim
    end  # end loop j
  end  # end loop i


  return nothing
end


"""
  Calculate the face integral in an entropy stable manner using Lax-Friedrich
  type dissipation.  
  This uses calcECFaceIntegral and calcLFEntropyPenaltyIntegral internally, 
  see those functions for details.
"""
function calcESLFFaceIntegral{Tdim, Tsol, Tres, Tmsh}(
                             params::AbstractParamType{Tdim}, 
                             sbpface::AbstractFace, 
                             iface::Interface,
                             qL::AbstractMatrix{Tsol}, 
                             qR::AbstractMatrix{Tsol}, 
                             aux_vars::AbstractMatrix{Tres}, 
                             nrm_face::AbstractMatrix{Tmsh},
                             functor::FluxType, 
                             resL::AbstractMatrix{Tres}, 
                             resR::AbstractMatrix{Tres})

  calcECFaceIntegral(params, sbpface, iface, qL, qR, aux_vars, nrm_face, 
                     functor, resL, resR)
  calcLFEntropyPenaltyIntegral(params, sbpface, iface, qL, qR, aux_vars, 
                               nrm_face, resL, resR)

  return nothing
end

"""
  Calculate the face integral in an entropy stable manner using approximate
  Lax-Wendroff type dissipation.  
  This uses calcECFaceIntegral and calcLWEntropyPenaltyIntegral internally, 
  see those functions for details.
"""
function calcESLWFaceIntegral{Tdim, Tsol, Tres, Tmsh}(
                             params::AbstractParamType{Tdim}, 
                             sbpface::AbstractFace, 
                             iface::Interface,
                             qL::AbstractMatrix{Tsol}, 
                             qR::AbstractMatrix{Tsol}, 
                             aux_vars::AbstractMatrix{Tres}, 
                             nrm_face::AbstractMatrix{Tmsh},  # dxidx or nrm
                             functor::FluxType, 
                             resL::AbstractMatrix{Tres}, 
                             resR::AbstractMatrix{Tres})

  calcECFaceIntegral(params, sbpface, iface, qL, qR, aux_vars, nrm_face, 
                     functor, resL, resR)
  calcLWEntropyPenaltyIntegral(params, sbpface, iface, qL, qR, aux_vars, 
                               nrm_face, resL, resR)

  return nothing
end

"""
  Calculate the face integral in an entropy stable manner using
  Lax-Wendroff type dissipation.  
  This uses calcECFaceIntegral and calcLW2EntropyPenaltyIntegral internally, 
  see those functions for details.
"""
function calcESLW2FaceIntegral{Tdim, Tsol, Tres, Tmsh}(
                             params::AbstractParamType{Tdim}, 
                             sbpface::AbstractFace, 
                             iface::Interface,
                             qL::AbstractMatrix{Tsol}, 
                             qR::AbstractMatrix{Tsol}, 
                             aux_vars::AbstractMatrix{Tres}, 
                             nrm_face::AbstractMatrix{Tmsh}, # dxidx or nrm
                             functor::FluxType, 
                             resL::AbstractMatrix{Tres}, 
                             resR::AbstractMatrix{Tres})

  calcECFaceIntegral(params, sbpface, iface, qL, qR, aux_vars, nrm_face, 
                     functor, resL, resR)
  calcLW2EntropyPenaltyIntegral(params, sbpface, iface, qL, qR, aux_vars, 
                                nrm_face, resL, resR)

  return nothing
end

#-----------------------------------------------------------------------------
# Internal functions that calculate the penalties

"""
  Calculate a term that provably dissipates (mathematical) entropy using a 
  Lax-Friedrich type of dissipation.  
  This
  requires data from the left and right element volume nodes, rather than
  face nodes for a regular face integral.

  Note that nrm_face must contain the scaled face normal vector in x-y space
  at the face nodes, and qL, qR, resL, and resR are the arrays for the
  entire element, not just the face.


  Aliasing restrictions: params.nrm2, params.A0, w_vals_stencil, w_vals2_stencil
"""
function calcLFEntropyPenaltyIntegral{Tdim, Tsol, Tres, Tmsh}(
             params::ParamType{Tdim, :conservative, Tsol, Tres, Tmsh},
             sbpface::AbstractFace, iface::Interface, 
             qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, 
             aux_vars::AbstractMatrix{Tres}, nrm_face::AbstractArray{Tmsh, 2},
             resL::AbstractMatrix{Tres}, resR::AbstractMatrix{Tres})

  numDofPerNode = size(qL, 1)

  # convert qL and qR to entropy variables (only the nodes that will be used)
  wL = params.w_vals_stencil
  wR = params.w_vals2_stencil
#  wL = zeros(Tsol, numDofPerNode, sbpface.stencilsize)
#  wR = zeros(wL)

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

  # convert to IR entropy variables

  # accumulate wL at the node
  wL_i = params.v_vals
  wR_i = params.v_vals2
  qL_i = params.q_vals
  qR_i = params.q_vals2

#  wL_i = zeros(Tsol, numDofPerNode)
#  wR_i = zeros(Tsol, numDofPerNode)
  # convert wL at the node back to qL
#  qL_i = zeros(Tsol, numDofPerNode)
#  qR_i = zeros(Tsol, numDofPerNode)
#  dir = params.nrm2
  A0 = params.A0
  fastzero!(A0)

  @simd for i=1:sbpface.numnodes  # loop over face nodes
    ni = sbpface.nbrperm[i, iface.orient]
    dir = ro_sview(nrm_face, :, i)
    fastzero!(wL_i)
    fastzero!(wR_i)
#    fastzero!(qL_i)
#    fastzero!(qR_i)

    # interpolate wL and wR to this node
    @simd for j=1:sbpface.stencilsize
      interpL = sbpface.interp[j, i]
      interpR = sbpface.interp[j, ni]

      @simd for k=1:numDofPerNode
        wL_i[k] += interpL*wL[k, j]
        wR_i[k] += interpR*wR[k, j]
#        qL_i[k] += interpL*qL[k, j]
#        qR_i[k] += interpR*qR[k, j]
      end
    end

    #TODO: write getLambdaMaxSimple and getIRA0 in terms of the entropy
    #      variables to avoid the conversion
    convertToConservativeFromIR_(params, wL_i, qL_i)
    convertToConservativeFromIR_(params, wR_i, qR_i)
    # get lambda * IRA0
    lambda_max = getLambdaMaxSimple(params, qL_i, qR_i, dir)
#    @assert lambda_max > 0
#    lambda_max *= sqrt(params.h)
    # poor mans entropy fix
#    lambda_max *= 0.1
#    lambda_max += 0.25
#    lambda_max *= 2
#    lambda_max = 3.0
    
    # compute average qL
    # also delta w (used later)
    @simd for j=1:numDofPerNode
      qL_i[j] = 0.5*(qL_i[j] + qR_i[j])
      wL_i[j] -= wR_i[j]
    end

    getIRA0(params, qL_i, A0)
#    for j=1:numDofPerNode
#      A0[j, j] = 1
#    end

    # wface[i] * lambda_max * A0 * delta w
    smallmatvec!(A0, wL_i, wR_i)
    fastscale!(wR_i, sbpface.wface[i]*lambda_max)

#    middle_term = scale(A0, sbpface.wface[i]*lambda_max)
#    println("middle_term = \n", middle_term)

    # interpolate back to volume nodes
    @simd for j=1:sbpface.stencilsize
      j_pL = sbpface.perm[j, iface.faceL]
      j_pR = sbpface.perm[j, iface.faceR]

      @simd for p=1:numDofPerNode
        res_old = resL[p, j_pL]  # DEBUGGING
        resL[p, j_pL] -= sbpface.interp[j, i]*wR_i[p]
        resR[p, j_pR] += sbpface.interp[j, ni]*wR_i[p]
      end
    end

  end  # end loop i

  return nothing
end

"""
  Calculate a term that provably dissipates (mathematical) entropy using a 
  an approximation to Lax-Wendroff type of dissipation.  
  This requires data from the left and right element volume nodes, rather than
  face nodes for a regular face integral.

  Note that nrm_face must contain the scaled normal vector in x-y space
  at the face nodes, and qL, qR, resL, and resR are the arrays for the
  entire element, not just the face.

  The approximation to Lax-Wendroff is the computation of

  for i=1:Tdim
    abs(ni*Y_i*S2_i*Lambda_i*Y_i.')
  end

  rather than computing the flux jacobian in the normal direction.

  Aliasing restrictions: from params the following fields are used:
    Y, S2, Lambda, res_vals1, res_vals2, res_vals3,  w_vals_stencil, 
    w_vals2_stencil, v_vals, v_vals2, q_vals, q_vals2
"""
function calcLWEntropyPenaltyIntegral{Tdim, Tsol, Tres, Tmsh}(
             params::ParamType{Tdim, :conservative, Tsol, Tres, Tmsh},
             sbpface::AbstractFace, iface::Interface, 
             qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, 
             aux_vars::AbstractMatrix{Tres}, nrm_face::AbstractArray{Tmsh, 2},
             resL::AbstractMatrix{Tres}, resR::AbstractMatrix{Tres})

#  println("----- entered calcEntropyLWEntropyPenaltyIntegral -----")

  numDofPerNode = size(qL, 1)

  # convert qL and qR to entropy variables (only the nodes that will be used)
  wL = params.w_vals_stencil
  wR = params.w_vals2_stencil

  Y = params.A0  # eigenvectors of flux jacobian
  S2 = params.S2  # diagonal scaling matrix squared
                       # S is defined s.t. (YS)*(YS).' = A0
  Lambda = params.Lambda  # diagonal matrix of eigenvalues
  tmp1 = params.res_vals1  # work vectors
  tmp2 = params.res_vals2
  tmp3 = params.res_vals3  # accumulate result vector

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

  # convert to IR entropy variables

  # accumulate wL at the node
  wL_i = params.v_vals
  wR_i = params.v_vals2
  qL_i = params.q_vals
  qR_i = params.q_vals2

  for i=1:sbpface.numnodes  # loop over face nodes
    ni = sbpface.nbrperm[i, iface.orient]
    fill!(wL_i, 0.0)
    fill!(wR_i, 0.0)
    fill!(tmp3, 0.0)
    # interpolate wL and wR to this node
    for j=1:sbpface.stencilsize
      interpL = sbpface.interp[j, i]
      interpR = sbpface.interp[j, ni]

      for k=1:numDofPerNode
        wL_i[k] += interpL*wL[k, j]
        wR_i[k] += interpR*wR[k, j]
      end
    end

    # need conservative variables for flux jacobian calculation
    convertToConservativeFromIR_(params, wL_i, qL_i)
    convertToConservativeFromIR_(params, wR_i, qR_i)

    for j=1:numDofPerNode
      # use flux jacobian at arithmetic average state
      qL_i[j] = 0.5*( qL_i[j] + qR_i[j])
      # put delta w into wL_i
      wL_i[j] -= wR_i[j]
    end


    # get the normal vector (scaled)

    for dim =1:Tdim
      nrm_dim = nrm_face[dim, i]

      # get the eigensystem in the current direction
      if dim == 1
        calcEvecsx(params, qL_i, Y)
        calcEvalsx(params, qL_i, Lambda)
        calcEScalingx(params, qL_i, S2)
      elseif dim == 2
        calcEvecsy(params, qL_i, Y)
        calcEvalsy(params, qL_i, Lambda)
        calcEScalingy(params, qL_i, S2)
      elseif dim == 3
        calcEvecsz(params, qL_i, Y)
        calcEvalsz(params, qL_i, Lambda)
        calcEScalingz(params, qL_i, S2)
      end

      # DEBUGGING: turn this into Lax-Friedrich
#      lambda_max = maximum(absvalue(Lambda))
#      fill!(Lambda, lambda_max)

      # compute the Lax-Wendroff term, returned in tmp2
      applyEntropyLWUpdate(Y, Lambda, S2, wL_i, absvalue(nrm_dim), tmp1, tmp2)
      # accumulate result
      for j=1:length(tmp3)
        tmp3[j] += tmp2[j]
      end
    end

    # scale by wface[i]
    for j=1:length(tmp3)
      tmp3[j] *= sbpface.wface[i]
    end

    # interpolate back to volume nodes
    for j=1:sbpface.stencilsize
      j_pL = sbpface.perm[j, iface.faceL]
      j_pR = sbpface.perm[j, iface.faceR]

      for p=1:numDofPerNode
        resL[p, j_pL] -= sbpface.interp[j, i]*tmp3[p]
        resR[p, j_pR] += sbpface.interp[j, ni]*tmp3[p]
      end
    end

  end  # end loop i

  return nothing
end

@inline function applyEntropyLWUpdate(Y::AbstractMatrix, 
           Lambda::AbstractVector, S2::AbstractVector, delta_v::AbstractVector, 
           ni::Number, tmp1::AbstractVector, tmp2::AbstractVector)
# this is the computation kernel Lax-Wendroff entropy dissipation
# the result is returned in tmp2

  # multiply delta_v by Y.'
  smallmatTvec!(Y, delta_v, tmp1)
  # multiply by diagonal terms, normal vector component
  for i=1:length(delta_v)
    tmp1[i] *= ni*S2[i]*absvalue(Lambda[i])
  end
  # multiply by Y
  smallmatvec!(Y, tmp1, tmp2)

  return nothing
end

"""
  Calculate a term that provably dissipates (mathematical) entropy using a 
  Lax-Wendroff type of dissipation.  
  This requires data from the left and right element volume nodes, rather than
  face nodes for a regular face integral.

  Note nrm_face must contain the scaled normal vector in x-y space
  at the face nodes, and qL, qR, resL, and resR are the arrays for the
  entire element, not just the face.

  Implementation Detail:
    Because the scaling does not exist in arbitrary directions for 3D, 
    the function projects q into n-t coordinates, computes the
    eigendecomposition there, and then rotates back

  Aliasing restrictions: from params the following fields are used:
    Y, S2, Lambda, res_vals1, res_vals2,  w_vals_stencil, 
    w_vals2_stencil, v_vals, v_vals2, q_vals, q_vals2, nrm2, P

"""
function calcLW2EntropyPenaltyIntegral{Tdim, Tsol, Tres, Tmsh}(
             params::ParamType{Tdim, :conservative, Tsol, Tres, Tmsh},
             sbpface::AbstractFace, iface::Interface, 
             qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, 
             aux_vars::AbstractMatrix{Tres}, nrm_face::AbstractArray{Tmsh, 2},
             resL::AbstractMatrix{Tres}, resR::AbstractMatrix{Tres})

#  println("----- entered calcLW2EntropyPenaltyIntegral -----")
  numDofPerNode = size(qL, 1)

  # convert qL and qR to entropy variables (only the nodes that will be used)
  wL = params.w_vals_stencil
  wR = params.w_vals2_stencil

  Y = params.A0  # eigenvectors of flux jacobian
  S2 = params.S2  # diagonal scaling matrix squared
                       # S is defined s.t. (YS)*(YS).' = A0
  Lambda = params.Lambda  # diagonal matrix of eigenvalues
  tmp1 = params.res_vals1  # work vectors
  tmp2 = params.res_vals2

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

  # convert to IR entropy variables

  # accumulate wL at the node
  wL_i = params.v_vals
  wR_i = params.v_vals2
  qL_i = params.q_vals
  qR_i = params.q_vals2
  nrm = params.nrm2
  P = params.P  # projection matrix

  for i=1:sbpface.numnodes  # loop over face nodes
    ni = sbpface.nbrperm[i, iface.orient]
    fill!(wL_i, 0.0)
    fill!(wR_i, 0.0)
    # interpolate wL and wR to this node
    for j=1:sbpface.stencilsize
      interpL = sbpface.interp[j, i]
      interpR = sbpface.interp[j, ni]

      for k=1:numDofPerNode
        wL_i[k] += interpL*wL[k, j]
        wR_i[k] += interpR*wR[k, j]
      end
    end

    # need conservative variables for flux jacobian calculation
    convertToConservativeFromIR_(params, wL_i, qL_i)
    convertToConservativeFromIR_(params, wR_i, qR_i)

    for j=1:numDofPerNode
      # use flux jacobian at arithmetic average state
      qL_i[j] = 0.5*( qL_i[j] + qR_i[j])
      # put delta w into wL_i
      wL_i[j] -= wR_i[j]
    end


    # get the normal vector (scaled)
    for dim=1:Tdim
      nrm[dim] = nrm_face[dim, i]
    end
    # normalize direction vector
    len_fac = calcLength(params, nrm)
    for dim=1:Tdim
      nrm[dim] = nrm[dim]/len_fac
    end

    # project q into n-t coordinate system
    getProjectionMatrix(params, nrm, P)
    projectToNT(params, P, qL_i, qR_i)  # qR_i is qprime

    # get eigensystem in the normal direction, which is equivalent to
    # the x direction now that q has been rotated
    calcEvecsx(params, qR_i, Y)
    calcEvalsx(params, qR_i, Lambda)
    calcEScalingx(params, qR_i, S2)

    calcEntropyFix(params, Lambda)
 
    #DEBUGGING: make this into LF
#    fill!(Lambda, maximum(abs(Lambda)))
#    lambda_max_scaled = len_fac*maximum(abs(Lambda))

    # compute LF term in n-t coordinates, then rotate back to x-y
    projectToNT(params, P, wL_i, tmp1)
    smallmatTvec!(Y, tmp1, tmp2)
    # multiply by diagonal Lambda and S2, also include the scalar
    # wface and len_fac components
    for j=1:length(tmp2)
      tmp2[j] *= len_fac*sbpface.wface[i]*absvalue(Lambda[j])*S2[j]
    end
    smallmatvec!(Y, tmp2, tmp1)
    projectToXY(params, P, tmp1, tmp2)

    # interpolate back to volume nodes
    for j=1:sbpface.stencilsize
      j_pL = sbpface.perm[j, iface.faceL]
      j_pR = sbpface.perm[j, iface.faceR]

      for p=1:numDofPerNode
        res_old = resL[p, j_pL]  # DEBUGGING
        resL[p, j_pL] -= sbpface.interp[j, i]*tmp2[p]
        resR[p, j_pR] += sbpface.interp[j, ni]*tmp2[p]
      end
    end

  end  # end loop i

  return nothing
end

"""
  This function modifies the eigenvalues of the euler flux jacobian such
  that if any value is zero, a little dissipation is still added.  The
  absolute values of the eigenvalues modified eigenvalues are calculated.

  Methods are available for 2 and 3 dimensions

  This function depends on the ordering of the eigenvalues produced by
  calcEvals.

  Inputs:
    params: ParamType, used to dispatch to 2 or 3D method

  Inputs/Outputs:
    Lambda: vector of eigenvalues to be modified

  Aliasing restrictions: none
"""
function calcEntropyFix(params::ParamType{2}, Lambda::AbstractVector)
  
  # entropy fix parameters
  sat_Vn = 0.025
  sat_Vl = 0.05


  # this is dependent on the ordering of the eigenvalues produced
  # by calcEvals
  lambda3 = Lambda[2]  # Un
  lambda4 = Lambda[3]  # Un + a
  lambda5 = Lambda[4]  # Un - a


  # if any eigenvalue is zero, introduce dissipation that is a small
  # fraction of the maximum eigenvalue
  rhoA = max(absvalue(lambda4), absvalue(lambda5))  # absvalue(Un) + a
  lambda3 = max( absvalue(lambda3), sat_Vl*rhoA)
  lambda4 = max( absvalue(lambda4), sat_Vn*rhoA)
  lambda5 = max( absvalue(lambda5), sat_Vn*rhoA)

  Lambda[1] = lambda3
  Lambda[2] = lambda3
  Lambda[3] = lambda4
  Lambda[4] = lambda5
  
  return nothing
end

function calcEntropyFix(params::ParamType{3}, Lambda::AbstractVector)
  
  # entropy fix parameters
  sat_Vn = 0.025
  sat_Vl = 0.05


  # this is dependent on the ordering of the eigenvalues produced
  # by calcEvals
  lambda3 = Lambda[3]  # Un
  lambda4 = Lambda[4]  # Un + a
  lambda5 = Lambda[5]  # Un - a


  # if any eigenvalue is zero, introduce dissipation that is a small
  # fraction of the maximum eigenvalue
  rhoA = max(absvalue(lambda4), absvalue(lambda5))  # absvalue(Un) + a
  lambda3 = max( absvalue(lambda3), sat_Vl*rhoA)
  lambda4 = max( absvalue(lambda4), sat_Vn*rhoA)
  lambda5 = max( absvalue(lambda5), sat_Vn*rhoA)

  Lambda[1] = lambda3
  Lambda[2] = lambda3
  Lambda[3] = lambda3
  Lambda[4] = lambda4
  Lambda[5] = lambda5
  
  return nothing
end


#-----------------------------------------------------------------------------
# do the functor song and dance


"""
  Entropy conservative term only
"""
type ECFaceIntegral <: FaceElementIntegralType
end

function call{Tsol, Tres, Tmsh, Tdim}(obj::ECFaceIntegral, 
              params::AbstractParamType{Tdim}, 
              sbpface::AbstractFace, iface::Interface,
              qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, 
              aux_vars::AbstractMatrix{Tres}, nrm_face::AbstractMatrix{Tmsh},
              functor::FluxType, 
              resL::AbstractMatrix{Tres}, resR::AbstractMatrix{Tres})


  calcECFaceIntegral(params, sbpface, iface, qL, qR, aux_vars, nrm_face, 
                      functor, resL, resR)

end


"""
  Entropy conservative integral + Lax-Friedrich penalty
"""
type ESLFFaceIntegral <: FaceElementIntegralType
end

@inline function call{Tsol, Tres, Tmsh, Tdim}(obj::ESLFFaceIntegral, 
              params::AbstractParamType{Tdim}, 
              sbpface::AbstractFace, iface::Interface,
              qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, 
              aux_vars::AbstractMatrix{Tres}, nrm_face::AbstractMatrix{Tmsh},
              functor::FluxType, 
              resL::AbstractMatrix{Tres}, resR::AbstractMatrix{Tres})


  calcESLFFaceIntegral(params, sbpface, iface, qL, qR, aux_vars, nrm_face, functor, resL, resR)

end

"""
  Lax-Friedrich entropy penalty term only
"""
type ELFPenaltyFaceIntegral <: FaceElementIntegralType
end

function call{Tsol, Tres, Tmsh, Tdim}(obj::ELFPenaltyFaceIntegral, 
              params::AbstractParamType{Tdim}, 
              sbpface::AbstractFace, iface::Interface,
              qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, 
              aux_vars::AbstractMatrix{Tres}, nrm_face::AbstractMatrix{Tmsh},
              functor::FluxType, 
              resL::AbstractMatrix{Tres}, resR::AbstractMatrix{Tres})


  calcLFEntropyPenaltyIntegral(params, sbpface, iface, qL, qR, aux_vars, nrm_face, resL, resR)

end

"""
  Entropy conservative integral + approximate Lax-Wendroff penalty
"""
type ESLWFaceIntegral <: FaceElementIntegralType
end

function call{Tsol, Tres, Tmsh, Tdim}(obj::ESLWFaceIntegral, 
              params::AbstractParamType{Tdim}, 
              sbpface::AbstractFace, iface::Interface,
              qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, 
              aux_vars::AbstractMatrix{Tres}, nrm_face::AbstractMatrix{Tmsh},
              functor::FluxType, 
              resL::AbstractMatrix{Tres}, resR::AbstractMatrix{Tres})


  calcESLWFaceIntegral(params, sbpface, iface, qL, qR, aux_vars, nrm_face, functor, resL, resR)

end

"""
  Approximate Lax-Wendroff entropy penalty term only
"""
type ELWPenaltyFaceIntegral <: FaceElementIntegralType
end

function call{Tsol, Tres, Tmsh, Tdim}(obj::ELWPenaltyFaceIntegral, 
              params::AbstractParamType{Tdim}, 
              sbpface::AbstractFace, iface::Interface,
              qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, 
              aux_vars::AbstractMatrix{Tres}, nrm_face::AbstractMatrix{Tmsh},
              functor::FluxType, 
              resL::AbstractMatrix{Tres}, resR::AbstractMatrix{Tres})


  calcLWEntropyPenaltyIntegral(params, sbpface, iface, qL, qR, aux_vars, nrm_face, resL, resR)

end

"""
  Entropy conservative integral + Lax-Wendroff penalty
"""
type ESLW2FaceIntegral <: FaceElementIntegralType
end

function call{Tsol, Tres, Tmsh, Tdim}(obj::ESLW2FaceIntegral, 
              params::AbstractParamType{Tdim}, 
              sbpface::AbstractFace, iface::Interface,
              qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, 
              aux_vars::AbstractMatrix{Tres}, nrm_face::AbstractMatrix{Tmsh},
              functor::FluxType, 
              resL::AbstractMatrix{Tres}, resR::AbstractMatrix{Tres})

  calcESLW2FaceIntegral(params, sbpface, iface, qL, qR, aux_vars, nrm_face, functor, resL, resR)

end

"""
  Lax-Wendroff entropy penalty term only
"""
type ELW2PenaltyFaceIntegral <: FaceElementIntegralType
end

function call{Tsol, Tres, Tmsh, Tdim}(obj::ELW2PenaltyFaceIntegral, 
              params::AbstractParamType{Tdim}, 
              sbpface::AbstractFace, iface::Interface,
              qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, 
              aux_vars::AbstractMatrix{Tres}, nrm_face::AbstractMatrix{Tmsh},
              functor::FluxType, 
              resL::AbstractMatrix{Tres}, resR::AbstractMatrix{Tres})


  calcLW2EntropyPenaltyIntegral(params, sbpface, iface, qL, qR, aux_vars, nrm_face, resL, resR)

end



global const FaceElementDict = Dict{ASCIIString, FaceElementIntegralType}(
"ECFaceIntegral" => ECFaceIntegral(),
"ESLFFaceIntegral" => ESLFFaceIntegral(),
"ELFPenaltyFaceIntegral" => ELFPenaltyFaceIntegral(),
"ESLWFaceIntegral" => ESLWFaceIntegral(),
"ELWPenaltyFaceIntegral" => ELWPenaltyFaceIntegral(),
"ESLW2FaceIntegral" => ESLW2FaceIntegral(),
"ELW2PenaltyFaceIntegral" => ELW2PenaltyFaceIntegral(),


)

"""
  Populates the field(s) of the EulerData object with
  [`FaceElementIntegralType`](@ref) functors as specified by the options
  dictionary

  **Inputs**

   * mesh: an AbstractMesh
   * sbp: an SBP operator
   * opts: the options dictionary

  **Inputs/Outputs**

   * eqn: the EulerData object
"""
function getFaceElementFunctors(mesh, sbp, eqn::AbstractEulerData, opts)

  eqn.face_element_integral_func = FaceElementDict[opts["FaceElementIntegral_name"]]
  return nothing
end
