# shock capturing using the SBP-BR2 discretization of the second derivative term

"""
  Applies [`SBPParabolic`](@ref) shock capturing.
"""
function applyShockCapturing(mesh::AbstractMesh, sbp::AbstractOperator,
                             eqn::EulerData, opts,
                             capture::SBPParabolicSC{Tsol, Tres},
                             shockmesh::ShockedElements) where {Tsol, Tres}


  computeGradW(mesh, sbp, eqn, opts, capture, shockmesh,
               capture.convert_entropy, capture.diffusion)

  computeVolumeTerm(mesh, sbp, eqn, opts, capture, shockmesh)

  computeFaceTerm(mesh, sbp, eqn, opts, capture, shockmesh, capture.diffusion,
                  capture.penalty)

  computeBoundaryTerm(mesh, sbp, eqn, opts, capture, shockmesh)


  return nothing
end

"""
  Computes:

  [ grad_x q  =  [ lambda_xx lambda_xy  [ Dx * w
    grad_y q]      lambda_yx lambda_yy]   Dy * w]

  and stores it in `capture.grad_w`.  Note that lambda = 0 for elements that
  do not have shocks in them, so the `grad_q` is set to zero there.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * capture: [`SBPParabolicSC`](@ref)
   * shockmesh: the `ShockedElements`
   * convert_entropy: function that converts conservative variables to
                      entropy variables.  Signature must be
      `convert_entropy(params::ParamType, q::AbstractVector, w::AbstractVector)`
   * diffusion: an [`AbstractDiffusion`](@ref)
"""
function computeGradW(mesh, sbp, eqn, opts, capture::SBPParabolicSC{Tsol, Tres},
                      shockmesh::ShockedElements, convert_entropy,
                      diffusion::AbstractDiffusion
                     ) where {Tsol, Tres}

  wxi = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  grad_w = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  @simd for i=1:shockmesh.numShock
    i_full = shockmesh.elnums_all[i]
    @simd for j=1:mesh.numNodesPerElement
      q_j = ro_sview(eqn.q, :, j, i_full)
      w_j = sview(capture.w_el, :, j, i)

      # convert to entropy variables
      convert_entropy(eqn.params, q_j, w_j)
    end

    # apply D operator
    w_i = ro_sview(capture.w_el, :, :, i)
    dxidx_i = ro_sview(mesh.dxidx, :, :, :, i_full)
    jac_i = ro_sview(mesh.jac, :, i_full)
    fill!(grad_w, 0)
    applyDx(sbp, w_i, dxidx_i, jac_i, wxi, grad_w)

    # apply diffusion tensor
    lambda_gradq_i = sview(capture.grad_w, :, :, :, i)
    applyDiffusionTensor(diffusion, w_i,  i, grad_w, lambda_gradq_i)
  end

  # the diffusion is zero in the neighboring elements, so convert to entropy
  # but zero out grad_w
  @simd for i=(shockmesh.numShock+1):shockmesh.numEl
    i_full = shockmesh.elnums_all[i]
    @simd for j=1:mesh.numNodesPerElement
      q_j = ro_sview(eqn.q, :, j, i_full)
      w_j = sview(capture.w_el, :, j, i)
      convert_entropy(eqn.params, q_j, w_j)
    end

    gradw_i = sview(capture.grad_w, :, :, :, i)
    fill!(gradw_i, 0)
  end

  return nothing
end


"""
  Computes the volume terms, using the intermediate variable calcualted by
  [`computeGradW`](@ref).  The residual is updated with

  [Qx Qy] * [Lambda] * [Dx * w
                        Dy * w]

  Note that this used Qx, not Qx^T, so the term has not been integrated by
  parts.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * capture: [`SBPParbaolicSC`](@ref)
   * shockmesh
"""
function computeVolumeTerm(mesh, sbp, eqn, opts,
                           capture::SBPParabolicSC{Tsol, Tres},
                           shockmesh::ShockedElements) where {Tsol, Tres}

  # computeGradW computes Lambda * D * w, so all that remains to do is
  # compute Qx * grad_q_x
  # Note that this term is not entropy stable by itself, because Qx was
  # not replaced by -Qx^T + Ex.  The entire discretization should be
  # entropy-stable however.
  work = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  op = SummationByParts.Subtract()
  for i=1:shockmesh.numShock
    i_full = shockmesh.elnums_all[i]

    gradq_i = ro_sview(capture.grad_w, :, :, :, i)
    dxidx_i = ro_sview(mesh.dxidx, :, :, :, i_full)
    res_i = sview(eqn.res, :, :, i_full)
    applyQx(sbp, gradq_i, dxidx_i, work, res_i)
    #applyQxTransposed(sbp, gradq_i, dxidx_i, work, res_i, op)
  end

  return nothing
end


"""
  Computes the interior face term for any [`SBPParabolic`](@ref) shock
  capturing scheme.

  **Inputs**

   * mesh
   * sbp
   * eqn: `eqn.res` is updated with the result
   * opts
   * shockmesh
   * diffusion: the [`AbstractDiffusion`](@ref) to use (must be consistent
                with the one passed to [`computeGradW`](@ref).
   * penalty: [`AbstractDiffusionPenalty`](@ref) specifying which scheme to use

"""
function computeFaceTerm(mesh, sbp, eqn, opts,
                      capture::SBPParabolicSC{Tsol, Tres},
                      shockmesh::ShockedElements, diffusion::AbstractDiffusion,
                      penalty::AbstractDiffusionPenalty) where {Tsol, Tres}

  delta_w = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  theta = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)

  t1L = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  t1R = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  t2L = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  t2R = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  op = SummationByParts.Subtract()

  for i=1:shockmesh.numInterfaces
    iface_red = shockmesh.ifaces[i].iface
    iface_idx = shockmesh.ifaces[i].idx_orig
    elnumL = shockmesh.elnums_all[iface_red.elementL]
    elnumR = shockmesh.elnums_all[iface_red.elementR]

    wL = ro_sview(capture.w_el, :, :, iface_red.elementL)
    wR = ro_sview(capture.w_el, :, :, iface_red.elementR)
    gradwL = ro_sview(capture.grad_w, :, :, :, iface_red.elementL)
    gradwR = ro_sview(capture.grad_w, :, :, :, iface_red.elementR)
    nrm_face = ro_sview(mesh.nrm_face, :, :, iface_idx)

    # get data needed for next steps
    nrm_face = ro_sview(mesh.nrm_face, :, :, iface_idx)
    dxidxL = ro_sview(mesh.dxidx, :, :, :, elnumL)
    dxidxR = ro_sview(mesh.dxidx, :, :, :, elnumR)
    jacL = ro_sview(mesh.jac, :, elnumL)
    jacR = ro_sview(mesh.jac, :, elnumR)
    alphas = ro_sview(capture.alpha, :, i)
    resL = sview(eqn.res, :, :, elnumL)
    resR = sview(eqn.res, :, :, elnumR)



    # compute delta w tilde and theta_bar = Dgk w_k + Dgn w_n
    getFaceVariables(capture, mesh.sbpface, iface_red, wL, wR, gradwL, gradwR,
                     nrm_face, delta_w, theta)

    # apply the penalty coefficient matrix
    applyPenalty(penalty, sbp, mesh.sbpface, diffusion, iface_red, delta_w, theta,
                 wL, wR, nrm_face, alphas, jacL, jacR, t1L, t1R, t2L, t2R)

    # apply Rgk^T, Rgn^T, Dgk^T, Dgn^T
    # need to apply R^T * t1, not R^T * B * t1, so
    # interiorFaceIntegrate won't work.  Use the reverse mode instead
    for j=1:mesh.numNodesPerFace
      for k=1:mesh.numDofPerNode
        t1L[k, j] = -t1L[k, j]  # the SAT has a minus sign in front of it
        t1R[k, j] = -t1R[k, j]
      end
    end
    interiorFaceInterpolate_rev!(mesh.sbpface, iface_red, resL, resR, t1L, t1R)

    # apply Dgk^T and Dgn^T
    applyDgkTranspose(capture, sbp, mesh.sbpface, iface_red, diffusion, t2L, t2R,
                      wL, wR, nrm_face, dxidxL, dxidxR, jacL, jacR, resL, resR,
                      op)

  end  # end loop i

  return nothing
end


"""
  Computes a Neumann boundary condition Lambda * grad w = 0 for all faces
  on the boundary.  This term is required to make the entire shock capturing
  scheme entropy stable.

  This term does not depend on what penalty scheme is used, so it does not
  take an `AbstractPenalty` object.  It uses `capture.grad_w`, so it does
  not need an `AbstractDiffusion` object either.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * capture: an [`SBPParabolic`](@ref) object
   * shockmesh
"""
function computeBoundaryTerm(mesh, sbp, eqn, opts,
                      capture::SBPParabolicSC{Tsol, Tres},
                      shockmesh::ShockedElements,
                      ) where {Tsol, Tres}

  # for shock capturing, apply the Neumann boundary condition
  # Lambda * grad_w = 0.  The resulting term is -R^T * B * Dgk * u

  temp_face = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  temp2_face = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  op = SummationByParts.Subtract()
  for i=1:shockmesh.numBoundaryFaces
    bndry_i = shockmesh.bndryfaces[i].bndry
    idx_orig = shockmesh.bndryfaces[i].idx_orig
    elnum_orig = shockmesh.elnums_all[bndry_i.element]


    nrm_i = ro_sview(mesh.nrm_bndry, :, :, idx_orig)
    res_i = sview(eqn.res, :, :, elnum_orig)

    # grad_w already has Lambda * [Dx; Dy]*u, now apply R and N
    fill!(temp2_face, 0)
    for d1=1:mesh.dim
      gradw_d = ro_sview(capture.grad_w, :, :, d1, bndry_i.element)
      fill!(temp_face, 0)
      boundaryFaceInterpolate!(mesh.sbpface, bndry_i.face, gradw_d, temp_face)

      @simd for j=1:mesh.numNodesPerFace
        @simd for k=1:mesh.numDofPerNode
          temp2_face[k, j] += nrm_i[d1, j]*temp_face[k, j]
        end
      end
    end  # end d1

    # apply R^T B
    boundaryFaceIntegrate!(mesh.sbpface, bndry_i.face, temp2_face, res_i, op)

  end  # end i

  return nothing
end





"""
  Compute delta_w = Rgk * wL - Rgn - wR and theta = Dgk * wL + Dgn * wR
  at the face.

  **Inputs**

   * capture: an [`SBPParabolic`](@ref) object
   * sbpface
   * iface_red: an `Interface` object from the `shockmesh`
   * wL: the entropy variables for the left element, `numDofPerNode` x
         `numNodesPerElement`
   * wR: the entropy variables for the right element, same size as `wL`
   * gradwL: Lambda * [Dx, Dy]*wL (as computed by [`computeGradW`](@ref),
             `numDofPerNode` x `numNodesPerElement` x `dim`
   * gradwR: same as `gradwR`, but for the right element
   * nrm_face: the (scaled) normal vectors for the face, `dim` x
               `numNodesPerFace`

  **Inputs/Outputs**

   * delta_w: `numDofPerNode` x `numNodesPerFace`, overwritten
   * theta: `numDofPerNode` x `numNodesPerFace`, overwritten
"""
function getFaceVariables(capture::SBPParabolicSC{Tsol, Tres},
                          sbpface::AbstractFace, iface_red::Interface,
                          wL::AbstractMatrix, wR::AbstractMatrix,
                          gradwL::Abstract3DArray, gradwR::Abstract3DArray,
                          nrm_face::AbstractMatrix,
                          delta_w::AbstractMatrix, theta::AbstractArray
                         ) where {Tsol, Tres}

  numDofPerNode, numNodesPerElement = size(wL)
  numNodesPerFace = size(delta_w, 2)
  dim = size(gradwL, 3)

  @unpack capture w_faceL w_faceR grad_faceL grad_faceR

  interiorFaceInterpolate!(sbpface, iface_red, wL, wR, w_faceL, w_faceR)

  @simd for j=1:numNodesPerFace
    @simd for k=1:numDofPerNode
      delta_w[k, j] = w_faceL[k, j] - w_faceR[k, j]
    end
  end

  fill!(theta, 0)
  for d=1:dim
    gradwL_d = ro_sview(gradwL, :, :, d)
    gradwR_d = ro_sview(gradwR, :, :, d)
    interiorFaceInterpolate!(sbpface, iface_red, gradwL_d, gradwR_d,
                             grad_faceL, grad_faceR)

    @simd for j=1:numNodesPerFace
      @simd for k=1:numDofPerNode
        theta[k, j] += nrm_face[d, j]*(grad_faceL[k, j] - grad_faceR[k, j])
      end
    end
  end  # end d

  return nothing
end




"""
  This function applies Dgk^T and Dgn^T.

  **Inputs**

   * capture: an [`SBPParabolic`](@ref) object
   * sbp
   * sbpface
   * iface: an `Interface` from `shockmesh`
   * diffusion: an [`AbstractDiffusion`](@ref)
   * t2L: array to multiply Dgk^T with, `numDofPerNode` x `numNodesPerFace`
   * t2R: array to multiply Dgn^T with, `numDofPerNode` x `numNodesPerFace`
   * wL: entropy variables for left element, `numDofPerNode` x
         `numNodesPerElement`
   * wR: entropy variables for the right element, same size as `wL`
   * nrm_face: (scaled) normal vectors at the face, `dim` x `numNodesPerFace`
   * dxidxL: (scaled) mapping jacobian at left element, `dim` x `dim` x 
             `numNodesPerElement`
   * dxidxR: (scaled) mapping jacobian at the right element, same size as
             `dxidxL`
   * jacL: mapping jacobian determinant for left element, `numNodesPerElement`
   * jacR: mapping jacobian determinant for the right element,
           `numNodesPerElement`
   * op: a `SummationByParts.UnaryFunctor`, determins if the output arrays are
         added to or subtracted from

  **Inputs/Outputs**

   * resL: array to update with Dgk^T * t2L, `numDofPerNode` x
           `numNodesPerElement`
   * resR: array to update with Dgn^T * t2R, same size as `resL`

"""
function applyDgkTranspose(capture::SBPParabolicSC{Tsol, Tres}, sbp,
                           sbpface, iface::Interface,
                           diffusion::AbstractDiffusion,
                           t2L::AbstractMatrix, t2R::AbstractMatrix,
                           wL::AbstractMatrix, wR::AbstractMatrix,
                           nrm_face::AbstractMatrix,
                           dxidxL::Abstract3DArray, dxidxR::Abstract3DArray,
                           jacL::AbstractVector, jacR::AbstractVector,
                           resL::AbstractMatrix, resR::AbstractMatrix,
                           op::SummationByParts.UnaryFunctor=SummationByParts.Add()) where {Tsol, Tres}

  dim, numNodesPerFace = size(nrm_face)
  numNodesPerElement = size(resL, 2)
  numDofPerNode = size(wL, 1)

  @unpack capture temp1L temp1R temp2L temp2R temp3L temp3R work
  fill!(temp2L, 0); fill!(temp2R, 0)

  # apply N and R^T
  @simd for d=1:dim
    @simd for j=1:numNodesPerFace
      @simd for k=1:numDofPerNode
        temp1L[k, j] =  nrm_face[d, j]*t2L[k, j]
        temp1R[k, j] = -nrm_face[d, j]*t2R[k, j]
      end
    end

    tmp2L = sview(temp2L, :, :, d); tmp2R = sview(temp2R, :, :, d)
    interiorFaceInterpolate_rev!(sbpface, iface, tmp2L, tmp2R, temp1L, temp1R)
  end


  # multiply by D^T Lambda
  applyDiffusionTensor(diffusion, wL, iface.elementL, temp2L, temp3L)
  applyDiffusionTensor(diffusion, wR, iface.elementR, temp2R, temp3R)

  # saving temp3 to a mesh-wide array and then applying Dx^T would save
  # a lot of flops.
  applyDxTransposed(sbp, temp3L, dxidxL, jacL, work, resL, op)
  applyDxTransposed(sbp, temp3R, dxidxR, jacR, work, resR, op)
  
  return nothing
end


"""
  Applies the penalty for the SBP-SAT generalization of the modified scheme
  of Bassi and Rebay (BR2).  More specifically, it applies the penalty matrix
  for both sides of the face at the same time:

  [res1L   = [T1 T3  [delta_w
   res2L]     T2 T4]  theta]
 
  and similarly for the right element.  Note that delta_w is negated for
  the right element, and this function performs the required negation.

  **Inputs**

   * penalty: the [`BR2Penalty`](@ref) object
   * sbp
   * sbpface
   * diffusion: the [`AbstractDiffusion`](@ref) object
   * iface: `Interface` object
   * delta_w: difference of entropy variables at the face, as compuated by
              [`getFaceVariables`](@ref)
   * theta: Dgk * wL + Dgn * wR, as compuated by `getFaceVariables`
   * wL: entropy variables for the left element, `numDofPerNode` x
         `numNodesPerElement`
   * wR: entropy variables for the right element, same size as `wR`
   * nrm_face: (scaled) normal vectors at the face, `dim` x `numNodesPerFace`
   * alphas: vector of length 2 containing alpha_gk and alpha_gn
   * jacL: mapping jacobian determinant for the left element,
           `numNodesPerElement`
   * jacR: mapping jacobian determinant for the right element, same size as
           `jacR`

  **Inputs/Outputs**
  
   * res1L: `numDofPerNode` x `numNodesPerElement`
   * res1R: same size as above
   * res2L: same size as above
   * res2R: same size as above

  All output arrays are overwritten
"""
function applyPenalty(penalty::BR2Penalty{Tsol, Tres}, sbp, sbpface,
                      diffusion::AbstractDiffusion, iface::Interface,
                      delta_w::AbstractMatrix{Tsol}, theta::AbstractMatrix{Tres},
                      wL::AbstractMatrix, wR::AbstractMatrix,
                      nrm_face::AbstractMatrix,
                      alphas::AbstractVector,
                      jacL::AbstractVector, jacR::AbstractVector,
                      res1L::AbstractMatrix, res1R::AbstractMatrix,
                      res2L::AbstractMatrix, res2R::AbstractMatrix
                     ) where {Tsol, Tres}

  fill!(res1L, 0); fill!(res1R, 0)
  fill!(res2L, 0); fill!(res2R, 0)
  numDofPerNode, numNodesPerFace = size(delta_w)
  numNodesPerElement = length(jacL)
  dim = size(nrm_face, 1)

  @unpack penalty delta_w_n qL qR t1L t1R t2L t2R
  fill!(qL, 0); fill!(qR, 0)

  #--------------------------
  # apply T1
  # multiply by normal vector, then R^T B
  #alpha_g = 1/(dim + 1)  # = 1/number of faces of a simplex
  @simd for d1=1:dim
    @simd for j=1:numNodesPerFace
      @simd for k=1:numDofPerNode
        delta_w_n[k, j] = 0.25*delta_w[k, j]*nrm_face[d1, j]
      end
    end
    
    qL_d = sview(qL, :, :, d1); qR_d = sview(qR, :, :, d1)
    interiorFaceIntegrate!(sbpface, iface, delta_w_n, qL_d, qR_d)

    # interiorFaceIntegrates -= the second output
    @simd for j=1:numNodesPerElement
      @simd for k=1:numDofPerNode
        qR_d[k, j] = -qR_d[k, j]
      end
    end
  end

  # apply Lambda matrix
  applyDiffusionTensor(diffusion, wL, iface.elementL, qL, t1L)
  applyDiffusionTensor(diffusion, wR, iface.elementR, qR, t1R)

  # apply inverse mass matrix, then apply B*Nx*R*t2L_x + B*Ny*R*t2L_y
  @simd for d1=1:dim
    @simd for j=1:numNodesPerElement
      facL = (1/alphas[1])*jacL[j]/sbp.w[j]
      facR = (1/alphas[2])*jacR[j]/sbp.w[j]
      @simd for k=1:numDofPerNode
        t1L[k, j, d1] *= facL
        t1R[k, j, d1] *= facR
      end
    end

    #TODO: make t2L 2D rather than 3D
    t1L_d = ro_sview(t1L, :, :, d1); t1R_d = ro_sview(t1R, :, :, d1)
    t2L_d = sview(t2L, :, :, d1);    t2R_d = sview(t2R, :, :, d1)
    interiorFaceInterpolate!(sbpface, iface, t1L_d, t1R_d, t2L_d, t2R_d)

    @simd for j=1:numNodesPerFace
      @simd for k=1:numDofPerNode
        val = sbpface.wface[j]*nrm_face[d1, j]*(t2L_d[k, j] + t2R_d[k, j])

        res1L[k, j] += val
        res1R[k, j] -= val  # delta_w is reversed for elementR
      end
    end

  end  # end d1

  #----------------------------
  # apply T2 and T3
  @simd for j=1:numNodesPerFace
    @simd for k=1:numDofPerNode
      val1 =  0.5*sbpface.wface[j]*theta[k, j]  # this one is the problem
      val2 = -0.5*sbpface.wface[j]*delta_w[k, j]
      res1L[k, j] += val1
      res1R[k, j] += val1
      res2L[k, j] += val2
      res2R[k, j] -= val2  # delta_w is reversed for elementR
    end
  end

  return nothing
end
