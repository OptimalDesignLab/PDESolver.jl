# shock capturing using the SBP-BR2 discretization of the second derivative term


function applyShockCapturing(mesh::AbstractMesh, sbp::AbstractOperator,
                             eqn::EulerData, opts,
                             capture::SBPParabolicSC{Tsol, Tres},
                             shockmesh::ShockedElements) where {Tsol, Tres}


  computeGradW(mesh, sbp, eqn, opts, capture, shockmesh,
               capture.convert_entropy, capture.diffusion)

  computeVolumeTerm(mesh, sbp, eqn, opts, capture, shockmesh)

  computeFaceTerm(mesh, sbp, eqn, opts, capture, shockmesh, capture.diffusion,
                  capture.penalty)

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
  for i=1:shockmesh.numShock
    i_full = shockmesh.elnums_all[i]
    for j=1:mesh.numNodesPerElement
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
  for i=(shockmesh.numShock+1):shockmesh.numEl
    i_full = shockmesh.elnums_all[i]
    for j=1:mesh.numNodesPerElement
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
  [`computeGradW`](@ref)
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
  for i=1:shockmesh.numShock
    i_full = shockmesh.elnums_all[i]

    gradq_i = ro_sview(capture.grad_w, :, :, :, i)
    dxidx_i = ro_sview(mesh.dxidx, :, :, :, i_full)
    res_i = sview(eqn.res, :, :, i_full)
    applyQx(sbp, gradq_i, dxidx_i, work, res_i)
  end

  return nothing
end


function computeFaceTerm(mesh, sbp, eqn, opts,
                      capture::SBPParabolicSC{Tsol, Tres},
                      shockmesh::ShockedElements, diffusion::AbstractDiffusion,
                      penalty::AbstractDiffusionPenalty) where {Tsol, Tres}

  delta_w = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  theta = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)

  t1 = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  t2 = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  op = SummationByParts.Subtract()

  for i=1:shockmesh.numInterfaces
    iface_red = shockmesh.ifaces[i].iface
    iface_idx = shockmesh.ifaces[i].idx_orig
    elnumL = shockmesh.elnums_all[iface_red.elementL]
    elnumR = shockmesh.elnums_all[iface_red.elementR]

    # compute delta w tilde and theta_bar = Dgk w_k + Dgn w_n
    wL = ro_sview(capture.w_el, :, :, iface_red.elementL)
    wR = ro_sview(capture.w_el, :, :, iface_red.elementR)
    gradwL = ro_sview(capture.grad_w, :, :, :, iface_red.elementL)
    gradwR = ro_sview(capture.grad_w, :, :, :, iface_red.elementR)
    nrm_face = ro_sview(mesh.nrm_face, :, :, iface_idx)

    getFaceVariables(capture, mesh.sbpface, iface_red, wL, wR, gradwL, gradwR,
                     nrm_face, delta_w, theta)


    # get data needed for next steps
    nrm_face = ro_sview(mesh.nrm_face, :, :, iface_idx)
    dxidxL = ro_sview(mesh.dxidx, :, :, :, elnumL)
    dxidxR = ro_sview(mesh.dxidx, :, :, :, elnumR)
    jacL = ro_sview(mesh.jac, :, elnumL)
    jacR = ro_sview(mesh.jac, :, elnumR)


    # apply the penalty coefficient matrix
    applyPenalty(penalty, sbp, mesh.sbpface, diffusion, iface_red, delta_w, theta,
                 wL, wR, nrm_face, jacL, jacR, t1, t2)

    # apply Rgk^T, Rgn^T, Dgk^T, Dgn^T
    resL = sview(eqn.res, :, :, elnumL)
    resR = sview(eqn.res, :, :, elnumR)

    # need to apply R^T * t1, not R^T * B * t1, so
    # interiorFaceIntegrate won't work.  Use the reverse mode instead
    interiorFaceInterpolate_rev!(mesh.sbpface, iface_red, resL, resR, t1, t1)

    # apply Dgk^T and Dgn^T
    applyDgkTranspose(capture, sbp, mesh.sbpface, iface_red, diffusion, t2,
                      wL, wR, nrm_face, dxidxL, dxidxR, jacL, jacR, resL, resR,
                      op)

  end  # end loop i

  return nothing
end


"""
  Compute delta_w and theta.
"""
function getFaceVariables(capture::SBPParabolicSC{Tsol, Tres},
                          sbpface::AbstractFace, iface_red::Interface,
                          wL::AbstractMatrix, wR::AbstractMatrix,
                          gradwL::Abstract3DArray, gradwR::Abstract3DArray,
                          nrm_face::AbstractMatrix,
                          delta_w::AbstractArray, theta::AbstractArray
                         ) where {Tsol, Tres}

  numDofPerNode, numNodesPerElement = size(wL)
  numNodesPerFace = size(delta_w, 2)
  dim = size(gradwL, 3)

  @unpack capture w_faceL w_faceR grad_faceL grad_faceR

  interiorFaceInterpolate!(sbpface, iface_red, wL, wR, w_faceL, w_faceR)

  for j=1:numNodesPerFace
    for k=1:numDofPerNode
      delta_w[k, j] = w_faceL[k, j] - w_faceR[k, j]
    end
  end

  fill!(theta, 0)
  for d=1:dim
    gradwL_d = ro_sview(gradwL, :, :, d)
    gradwR_d = ro_sview(gradwR, :, :, d)
    interiorFaceInterpolate!(sbpface, iface_red, gradwL_d, gradwR_d,
                             grad_faceL, grad_faceR)

    for j=1:numNodesPerFace
      for k=1:numDofPerNode
        theta[k, j] += nrm_face[d, j]*(grad_faceL[k, j] - grad_faceR[k, j])
      end
    end
  end  # end d

  return nothing
end




"""
"""
function applyDgkTranspose(capture::SBPParabolicSC{Tsol, Tres}, sbp,
                           sbpface, iface::Interface,
                           diffusion::AbstractDiffusion,
                           t2::AbstractMatrix,
                           wL::AbstractMatrix, wR::AbstractMatrix,
                           nrm_face::AbstractMatrix,
                           dxidxL::Abstract3DArray, dxidxR::Abstract3DArray,
                           jacL::AbstractVector, jacR::AbstractVector,
                           resL::AbstractMatrix, resR::AbstractMatrix,
                           op::SummationByParts.UnaryFunctor=SummationByParts.Add()) where {Tsol, Tres}

  dim, numNodesPerFace = size(nrm_face)
  numNodesPerElement = size(resL, 2)
  numDofPerNode = size(wL, 1)

  @unpack capture temp1 temp2L temp2R temp3L temp3R work
  fill!(temp2L, 0); fill!(temp2R, 0)

  # apply N and R^T
  for d=1:dim
    for j=1:numNodesPerFace
      for k=1:numDofPerNode
        temp1[k, j] = nrm_face[d, j]*t2[k, j]
      end
    end

    tmp2L = sview(temp2L, :, :, d); tmp2R = sview(temp2R, :, :, d)
    interiorFaceInterpolate_rev!(sbpface, iface, tmp2L, tmp2R, temp1, temp1)

    # the normal vector for elementR is negated.  Could have done this before,
    # but that would require 2 calls 2 interiorFaceInterpolate_rev
    # TODO: it would be better to add a function to SBP that gets the
    # inverse operation for a UnaryFunctor
    for j=1:numNodesPerElement
      for k=1:numDofPerNode
        tmp2R[k, j] = -tmp2R[k, j]
      end
    end
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
"""
function applyPenalty(penalty::BR2Penalty{Tsol, Tres}, sbp, sbpface,
                      diffusion::AbstractDiffusion, iface::Interface,
                      delta_w::AbstractMatrix{Tsol}, theta::AbstractMatrix{Tres},
                      wL::AbstractMatrix, wR::AbstractMatrix,
                      nrm_face::AbstractMatrix,
                      jacL::AbstractVector, jacR::AbstractVector,
                      res1::AbstractMatrix, res2::AbstractMatrix) where {Tsol, Tres}

  fill!(res1, 0); fill!(res2, 0)
  numDofPerNode, numNodesPerFace = size(delta_w)
  numNodesPerElement = length(jacL)
  dim = size(nrm_face, 1)

  @unpack penalty delta_w_n qL qR t1L t1R t2L t2R
  fill!(qL, 0); fill!(qR, 0)

  # multiply by normal vector, then R^T B
  alpha_g = 1/(dim + 1)  # = 1/number of faces of a simplex
  for d1=1:dim
    for j=1:numNodesPerFace
      for k=1:numDofPerNode
        delta_w_n[k, j] = alpha_g*delta_w[k, j]*nrm_face[d1, j]
      end
    end
    
    qL_d = sview(qL, :, :, d1); qR_d = sview(qR, :, :, d1)
    interiorFaceIntegrate!(sbpface, iface, delta_w_n, qL_d, qR_d)
    scale!(qR_d, -1)  # interiorFaceIntegrates -= the second output
  end

  # apply Lambda matrix
  applyDiffusionTensor(diffusion, wL, iface.elementL, qL, t1L)
  applyDiffusionTensor(diffusion, wR, iface.elementR, qR, t1R)

  # apply inverse mass matrix, then apply B*Nx*R*t2L_x + B*Ny*R*t2L_y
  for d1=1:dim
    for j=1:numNodesPerElement
      facL = jacL[j]/sbp.w[j]
      facR = jacR[j]/sbp.w[j]
      for k=1:numDofPerNode
        t1L[k, j, d1] *= facL
        t1R[k, j, d1] *= facR
      end
    end

    t1L_d = ro_sview(t1L, :, :, d1); t1R_d = ro_sview(t1R, :, :, d1)
    t2L_d = sview(t2L, :, :, d1);    t2R_d = sview(t2R, :, :, d1)
    interiorFaceInterpolate!(sbpface, iface, t1L_d, t1R_d, t2L_d, t2R_d)

    for j=1:numNodesPerFace
      for k=1:numDofPerNode
        res1[k, j] += sbpface.wface[j]*nrm_face[d1, j]*(t2L_d[k, j] + t2R_d[k, j])
#        res1[k, j] += sbpface.wface[j]*(t2L_d[k, j] + t2R_d[k, j])
      end
    end

  end  # end d1

  # apply T2 and T3
  for j=1:numNodesPerFace
    for k=1:numDofPerNode
      res1[k, j] += -0.5*sbpface.wface[j]*theta[k, j]
      res2[k, j] +=  0.5*sbpface.wface[j]*delta_w[k, j]
    end
  end

  return nothing
end
