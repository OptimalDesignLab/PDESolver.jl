# shock capturing using any the SBP-SAT discretization of the second derivative
# term

"""
  Applies [`SBPParabolic`](@ref) shock capturing.
"""
function calcShockCapturing(mesh::AbstractMesh, sbp::AbstractOperator,
                             eqn::EulerData, opts,
                             capture::SBPParabolicSC{Tsol, Tres},
                             shockmesh::ShockedElements) where {Tsol, Tres}

  println("initially, residual norm = ", calcNorm(eqn, eqn.res))

  computeGradW(mesh, sbp, eqn, opts, capture, shockmesh,
               capture.entropy_vars, capture.diffusion)

  computeVolumeTerm(mesh, sbp, eqn, opts, capture, shockmesh)

  println("after volume term, residual norm = ", calcNorm(eqn, eqn.res))

  computeFaceTerm(mesh, sbp, eqn, opts, capture, shockmesh, capture.diffusion,
                  capture.penalty)

  #println("after face term, residual norm = ", calcNorm(eqn, eqn.res))
  if shockmesh.isNeumann
    computeNeumannBoundaryTerm(mesh, sbp, eqn, opts, capture, shockmesh)
  else
    computeDirichletBoundaryTerm(mesh, sbp, eqn, opts, capture, shockmesh)
  end

  println("after boundary term, residual norm = ", calcNorm(eqn, eqn.res))
  computeSharedFaceTerm(mesh, sbp, eqn, opts, capture, shockmesh,
                              capture.diffusion, capture.penalty)



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
   * entropy_vars: an [`AbstractVariables`](@ref) object
   * diffusion: an [`AbstractDiffusion`](@ref)
"""
function computeGradW(mesh, sbp, eqn, opts, capture::SBPParabolicSC{Tsol, Tres},
                      shockmesh::ShockedElements,
                      entropy_vars::AbstractVariables,
                      diffusion::AbstractDiffusion,
                     ) where {Tsol, Tres}

  wxi = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  grad_w = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)

  # do local elements
  @simd for i in shockmesh.local_els
    i_full = shockmesh.elnums_all[i]
    @simd for j=1:mesh.numNodesPerElement
      q_j = ro_sview(eqn.q, :, j, i_full)
      w_j = sview(capture.w_el, :, j, i)

      # convert to entropy variables
      convertToEntropy(entropy_vars, eqn.params, q_j, w_j)
    end

    # apply D operator
    q_i = ro_sview(eqn.q, :, :, i_full)
    w_i = ro_sview(capture.w_el, :, :, i)
    coords = ro_sview(mesh.coords, :, :, i_full)
    dxidx_i = ro_sview(mesh.dxidx, :, :, :, i_full)
    jac_i = ro_sview(mesh.jac, :, i_full)
    fill!(grad_w, 0)
    applyDx(sbp, w_i, dxidx_i, jac_i, wxi, grad_w)

    # apply diffusion tensor
    lambda_gradq_i = sview(capture.grad_w, :, :, :, i)

    applyDiffusionTensor(diffusion, sbp, eqn.params, q_i, w_i, coords, dxidx_i,
                         jac_i, i, grad_w, lambda_gradq_i)
  end

  # the diffusion is zero in the neighboring elements, so convert to entropy
  # but zero out grad_w
  @simd for i in shockmesh.neighbor_els
    i_full = shockmesh.elnums_all[i]
    @simd for j=1:mesh.numNodesPerElement
      q_j = ro_sview(eqn.q, :, j, i_full)
      w_j = sview(capture.w_el, :, j, i)
      convertToEntropy(entropy_vars, eqn.params, q_j, w_j)
    end

    gradw_i = sview(capture.grad_w, :, :, :, i)
    fill!(gradw_i, 0)
  end

  # do shared elements
  for peer=1:shockmesh.npeers
    peer_full = shockmesh.peer_indices[peer]
    data = eqn.shared_data[peer_full]
    metrics = mesh.remote_metrics[peer_full]

    for i in shockmesh.shared_els[peer]
      i_full = getSharedElementIndex(shockmesh, mesh, peer, i)

      # compute entropy variables
      for j=1:mesh.numNodesPerElement
        q_j = ro_sview(data.q_recv, :, j, i_full)
        w_j = sview(capture.w_el, :, j, i)
        convertToEntropy(entropy_vars, eqn.params, q_j, w_j)
      end

      lambda_gradq_i = sview(capture.grad_w, :, :, :, i)

      # some of these have ee = 0, but we don't have the sensor so we can't
      # check
      q_i = ro_sview(data.q_recv, :, :, i_full)
      w_i = ro_sview(capture.w_el, :, :, i)
      coords_i = ro_sview(metrics.coords, :, :, i_full)
      dxidx_i = ro_sview(metrics.dxidx, :, :, :, i_full)
      jac_i = ro_sview(metrics.jac, :, i_full)
      fill!(grad_w, 0)

      applyDx(sbp, w_i, dxidx_i, jac_i, wxi, grad_w)
      applyDiffusionTensor(diffusion, sbp, eqn.params, q_i,  w_i, coords_i,
                           dxidx_i, jac_i, i, grad_w, lambda_gradq_i)

    end  # end i
  end  # end peer

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

    # get data needed for next steps
    qL = ro_sview(eqn.q, :, :, elnumL)
    qR = ro_sview(eqn.q, :, :, elnumR)
    wL = ro_sview(capture.w_el, :, :, iface_red.elementL)
    wR = ro_sview(capture.w_el, :, :, iface_red.elementR)
    gradwL = ro_sview(capture.grad_w, :, :, :, iface_red.elementL)
    gradwR = ro_sview(capture.grad_w, :, :, :, iface_red.elementR)

    coordsL = ro_sview(mesh.coords, :, :, elnumL)
    coordsR = ro_sview(mesh.coords, :, :, elnumR)
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
    applyPenalty(penalty, sbp, eqn.params, mesh.sbpface, diffusion,
                 iface_red, delta_w, theta,
                 qL, qR, wL, wR, coordsL, coordsR, nrm_face, alphas, dxidxL,
                 dxidxR, jacL, jacR, t1L, t1R, t2L, t2R)

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
    # TODO: there is an allocation here caused by the number of arguments
    #       to the function, however the allocation is context dependent.
    #       Pulling this function out of PDESolver removes the allocation.
    #       The allocation can be avoided by packing dxidxL, jacL, and resL
    #       into a tuple (and similarly for the corresponding R arguments),
    #       however the allocation does not appear to cause a significant
    #       performance problem.
    applyDgkTranspose(capture, sbp, eqn.params, mesh.sbpface, iface_red,
                      diffusion, t2L, t2R,
                      qL, qR, wL, wR, coordsL, coordsR, nrm_face, dxidxL,
                      dxidxR, jacL, jacR, resL, resR, op)

  end  # end loop i

  return nothing
end


"""
  Does the same thing as `compuateFaceTerm`, but for the shared faces, updating
  the residual on the local element only
"""
function computeSharedFaceTerm(mesh, sbp, eqn, opts,
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


  # don't care about resR for shared faces
  resR = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)
  for _peer=1:shockmesh.npeers
    # finish parallel communication for this set of faces
    peer = receiveAlpha(capture, shockmesh)

    peer_full = shockmesh.peer_indices[peer]
    metrics = mesh.remote_metrics[peer_full]
    data = eqn.shared_data[peer_full]

    for i=1:shockmesh.numSharedInterfaces[peer]
      iface_red = shockmesh.shared_interfaces[peer][i].iface
      iface_idx = shockmesh.shared_interfaces[peer][i].idx_orig
      elnumL = shockmesh.elnums_all[iface_red.elementL]
      elnumR = getSharedElementIndex(shockmesh, mesh, peer, iface_red.elementR)

      # get data needed for next steps
      qL = ro_sview(eqn.q, :, :, elnumL)
      qR = ro_sview(data.q_recv, :, :, elnumR)
      wL = ro_sview(capture.w_el, :, :, iface_red.elementL)
      wR = ro_sview(capture.w_el, :, :, iface_red.elementR)
      gradwL = ro_sview(capture.grad_w, :, :, :, iface_red.elementL)
      gradwR = ro_sview(capture.grad_w, :, :, :, iface_red.elementR)

      coordsL = ro_sview(mesh.coords   , :, :, elnumL)
      coordsR = ro_sview(metrics.coords, :, :, elnumR)
      nrm_face = ro_sview(mesh.nrm_sharedface[peer_full], :, :, iface_idx)
      dxidxL = ro_sview(mesh.dxidx,    :, :, :, elnumL)
      dxidxR = ro_sview(metrics.dxidx, :, :, :, elnumR)
      jacL = ro_sview(mesh.jac,    :, elnumL)
      jacR = ro_sview(metrics.jac, :, elnumR)
      alphas = ro_sview(capture.alpha_parallel[peer], :, i)
      resL = sview(eqn.res, :, :, elnumL)

      # compute delta w tilde and theta_bar = Dgk w_k + Dgn w_n
      getFaceVariables(capture, mesh.sbpface, iface_red, wL, wR, gradwL, gradwR,
                       nrm_face, delta_w, theta)

      # apply the penalty coefficient matrix
      applyPenalty(penalty, sbp, eqn.params, mesh.sbpface, diffusion,
                   iface_red, delta_w,
                   theta, qL, qR,  wL, wR, coordsL, coordsR, nrm_face, alphas,
                   dxidxL, dxidxR, jacL, jacR, t1L, t1R, t2L, t2R)

      # apply Rgk^T, Rgn^T, Dgk^T, Dgn^T
      # need to apply R^T * t1, not R^T * B * t1, so
      # interiorFaceIntegrate won't work.  Use the reverse mode instead
      for j=1:mesh.numNodesPerFace
        for k=1:mesh.numDofPerNode
          t1L[k, j] = -t1L[k, j]  # the SAT has a minus sign in front of it
          t1R[k, j] = -t1R[k, j]
        end
      end
      #TODO: this could be a boundaryIntegrate now that resR doesn't matter
      interiorFaceInterpolate_rev!(mesh.sbpface, iface_red, resL, resR,
                                   t1L, t1R)

      # apply Dgk^T and Dgn^T
      applyDgkTranspose(capture, sbp, eqn.params, mesh.sbpface, iface_red,
                        diffusion,
                        t2L, qL, wL, coordsL, nrm_face, dxidxL, jacL, resL, op)
    end  # end i
  end  # end peer

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
function computeNeumannBoundaryTerm(mesh::AbstractMesh{Tmsh}, sbp, eqn, opts,
                      capture::SBPParabolicSC{Tsol, Tres},
                      shockmesh::ShockedElements,
                      ) where {Tsol, Tres, Tmsh}

  # for shock capturing, apply the Neumann boundary condition
  # Lambda * grad_w = 0.  The resulting term is -R^T * B * Dgk * u

  @assert mesh.coord_order == 1

  temp_face = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  temp2_face = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  nrm_i = zeros(Tmsh, mesh.dim, mesh.numNodesPerFace)
  op = SummationByParts.Subtract()
  for i=1:shockmesh.numBoundaryFaces
    bndry_i = shockmesh.bndryfaces[i].bndry
    idx_orig = shockmesh.bndryfaces[i].idx_orig
    elnum_orig = shockmesh.elnums_all[bndry_i.element]


    if i < shockmesh.bndry_offsets[end]
      #nrm_i = ro_sview(mesh.nrm_bndry, :, :, idx_orig)
      for j=1:mesh.numNodesPerFace
        for d=1:mesh.dim
          nrm_i[d, j] = mesh.nrm_bndry[d, j, idx_orig]
        end
      end
    else
      #nrm_i = ro_sview(mesh.nrm_face, :, :, idx_orig)
      #TODO: nbrperm
      fac = shockmesh.bndryfaces[i].fac
      for j=1:mesh.numNodesPerFace
        for d=1:mesh.dim
          nrm_i[d, j] = fac*mesh.nrm_face[d, j, idx_orig]
        end
      end
    end
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
  Computes a Dirichlet BC that is consistent with the inviscid one
"""
function computeDirichletBoundaryTerm(mesh, sbp, eqn, opts,
                      capture::SBPParabolicSC{Tsol, Tres},
                      shockmesh::ShockedElements,
                      ) where {Tsol, Tres}

  for i=1:mesh.numBC
    bc_range = (shockmesh.bndry_offsets[i]:(shockmesh.bndry_offsets[i+1]-1))
    bc_func = mesh.bndry_funcs[i]

    calcBoundaryFlux(mesh, sbp, eqn, opts, shockmesh, capture,
                     capture.penalty, capture.diffusion, 
                     capture.entropy_vars, bc_func, bc_range)
  end

  return nothing
end

function calcBoundaryFlux(mesh::AbstractMesh, sbp, eqn::EulerData, opts,
                      shockmesh::ShockedElements,
                      capture::SBPParabolicSC{Tsol, Tres},
                      penalty::AbstractDiffusionPenalty,
                      diffusion::AbstractDiffusion,
                      entropy_vars::AbstractVariables, bc_func::BCType,
                      bc_range::AbstractVector,
                      ) where {Tsol, Tres}


  delta_w = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  res1 = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  op = SummationByParts.Subtract()

  for i in bc_range
    bndry_i = shockmesh.bndryfaces[i].bndry
    idx_orig = shockmesh.bndryfaces[i].idx_orig
    elnum_orig = shockmesh.elnums_all[bndry_i.element]

    w_i = ro_sview(capture.w_el, :, :, bndry_i.element)
    q_i = ro_sview(eqn.q, :, :, elnum_orig)
    coords_i = ro_sview(mesh.coords_bndry, :, :, idx_orig)
    nrm_i = ro_sview(mesh.nrm_bndry, :, :, idx_orig)
    alpha = capture.alpha_b[i]
    dxidxL = ro_sview(mesh.dxidx, :, :, :, elnum_orig)
    jacL = ro_sview(mesh.jac, :, elnum_orig)
    res_i = sview(eqn.res, :, :, elnum_orig)

    computeDirichletDelta(capture, eqn.params, mesh.sbpface, bndry_i,
                          bc_func, entropy_vars, w_i, q_i, coords_i,
                          nrm_i, delta_w)

    applyDirichletPenalty(penalty, sbp, eqn.params, mesh.sbpface, diffusion,
                          bndry_i, delta_w, q_i, w_i, coords_i, nrm_i, alpha,
                          dxidxL, jacL, res1)

    # apply R^T to T_D * delta_w
    scale!(res1, -1)  # SAT has a - sign in front of it

    boundaryFaceInterpolate_rev!(mesh.sbpface, bndry_i.face, res_i, res1)

    # apply B and then Dgk^T to delta_w
    for j=1:mesh.numNodesPerFace
      for i=1:mesh.numDofPerNode
        delta_w[i, j] *= -mesh.sbpface.wface[j]
      end
    end

    applyDgkTranspose(capture, sbp, eqn.params, mesh.sbpface, bndry_i,
                      diffusion, delta_w, q_i, w_i, coords_i, nrm_i, dxidxL,
                      jacL, res_i, op)

  end

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
  Computes the difference between the solution at the face and the
  prescribed Dirichlet boundary condition.  Specifically, computes the
  difference between the entropy variables interpolated to the face
  and the Dirichlet state computed from the interpolated conservative
  variables (for consistency with the inviscid BC).

  This only works for boundary conditions that define `getDirichletState`.

  **Inputs**

   * capture
   * params
   * sbpface
   * bndry: `Boundary` object
   * bndry_node: `BoundaryNode` object
   * func: `BCType` for this boundary condition
   * entropy_vars: an [`AbstractVariables`](@ref) 
   * wL: entropy variables for the element, `numDofPerNode` x
         `numNodesPerElement`
   * qL: conservative variables for the element, same size as `wL`
   * coords: `dim` x `numNodesPerFace` array containing the xyz coordinates
             of the face nodes
   * nrm_face: `dim` x `numNodesPerFace` containing the normal vector at the
               face nodes (can be scaled or unscaled, `getDirichletState`
               should not care).
   
  **Inputs/Outputs**

   * delta_w: array to be overwritten with the result. `numDofPerNode` x
              `numNodesPerFace`
"""
function computeDirichletDelta(capture::SBPParabolicSC{Tsol, Tres},
                               params::ParamType,
                               sbpface::AbstractFace, bndry::Boundary,
                               func::BCType, entropy_vars::AbstractVariables,
                               wL::AbstractMatrix, qL::AbstractMatrix,
                               coords::AbstractMatrix, nrm_face::AbstractMatrix,
                               delta_w::AbstractMatrix
                              ) where {Tsol, Tres}

  numDofPerNode, numNodesPerFace = size(delta_w)
  numNodesPerElement = size(wL, 2)

  #TODO: maybe interpolating the conservative variables and then converting
  #      would be more consistent, but Eq. 31 requires the same interpolation
  #      Rgk * w for both the interface and the boundary term.
  #      For the boundary state qg we have more leeway because it is zero
  #      for the energy stability analysis (-> entropy stability when u = w)
  wface = zeros(Tsol, numDofPerNode, numNodesPerFace)
  boundaryFaceInterpolate!(sbpface, bndry.face, wL, wface)

  qface = zeros(Tsol, numDofPerNode, numNodesPerFace)
  qgface = zeros(Tres, numDofPerNode)
  wgface = zeros(Tres, numDofPerNode)
  aux_vars = Tres[]

  # the inviscid boundary conditions interpolate q to the face and then
  # compute qg, so do that here as well for consistency
  boundaryFaceInterpolate!(sbpface, bndry.face, qL, qface)
  for i=1:numNodesPerFace
    q_i = sview(qface, :, i)
    coords_i = sview(coords, :, i)
    nrm_i = sview(nrm_face, :, i)

    getDirichletState(func, params, q_i, aux_vars, coords_i, nrm_i, qgface)

    convertToEntropy(entropy_vars, params, qgface, wgface)
#    fill!(wgface, 0)  #TODO: undo this
#    println("wgface = ", wgface)

    # compute the delta
    for j=1:numDofPerNode
      delta_w[j, i] = wface[j, i] - wgface[j]  
    end
  end

  return nothing
end






"""
  This function applies Dgk^T * t2L and Dgn^T * t2R.

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
                           params::ParamType,
                           sbpface, iface::Interface,
                           diffusion::AbstractDiffusion,
                           t2L::AbstractMatrix, t2R::AbstractMatrix,
                           qL::AbstractMatrix, qR::AbstractMatrix,
                           wL::AbstractMatrix, wR::AbstractMatrix,
                           coordsL::AbstractMatrix, coordsR::AbstractMatrix,
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
  applyDiffusionTensor(diffusion, sbp, params, qL, wL, coordsL, dxidxL, jacL,
                       iface.elementL, temp2L, temp3L)
  applyDiffusionTensor(diffusion, sbp, params, qR, wR,  coordsR, dxidxR, jacR,
                       iface.elementR, temp2R, temp3R)

  # saving temp3 to a mesh-wide array and then applying Dx^T would save
  # a lot of flops.
  applyDxTransposed(sbp, temp3L, dxidxL, jacL, work, resL, op)
  applyDxTransposed(sbp, temp3R, dxidxR, jacR, work, resR, op)

  return nothing
end

"""
  Method to apply Dgk^T * t2L only.
"""
function applyDgkTranspose(capture::SBPParabolicSC{Tsol, Tres}, sbp,
                           params::ParamType,
                           sbpface, bndry::Union{Interface, Boundary},
                           diffusion::AbstractDiffusion,
                           t2L::AbstractMatrix,
                           qL::AbstractMatrix,
                           wL::AbstractMatrix,
                           coordsL::AbstractMatrix,
                           nrm_face::AbstractMatrix,
                           dxidxL::Abstract3DArray,
                           jacL::AbstractVector,
                           resL::AbstractMatrix,
                           op::SummationByParts.UnaryFunctor=SummationByParts.Add()) where {Tsol, Tres}

  dim, numNodesPerFace = size(nrm_face)
  numNodesPerElement = size(resL, 2)
  numDofPerNode = size(wL, 1)

  @unpack capture temp1L temp2L temp3L work
  fill!(temp2L, 0)

  # apply N and R^T
  @simd for d=1:dim
    @simd for j=1:numNodesPerFace
      @simd for k=1:numDofPerNode
        temp1L[k, j] =  nrm_face[d, j]*t2L[k, j]
      end
    end

    tmp2L = sview(temp2L, :, :, d);
    boundaryFaceInterpolate_rev!(sbpface, getFaceL(bndry), tmp2L, temp1L)
  end

  # multiply by D^T Lambda
  applyDiffusionTensor(diffusion, sbp, params, qL, wL, coordsL, dxidxL, jacL,
                       getElementL(bndry), temp2L, temp3L)

  applyDxTransposed(sbp, temp3L, dxidxL, jacL, work, resL, op)

  return nothing
end



#------------------------------------------------------------------------------
# helper functions

"""
  Computes the alpha parameter.  Starts parallel communication so alpha
  for remote elements can be computed later.

  **Inputs**

   * caputure: `capture.alpha` is overwritten.
   * mesh
   * shockmesh
"""
function computeAlpha(capture::SBPParabolicSC, mesh::AbstractMesh,
                      shockmesh::ShockedElements)

  # I think the alpha_gk parameter should be 0 for faces on the Neumann boundary
  # Count the number of faces each element has in the mesh to compute alpha_gk
  # such that it sums to 1.
  if size(capture.alpha, 2) < shockmesh.numInterfaces
    capture.alpha = zeros(2, shockmesh.numInterfaces)
  else
    fill!(capture.alpha, 0)
  end

  if !shockmesh.isNeumann
    if length(capture.alpha_b) < shockmesh.numBoundaryFaces
      capture.alpha_b = zeros(shockmesh.numBoundaryFaces)
    else
      fill!(capture.alpha_b, 0)
    end
  end

  oldsize = length(capture.alpha_parallel)
  resize!(capture.alpha_parallel, shockmesh.npeers)
  for peer=1:shockmesh.npeers
    len_i = shockmesh.numSharedInterfaces[peer]
    # don't try to get the size of alpha_parallel[peer] if it is undefined
    if peer > oldsize || size(capture.alpha_parallel[peer], 2) < len_i
      capture.alpha_parallel[peer] = zeros(2, len_i)
    else
      fill!(capture.alpha_parallel[peer], 0)
    end
  end

  # count number of faces of each element not on boundary
  el_counts = zeros(UInt8, shockmesh.numEl)
  for i=1:shockmesh.numInterfaces
    iface_i = shockmesh.ifaces[i].iface
    elnumL = shockmesh.elnums_all[iface_i.elementL]
    elnumR = shockmesh.elnums_all[iface_i.elementR]
    el_counts[iface_i.elementL] += 1
    el_counts[iface_i.elementR] += 1
  end

  if !shockmesh.isNeumann
    for i=1:shockmesh.numBoundaryFaces
      bndry_i = shockmesh.bndryfaces[i].bndry
      el_counts[bndry_i.element] += 1
    end
  end

  for peer=1:shockmesh.npeers
    for i=1:shockmesh.numSharedInterfaces[peer]
      iface_i = shockmesh.shared_interfaces[peer][i].iface
      el_counts[iface_i.elementL] += 1
      # only do local elements, because we can't see all faces of remote
      # elements
    end
  end

  # start parallel communications
  sendAlphas(capture.alpha_comm, shockmesh, el_counts)

  # compute alpha from face counts
  # for neighbor elements alpha doesn't sum to 1, but thats ok because
  # lambda = 0 there, so alpha multiplies zero.
  for i=1:shockmesh.numInterfaces
    iface_i = shockmesh.ifaces[i].iface
    capture.alpha[1, i] = 1/el_counts[iface_i.elementL]
    capture.alpha[2, i] = 1/el_counts[iface_i.elementR]
  end

  if !shockmesh.isNeumann
    for i=1:shockmesh.numBoundaryFaces
      bndry_i = shockmesh.bndryfaces[i].bndry
      capture.alpha_b[i] = 1/el_counts[bndry_i.element]
    end
  end

  for peer=1:shockmesh.npeers
    for i=1:shockmesh.numSharedInterfaces[peer]
      iface_i = shockmesh.shared_interfaces[peer][i].iface
      capture.alpha_parallel[peer][1, i] = 1/el_counts[iface_i.elementL]
    end
  end

  return nothing
end



"""
  Send the number of faces not on a Neumann boundary each element has to
  other processes.

  **Inputs**

   * comm: `AlphaComm` object
   * shockmesh: `ShockedElements`
   * el_counts: array containing the number of faces each local element has
                that are not on a Neumann boundary
               
"""
function sendAlphas(comm::AlphaComm, shockmesh::ShockedElements,
                    el_counts::AbstractVector{I}) where {I <: Integer}

  # extract the required values into the send buffers

  # allocateArrays checks that previous communications are complete

  # We only need to send number of faces per element, but that would require
  # uniquely identifying the elements on either side of the partition boundary.
  # Instead, send the same value for all faces of the element and let the
  # the receiver sort it out

  # post receives
  for peer=1:shockmesh.npeers
    Irecv!(comm, peer)
  end

  # send
  for peer=1:shockmesh.npeers
    println("Sending alphas to process ", comm.peer_parts_red[peer]); flush(STDOUT)
    # pack the buffer
    for j=1:shockmesh.numSharedInterfaces[peer]
      iface_i = shockmesh.shared_interfaces[peer][j].iface
      comm.send_bufs[peer][j] = el_counts[iface_i.elementL]
    end

    Isend(comm, peer)
  end

  return nothing
end


"""
  Receives the parallel communication started by `sendAlphas`.  Specifically,
  it receives the communication from one peer process and returns the index
  of that process.  This function should be called exactly the same number of
  times that there are peer process.

  On exit, `capture.alpha` will be populated with the correct alpha values
  for the element shared with the peer process.

  **Inputs**

   * capture: `SBPParabolicSC` shock capturing object
   * shockmesh: `ShockedElements`

  **Outputs**

   * idx: the index of the peer that communication was received from.
"""
function receiveAlpha(capture::SBPParabolicSC, shockmesh::ShockedElements)

  peer = Waitany!(capture.alpha_comm)
  comm = capture.alpha_comm

  # unpack buffer and calculate alpha
  for j=1:shockmesh.numSharedInterfaces[peer]
    iface_i = shockmesh.shared_interfaces[peer][j]
    capture.alpha_parallel[peer][2, j] = 1/comm.recv_bufs[peer][j]
    # the buffer was packed with el_counts for the parent element of each
    # face, so all thats required it so take 1/the value
  end

  return peer
end



