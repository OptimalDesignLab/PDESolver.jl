# shock capturing using any the SBP-SAT discretization of the second derivative
# term

import SummationByParts: UnaryFunctor, Add, Subtract

"""
  Applies [`SBPParabolic`](@ref) shock capturing.
"""
function calcShockCapturing(mesh::AbstractMesh, sbp::AbstractOperator,
                             eqn::EulerData, opts,
                             capture::SBPParabolicReducedSC{Tsol, Tres},
                             shockmesh::ShockedElements) where {Tsol, Tres}

  computeVolumeTerm(mesh, sbp, eqn, opts, capture, capture.entropy_vars,
                    capture.diffusion, shockmesh)

  finishConvertEntropy(mesh, sbp, eqn, opts, capture, capture.entropy_vars,
                       shockmesh)


  computeFaceTerm(mesh, sbp, eqn, opts, capture, shockmesh, capture.sensor_const,
                  capture.penalty)

  #println("after face term, residual norm = ", calcNorm(eqn, eqn.res))
#=
  if shockmesh.isNeumann
    computeNeumannBoundaryTerm(mesh, sbp, eqn, opts, capture, shockmesh)
  else
    computeDirichletBoundaryTerm(mesh, sbp, eqn, opts, capture, shockmesh)
  end
=#
#  computeSharedFaceTerm(mesh, sbp, eqn, opts, capture, shockmesh,
#                              capture.diffusion, capture.penalty)



  return nothing
end

function calcShockCapturing_diff(mesh::AbstractMesh, sbp::AbstractOperator,
                             eqn::EulerData, opts,
                             capture::SBPParabolicReducedSC{Tsol, Tres},
                             shockmesh::ShockedElements,
                             assem::AssembleElementData) where {Tsol, Tres}

  println("\nDoing volume term")
  computeVolumeTerm_diff(mesh, sbp, eqn, opts, capture,
                               capture.diffusion, capture.entropy_vars,
                               shockmesh, assem)
  finishConvertEntropy(mesh, sbp, eqn, opts, capture, capture.entropy_vars,
                       shockmesh)

  println("\nDoing face term")
  computeFaceTerm_diff(mesh, sbp, eqn, opts, capture, shockmesh,
                       capture.sensor_const, capture.entropy_vars,
                       capture.penalty, assem)
#=
  if shockmesh.isNeumann
    println("computing Neumann boundary condition")
    computeNeumannBoundaryTerm_diff(mesh, sbp, eqn, opts, capture, shockmesh,
                              capture.diffusion, capture.entropy_vars, assem)
  else
    println("computing Dirichlet boundary condition")
    computeDirichletBoundaryTerm_diff(mesh, sbp, eqn, opts, capture, shockmesh,
                                      assem)
  end
=#
  #@time computeSharedFaceTerm_diff(mesh, sbp, eqn, opts, capture, shockmesh,
  #                                 capture.diffusion, capture.penalty)



  return nothing
end


"""
  Compute the entropy variables for every element not done by
  `computeVolumeTerm`.
"""
@noinline function finishConvertEntropy(mesh, sbp, eqn, opts,
                              capture::SBPParabolicReducedSC,
                              entropy_vars::AbstractVariables,
                              shockmesh::ShockedElements)

  for i in shockmesh.neighbor_els
    i_full = shockmesh.elnums_all[i]
    for j=1:mesh.numNodesPerElement
      q_j = ro_sview(eqn.q, :, j, i_full)
      w_j = sview(capture.w_el, :, j, i)
      convertToEntropy(entropy_vars, eqn.params, q_j, w_j)
    end
  end

  return nothing
end


"""
  Does the volume term Q^T * Lambda D * w.  Also saves the entropy variables
  for all shocked elements in `capture.w_el`
"""
@noinline function computeVolumeTerm(mesh, sbp, eqn, opts,
                      capture::SBPParabolicReducedSC{Tsol, Tres},
                      entropy_vars::AbstractVariables,
                      diffusion::AbstractDiffusion,
                      shockmesh::ShockedElements
                     ) where {Tsol, Tres}

  @assert eqn.params.use_Minv != 1

  work = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  grad_w = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  lambda_gradw = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement,
                             mesh.dim)
  op = SummationByParts.Subtract()

  @simd for i=1:shockmesh.numShock
    i_full = shockmesh.elnums_all[i]
    w_i = sview(capture.w_el, :, :, i)  # save for use in face integrals
    @simd for j=1:mesh.numNodesPerElement
      q_j = ro_sview(eqn.q, :, j, i)
      w_j = sview(w_i, :, j)

      # convert to entropy variables
      convertToEntropy(entropy_vars, eqn.params, q_j, w_j)
    end

    # apply D operator
    q_i = ro_sview(eqn.q, :, :, i_full)
    coords = ro_sview(mesh.coords, :, :, i_full)
    dxidx_i = ro_sview(mesh.dxidx, :, :, :, i_full)
    jac_i = ro_sview(mesh.jac, :, i_full)

    fill!(grad_w, 0)
    applyDx(sbp, w_i, dxidx_i, jac_i, work, grad_w)

    # apply diffusion tensor
    applyDiffusionTensor(diffusion, sbp, eqn.params, q_i, w_i, coords, dxidx_i,
                         jac_i, i, grad_w, lambda_gradw)

    # apply Q^T
    res_i = sview(eqn.res, :, :, i)
    applyQxTransposed(sbp, lambda_gradw, dxidx_i, work, res_i, op)

  end  # end i

  return nothing
end


@noinline function computeVolumeTerm_diff(mesh, sbp, eqn, opts,
                           capture::SBPParabolicReducedSC{Tsol, Tres},
                           diffusion::AbstractDiffusion,
                           entropy_vars::AbstractVariables,
                           shockmesh::ShockedElements,
                           assembler::AssembleElementData) where {Tsol, Tres}

  @assert eqn.params.use_Minv != 1

  work = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  Dx = zeros(mesh.numNodesPerElement, mesh.numNodesPerElement, mesh.dim)
  w_dot = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerElement)
  t1_dot = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode, mesh.dim,
                 mesh.numNodesPerElement, mesh.numNodesPerElement)
  t2_dot = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode, mesh.dim,
                 mesh.numNodesPerElement, mesh.numNodesPerElement)
  gradq_i = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  # if dim was the last dimension of t1_dot, we could re-use it, but instead
  # we have to use a 3rd array
  res_jac = eqn.params.calc_volume_integrals_data.res_jac
  op = SummationByParts.Subtract()

  for i=1:shockmesh.numShock
    i_full = shockmesh.elnums_all[i]

    q_i = ro_sview(eqn.q, :, :, i_full)
    w_i = sview(capture.w_el, :, :, i)
    coords_i = ro_sview(mesh.coords, :, :, i_full)
    dxidx_i = ro_sview(mesh.dxidx, :, :, :, i_full)
    jac_i = ro_sview(mesh.jac, :, i_full)

    # to compute the jacobian, dq/dq = I, so dw/dq = the dw/dq evaluated at
    # each node
    for j=1:mesh.numNodesPerElement
      q_j = ro_sview(eqn.q, :, j, i_full)
      w_j = sview(w_i, :, j)
      convertToEntropy(entropy_vars, eqn.params, q_j, w_j)
      w_dot_j = sview(w_dot, :, :, j)
      getA0inv(entropy_vars, eqn.params, q_j, w_dot_j)
    end

    # apply Dx * w_dot
    calcDx(sbp, dxidx_i, jac_i, Dx)
    applyOperatorJac(Dx, w_dot, t1_dot)

    # need Dx * w for applyDiffusionTensor
    for d=1:mesh.dim
      Dx_d = ro_sview(Dx, :, :, d)
      gradq_d = sview(gradq_i, :, :, d)
      smallmatmatT!(w_i, Dx_d, gradq_d)
    end

    # apply the diffusion tensor to all nodes of the element
    applyDiffusionTensor_diff(diffusion, sbp, eqn.params, q_i, w_i, coords_i,
                              dxidx_i, jac_i, i, gradq_i, t1_dot, t2_dot)

    # apply Qx^T
    calcQxTransposed(sbp, dxidx_i, Dx)
    applyOperatorJac(Dx, t2_dot, res_jac, true, op)

    assembleElement(assembler, mesh, i_full, res_jac)
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
@noinline function computeFaceTerm(mesh, sbp, eqn, opts,
                      capture::SBPParabolicReducedSC{Tsol, Tres},
                      shockmesh::ShockedElements, sensor::ShockSensorHApprox,
                      penalty::AbstractDiffusionPenalty) where {Tsol, Tres}

  setElementCounts(sensor, shockmesh.numShock, shockmesh.numEl)

  delta_w = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)

  t1L = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  t1R = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
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

#    coordsL = ro_sview(mesh.coords, :, :, elnumL)
#    coordsR = ro_sview(mesh.coords, :, :, elnumR)
    nrm_face = ro_sview(mesh.nrm_face, :, :, iface_idx)
#    dxidxL = ro_sview(mesh.dxidx, :, :, :, elnumL)
#    dxidxR = ro_sview(mesh.dxidx, :, :, :, elnumR)
    jacL = ro_sview(mesh.jac, :, elnumL)
    jacR = ro_sview(mesh.jac, :, elnumR)
#    alphas = ro_sview(capture.alpha, :, i)
    resL = sview(eqn.res, :, :, elnumL)
    resR = sview(eqn.res, :, :, elnumR)

    # compute delta w tilde
#    getFaceVariables(capture, mesh.sbpface, iface_red, wL, wR, delta_w)

    # apply the penalty coefficient matrix
#    applyReducedPenalty(penalty, sbp, eqn.params, mesh.sbpface, diffusion,
#                 iface_red, delta_w,
#                 qL, qR, wL, wR, coordsL, coordsR, nrm_face, alphas, dxidxL,
#                 dxidxR, jacL, jacR, t1L, t1R)
    applyReducedPenalty(penalty, sbp, eqn.params, mesh.sbpface, sensor,
                        iface_red, wL, wR, nrm_face, jacL, jacR, resL, resR,
                        op)

#=
    # apply Rgk^T, Rgn^T,
    # need to apply R^T * t1, not R^T * B * t1, so
    # interiorFaceIntegrate won't work.  Use the reverse mode instead
    for j=1:mesh.numNodesPerFace
      for k=1:mesh.numDofPerNode
        t1L[k, j] = -t1L[k, j]  # the SAT has a minus sign in front of it
        t1R[k, j] = -t1R[k, j]
      end
    end
    interiorFaceInterpolate_rev!(mesh.sbpface, iface_red, resL, resR, t1L, t1R)
=#
  end  # end loop i

  return nothing
end


@noinline function computeFaceTerm_diff(mesh, sbp, eqn, opts,
                      capture::SBPParabolicReducedSC{Tsol, Tres},
                      shockmesh::ShockedElements, sensor::ShockSensorHApprox,
                      entropy_vars::AbstractVariables,
                      penalty::AbstractDiffusionPenalty,
                      assem::AssembleElementData) where {Tsol, Tres}

  setElementCounts(sensor, shockmesh.numShock, shockmesh.numEl)

  delta_w = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  delta_w_dotL = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode,
                             mesh.numNodesPerFace, mesh.numNodesPerElement)
  delta_w_dotR = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode,
                             mesh.numNodesPerFace, mesh.numNodesPerElement)

#  res1L_dotL = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode,
#                           mesh.numNodesPerFace, mesh.numNodesPerElement)
#  res1L_dotR = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode,
#                           mesh.numNodesPerFace, mesh.numNodesPerElement)
#
#  res1R_dotL = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode,
#                           mesh.numNodesPerFace, mesh.numNodesPerElement)
#  res1R_dotR = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode,
#                           mesh.numNodesPerFace, mesh.numNodesPerElement)


  data = eqn.params.calc_face_integrals_data
  @unpack data res_jacLL res_jacLR res_jacRL res_jacRR

  fill!(res_jacLL, 0); fill!(res_jacRL, 0)
  fill!(res_jacLR, 0); fill!(res_jacRR, 0)

  #t1L = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  #t1R = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  t2L = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  t2R = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
#  subtract = SummationByParts.Subtract()
#  add = SummationByParts.Add()

  op = Subtract()

  for i=1:shockmesh.numInterfaces
    iface_red = shockmesh.ifaces[i].iface
    iface_idx = shockmesh.ifaces[i].idx_orig
    elnumL = shockmesh.elnums_all[iface_red.elementL]
    elnumR = shockmesh.elnums_all[iface_red.elementR]
    iface_full = replace_interface(iface_red, elnumL, elnumR)

    # get data needed for next steps
#    wL = ro_sview(capture.w_el, :, :, iface_red.elementL)
#    wR = ro_sview(capture.w_el, :, :, iface_red.elementR)
    qL = ro_sview(eqn.q, :, :, elnumL)
    qR = ro_sview(eqn.q, :, :, elnumR)

#    coordsL = ro_sview(mesh.coords, :, :, elnumL)
#    coordsR = ro_sview(mesh.coords, :, :, elnumR)
    nrm_face = ro_sview(mesh.nrm_face, :, :, iface_idx)
#    dxidxL = ro_sview(mesh.dxidx, :, :, :, elnumL)
#    dxidxR = ro_sview(mesh.dxidx, :, :, :, elnumR)
    jacL = ro_sview(mesh.jac, :, elnumL)
    jacR = ro_sview(mesh.jac, :, elnumR)
#    alphas = ro_sview(capture.alpha, :, i)

    #TODO: stop doing this
    fill!(res_jacLL, 0); fill!(res_jacLR, 0)
    fill!(res_jacRL, 0); fill!(res_jacRR, 0)

#    has_imag = maximum(imag(qL)) > 1e-13 || maximum(imag(qR)) > 1e-13

    # compute delta w tilde
#    getFaceVariables_diff(eqn.params, capture, entropy_vars,
#                          sbp, mesh.sbpface,
#                          iface_red, wL, wR, qL, qR,
#                          delta_w, delta_w_dotL, delta_w_dotR)

    # apply the penalty coefficient matrix
    applyReducedPenalty_diff(penalty, sbp, eqn.params,  mesh.sbpface,
                      sensor, entropy_vars, iface_red, qL, qR,  nrm_face,
                      jacL, jacR,
                      res_jacLL, res_jacLR, res_jacRL, res_jacRR, op)


    # apply Rgk^T, Rgn^T,
#    fill!(res_jacLL, 0); fill!(res_jacRL, 0)
#    fill!(res_jacLR, 0); fill!(res_jacRR, 0)
#    interiorFaceIntegrate_jac!(mesh.sbpface, iface_red, res1L_dotL, res1R_dotL,
#                               res_jacLL, res_jacRL, subtract, subtract, false)
#    interiorFaceIntegrate_jac!(mesh.sbpface, iface_red, res1L_dotR, res1R_dotR,
#                               res_jacLR, res_jacRR, subtract, subtract, false)

#    assembleInterfaceVisc(assem, mesh.sbpface, mesh, iface_full, res_jacLL,
#                          res_jacLR, res_jacRL, res_jacRR)
    assembleInterface(assem, mesh.sbpface, mesh, iface_full, res_jacLL,
                          res_jacLR, res_jacRL, res_jacRR)

  end  # end loop i

  return nothing
end


#=
"""
  Does the same thing as `compuateFaceTerm`, but for the shared faces, updating
  the residual on the local element only
"""
function computeSharedFaceTerm(mesh, sbp, eqn, opts,
                      capture::SBPParabolicReducedSC{Tsol, Tres},
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
=#

#=
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
                      capture::SBPParabolicReducedSC{Tsol, Tres},
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
                      capture::SBPParabolicReducedSC{Tsol, Tres},
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
                      capture::SBPParabolicReducedSC{Tsol, Tres},
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
=#

#=
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

  **Inputs/Outputs**

   * delta_w: `numDofPerNode` x `numNodesPerFace`, overwritten
"""
function getFaceVariables(capture::SBPParabolicReducedSC{Tsol, Tres},
                          sbpface::AbstractFace, iface_red::Interface,
                          wL::AbstractMatrix, wR::AbstractMatrix,
                          delta_w::AbstractMatrix
                         ) where {Tsol, Tres}

  numDofPerNode, numNodesPerElement = size(wL)
  numNodesPerFace = size(delta_w, 2)

  @unpack capture w_faceL w_faceR

  interiorFaceInterpolate!(sbpface, iface_red, wL, wR, w_faceL, w_faceR)


  @simd for j=1:numNodesPerFace
    @simd for k=1:numDofPerNode
      delta_w[k, j] = w_faceL[k, j] - w_faceR[k, j]
    end
  end

  return nothing
end
=#


#=
"""
  Differentiated version of [`getFaceVariables`](@ref).
"""
function getFaceVariables_diff(params::ParamType,
                          capture::SBPParabolicReducedSC{Tsol, Tres},
                          entropy_vars::AbstractVariables,
                          sbp::AbstractOperator,
                          sbpface::AbstractFace, iface_red::Interface,
                          wL::AbstractMatrix, wR::AbstractMatrix,
                          qL::AbstractMatrix, qR::AbstractMatrix,
                          delta_w::AbstractMatrix,
                          delta_w_dotL::Abstract4DArray, delta_w_dotR::Abstract4DArray,
                         ) where {Tsol, Tres}

  numDofPerNode, numNodesPerElement = size(wL)
  numNodesPerFace = size(delta_w, 2)

  getFaceVariables(capture, sbpface, iface_red, wL, wR, delta_w)

  @unpack capture wL_dot wR_dot

  # derivative of delta_w wrt qL and qR
  for i=1:numNodesPerElement
    wL_dot_i = sview(wL_dot, :, :, i)
    wR_dot_i = sview(wR_dot, :, :, i)
    qL_i = ro_sview(qL, :, i)
    qR_i = ro_sview(qR, :, i)
    getA0inv(entropy_vars, params, qL_i, wL_dot_i)
    getA0inv(entropy_vars, params, qR_i, wR_dot_i)
  end

  fill!(delta_w_dotL, 0); fill!(delta_w_dotR, 0)
  interiorFaceInterpolate_jac!(sbpface, iface_red, wL_dot, wR_dot,
                              delta_w_dotL, delta_w_dotR)

  @simd for i=1:length(delta_w_dotR)
    delta_w_dotR[i] = -delta_w_dotR[i]  # wR is negated in wL - wR
  end

  return nothing
end
=#


#=
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
function computeDirichletDelta(capture::SBPParabolicReducedSC{Tsol, Tres},
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
=#



#------------------------------------------------------------------------------
# helper functions

#TODO: see if this can be shared with the SBPParabolicSC
 
"""
  Computes the alpha parameter.  Starts parallel communication so alpha
  for remote elements can be computed later.

  **Inputs**

   * caputure: `capture.alpha` is overwritten.
   * mesh
   * shockmesh
"""
function computeAlpha(capture::SBPParabolicReducedSC, mesh::AbstractMesh,
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
  Receives the parallel communication started by `sendAlphas`.  Specifically,
  it receives the communication from one peer process and returns the index
  of that process.  This function should be called exactly the same number of
  times that there are peer process.

  On exit, `capture.alpha` will be populated with the correct alpha values
  for the element shared with the peer process.

  **Inputs**

   * capture: `SBPParabolicReducedSC` shock capturing object
   * shockmesh: `ShockedElements`

  **Outputs**

   * idx: the index of the peer that communication was received from.
"""
function receiveAlpha(capture::SBPParabolicReducedSC, shockmesh::ShockedElements)

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


#------------------------------------------------------------------------------
# Penalty functions

#=
"""
  Applies a reduced form of the penalty for the SBP-SAT generalization of the
  modified scheme of Bassi and Rebay (BR2).
"""
function applyReducedPenalty(penalty::BR2Penalty{Tsol, Tres}, sbp,
                      params::AbstractParamType, sbpface,
                      diffusion::AbstractDiffusion, iface::Interface,
                      delta_w::AbstractMatrix{Tsol},
                      qL_el::AbstractMatrix, qR_el::AbstractMatrix,
                      wL::AbstractMatrix, wR::AbstractMatrix,
                      coordsL::AbstractMatrix, coordsR::AbstractMatrix,
                      nrm_face::AbstractMatrix,
                      alphas::AbstractVector,
                      dxidxL::Abstract3DArray, dxidxR::Abstract3DArray,
                      jacL::AbstractVector, jacR::AbstractVector,
                      res1L::AbstractMatrix, res1R::AbstractMatrix,
                     ) where {Tsol, Tres}

  fill!(res1L, 0); fill!(res1R, 0)
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
  applyDiffusionTensor(diffusion, sbp, params, qL_el, wL, coordsL, dxidxL,
                       jacL, iface.elementL, sbpface, iface.faceL,
                       qL, t1L)
  applyDiffusionTensor(diffusion, sbp, params, qR_el, wR, coordsR, dxidxR,
                       jacR, iface.elementR, sbpface, iface.faceR,
                       qR, t1R)

  # apply inverse mass matrix, then apply B*Nx*R*t2L_x + B*Ny*R*t2L_y
  @simd for d1=1:dim
    @simd for j=1:numNodesPerElement
      facL = (1/alphas[1])*jacL[j]/sbp.w[j]
      facR = (1/alphas[2])*jacR[j]/sbp.w[j]
      #facL = (3)*jacL[j]/sbp.w[j]
      #facR = (3)*jacR[j]/sbp.w[j]

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

  return nothing
end
=#

function applyReducedPenalty(penalty::BR2Penalty{Tsol, Tres}, sbp,
                      params::AbstractParamType{Tdim}, sbpface::SparseFace,
                      sensor::ShockSensorHApprox, iface::Interface,
                      wL::AbstractMatrix, wR::AbstractMatrix,
                      nrm_face::AbstractMatrix,
                      jacL::AbstractVector{Tmsh}, jacR::AbstractVector,
                      res1L::AbstractMatrix, res1R::AbstractMatrix,
                      op::UnaryFunctor=Add()
                     ) where {Tsol, Tres, Tmsh, Tdim}

  numDofPerNode = size(wL, 1)
  numNodesPerFace = size(nrm_face, 2)
  numNodesPerElement = length(jacL)

  eps_L = getShockSensor(params, sbp, sensor, iface.elementL, jacL)
  eps_R = getShockSensor(params, sbp, sensor, iface.elementR, jacR)

  alphaL = eps_L/4
  alphaR = eps_R/4
  for i=1:numNodesPerFace
    # figure out the volume node index
    # This is basically what interiorFaceInterpolate does
    nperm = sbpface.nbrperm[i, iface.orient]
    iL = sbpface.perm[i, iface.faceL]; iR = sbpface.perm[nperm, iface.faceR]

    # the term to compute is:
    # [ B * Nx * Rgk * H^-1 * Rgk^T * Nx * B * delta_w
    #   B * Ny * Rgk * H^-1 * Rgk^T * Ny * B * delta_w] 
    # and this has to be done for both elements kappa and nu.
    # This term is applied to all mesh.numDofPerNode components of delta_w

    B2_i = sbpface.wface[i]*sbpface.wface[i]
    HinvL = jacL[iL]/sbp.w[iL]; HinvR = jacR[iR]/sbp.w[iR]
    facL = zero(Tres); facR = zero(Tres)
    @simd for d=1:Tdim
      N2_d = nrm_face[d, i]*nrm_face[d, i]
      facL += B2_i*N2_d*HinvL
      facR += B2_i*N2_d*HinvR
    end
    fac = op(alphaL*facL + alphaR*facR)

    @simd for j=1:numDofPerNode
      delta_w_j = wL[j, iL] - wR[j, iR]
      res1L[j, iL] += fac*delta_w_j
      res1R[j, iR] -= fac*delta_w_j  # delta_w is reversed for elementR
    end  # end j
  end  # end i

  return nothing
end


function applyReducedPenalty_diff(penalty::BR2Penalty{Tsol, Tres}, sbp,
                      params::AbstractParamType{Tdim}, sbpface::SparseFace,
                      sensor::ShockSensorHApprox,
                      entropy_vars::AbstractVariables, iface::Interface,
                      qL_el::AbstractMatrix, qR_el::AbstractMatrix,
                      nrm_face::AbstractMatrix{Tmsh},
                      jacL::AbstractVector, jacR::AbstractVector,
                      res1L_dotL::Abstract4DArray, res1L_dotR::Abstract4DArray,
                      res1R_dotL::Abstract4DArray, res1R_dotR::Abstract4DArray,
                      op::UnaryFunctor=Add()
                     ) where {Tsol, Tres, Tmsh, Tdim}
# this function does not take in delta_w_dot to save on computation

  numDofPerNode = size(qL_el, 1)
  numNodesPerFace = size(nrm_face, 2)
  numNodesPerElement = length(jacL)

  eps_L = getShockSensor(params, sbp, sensor, iface.elementL, jacL)
  eps_R = getShockSensor(params, sbp, sensor, iface.elementR, jacR)


  alphaL = eps_L/4
  alphaR = eps_R/4

  A0invL = zeros(Tsol, numDofPerNode, numDofPerNode)
  A0invR = zeros(Tsol, numDofPerNode, numDofPerNode)
  for i=1:numNodesPerFace
    # figure out the volume node index
    # This is basically what interiorFaceInterpolate does
    nperm = sbpface.nbrperm[i, iface.orient]
    iL = sbpface.perm[i, iface.faceL]; iR = sbpface.perm[nperm, iface.faceR]

    # the term to compute is:
    # [ B * Nx * Rgk * H^-1 * Rgk^T * Nx * B * delta_w
    #   B * Ny * Rgk * H^-1 * Rgk^T * Ny * B * delta_w] 
    # and this has to be done for both elements kappa and nu.
    # This term is applied to all mesh.numDofPerNode components of delta_w

    B2_i = sbpface.wface[i]*sbpface.wface[i]
    HinvL = jacL[iL]/sbp.w[iL]; HinvR = jacR[iR]/sbp.w[iR]
    facL = zero(Tres); facR = zero(Tres)
    @simd for d=1:Tdim
      N2_d = nrm_face[d, i]*nrm_face[d, i]
      facL += B2_i*N2_d*HinvL
      facR += B2_i*N2_d*HinvR
    end
    fac = op(alphaL*facL + alphaR*facR)

    getA0inv(entropy_vars, params, ro_sview(qL_el, :, iL), A0invL)
    getA0inv(entropy_vars, params, ro_sview(qR_el, :, iR), A0invR)
    @simd for j=1:numDofPerNode
      @simd for k=1:numDofPerNode
        res1L_dotL[k, j, iL, iL] =  fac*A0invL[k, j]
        res1L_dotR[k, j, iL, iR] = -fac*A0invR[k, j]
        res1R_dotL[k, j, iR, iL] = -fac*A0invL[k, j]
        res1R_dotR[k, j, iR, iR] =  fac*A0invR[k, j]
      end  # end k
    end  # end j
  end  # end i

  return nothing
end




#=
"""
  Differentiated version of the reduced BR2 penalty
"""
function applyReducedPenalty_diff(penalty::BR2Penalty{Tsol, Tres}, sbp,
                      params::AbstractParamType, sbpface,
                      diffusion::AbstractDiffusion, iface::Interface,
                      delta_w::AbstractMatrix{Tsol},
                      delta_w_dotL::Abstract4DArray, delta_w_dotR::Abstract4DArray,
                      qL_el::AbstractMatrix, qR_el::AbstractMatrix,
                      wL::AbstractMatrix, wR::AbstractMatrix,
                      coordsL::AbstractMatrix, coordsR::AbstractMatrix,
                      nrm_face::AbstractMatrix,
                      alphas::AbstractVector,
                      dxidxL::Abstract3DArray, dxidxR::Abstract3DArray,
                      jacL::AbstractVector, jacR::AbstractVector,
                      res1L_dotL::Abstract4DArray, res1L_dotR::Abstract4DArray,
                      res1R_dotL::Abstract4DArray, res1R_dotR::Abstract4DArray,
                     ) where {Tsol, Tres}

  numDofPerNode, numNodesPerFace = size(delta_w)
  numNodesPerElement = length(jacL)
  dim = size(nrm_face, 1)
  add = SummationByParts.Add()

  # nodes in the stencil of R
  nodesL = sview(sbpface.perm, :, iface.faceL)
  nodesR = sview(sbpface.perm, :, iface.faceR)

  # the first L/R indicates if this is element kappa or nu
  # the second L/R indices if the derivative is wrt qL or qR
  @unpack penalty delta_w_n qL qR t1_dotL t1_dotR t2LL t2LR t2RL t2RR
  @unpack penalty t3LL t3LR t3RL t3RR t4RL t4RR
  t4LL = t1_dotL; t4LR = t1_dotR

  fill!(qL, 0); fill!(qR, 0)

  #--------------------------
  # T1
  
  # compute primal quantity needed later
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

  # compute derivative quantity
  # apply B and [Nx, Ny]
  @simd for q=1:numNodesPerElement
    @simd for p=1:numNodesPerFace
      @simd for d=1:dim
        fac = 0.25*sbpface.wface[p]*nrm_face[d, p]
        @simd for j=1:numDofPerNode
          @simd for i=1:numDofPerNode
            t1_dotL[i, j, d, p, q] = fac*delta_w_dotL[i, j, p, q]
            t1_dotR[i, j, d, p, q] = fac*delta_w_dotR[i, j, p, q]
          end
        end
      end
    end
  end

  # interpolate back to volume nodes
  fill!(t2LL, 0); fill!(t2RL, 0); fill!(t2LR, 0); fill!(t2RR, 0)
  interiorFaceIntegrate_jac!(sbpface, iface, t1_dotL, t1_dotL, t2LL, t2RL, add,
                             add, false)
  interiorFaceIntegrate_jac!(sbpface, iface, t1_dotR, t1_dotR, t2LR, t2RR, add,
                             add, false)

  # apply diffusion tensor
  applyDiffusionTensor_diff(diffusion, sbp, params, qL_el, wL, coordsL, dxidxL,
                            jacL, iface.elementL, sbpface,
                            iface.faceL, qL, t2LL, t2LR, t3LL, t3LR)
  # reverse RR and RL here because Lambda is a function of qR, not qL
  applyDiffusionTensor_diff(diffusion, sbp, params, qR_el, wR, coordsR, dxidxR,
                            jacR, iface.elementR, sbpface,
                            iface.faceR, qR, t2RR, t2RL, t3RR, t3RL)

  # apply inverse mass matrix (only nodes needed for interpolation later)
  @simd for q=1:numNodesPerElement
    # elementL
    @simd for p in nodesL
      facL = (1/alphas[1])*jacL[p]/sbp.w[p]
      @simd for d=1:dim
        @simd for j=1:numDofPerNode
          @simd for i=1:numDofPerNode
            t3LL[i, j, d, p, q] *= facL
            t3LR[i, j, d, p, q] *= facL
          end
        end
      end
    end
    # elementR
    @simd for p in nodesR
      facR = (1/alphas[2])*jacR[p]/sbp.w[p]
      @simd for d=1:dim
        @simd for j=1:numDofPerNode
          @simd for i=1:numDofPerNode
            t3RL[i, j, d, p, q] *= facR
            t3RR[i, j, d, p, q] *= facR
          end
        end
      end
    end
  end  # end q


  # interpolate back to the face
  fill!(t4LL, 0); fill!(t4RL, 0); fill!(t4LR, 0); fill!(t4RR, 0)
  interiorFaceInterpolate_jac!(sbpface, iface, t3LL, t3RL, t4LL, t4RL)
  interiorFaceInterpolate_jac!(sbpface, iface, t3LR, t3RR, t4LR, t4RR)

  # apply B * [Nx, Ny] and sum
  @simd for q=1:numNodesPerElement
    @simd for p=1:numNodesPerFace
      @simd for j=1:numDofPerNode
        @simd for i=1:numDofPerNode
          res1L_dotL[i, j, p, q] = 0
          res1L_dotR[i, j, p, q] = 0
          res1R_dotL[i, j, p, q] = 0
          res1R_dotR[i, j, p, q] = 0
        end
      end

      @simd for d=1:dim
        fac = sbpface.wface[p]*nrm_face[d, p]
        @simd for j=1:numDofPerNode
          @simd for i=1:numDofPerNode
            valL = fac*(t4LL[i, j, d, p, q] + t4RL[i, j, d, p, q])
            valR = fac*(t4LR[i, j, d, p, q] + t4RR[i, j, d, p, q])

            res1L_dotL[i, j, p, q] +=  valL
            res1L_dotR[i, j, p, q] +=  valR
            res1R_dotL[i, j, p, q] += -valL  # delta_w is reversed for elementR
            res1R_dotR[i, j, p, q] += -valR
          end
        end
      end
    end
  end

  return false
end

=#
