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

  # set up the approximate shock sensor
  setShockedElements(capture.sensor_const, mesh, sbp, eqn, opts,
                     getShockSensor(capture), shockmesh)

  computeFaceTerm(mesh, sbp, eqn, opts, capture, shockmesh, capture.sensor_const,
                  capture.penalty)

  computeSharedFaceTerm(mesh, sbp, eqn, opts, capture, shockmesh,
                        capture.sensor_const, capture.penalty)


  return nothing
end

function calcShockCapturing_diff(mesh::AbstractMesh, sbp::AbstractOperator,
                             eqn::EulerData, opts,
                             capture::SBPParabolicReducedSC{Tsol, Tres},
                             shockmesh::ShockedElements,
                             assem::AssembleElementData) where {Tsol, Tres}

  computeVolumeTerm_diff(mesh, sbp, eqn, opts, capture,
                               capture.diffusion, capture.entropy_vars,
                               shockmesh, assem)
  finishConvertEntropy(mesh, sbp, eqn, opts, capture, capture.entropy_vars,
                       shockmesh)

  # set up the approximate shock sensor
  setShockedElements(capture.sensor_const, mesh, sbp, eqn, opts,
                     getShockSensor(capture), shockmesh)

  computeFaceTerm_diff(mesh, sbp, eqn, opts, capture, shockmesh,
                       capture.sensor_const, capture.entropy_vars,
                       capture.penalty, assem)

  computeSharedFaceTerm_diff(mesh, sbp, eqn, opts, capture, shockmesh,
                                   capture.sensor_const, capture.entropy_vars,
                                   capture.penalty, assem)



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
    for j=1:mesh.numNodesPerElement #j in mesh.sbpface.perm  # only do face nodes
      q_j = ro_sview(eqn.q, :, j, i_full)
      w_j = sview(capture.w_el, :, j, i)
      convertToEntropy(entropy_vars, eqn.params, q_j, w_j)
    end
  end

  for peer=1:shockmesh.npeers
    peer_full = shockmesh.peer_indices[peer]
    data = eqn.shared_data[peer_full]

    for i in shockmesh.shared_els[peer]
      i_full = getSharedElementIndex(shockmesh, mesh, peer, i)

      # compute entropy variables
      for j=1:mesh.numNodesPerElement # in mesh.sbpface.perm
        q_j = ro_sview(data.q_recv, :, j, i_full)
        w_j = sview(capture.w_el, :, j, i)
        convertToEntropy(entropy_vars, eqn.params, q_j, w_j)
      end
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

  op = SummationByParts.Subtract()

  for i=1:shockmesh.numInterfaces
    iface_red = shockmesh.ifaces[i].iface
    iface_idx = shockmesh.ifaces[i].idx_orig
    elnumL = shockmesh.elnums_all[iface_red.elementL]
    elnumR = shockmesh.elnums_all[iface_red.elementR]

    # get data needed for next steps
    wL = ro_sview(capture.w_el, :, :, iface_red.elementL)
    wR = ro_sview(capture.w_el, :, :, iface_red.elementR)

    nrm_face = ro_sview(mesh.nrm_face, :, :, iface_idx)
    jacL = ro_sview(mesh.jac, :, elnumL)
    jacR = ro_sview(mesh.jac, :, elnumR)
    resL = sview(eqn.res, :, :, elnumL)
    resR = sview(eqn.res, :, :, elnumR)

    applyReducedPenalty(penalty, sbp, eqn.params, mesh.sbpface, sensor,
                        iface_red, wL, wR, nrm_face, jacL, jacR, resL, resR,
                        op)
  end  # end loop i

  return nothing
end


@noinline function computeFaceTerm_diff(mesh, sbp, eqn, opts,
                      capture::SBPParabolicReducedSC{Tsol, Tres},
                      shockmesh::ShockedElements, sensor::ShockSensorHApprox,
                      entropy_vars::AbstractVariables,
                      penalty::AbstractDiffusionPenalty,
                      assem::AssembleElementData) where {Tsol, Tres}

#  println("\n\nDoing shock capturing face integrals")

  data = eqn.params.calc_face_integrals_data
  @unpack data res_jacLL res_jacLR res_jacRL res_jacRR

  op = Subtract()

  for i=1:shockmesh.numInterfaces
    iface_red = shockmesh.ifaces[i].iface
    iface_idx = shockmesh.ifaces[i].idx_orig
    elnumL = shockmesh.elnums_all[iface_red.elementL]
    elnumR = shockmesh.elnums_all[iface_red.elementR]
    iface_full = replace_interface(iface_red, elnumL, elnumR)

    # get data needed for next steps
    qL = ro_sview(eqn.q, :, :, elnumL)
    qR = ro_sview(eqn.q, :, :, elnumR)

    nrm_face = ro_sview(mesh.nrm_face, :, :, iface_idx)
    jacL = ro_sview(mesh.jac, :, elnumL)
    jacR = ro_sview(mesh.jac, :, elnumR)

    # apply the penalty coefficient matrix
    applyReducedPenalty_diff(penalty, sbp, eqn.params,  mesh.sbpface,
                      sensor, entropy_vars, iface_red, qL, qR,  nrm_face,
                      jacL, jacR,
                      res_jacLL, res_jacLR, res_jacRL, res_jacRR, op)

#    println("assembling shock capturing interface ", iface_full)
    assembleInterface(assem, mesh.sbpface, mesh, iface_full, res_jacLL,
                          res_jacLR, res_jacRL, res_jacRR)
  end  # end loop i

  return nothing
end



"""
  Does the same thing as `compuateFaceTerm`, but for the shared faces, updating
  the residual on the local element only
"""
@noinline function computeSharedFaceTerm(mesh, sbp, eqn, opts,
                      capture::SBPParabolicReducedSC{Tsol, Tres},
                      shockmesh::ShockedElements, sensor::ShockSensorHApprox,
                      penalty::AbstractDiffusionPenalty) where {Tsol, Tres}


  op = SummationByParts.Subtract()

  # don't care about resR for shared faces
  resR = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)
  for peer=1:shockmesh.npeers

    peer_full = shockmesh.peer_indices[peer]
    metrics = mesh.remote_metrics[peer_full]
    data = eqn.shared_data[peer_full]

    for i=1:shockmesh.numSharedInterfaces[peer]
      iface_red = shockmesh.shared_interfaces[peer][i].iface
      iface_idx = shockmesh.shared_interfaces[peer][i].idx_orig
      elnumL = shockmesh.elnums_all[iface_red.elementL]
      elnumR = getSharedElementIndex(shockmesh, mesh, peer, iface_red.elementR)

      # get data needed for next steps
      wL = ro_sview(capture.w_el, :, :, iface_red.elementL)
      wR = ro_sview(capture.w_el, :, :, iface_red.elementR)

      nrm_face = ro_sview(mesh.nrm_sharedface[peer_full], :, :, iface_idx)
      jacL = ro_sview(mesh.jac,    :, elnumL)
      jacR = ro_sview(metrics.jac, :, elnumR)
      resL = sview(eqn.res, :, :, elnumL)

      resL_orig = eqn.res[:, :, elnumL]

      applyReducedPenalty(penalty, sbp, eqn.params, mesh.sbpface, sensor,
                        iface_red, wL, wR, nrm_face, jacL, jacR, resL, resR,
                        op)
    end  # end i
  end  # end peer

  return nothing
end


@noinline function computeSharedFaceTerm_diff(mesh, sbp, eqn, opts,
                      capture::SBPParabolicReducedSC{Tsol, Tres},
                      shockmesh::ShockedElements, sensor::ShockSensorHApprox,
                      entropy_vars::AbstractVariables,
                      penalty::AbstractDiffusionPenalty,
                      assem::AssembleElementData
                     ) where {Tsol, Tres}

  data = eqn.params.calc_face_integrals_data
  @unpack data res_jacLL res_jacLR res_jacRL res_jacRR

  # do this once at the beginning only, because applyReducedPenalty_diff
  # overwrites the relevent entries
  fill!(res_jacLL, 0); fill!(res_jacRL, 0)
  fill!(res_jacLR, 0); fill!(res_jacRR, 0)

  op = SummationByParts.Subtract()

  # don't care about resR for shared faces
  resR = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)
  for peer=1:shockmesh.npeers

    peer_full = shockmesh.peer_indices[peer]
    metrics = mesh.remote_metrics[peer_full]
    data = eqn.shared_data[peer_full]

    for i=1:shockmesh.numSharedInterfaces[peer]
      iface_red = shockmesh.shared_interfaces[peer][i].iface
      iface_idx = shockmesh.shared_interfaces[peer][i].idx_orig
      elnumL = shockmesh.elnums_all[iface_red.elementL]
      elnumR = getSharedElementIndex(shockmesh, mesh, peer, iface_red.elementR)
      iface_full = replace_interface(iface_red, elnumL,
                        mesh.shared_interfaces[peer_full][iface_idx].elementR)

      # get data needed for next steps
      qL = ro_sview(eqn.q, :, :, elnumL)
      qR = ro_sview(data.q_recv, :, :, elnumR)

      nrm_face = ro_sview(mesh.nrm_sharedface[peer_full], :, :, iface_idx)
      jacL = ro_sview(mesh.jac,    :, elnumL)
      jacR = ro_sview(metrics.jac, :, elnumR)


      # apply the penalty coefficient matrix
      applyReducedPenalty_diff(penalty, sbp, eqn.params,  mesh.sbpface,
                        sensor, entropy_vars, iface_red, qL, qR,  nrm_face,
                        jacL, jacR,
                        res_jacLL, res_jacLR, res_jacRL, res_jacRR, op)

      assembleSharedFace(assem, mesh.sbpface, mesh, iface_full, res_jacLL, res_jacLR)
    end  # end i
  end  # end peer

  return nothing
end




#------------------------------------------------------------------------------
# Penalty functions

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

#  println("epsL = ", real(eps_L), ", eps_R = ", real(eps_R))

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

