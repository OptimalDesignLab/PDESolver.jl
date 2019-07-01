# shock capturing using a reduced from of the SBP-SAT 2nd derivative penalty
# Only works for diagonal E

import SummationByParts: UnaryFunctor, Add, Subtract


#------------------------------------------------------------------------------
# main functions

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

#  computeFaceTerm(mesh, sbp, eqn, opts, capture, shockmesh, capture.sensor_const,
#                  capture.penalty)

  computeFaceTerm(mesh, sbp, eqn, opts, capture, capture.entropy_vars,
                  shockmesh, capture.sensor_const, capture.penalty)



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


function calcShockCapturing_revq(mesh::AbstractMesh, sbp::AbstractOperator,
                             eqn::EulerData, opts,
                             capture::SBPParabolicReducedSC{Tsol, Tres},
                             shockmesh::ShockedElements) where {Tsol, Tres}

  computeVolumeTerm_revq(mesh, sbp, eqn, opts, capture, capture.entropy_vars,
                    capture.diffusion, shockmesh)

  #TODO: unneeded?
  finishConvertEntropy(mesh, sbp, eqn, opts, capture, capture.entropy_vars,
                       shockmesh)

  # set up the approximate shock sensor
  setShockedElements(capture.sensor_const, mesh, sbp, eqn, opts,
                     getShockSensor(capture), shockmesh)

  computeFaceTerm_revq(mesh, sbp, eqn, opts, capture, shockmesh,
                       capture.sensor_const, capture.entropy_vars,
                       capture.penalty)

  computeSharedFaceTerm_revq(mesh, sbp, eqn, opts, capture, shockmesh,
                             capture.sensor_const, capture.entropy_vars,
                             capture.penalty)


  return nothing
end


function calcShockCapturing_revm(mesh::AbstractMesh, sbp::AbstractOperator,
                             eqn::EulerData, opts,
                             capture::SBPParabolicReducedSC{Tsol, Tres},
                             shockmesh::ShockedElements) where {Tsol, Tres}

  computeVolumeTerm_revm(mesh, sbp, eqn, opts, capture, capture.entropy_vars,
                    capture.diffusion, shockmesh)

  finishConvertEntropy(mesh, sbp, eqn, opts, capture, capture.entropy_vars,
                       shockmesh)

  # set up the approximate shock sensor
  setShockedElements(capture.sensor_const, mesh, sbp, eqn, opts,
                     getShockSensor(capture), shockmesh)

  computeFaceTerm_revm(mesh, sbp, eqn, opts, capture, shockmesh,
                       capture.sensor_const, capture.penalty)

  computeSharedFaceTerm_revm(mesh, sbp, eqn, opts, capture, shockmesh,
                             capture.sensor_const, capture.penalty)

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


#-----------------------------------------------------------------------------
# Volume terms

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
      q_j = ro_sview(eqn.q, :, j, i_full)
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
    res_i = sview(eqn.res, :, :, i_full)
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

@noinline function computeVolumeTerm_revq(mesh, sbp, eqn, opts,
                      capture::SBPParabolicReducedSC{Tsol, Tres},
                      entropy_vars::AbstractVariables,
                      diffusion::AbstractDiffusion,
                      shockmesh::ShockedElements,
                     ) where {Tsol, Tres}

  @assert eqn.params.use_Minv != 1

  work = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  grad_w = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  lambda_gradw = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement,
                             mesh.dim)
  #w_i = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerElement)
  op = SummationByParts.Subtract()

  w_bar_i = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerElement)
  grad_w_bar = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  lambda_gradw_bar = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement,
                             mesh.dim)

  A0inv = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode)

  #TODO: verify usage of i_full
  # do local elements
  @simd for i=1:shockmesh.numShock
    i_full = shockmesh.elnums_all[i]
    w_i = sview(capture.w_el, :, :, i)
    @simd for j=1:mesh.numNodesPerElement
      q_j = ro_sview(eqn.q, :, j,i_full)  
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
#    applyDiffusionTensor(diffusion, sbp, eqn.params, q_i, w_i, coords, dxidx_i,
#                         jac_i, i, grad_w, lambda_gradw)

    # apply Q^T
#    res_i = sview(eqn.res, :, :, i)
#    applyQxTransposed(sbp, lambda_gradw, dxidx_i, work, res_i, op)


    # reverse sweep
    res_bar_i = ro_sview(eqn.res_bar, :, :, i_full)
    q_bar_i = sview(eqn.q_bar, :, :, i_full)
    fill!(w_bar_i, 0); fill!(lambda_gradw_bar, 0); fill!(grad_w_bar, 0)

    applyQxTransposed_revq(sbp, lambda_gradw_bar, dxidx_i, work, res_bar_i, op)

    applyDiffusionTensor_revq(diffusion, sbp, eqn.params, q_i, q_bar_i, w_i, w_bar_i,
                         coords, dxidx_i, jac_i, i, grad_w, grad_w_bar,
                         lambda_gradw, lambda_gradw_bar)

    applyDx_revq(sbp, w_bar_i, dxidx_i, jac_i, work, grad_w_bar)

    # convertToEntropy reverse mode
    for j=1:mesh.numNodesPerElement
      q_j = ro_sview(eqn.q, :, j, i_full)
      q_bar_j = sview(eqn.q_bar, :, j, i_full)
      w_bar_j = ro_sview(w_bar_i, :, j)

      getA0inv(entropy_vars, eqn.params, q_j, A0inv)
      smallmatvec_kernel!(A0inv, w_bar_j, q_bar_j, 1, 1)
    end

  end  # end i

  return nothing
end

@noinline function computeVolumeTerm_revm(mesh, sbp, eqn, opts,
                      capture::SBPParabolicReducedSC{Tsol, Tres},
                      entropy_vars::AbstractVariables,
                      diffusion::AbstractDiffusion,
                      shockmesh::ShockedElements,
                     ) where {Tsol, Tres}

  @assert eqn.params.use_Minv != 1

  work = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  grad_w = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  lambda_gradw = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement,
                             mesh.dim)
  op = SummationByParts.Subtract()

  w_bar_i = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerElement)
  grad_w_bar = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  lambda_gradw_bar = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement,
                             mesh.dim)

  # do local elements
  @simd for i=1:shockmesh.numShock
    i_full = shockmesh.elnums_all[i]
    w_i = sview(capture.w_el, :, :, i)
    @simd for j=1:mesh.numNodesPerElement
      q_j = ro_sview(eqn.q, :, j, i_full)
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
#    res_i = sview(eqn.res, :, :, i)
#    applyQxTransposed(sbp, lambda_gradw, dxidx_i, work, res_i, op)


    # reverse sweep
    res_bar_i = ro_sview(eqn.res_bar, :, :, i_full)
    dxidx_bar_i = sview(mesh.dxidx_bar, :, :, :, i_full)
    jac_bar_i = sview(mesh.jac_bar, :, i_full)
    fill!(w_bar_i, 0); fill!(lambda_gradw_bar, 0); fill!(grad_w_bar, 0)

    applyQxTransposed_revm(sbp, lambda_gradw, lambda_gradw_bar, dxidx_i,
                           dxidx_bar_i, work, res_bar_i, op)

    applyDiffusionTensor_revm(diffusion, sbp, eqn.params, q_i, w_i,
                         coords, dxidx_i, dxidx_bar_i, jac_i, jac_bar_i,
                         i, grad_w, grad_w_bar, lambda_gradw, lambda_gradw_bar)

    applyDx_revm(sbp, w_i, w_bar_i, dxidx_i, dxidx_bar_i, jac_i, jac_bar_i,
                 work, grad_w_bar)

    # no metrics involved in convert to entropy variables, don't need to 
    # compute reverse mode.
  end  # end i

  return nothing
end





#------------------------------------------------------------------------------
# Face terms
#=
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
=#
@noinline function computeFaceTerm(mesh, sbp, eqn, opts,
                      capture::SBPParabolicReducedSC{Tsol, Tres},
                      entropy_vars::AbstractVariables,
                      shockmesh::ShockedElements, sensor::ShockSensorHApprox,
                      penalty::AbstractDiffusionPenalty) where {Tsol, Tres}

  op = SummationByParts.Subtract()

  wL = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerElement)
  wR = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerElement)

  for i=1:mesh.numInterfaces
    iface = mesh.interfaces[i]
    elnumL = iface.elementL
    elnumR = iface.elementR

    for j=1:mesh.numNodesPerElement
      qL_j = ro_sview(eqn.q, :, j, elnumL)
      qR_j = ro_sview(eqn.q, :, j, elnumR)
      wL_j = sview(wL, :, j)
      wR_j = sview(wR, :, j)
      convertToEntropy(entropy_vars, eqn.params, qL_j, wL_j)
      convertToEntropy(entropy_vars, eqn.params, qR_j, wR_j)
    end

    nrm_face = ro_sview(mesh.nrm_face, :, :, i)
    jacL = ro_sview(mesh.jac, :, elnumL)
    jacR = ro_sview(mesh.jac, :, elnumR)
    resL = sview(eqn.res, :, :, elnumL)
    resR = sview(eqn.res, :, :, elnumR)

    applyReducedPenalty(penalty, sbp, eqn.params, mesh.sbpface, sensor,
                        iface, wL, wR, nrm_face, jacL, jacR, resL, resR,
                        op)
  end  # end loop i

  return nothing
end



#=
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
=#

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

  for i=1:mesh.numInterfaces
    iface = mesh.interfaces[i]
    elnumL = iface.elementL
    elnumR = iface.elementR

    # get data needed for next steps
    qL = ro_sview(eqn.q, :, :, elnumL)
    qR = ro_sview(eqn.q, :, :, elnumR)

    nrm_face = ro_sview(mesh.nrm_face, :, :, i)
    jacL = ro_sview(mesh.jac, :, elnumL)
    jacR = ro_sview(mesh.jac, :, elnumR)

    # apply the penalty coefficient matrix
    applyReducedPenalty_diff(penalty, sbp, eqn.params,  mesh.sbpface,
                      sensor, entropy_vars, iface, qL, qR,  nrm_face,
                      jacL, jacR,
                      res_jacLL, res_jacLR, res_jacRL, res_jacRR, op)

#    println("assembling shock capturing interface ", iface_full)
    assembleInterface(assem, mesh.sbpface, mesh, iface, res_jacLL,
                      res_jacLR, res_jacRL, res_jacRR)
  end  # end loop i

  return nothing
end



@noinline function computeFaceTerm_revq(mesh, sbp, eqn, opts,
                      capture::SBPParabolicReducedSC{Tsol, Tres},
                      shockmesh::ShockedElements, sensor::ShockSensorHApprox,
                      entropy_vars::AbstractVariables,
                      penalty::AbstractDiffusionPenalty) where {Tsol, Tres}

  op = SummationByParts.Subtract()

  A0invL = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode)
  A0invR = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode)
  wL_bar = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)
  wR_bar = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)
  for i=1:shockmesh.numInterfaces
    iface_red = shockmesh.ifaces[i].iface
    iface_idx = shockmesh.ifaces[i].idx_orig
    elnumL = shockmesh.elnums_all[iface_red.elementL]
    elnumR = shockmesh.elnums_all[iface_red.elementR]

    # get data needed for next steps
    #wL = ro_sview(capture.w_el, :, :, iface_red.elementL)
    #wR = ro_sview(capture.w_el, :, :, iface_red.elementR)

    nrm_face = ro_sview(mesh.nrm_face, :, :, iface_idx)
    jacL = ro_sview(mesh.jac, :, elnumL)
    jacR = ro_sview(mesh.jac, :, elnumR)
    resL_bar = sview(eqn.res_bar, :, :, elnumL)
    resR_bar = sview(eqn.res_bar, :, :, elnumR)

    # applyReducedPenalty is linear wrt q and symmetric, so we don't need to
    # write a separate applyReducedPenalty_revq, just call it with w and res
    # reversed and replaced with the _bar version
    fill!(wL_bar, 0); fill!(wR_bar, 0)
    applyReducedPenalty(penalty, sbp, eqn.params, mesh.sbpface, sensor,
                        iface_red, resL_bar, resR_bar, nrm_face, jacL, jacR,
                        wL_bar, wR_bar, op)

    # convertToEntropy reverse mode
    for j=1:mesh.numNodesPerFace
      permL = mesh.sbpface.perm[j, iface_red.faceL]
      permR = mesh.sbpface.perm[j, iface_red.faceR]

      qL_j = ro_sview(eqn.q, :, permL, elnumL)
      qR_j = ro_sview(eqn.q, :, permR, elnumR)
      qL_bar_j = sview(eqn.q_bar, :, permL, elnumL)
      qR_bar_j = sview(eqn.q_bar, :, permR, elnumR)
      wL_bar_j = ro_sview(wL_bar, :, permL)
      wR_bar_j = ro_sview(wR_bar, :, permR)

      getA0inv(entropy_vars, eqn.params, qL_j, A0invL)
      getA0inv(entropy_vars, eqn.params, qR_j, A0invR)

      smallmatvec_kernel!(A0invL, wL_bar_j, qL_bar_j, 1, 1)
      smallmatvec_kernel!(A0invR, wR_bar_j, qR_bar_j, 1, 1)
    end
  end  # end loop i

  return nothing
end


@noinline function computeFaceTerm_revm(mesh, sbp, eqn, opts,
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
    nrm_face_bar = sview(mesh.nrm_face_bar, :, :, iface_idx)
    jacL = ro_sview(mesh.jac, :, elnumL)
    jacL_bar = sview(mesh.jac_bar, :, elnumL)
    jacR = ro_sview(mesh.jac, :, elnumR)
    jacR_bar = sview(mesh.jac_bar, :, elnumR)
    resL_bar = sview(eqn.res_bar, :, :, elnumL)
    resR_bar = sview(eqn.res_bar, :, :, elnumR)

    applyReducedPenalty_revm(penalty, sbp, eqn.params, mesh.sbpface, sensor,
                        iface_red, wL, wR, nrm_face, nrm_face_bar, 
                        jacL, jacL_bar, jacR, jacR_bar, resL_bar, resR_bar, op)
  end  # end loop i

  return nothing
end


#------------------------------------------------------------------------------
# Shared face terms

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


@noinline function computeSharedFaceTerm_revq(mesh, sbp, eqn, opts,
                      capture::SBPParabolicReducedSC{Tsol, Tres},
                      shockmesh::ShockedElements, sensor::ShockSensorHApprox,
                      entropy_vars::AbstractVariables,
                      penalty::AbstractDiffusionPenalty) where {Tsol, Tres}


  op = SummationByParts.Subtract()

  # don't care about resR for shared faces
  A0invL = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode)
  wL_bar = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)
  wR_bar = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)
  for peer=1:shockmesh.npeers

    peer_full = shockmesh.peer_indices[peer]
    metrics = mesh.remote_metrics[peer_full]
    data = eqn.shared_data[peer_full]
    data_bar = eqn.shared_data_bar[peer_full]

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
      resL_bar = ro_sview(eqn.res_bar, :, :, elnumL)
      resR_bar = ro_sview(data_bar.q_recv, :, :, elnumR)

      fill!(wL_bar, 0); fill!(wR_bar, 0)
      applyReducedPenalty(penalty, sbp, eqn.params, mesh.sbpface, sensor,
                        iface_red, resL_bar, resR_bar, nrm_face, jacL, jacR,
                        wL_bar, wR_bar, op)

      # convertToEntropy reverse mode
      for j=1:mesh.numNodesPerFace
        permL = mesh.sbpface.perm[j, iface_red.faceL]

        qL_j = ro_sview(eqn.q, :, permL, elnumL)
        qL_bar_j = sview(eqn.q_bar, :, permL, elnumL)
        wL_bar_j = ro_sview(wL_bar, :, permL)

        getA0inv(entropy_vars, eqn.params, qL_j, A0invL)
        smallmatvec_kernel!(A0invL, wL_bar_j, qL_bar_j, 1, 1)
      end
    end  # end i
  end  # end peer

  return nothing
end


@noinline function computeSharedFaceTerm_revm(mesh, sbp, eqn, opts,
                      capture::SBPParabolicReducedSC{Tsol, Tres},
                      shockmesh::ShockedElements, sensor::ShockSensorHApprox,
                      penalty::AbstractDiffusionPenalty) where {Tsol, Tres}


  op = SummationByParts.Subtract()

  # the primal function does (metricsL, metricsR) -> resL, so the
  # reverse mode is resL -> (metricsL, metricsR)
  resR_bar = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)
  for peer=1:shockmesh.npeers
    peer_full = shockmesh.peer_indices[peer]
    metrics = mesh.remote_metrics[peer_full]
    metrics_bar = mesh.remote_metrics_bar[peer_full]
    data = eqn.shared_data[peer_full]
    #data_bar = eqn.shared_data_bar[peer_full]

    for i=1:shockmesh.numSharedInterfaces[peer]
      iface_red = shockmesh.shared_interfaces[peer][i].iface
      iface_idx = shockmesh.shared_interfaces[peer][i].idx_orig
      elnumL = shockmesh.elnums_all[iface_red.elementL]
      elnumR = getSharedElementIndex(shockmesh, mesh, peer, iface_red.elementR)

      # get data needed for next steps
      wL = ro_sview(capture.w_el, :, :, iface_red.elementL)
      wR = ro_sview(capture.w_el, :, :, iface_red.elementR)

      nrm_face = ro_sview(mesh.nrm_sharedface[peer_full], :, :, iface_idx)
      nrm_face_bar = sview(mesh.nrm_sharedface_bar[peer_full], :, :, iface_idx)
      jacL     = ro_sview(mesh.jac, :, elnumL)
      jacL_bar = sview(mesh.jac_bar, :, elnumL) 
      jacR     = ro_sview(metrics.jac, :, elnumR)
      jacR_bar = sview(metrics_bar.jac, :, elnumR)
      resL_bar = ro_sview(eqn.res_bar, :, :, elnumL)
#      resR_bar = ro_sview(data_bar.q_recv, :, :, elnumR)

      applyReducedPenalty_revm(penalty, sbp, eqn.params, mesh.sbpface, sensor,
                        iface_red, wL, wR, nrm_face, nrm_face_bar, 
                        jacL, jacL_bar, jacR, jacR_bar, resL_bar, resR_bar, op)
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


function applyReducedPenalty_revm(penalty::BR2Penalty{Tsol, Tres}, sbp,
                      params::AbstractParamType{Tdim}, sbpface::SparseFace,
                      sensor::ShockSensorHApprox, iface::Interface,
                      wL::AbstractMatrix, wR::AbstractMatrix,
                      nrm_face::AbstractMatrix,
                      nrm_face_bar::AbstractMatrix,
                      jacL::AbstractVector{Tmsh}, jacL_bar::AbstractVector,
                      jacR::AbstractVector, jacR_bar::AbstractVector,
                      resL_bar::AbstractMatrix, resR_bar::AbstractMatrix,
                      op::UnaryFunctor=Add()
                     ) where {Tsol, Tres, Tmsh, Tdim}

  numDofPerNode = size(wL, 1)
  numNodesPerFace = size(nrm_face, 2)
  numNodesPerElement = length(jacL)

  eps_L = getShockSensor(params, sbp, sensor, iface.elementL, jacL)
  eps_R = getShockSensor(params, sbp, sensor, iface.elementR, jacR)


  alphaL = eps_L/4; alphaL_bar = zero(Tres)
  alphaR = eps_R/4; alphaR_bar = zero(Tres)
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
    #fac = op(alphaL*facL + alphaR*facR)

    #-------------------------------------------
    # reverse sweep

    fac_bar = zero(Tres)
    @simd for j=1:numDofPerNode
      delta_w_j = wL[j, iL] - wR[j, iR]
      #res1L[j, iL] += fac*delta_w_j
      #res1R[j, iR] -= fac*delta_w_j  # delta_w is reversed for elementR
      fac_bar += delta_w_j*(resL_bar[j, iL] - resR_bar[j, iR])
    end  # end j

    alphaL_bar += op(facL*fac_bar); alphaR_bar += op(facR*fac_bar)
    facL_bar    = op(alphaL*fac_bar); facR_bar  = op(alphaR*fac_bar)

    HinvL_bar = zero(Tres); HinvR_bar = zero(Tres)
    @simd for d=1:Tdim
      n_d = nrm_face[d, i]
      N2_d = n_d*n_d
      nrm_face_bar[d, i] += 2*B2_i*n_d*HinvL*facL_bar
      nrm_face_bar[d, i] += 2*B2_i*n_d*HinvR*facR_bar
      HinvL_bar += B2_i*N2_d*facL_bar
      HinvR_bar += B2_i*N2_d*facR_bar
    end  # end d

    jacL_bar[iL] += HinvL_bar/sbp.w[iL]
    jacR_bar[iR] += HinvR_bar/sbp.w[iR]

  end  # end i

  epsL_bar = alphaL_bar/4
  epsR_bar = alphaR_bar/4

  getShockSensor_revm(params, sbp, sensor, iface.elementL, jacL, jacL_bar,
                      epsL_bar)
  getShockSensor_revm(params, sbp, sensor, iface.elementR, jacR, jacR_bar,
                      epsR_bar)

  return nothing
end


