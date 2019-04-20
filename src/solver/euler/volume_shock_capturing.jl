# shock capturing using Q^T * epsilon * D * w

function calcShockCapturing(mesh::AbstractMesh, sbp::AbstractOperator,
                             eqn::EulerData, opts,
                             sensor::AbstractShockSensor,
                             capture::VolumeShockCapturing)

  computeVolumeTerm(mesh, sbp, eqn, opts, capture, capture.entropy_vars,
                    capture.diffusion)

  return nothing
end


@noinline function computeVolumeTerm(mesh, sbp, eqn, opts,
                      capture::VolumeShockCapturing{Tsol, Tres},
                      entropy_vars::AbstractVariables,
                      diffusion::AbstractDiffusion,
                     ) where {Tsol, Tres}

  @assert eqn.params.use_Minv != 1

  work = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  grad_w = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  lambda_gradw = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement,
                             mesh.dim)
  w_i = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerElement)
  op = SummationByParts.Subtract()

  # do local elements
  @simd for i=1:mesh.numEl
    i_full = i
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


#------------------------------------------------------------------------------
# Differentiated version

function calcShockCapturing_diff(mesh::AbstractMesh, sbp::AbstractOperator,
                                 eqn::EulerData, opts,
                                 sensor::AbstractShockSensor,
                                 capture::VolumeShockCapturing,
                                 assem::AssembleElementData)

  computeVolumeTerm_diff(mesh, sbp, eqn, opts, capture, capture.diffusion,
                         capture.entropy_vars, assem)

  return nothing
end

@noinline function computeVolumeTerm_diff(mesh, sbp, eqn, opts,
                           capture::VolumeShockCapturing{Tsol, Tres},
                           diffusion::AbstractDiffusion,
                           entropy_vars::AbstractVariables,
                           assembler::AssembleElementData) where {Tsol, Tres}

  @assert eqn.params.use_Minv != 1

  w_i = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerElement)
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

  for i=1:mesh.numEl
    i_full = i  # in case a list of elements is introduced later

    q_i = ro_sview(eqn.q, :, :, i_full)
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

#------------------------------------------------------------------------------
# revq

function calcShockCapturing_revq(mesh::AbstractMesh, sbp::AbstractOperator,
                             eqn::EulerData, opts,
                             sensor::AbstractShockSensor,
                             capture::VolumeShockCapturing)

  computeVolumeTerm_revq(mesh, sbp, eqn, opts, capture, capture.entropy_vars,
                    capture.diffusion)

  return nothing
end


@noinline function computeVolumeTerm_revq(mesh, sbp, eqn, opts,
                      capture::VolumeShockCapturing{Tsol, Tres},
                      entropy_vars::AbstractVariables,
                      diffusion::AbstractDiffusion,
                     ) where {Tsol, Tres}

  @assert eqn.params.use_Minv != 1

  work = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  grad_w = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  lambda_gradw = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement,
                             mesh.dim)
  w_i = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerElement)
  op = SummationByParts.Subtract()

  w_bar_i = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerElement)
  grad_w_bar = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  lambda_gradw_bar = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement,
                             mesh.dim)

  A0inv = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode)

  # do local elements
  @simd for i=1:mesh.numEl
    i_full = i
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
      q_j = ro_sview(eqn.q, :, j, i)
      q_bar_j = sview(eqn.q_bar, :, j, i)
      w_bar_j = ro_sview(w_bar_i, :, j)

      getA0inv(entropy_vars, eqn.params, q_j, A0inv)
      smallmatvec_kernel!(A0inv, w_bar_j, q_bar_j, 1, 1)
    end

  end  # end i

  return nothing
end

#------------------------------------------------------------------------------
# revm

function calcShockCapturing_revm(mesh::AbstractMesh, sbp::AbstractOperator,
                             eqn::EulerData, opts,
                             sensor::AbstractShockSensor,
                             capture::VolumeShockCapturing)

  computeVolumeTerm_revm(mesh, sbp, eqn, opts, capture, capture.entropy_vars,
                    capture.diffusion)

  return nothing
end


@noinline function computeVolumeTerm_revm(mesh, sbp, eqn, opts,
                      capture::VolumeShockCapturing{Tsol, Tres},
                      entropy_vars::AbstractVariables,
                      diffusion::AbstractDiffusion,
                     ) where {Tsol, Tres}

  @assert eqn.params.use_Minv != 1

  work = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  grad_w = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  lambda_gradw = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement,
                             mesh.dim)
  w_i = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerElement)
  op = SummationByParts.Subtract()

  w_bar_i = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerElement)
  grad_w_bar = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  lambda_gradw_bar = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement,
                             mesh.dim)

  # do local elements
  @simd for i=1:mesh.numEl
    i_full = i
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


