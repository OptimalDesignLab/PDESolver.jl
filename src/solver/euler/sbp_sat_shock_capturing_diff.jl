# differentiated version of sbp_sat_shock_capturing.jl

function calcShockCapturing_diff(mesh::AbstractMesh, sbp::AbstractOperator,
                             eqn::EulerData, opts,
                             capture::SBPParabolicSC{Tsol, Tres},
                             shockmesh::ShockedElements,
                             assem::AssembleElementData) where {Tsol, Tres}

  println("\nEntered calcShockCapturing_diff")
  computeGradW(mesh, sbp, eqn, opts, capture, shockmesh,
               capture.entropy_vars, capture.diffusion)

  computeVolumeTerm_diff(mesh, sbp, eqn, opts, capture,
                               capture.diffusion, capture.entropy_vars,
                               shockmesh, assem)

  computeFaceTerm_diff(mesh, sbp, eqn, opts, capture, shockmesh,
                       capture.diffusion, capture.entropy_vars,
                       capture.penalty, assem)

  if shockmesh.isNeumann
    println("computing Neumann boundary condition")
    computeNeumannBoundaryTerm_diff(mesh, sbp, eqn, opts, capture, shockmesh,
                              capture.diffusion, capture.entropy_vars, assem)
  else
    println("computing Dirichlet boundary condition")
    computeDirichletBoundaryTerm_diff(mesh, sbp, eqn, opts, capture, shockmesh,
                                      assem)
  end

  #@time computeSharedFaceTerm_diff(mesh, sbp, eqn, opts, capture, shockmesh,
  #                                 capture.diffusion, capture.penalty)



  return nothing
end


function computeVolumeTerm_diff(mesh, sbp, eqn, opts,
                           capture::SBPParabolicSC{Tsol, Tres},
                           diffusion::AbstractDiffusion,
                           entropy_vars::AbstractVariables,
                           shockmesh::ShockedElements,
                           assembler::AssembleElementData) where {Tsol, Tres}

  # computeGradW computes Lambda * D * w, so all that remains to do is
  # compute Qx * grad_q_x
  # Note that this term is not entropy stable by itself, because Qx was
  # not replaced by -Qx^T + Ex.  The entire discretization should be
  # entropy-stable however.
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
    w_i = ro_sview(capture.w_el, :, :, i)
    coords_i = ro_sview(mesh.coords, :, :, i_full)
    dxidx_i = ro_sview(mesh.dxidx, :, :, :, i_full)
    jac_i = ro_sview(mesh.jac, :, i_full)

    # to compute the jacobian, dq/dq = I, so dw/dq = the dw/dq evaluated at
    # each node
    for j=1:mesh.numNodesPerElement
      q_j = ro_sview(eqn.q, :, j, i_full)
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

    # apply Qx and sum
    if eqn.params.use_Minv != 1
      @simd for d=1:mesh.dim
        @simd for j=1:mesh.numNodesPerElement
          @simd for i=1:mesh.numNodesPerElement
            Dx[i, j, d] *= sbp.w[i]/jac_i[i]  # Dx is now Qx
          end
        end
      end
    end
    applyOperatorJac(Dx, t2_dot, res_jac)

#    calcQxTransposed(sbp, dxidx_i, Dx)
#    applyOperatorJac(Dx, t2_dot, res_jac, true, op)

    assembleElement(assembler, mesh, i_full, res_jac)
  end

  return nothing
end


function computeFaceTerm_diff(mesh, sbp, eqn, opts,
                      capture::SBPParabolicSC{Tsol, Tres},
                      shockmesh::ShockedElements, diffusion::AbstractDiffusion,
                      entropy_vars::AbstractVariables,
                      penalty::AbstractDiffusionPenalty,
                      assem::AssembleElementData) where {Tsol, Tres}

  delta_w = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  delta_w_dotL = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode,
                             mesh.numNodesPerFace, mesh.numNodesPerElement)
  delta_w_dotR = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode,
                             mesh.numNodesPerFace, mesh.numNodesPerElement)

  theta = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  theta_dotL = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode,
                           mesh.numNodesPerFace, mesh.numNodesPerElement)
  theta_dotR = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode,
                           mesh.numNodesPerFace, mesh.numNodesPerElement)

  res1L_dotL = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode,
                           mesh.numNodesPerFace, mesh.numNodesPerElement)
  res1L_dotR = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode,
                           mesh.numNodesPerFace, mesh.numNodesPerElement)

  res1R_dotL = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode,
                           mesh.numNodesPerFace, mesh.numNodesPerElement)
  res1R_dotR = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode,
                           mesh.numNodesPerFace, mesh.numNodesPerElement)

  res2L_dotL = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode,
                           mesh.numNodesPerFace, mesh.numNodesPerElement)
  res2L_dotR = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode,
                           mesh.numNodesPerFace, mesh.numNodesPerElement)

  res2R_dotL = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode,
                           mesh.numNodesPerFace, mesh.numNodesPerElement)
  res2R_dotR = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode,
                           mesh.numNodesPerFace, mesh.numNodesPerElement)


  data = eqn.params.calc_face_integrals_data
  @unpack data res_jacLL res_jacLR res_jacRL res_jacRR

  #t1L = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  #t1R = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  t2L = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  t2R = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  subtract = SummationByParts.Subtract()
  add = SummationByParts.Add()

  for i=1:shockmesh.numInterfaces
    iface_red = shockmesh.ifaces[i].iface
    iface_idx = shockmesh.ifaces[i].idx_orig
    elnumL = shockmesh.elnums_all[iface_red.elementL]
    elnumR = shockmesh.elnums_all[iface_red.elementR]

    # get data needed for next steps
    wL = ro_sview(capture.w_el, :, :, iface_red.elementL)
    wR = ro_sview(capture.w_el, :, :, iface_red.elementR)
    qL = ro_sview(eqn.q, :, :, elnumL)
    qR = ro_sview(eqn.q, :, :, elnumR)
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

    fill!(res_jacLL, 0); fill!(res_jacLR, 0)
    fill!(res_jacRL, 0); fill!(res_jacRR, 0)

    # compute delta w tilde and theta_bar = Dgk w_k + Dgn w_n
    getFaceVariables_diff(eqn.params, capture, diffusion, entropy_vars,
                          sbp, mesh.sbpface,
                          iface_red, wL, wR, qL, qR, gradwL, gradwR,
                          coordsL, coordsR,
                          nrm_face, dxidxL, dxidxR, jacL, jacR,
                          delta_w, delta_w_dotL, delta_w_dotR,
                          theta, theta_dotL, theta_dotR)

    # apply the penalty coefficient matrix
    has_T4 = applyPenalty_diff(penalty, sbp, eqn.params,  mesh.sbpface,
                      diffusion,
                      iface_red, delta_w, delta_w_dotL, delta_w_dotR,
                      theta, theta_dotL, theta_dotR,
                      qL, qR, wL, wR, coordsL, coordsR,  nrm_face, alphas,
                      dxidxL, dxidxR, jacL, jacR,
                      res1L_dotL, res1L_dotR, res1R_dotL, res1R_dotR,
                      t2L, t2R,
                      res2L_dotL, res2L_dotR, res2R_dotL, res2R_dotR)


    # apply Rgk^T, Rgn^T, Dgk^T, Dgn^T
    fill!(res_jacLL, 0); fill!(res_jacRL, 0)
    fill!(res_jacLR, 0); fill!(res_jacRR, 0)
    interiorFaceIntegrate_jac!(mesh.sbpface, iface_red, res1L_dotL, res1R_dotL,
                               res_jacLL, res_jacRL, subtract, subtract, false)
    interiorFaceIntegrate_jac!(mesh.sbpface, iface_red, res1L_dotR, res1R_dotR,
                               res_jacLR, res_jacRR, subtract, subtract, false)

    # apply Dgk^T and Dgn^T
    applyDgkTranspose_diff(capture, sbp, eqn.params, mesh.sbpface, iface_red,
                      diffusion,
                      t2L, res2L_dotL, res2L_dotR,
                      t2R, res2R_dotL, res2R_dotR,
                      qL, qR, wL, wR, coordsL, coordsR,
                      nrm_face, dxidxL, dxidxR, jacL, jacR,
                      res_jacLL, res_jacLR, res_jacRL, res_jacRR,
                      subtract)

    iface_full = replace_interface(iface_red, elnumL, elnumR)
    if !has_T4
      assembleInterfaceVisc(assem, mesh.sbpface, mesh, iface_full, res_jacLL,
                            res_jacLR, res_jacRL, res_jacRR)
    else
      #TODO: the sparsity pattern of the Jacobian may be incorrect for this case
      assembleInterfaceFull(assem, mesh, iface_full, res_jacLL,
                            res_jacLR, res_jacRL, res_jacRR)
    end


  end  # end loop i

  return nothing
end


function computeNeumannBoundaryTerm_diff(mesh::AbstractMesh{Tmsh}, sbp, eqn,
                      opts,
                      capture::SBPParabolicSC{Tsol, Tres},
                      shockmesh::ShockedElements,
                      diffusion::AbstractDiffusion,
                      entropy_vars::AbstractVariables,
                      assem::AssembleElementData
                      ) where {Tsol, Tres, Tmsh}

  @assert eqn.params.use_Minv != 1
  @assert mesh.coord_order == 1  # because of the normal vector

  # for shock capturing, apply the Neumann boundary condition
  # Lambda * grad_w = 0.  The resulting term is -R^T * B * Dgk * u

  op = SummationByParts.Subtract()

  @unpack capture w_dot Dx t1 t1_dot t2_dot t3_dot t4_dot

  res_jac = eqn.params.calc_face_integrals_data.res_jacLL
  nrm_i = zeros(Tmsh, mesh.dim, mesh.numNodesPerFace)

  for i=1:shockmesh.numBoundaryFaces
    bndry_i = shockmesh.bndryfaces[i].bndry
    idx_orig = shockmesh.bndryfaces[i].idx_orig
    elnum_orig = shockmesh.elnums_all[bndry_i.element]

    q_i = ro_sview(eqn.q, :, :, elnum_orig)
    w_i = ro_sview(capture.w_el, :, :, bndry_i.element)
    coords_i = ro_sview(mesh.coords, :, :, elnum_orig)
    dxidx_i = ro_sview(mesh.dxidx, :, :, :, elnum_orig)
    jac_i = ro_sview(mesh.jac, :, elnum_orig)
    if i < shockmesh.bndry_offsets[end]
      #nrm_i = ro_sview(mesh.nrm_bndry, :, :, idx_orig)
      for j=1:mesh.numNodesPerFace
        for d=1:mesh.dim
          nrm_i[d, j] = mesh.nrm_bndry[d, j, idx_orig]
        end
      end

    else
      #nrm_i = ro_sview(mesh.nrm_face, :, :, idx_orig)
      fac = shockmesh.bndryfaces[i].fac
      for j=1:mesh.numNodesPerFace
        for d=1:mesh.dim
          nrm_i[d, j] = fac*mesh.nrm_face[d, j, idx_orig]
        end
      end
    end

    #res_i = sview(eqn.res, :, :, elnum_orig)
    # compute dw/dq
    for j=1:mesh.numNodesPerElement
      q_j = ro_sview(eqn.q, :, j, elnum_orig)
      w_dot_j = sview(w_dot, :, :, j)
      getA0inv(entropy_vars, eqn.params, q_j, w_dot_j)
    end

    # multiply by Dgk
    calcDx(sbp, dxidx_i, jac_i, Dx)
    applyOperatorJac(Dx, mesh.sbpface, bndry_i.face, w_dot, t1_dot)

    # need t1 = [Dx, Dy] * w
    for d=1:mesh.dim
      Dx_d = ro_sview(Dx, :, :, d)
      t1_d = sview(t1, :, :, d)
      smallmatmatT!(w_i, Dx_d, t1_d)
    end

    applyDiffusionTensor_diff(diffusion, sbp, eqn.params, q_i, w_i, coords_i,
                              dxidx_i, jac_i, bndry_i.element, mesh.sbpface,
                              bndry_i.face, t1, t1_dot, t2_dot)

    # apply N * B * R and do the sum over dimensions
    fill!(t3_dot, 0)
    boundaryFaceInterpolate_jac!(mesh.sbpface, bndry_i.face, t2_dot, t3_dot)

    fill!(t4_dot, 0)  #TODO: tile this
    @simd for q=1:mesh.numNodesPerElement
      @simd for p=1:mesh.numNodesPerFace
        @simd for d=1:mesh.dim
          fac = nrm_i[d, p] #*mesh.sbpface.wface[p]
          @simd for j=1:mesh.numDofPerNode
            @simd for k=1:mesh.numDofPerNode
              t4_dot[k, j, p, q] += fac*t3_dot[k, j, d, p, q]
            end
          end
        end
      end
    end

    # apply R^T
    fill!(res_jac, 0)
    boundaryFaceIntegrate_jac!(mesh.sbpface, bndry_i.face, t4_dot, res_jac, op)
    bndry_orig = replace_boundary(bndry_i, elnum_orig)
    assembleBoundaryFull(assem, mesh, bndry_orig, res_jac)
  end  # end i

  return nothing
end


"""
  Differentiated Dirichlet boundary term
"""
function computeDirichletBoundaryTerm_diff(mesh, sbp, eqn, opts,
                      capture::SBPParabolicSC{Tsol, Tres},
                      shockmesh::ShockedElements,
                      assem::AssembleElementData
                      ) where {Tsol, Tres}

  for i=1:mesh.numBC
    bc_range = (shockmesh.bndry_offsets[i]:(shockmesh.bndry_offsets[i+1]-1))
    bc_func = mesh.bndry_funcs[i]

    calcBoundaryFlux_diff(mesh, sbp, eqn, opts, shockmesh, capture,
                     capture.penalty, capture.diffusion, 
                     capture.entropy_vars, bc_func, bc_range, assem)
  end

  return nothing
end


function calcBoundaryFlux_diff(mesh::AbstractMesh, sbp, eqn::EulerData, opts,
                      shockmesh::ShockedElements,
                      capture::SBPParabolicSC{Tsol, Tres},
                      penalty::AbstractDiffusionPenalty,
                      diffusion::AbstractDiffusion,
                      entropy_vars::AbstractVariables, bc_func::BCType,
                      bc_range::AbstractVector,
                      assem::AssembleElementData,
                      ) where {Tsol, Tres}


  delta_w = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  delta_w_dot = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode,
                            mesh.numNodesPerFace, mesh.numNodesPerElement)
  delta_w_dotR = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode,
                             mesh.numNodesPerFace, mesh.numNodesPerElement)
  res1 = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  res1_dot = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode,
                         mesh.numNodesPerFace, mesh.numNodesPerElement)

  res_dot = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode, 
                        mesh.numNodesPerElement, mesh.numNodesPerElement)
  res_dotR = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode, 
                        mesh.numNodesPerElement, mesh.numNodesPerElement)

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

    computeDirichletDelta_diff(capture, eqn.params, eqn.params_complex,
                          mesh.sbpface, bndry_i,
                          bc_func, entropy_vars, w_i, q_i, coords_i,
                          nrm_i, delta_w, delta_w_dot)

    applyDirichletPenalty_diff(penalty, sbp, mesh.sbpface, diffusion,
                          bndry_i, delta_w, delta_w_dot, w_i, nrm_i, alpha,
                          jacL, res1_dot)

    # apply R^T to T_D * delta_w
    fill!(res_dot, 0)
    boundaryFaceIntegrate_jac!(mesh.sbpface, bndry_i.face, res1_dot, res_dot,
                               op, include_quadrature=false)

    # apply B and then Dgk^T to delta_w
    for q=1:mesh.numNodesPerElement
      for p=1:mesh.numNodesPerFace
        for j=1:mesh.numDofPerNode
          for i=1:mesh.numDofPerNode
            delta_w_dot[i, j, p, q] *= -mesh.sbpface.wface[p]
          end
        end
      end
    end

    for p=1:mesh.numNodesPerFace
      for i=1:mesh.numDofPerNode
        delta_w[i,p] *= -mesh.sbpface.wface[p]
      end
    end

    # the dotR variables aren't needed for boundary terms
    applyDgkTranspose_diff(capture, sbp, eqn.params, mesh.sbpface, bndry_i,
                      diffusion,
                      delta_w, delta_w_dot, delta_w_dotR, q_i, w_i, coords_i,
                      nrm_i, dxidxL, jacL, res_dot, res_dotR, op)

    assembleElement(assem, mesh, elnum_orig, res_dot)
  end

  return nothing
end



"""
  Differentiated version of [`getFaceVariables`](@ref).
"""
function getFaceVariables_diff(params::ParamType,
                          capture::SBPParabolicSC{Tsol, Tres},
                          diffusion::AbstractDiffusion,
                          entropy_vars::AbstractVariables,
                          sbp::AbstractOperator,
                          sbpface::AbstractFace, iface_red::Interface,
                          wL::AbstractMatrix, wR::AbstractMatrix,
                          qL::AbstractMatrix, qR::AbstractMatrix,
                          gradwL::Abstract3DArray, gradwR::Abstract3DArray,
                          coordsL::AbstractMatrix, coordsR::AbstractMatrix,
                          nrm_face::AbstractMatrix,
                          dxidxL::Abstract3DArray, dxidxR::Abstract3DArray,
                          jacL::AbstractVector, jacR::AbstractVector,
                          delta_w::AbstractMatrix,
                          delta_w_dotL::Abstract4DArray, delta_w_dotR::Abstract4DArray,
                          theta::AbstractArray,
                          theta_dotL::Abstract4DArray, theta_dotR::Abstract4DArray
                         ) where {Tsol, Tres}

  numDofPerNode, numNodesPerElement = size(wL)
  numNodesPerFace = size(delta_w, 2)
  dim = size(gradwL, 3)

  getFaceVariables(capture, sbpface, iface_red, wL, wR, gradwL, gradwR,
                   nrm_face, delta_w, theta)

  @unpack capture wL_dot wR_dot Dx t1 t1_dot t2L_dot t2R_dot t3L_dot t3R_dot

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

  # derivative of theta
  # do volume operations for wL (derivative of wgradL)
  calcDx(sbp, dxidxL, jacL, Dx)
  applyOperatorJac(Dx, sbpface, iface_red.faceL, wL_dot, t1_dot)
  for d=1:dim
    Dx_d = ro_sview(Dx, :, :, d)
    t1_d = sview(t1, :, :, d)
    smallmatmatT!(wL, Dx_d, t1_d)
  end
  applyDiffusionTensor_diff(diffusion, sbp, params, qL, wL, coordsL, dxidxL,
                            jacL, iface_red.elementL, sbpface,
                            iface_red.faceL, t1, t1_dot, t2L_dot)

  # do volume operations for wR (derivative of wgradR)
  calcDx(sbp, dxidxR, jacR, Dx)
  applyOperatorJac(Dx, sbpface, iface_red.faceR, wR_dot, t1_dot)
  for d=1:dim
    Dx_d = ro_sview(Dx, :, :, d)
    t1_d = sview(t1, :, :, d)
    smallmatmatT!(wR, Dx_d, t1_d)
  end
  applyDiffusionTensor_diff(diffusion, sbp, params, qR, wR, coordsR, dxidxR,
                            jacR, iface_red.elementR, sbpface,
                            iface_red.faceR, t1, t1_dot, t2R_dot)

  # interpolate to face and multiply by normal vector
  #TODO: it would be interesting to supply a diagonal matrix to the 5D version
  #      of interiorFaceInterpolate_jac so it can do the reduction and therefore
  #      avoid a second temporary array
  fill!(t3L_dot, 0); fill!(t3R_dot, 0)
  interiorFaceInterpolate_jac!(sbpface, iface_red, t2L_dot, t2R_dot, t3L_dot, t3R_dot)
  fill!(theta_dotL, 0); fill!(theta_dotR, 0)  # TODO: tile this
  @simd for q=1:numNodesPerElement
    @simd for p=1:numNodesPerFace
      @simd for d=1:dim
        @simd for j=1:numDofPerNode
          @simd for i=1:numDofPerNode
            theta_dotL[i, j, p, q] += nrm_face[d, p]*t3L_dot[i, j, d, p, q]
            theta_dotR[i, j, p, q] -= nrm_face[d, p]*t3R_dot[i, j, d, p, q]
          end
        end
      end
    end
  end

  return nothing
end



function computeDirichletDelta_diff(capture::SBPParabolicSC{Tsol, Tres},
                               params::ParamType, params_c::ParamType,
                               sbpface::AbstractFace, bndry::Boundary,
                               func::BCType, entropy_vars::AbstractVariables,
                               wL::AbstractMatrix, qL::AbstractMatrix,
                               coords::AbstractMatrix, nrm_face::AbstractMatrix,
                               delta_w::AbstractMatrix,
                               delta_w_dot::Abstract4DArray
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

  # jacobian
  wL_dot = zeros(Tsol, numDofPerNode, numDofPerNode, numNodesPerElement)
#  wface_dot = zeros(Tsol, numDofPerNode, numDofPerNode, numNodesPerFace,
#                    numNodesPerElement)
  for i=1:numNodesPerElement
    wL_dot_i = sview(wL_dot, :, :, i)
    qL_i = ro_sview(qL, :, i)
    getA0inv(entropy_vars, params, qL_i, wL_dot_i)
  end
#  println("qL = \n", real(qL))
#  println("wL_dot = \n", real(wL_dot))
#  println("perm = ", sbpface.perm[:, bndry.face])
  fill!(delta_w_dot, 0)
  boundaryFaceInterpolate_jac!(sbpface, bndry.face, wL_dot, delta_w_dot)

  qface = zeros(Tsol, numDofPerNode, numNodesPerFace)
  qgface = zeros(Tres, numDofPerNode, numNodesPerFace)
  wgface = zeros(Tres, numDofPerNode)
  aux_vars = Tres[]

  # the inviscid boundary conditions interpolate q to the face and then
  # compute qg, so do that here as well for consistency
  boundaryFaceInterpolate!(sbpface, bndry.face, qL, qface)
  for i=1:numNodesPerFace
    q_i = sview(qface, :, i)
    qgface_i = sview(qgface, :, i)  # need this later
    coords_i = sview(coords, :, i)
    nrm_i = sview(nrm_face, :, i)

    getDirichletState(func, params, q_i, aux_vars, coords_i, nrm_i, qgface_i)

    convertToEntropy(entropy_vars, params, qgface_i, wgface)

    # compute the delta
    for j=1:numDofPerNode
      delta_w[j, i] = wface[j, i] - wgface[j]
    end
  end

  # jacobian
  # Compute dqg/qface first, then R^T * dqg/qface to get dqg/dqvolume

  op = SummationByParts.Subtract()
  qg_dot = zeros(Tsol, numDofPerNode, numDofPerNode)
  A0inv = zeros(Tsol, numDofPerNode, numDofPerNode)
  wg_dotface = zeros(Tsol, numDofPerNode, numDofPerNode, numNodesPerFace)
  q_c = zeros(Complex128, numDofPerNode)
  qg_c = zeros(Complex128, numDofPerNode)
  h = 1e-20
  pert = Complex128(0, h)
  for i=1:numNodesPerFace
    for j=1:numDofPerNode
      q_c[j] = qface[j, i]
    end
    coords_i = sview(coords, :, i)
    nrm_i = sview(nrm_face, :, i)

    # compute jacobian of qgface wrt qface
    for j=1:numDofPerNode
      q_c[j] += pert

      getDirichletState(func, params_c, q_c, aux_vars, coords_i, nrm_i, qg_c)
      for k=1:numDofPerNode
        qg_dot[k, j] = imag(qg_c[k])/h
      end

      q_c[j] -= pert
    end

    # jacobian of convertToEntropy
    getA0inv(entropy_vars, params, sview(qgface, :, i), A0inv)
    wg_dot_i = sview(wg_dotface, :, :, i)
    smallmatmat!(A0inv, qg_dot, wg_dot_i)
  end

  # apply R^T to wg_dotface to get d wg/dqL, and subtract that from delta_w_dot
  # interpolation is linear, so apply the regular operation repatedly
  # It would be better if SBP had a jacobian operation for this
  faceToVolume_jac!(sbpface, bndry.face, wg_dotface, delta_w_dot, op)

  return nothing
end




function applyDgkTranspose_diff(capture::SBPParabolicSC{Tsol, Tres}, sbp,
                       params::ParamType, sbpface, iface::Interface,
                       diffusion::AbstractDiffusion,
                       t2L::AbstractMatrix, t2L_dotL::Abstract4DArray,
                       t2L_dotR::Abstract4DArray,
                       t2R::AbstractMatrix, t2R_dotL::Abstract4DArray,
                       t2R_dotR::Abstract4DArray,
                       qL::AbstractMatrix, qR::AbstractMatrix,
                       wL::AbstractMatrix, wR::AbstractMatrix,
                       coordsL::AbstractMatrix, coordsR::AbstractMatrix,
                       nrm_face::AbstractMatrix,
                       dxidxL::Abstract3DArray, dxidxR::Abstract3DArray,
                       jacL::AbstractVector, jacR::AbstractVector,
                       resL_dotL::Abstract4DArray, resL_dotR::Abstract4DArray,
                       resR_dotL::Abstract4DArray, resR_dotR::Abstract4DArray,
                       op::SummationByParts.UnaryFunctor=SummationByParts.Add()
                      ) where {Tsol, Tres}

  dim, numNodesPerFace = size(nrm_face)
  numNodesPerElement = size(dxidxL, 3)
  numDofPerNode = size(wL, 1)
  add = SummationByParts.Add()
  subtract = SummationByParts.Subtract()

  @unpack capture temp1L temp1R temp2L temp2R
  fill!(temp2L, 0); fill!(temp2R, 0)
  @unpack capture t3L_dotL t3L_dotR t3R_dotL t3R_dotR Dx
  @unpack capture t4L_dotL t4L_dotR t4R_dotL t4R_dotR
  @unpack capture t5L_dotL t5L_dotR t5R_dotL t5R_dotR

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

  # multiply by N
  @simd for q=1:numNodesPerElement
    @simd for p=1:numNodesPerFace
      @simd for d=1:dim
        @simd for j=1:numDofPerNode
          @simd for i=1:numDofPerNode
            # normally t3R would have a - sign because the normal vector
            # is reversed for element nu, but interiorFaceIntegrate -= the
            # second argument, so don't add the minus sign here
            t3L_dotL[i, j, d, p, q] = nrm_face[d, p]*t2L_dotL[i, j, p, q]
            t3L_dotR[i, j, d, p, q] = nrm_face[d, p]*t2L_dotR[i, j, p, q]
            t3R_dotL[i, j, d, p, q] = nrm_face[d, p]*t2R_dotL[i, j, p, q]
            t3R_dotR[i, j, d, p, q] = nrm_face[d, p]*t2R_dotR[i, j, p, q]
          end
        end
      end
    end
  end

  # apply R^T
  fill!(t4L_dotL, 0); fill!(t4L_dotR, 0); fill!(t4R_dotL, 0); fill!(t4R_dotR, 0)
  interiorFaceIntegrate_jac!(sbpface, iface, t3L_dotL, t3R_dotL,
                             t4L_dotL, t4R_dotL, add, subtract, false)
  interiorFaceIntegrate_jac!(sbpface, iface, t3L_dotR, t3R_dotR,
                             t4L_dotR, t4R_dotR, add, subtract, false)

  # multiply by D^T Lambda
  applyDiffusionTensor_diff(diffusion, sbp, params, qL, wL, coordsL, dxidxL,
                            jacL, iface.elementL, temp2L,
                            t4L_dotL, t4L_dotR, t5L_dotL, t5L_dotR)
  # reverse L and R here so the d/dq term will be added to t5R_dotR
  applyDiffusionTensor_diff(diffusion, sbp, params, qR,  wR, coordsR, dxidxR,
                            jacR,iface.elementR, temp2R,
                            t4R_dotR, t4R_dotL, t5R_dotR, t5R_dotL)

  calcDxTransposed(sbp, dxidxL, jacL, Dx)
  applyOperatorJac(Dx, t5L_dotL, resL_dotL, false, op)
  applyOperatorJac(Dx, t5L_dotR, resL_dotR, false, op)

  calcDxTransposed(sbp, dxidxR, jacR, Dx)
  applyOperatorJac(Dx, t5R_dotL, resR_dotL, false, op)
  applyOperatorJac(Dx, t5R_dotR, resR_dotR, false, op)

  return nothing
end

"""
  Jacobian of Dgk * t2L only
"""
function applyDgkTranspose_diff(capture::SBPParabolicSC{Tsol, Tres}, sbp,
                       params::ParamType,
                       sbpface, iface::Union{Interface, Boundary},
                       diffusion::AbstractDiffusion,
                       t2L::AbstractMatrix, t2L_dotL::Abstract4DArray,
                       t2L_dotR::Abstract4DArray,
                       qL::AbstractMatrix,
                       wL::AbstractMatrix,
                       coordsL::AbstractMatrix,
                       nrm_face::AbstractMatrix,
                       dxidxL::Abstract3DArray,
                       jacL::AbstractVector,
                       resL_dotL::Abstract4DArray, resL_dotR::Abstract4DArray,
                       op::SummationByParts.UnaryFunctor=SummationByParts.Add()
                      ) where {Tsol, Tres}

  dim, numNodesPerFace = size(nrm_face)
  numNodesPerElement = size(dxidxL, 3)
  numDofPerNode = size(wL, 1)
  add = SummationByParts.Add()
  subtract = SummationByParts.Subtract()

  @unpack capture temp1L temp2L
  fill!(temp2L, 0);
  @unpack capture t3L_dotL t3L_dotR Dx
  @unpack capture t4L_dotL t4L_dotR
  @unpack capture t5L_dotL t5L_dotR

  # apply N and R^T
  @simd for d=1:dim
    @simd for j=1:numNodesPerFace
      @simd for k=1:numDofPerNode
        temp1L[k, j] =  nrm_face[d, j]*t2L[k, j]
      end
    end

    tmp2L = sview(temp2L, :, :, d)
    boundaryFaceInterpolate_rev!(sbpface, getFaceL(iface), tmp2L, temp1L)
  end

  # multiply by N
  @simd for q=1:numNodesPerElement
    @simd for p=1:numNodesPerFace
      @simd for d=1:dim
        @simd for j=1:numDofPerNode
          @simd for i=1:numDofPerNode
            t3L_dotL[i, j, d, p, q] = nrm_face[d, p]*t2L_dotL[i, j, p, q]
            t3L_dotR[i, j, d, p, q] = nrm_face[d, p]*t2L_dotR[i, j, p, q]
          end
        end
      end
    end
  end

  # apply R^T
  fill!(t4L_dotL, 0); fill!(t4L_dotR, 0)
  #TODO: get rid of keyword arguments
  boundaryFaceIntegrate_jac!(sbpface, getFaceL(iface), t3L_dotL, t4L_dotL,
                             include_quadrature=false)
  boundaryFaceIntegrate_jac!(sbpface, getFaceL(iface), t3L_dotR, t4L_dotR,
                             include_quadrature=false)


  # multiply by D^T Lambda
  applyDiffusionTensor_diff(diffusion, sbp, params, qL, wL, coordsL, dxidxL,
                            jacL, getElementL(iface), temp2L,
                            t4L_dotL, t4L_dotR, t5L_dotL, t5L_dotR)

  calcDxTransposed(sbp, dxidxL, jacL, Dx)
  applyOperatorJac(Dx, t5L_dotL, resL_dotL, false, op)
  applyOperatorJac(Dx, t5L_dotR, resL_dotR, false, op)

  return nothing
end


