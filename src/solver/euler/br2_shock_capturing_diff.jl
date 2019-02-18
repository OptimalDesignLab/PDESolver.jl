# differentiated version of br2_shock_capturing.jl

function calcShockCapturing_diff(mesh::AbstractMesh, sbp::AbstractOperator,
                             eqn::EulerData, opts,
                             capture::SBPParabolicSC{Tsol, Tres},
                             shockmesh::ShockedElements,
                             assem::AssembleElementData) where {Tsol, Tres}


  println("entered calcShockCapturing_diff")
  computeGradW(mesh, sbp, eqn, opts, capture, shockmesh,
               capture.convert_entropy, capture.diffusion)

  computeVolumeTerm_diff(mesh, sbp, eqn, opts, capture,
                               capture.diffusion, shockmesh, assem)

  computeFaceTerm_diff(mesh, sbp, eqn, opts, capture, shockmesh, capture.diffusion,
                  capture.penalty, assem)

  computeBoundaryTerm_diff(mesh, sbp, eqn, opts, capture, shockmesh,
                            capture.diffusion, assem)

  #@time computeSharedFaceTerm(mesh, sbp, eqn, opts, capture, shockmesh,
  #                            capture.diffusion, capture.penalty)



  return nothing
end


function computeVolumeTerm_diff(mesh, sbp, eqn, opts,
                           capture::SBPParabolicSC{Tsol, Tres},
                           diffusion::AbstractDiffusion,
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

  for i=1:shockmesh.numShock
    i_full = shockmesh.elnums_all[i]

    w_i = ro_sview(capture.w_el, :, :, i)
    dxidx_i = ro_sview(mesh.dxidx, :, :, :, i_full)
    jac_i = ro_sview(mesh.jac, :, i_full)

    # to compute the jacobian, dq/dq = I, so dw/dq = the dw/dq evaluated at
    # each node
    for j=1:mesh.numNodesPerElement
      q_j = ro_sview(eqn.q, :, j, i_full)
      w_dot_j = sview(w_dot, :, :, j)
      getIRA0inv(eqn.params, q_j, w_dot_j)
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
    applyDiffusionTensor_diff(diffusion, w_i, i, gradq_i, t1_dot, t2_dot)

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

    assembleElement(assembler, mesh, i_full, res_jac)
  end

  return nothing
end


function computeFaceTerm_diff(mesh, sbp, eqn, opts,
                      capture::SBPParabolicSC{Tsol, Tres},
                      shockmesh::ShockedElements, diffusion::AbstractDiffusion,
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
    getFaceVariables_diff(eqn.params, capture, diffusion, sbp, mesh.sbpface,
                          iface_red, wL, wR, qL, qR, gradwL, gradwR,
                          nrm_face, dxidxL, dxidxR, jacL, jacR,
                          delta_w, delta_w_dotL, delta_w_dotR,
                          theta, theta_dotL, theta_dotR)

    # apply the penalty coefficient matrix
    has_T4 = applyPenalty_diff(penalty, sbp, mesh.sbpface, diffusion,
                      iface_red, delta_w, delta_w_dotL, delta_w_dotR,
                      theta, theta_dotL, theta_dotR,
                      wL, wR, nrm_face, alphas, jacL, jacR,
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
    applyDgkTranspose_diff(capture, sbp, mesh.sbpface, iface_red, diffusion,
                      t2L, res2L_dotL, res2L_dotR,
                      t2R, res2R_dotL, res2R_dotR,
                      wL, wR, nrm_face, dxidxL, dxidxR, jacL, jacR,
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


function computeBoundaryTerm_diff(mesh, sbp, eqn, opts,
                      capture::SBPParabolicSC{Tsol, Tres},
                      shockmesh::ShockedElements,
                      diffusion::AbstractDiffusion,
                      assem::AssembleElementData
                      ) where {Tsol, Tres}

  @assert eqn.params.use_Minv != 1

  # for shock capturing, apply the Neumann boundary condition
  # Lambda * grad_w = 0.  The resulting term is -R^T * B * Dgk * u

  op = SummationByParts.Subtract()

  @unpack capture w_dot Dx t1 t1_dot t2_dot t3_dot t4_dot

  res_jac = eqn.params.calc_face_integrals_data.res_jacLL

  for i=1:shockmesh.numBoundaryFaces
    bndry_i = shockmesh.bndryfaces[i].bndry
    idx_orig = shockmesh.bndryfaces[i].idx_orig
    elnum_orig = shockmesh.elnums_all[bndry_i.element]

    w_i = ro_sview(capture.w_el, :, :, bndry_i.element)
    dxidx_i = ro_sview(mesh.dxidx, :, :, :, elnum_orig)
    jac_i = ro_sview(mesh.jac, :, elnum_orig)
    nrm_i = ro_sview(mesh.nrm_bndry, :, :, idx_orig)
    #res_i = sview(eqn.res, :, :, elnum_orig)
    # compute dw/dq
    for j=1:mesh.numNodesPerElement
      q_j = ro_sview(eqn.q, :, j, elnum_orig)
      w_dot_j = sview(w_dot, :, :, j)
      getIRA0inv(eqn.params, q_j, w_dot_j)
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

    applyDiffusionTensor_diff(diffusion, w_i, bndry_i.element, mesh.sbpface,
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
  Differentiated version of [`getFaceVariables`](@ref).
"""
function getFaceVariables_diff(params::ParamType,
                          capture::SBPParabolicSC{Tsol, Tres},
                          diffusion::AbstractDiffusion,
                          sbp::AbstractOperator,
                          sbpface::AbstractFace, iface_red::Interface,
                          wL::AbstractMatrix, wR::AbstractMatrix,
                          qL::AbstractMatrix, qR::AbstractMatrix,
                          gradwL::Abstract3DArray, gradwR::Abstract3DArray,
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
    getIRA0inv(params, qL_i, wL_dot_i)
    getIRA0inv(params, qR_i, wR_dot_i)
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
  applyDiffusionTensor_diff(diffusion, wL, iface_red.elementL, sbpface,
                            iface_red.faceL, t1, t1_dot, t2L_dot)

  # do volume operations for wR (derivative of wgradR)
  calcDx(sbp, dxidxR, jacR, Dx)
  applyOperatorJac(Dx, sbpface, iface_red.faceR, wR_dot, t1_dot)
  for d=1:dim
    Dx_d = ro_sview(Dx, :, :, d)
    t1_d = sview(t1, :, :, d)
    smallmatmatT!(wR, Dx_d, t1_d)
  end
  applyDiffusionTensor_diff(diffusion, wR, iface_red.elementR, sbpface,
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



function applyDgkTranspose_diff(capture::SBPParabolicSC{Tsol, Tres}, sbp,
                       sbpface, iface::Interface,
                       diffusion::AbstractDiffusion,
                       t2L::AbstractMatrix, t2L_dotL::Abstract4DArray,
                       t2L_dotR::Abstract4DArray,
                       t2R::AbstractMatrix, t2R_dotL::Abstract4DArray,
                       t2R_dotR::Abstract4DArray,
                       wL::AbstractMatrix, wR::AbstractMatrix,
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
  applyDiffusionTensor_diff(diffusion, wL, iface.elementL, temp2L,
                            t4L_dotL, t4L_dotR, t5L_dotL, t5L_dotR)
  # reverse L and R here so the d/dq term will be added to t5R_dotR
  applyDiffusionTensor_diff(diffusion, wR, iface.elementR, temp2R,
                            t4R_dotR, t4R_dotL, t5R_dotR, t5R_dotL)

  calcDxTransposed(sbp, dxidxL, jacL, Dx)
  applyOperatorJac(Dx, t5L_dotL, resL_dotL, false, op)
  applyOperatorJac(Dx, t5L_dotR, resL_dotR, false, op)

  calcDxTransposed(sbp, dxidxR, jacR, Dx)
  applyOperatorJac(Dx, t5R_dotL, resR_dotL, false, op)
  applyOperatorJac(Dx, t5R_dotR, resR_dotR, false, op)

  return nothing
end


function applyPenalty_diff(penalty::BR2Penalty{Tsol, Tres}, sbp, sbpface,
                      diffusion::AbstractDiffusion, iface::Interface,
                      delta_w::AbstractMatrix{Tsol},
                      delta_w_dotL::Abstract4DArray, delta_w_dotR::Abstract4DArray,
                      theta::AbstractMatrix{Tres},
                      theta_dotL::Abstract4DArray, theta_dotR::Abstract4DArray,
                      wL::AbstractMatrix, wR::AbstractMatrix,
                      nrm_face::AbstractMatrix,
                      alphas::AbstractVector,
                      jacL::AbstractVector, jacR::AbstractVector,
                      res1L_dotL::Abstract4DArray, res1L_dotR::Abstract4DArray,
                      res1R_dotL::Abstract4DArray, res1R_dotR::Abstract4DArray,
                      res2L::AbstractMatrix, res2R::AbstractMatrix,
                      res2L_dotL::Abstract4DArray, res2L_dotR::Abstract4DArray,
                      res2R_dotL::Abstract4DArray, res2R_dotR::Abstract4DArray,
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
  applyDiffusionTensor_diff(diffusion, wL, iface.elementL, sbpface,
                            iface.faceL, qL, t2LL, t2LR, t3LL, t3LR)
  # reverse RR and RL here because Lambda is a function of qR, not qL
  applyDiffusionTensor_diff(diffusion, wR, iface.elementR, sbpface,
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
  fill!(res1L_dotL, 0); fill!(res1L_dotR, 0)  #TODO: tile this
  fill!(res1R_dotL, 0); fill!(res1R_dotR, 0)
  @simd for q=1:numNodesPerElement
    @simd for p=1:numNodesPerFace
      @simd for d=1:dim
        fac = sbpface.wface[p]*nrm_face[d, p]
        @simd for j=1:numDofPerNode
          @simd for i=1:numDofPerNode
            valL = fac*(t4LL[i, j, d, p, q] + t4RL[i, j, d, p, q])
            valR = fac*(t4LR[i, j, d, p, q] + t4RR[i, j, d, p, q])
            #println("t2L_d_dot = ", real(t4LL[i, j, d, p, q]))
            #println("t2R_d_dot = ", real(t4RL[i, j, d, p, q]))
            #println("val = ", real(valL))
            res1L_dotL[i, j, p, q] +=  valL
            res1L_dotR[i, j, p, q] +=  valR
            res1R_dotL[i, j, p, q] += -valL  # delta_w is reversed for elementR
            res1R_dotR[i, j, p, q] += -valR
          end
        end
      end
    end
  end

  # T2 and T3
  fill!(res2L_dotL, 0); fill!(res2L_dotR, 0)
  fill!(res2R_dotL, 0); fill!(res2R_dotR, 0)
  @simd for q=1:numNodesPerElement
    @simd for p=1:numNodesPerFace
      fac1 =  0.5*sbpface.wface[p]
      fac2 = -0.5*sbpface.wface[p]
      @simd for j=1:numDofPerNode
        @simd for i=1:numDofPerNode
          res1L_dotL[i, j, p, q] += fac1*theta_dotL[i, j, p, q]
          res1L_dotR[i, j, p, q] += fac1*theta_dotR[i, j, p, q]
          res1R_dotL[i, j, p, q] += fac1*theta_dotL[i, j, p, q]
          res1R_dotR[i, j, p, q] += fac1*theta_dotR[i, j, p, q]

          res2L_dotL[i, j, p, q] += fac2*delta_w_dotL[i, j, p, q]
          res2L_dotR[i, j, p, q] += fac2*delta_w_dotR[i, j, p, q]
          res2R_dotL[i, j, p, q] -= fac2*delta_w_dotL[i, j, p, q]
          res2R_dotR[i, j, p, q] -= fac2*delta_w_dotR[i, j, p, q]
        end
      end
    end
  end


  # also need res2 as output
  @simd for j=1:numNodesPerFace
    @simd for k=1:numDofPerNode
      val2 = -0.5*sbpface.wface[j]*delta_w[k, j]
      res2L[k, j] =  val2
      res2R[k, j] = -val2  # delta_w is reversed for elementR
    end
  end


  return false
end


