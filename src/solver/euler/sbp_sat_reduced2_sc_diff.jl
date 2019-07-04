# differentiated version of sbp_sat_shock_capturing.jl

function calcShockCapturing_diff(mesh::AbstractMesh, sbp::AbstractOperator,
                             eqn::EulerData, opts,
                             capture::SBPParabolicReduced2SC{Tsol, Tres},
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

  computeDirichletBoundaryTerm_diff(mesh, sbp, eqn, opts, capture, shockmesh,
                                    assem)
 
#  if shockmesh.isNeumann
#    println("computing Neumann boundary condition")
#    computeNeumannBoundaryTerm_diff(mesh, sbp, eqn, opts, capture, shockmesh,
#                              capture.diffusion, capture.entropy_vars, assem)
#  else
#    println("computing Dirichlet boundary condition")
#    computeDirichletBoundaryTerm_diff(mesh, sbp, eqn, opts, capture, shockmesh,
#                                      assem)
#  end

  #@time computeSharedFaceTerm_diff(mesh, sbp, eqn, opts, capture, shockmesh,
  #                                 capture.diffusion, capture.penalty)



  return nothing
end


function computeVolumeTerm_diff(mesh, sbp, eqn, opts,
                           capture::SBPParabolicReduced2SC{Tsol, Tres},
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
#=
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
=#
    calcQxTransposed(sbp, dxidx_i, Dx)
    applyOperatorJac(Dx, t2_dot, res_jac, true, op)

    assembleElement(assembler, mesh, i_full, res_jac)
  end

  return nothing
end


function computeFaceTerm_diff(mesh, sbp, eqn, opts,
                      capture::SBPParabolicReduced2SC{Tsol, Tres},
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
#=
  res2L_dotL = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode,
                           mesh.numNodesPerFace, mesh.numNodesPerElement)
  res2L_dotR = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode,
                           mesh.numNodesPerFace, mesh.numNodesPerElement)

  res2R_dotL = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode,
                           mesh.numNodesPerFace, mesh.numNodesPerElement)
  res2R_dotR = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode,
                           mesh.numNodesPerFace, mesh.numNodesPerElement)
=#

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
    resL = sview(eqn.res, :, :, elnumL)
    resR = sview(eqn.res, :, :, elnumR)

    # compute delta w tilde and theta_bar = Dgk w_k + Dgn w_n
    getFaceVariables_diff(eqn.params, capture, diffusion, entropy_vars,
                          sbp, mesh.sbpface,
                          iface_red, wL, wR, qL, qR, gradwL, gradwR,
                          coordsL, coordsR,
                          nrm_face, dxidxL, dxidxR, jacL, jacR,
                          delta_w, delta_w_dotL, delta_w_dotR)

    # apply the penalty coefficient matrix
    has_T4 = applyReduced2Penalty_diff(penalty, sbp, eqn.params,  mesh.sbpface,
                      diffusion,
                      iface_red, delta_w, delta_w_dotL, delta_w_dotR,
                      theta, theta_dotL, theta_dotR,
                      qL, qR, wL, wR, coordsL, coordsR,  nrm_face,
                      dxidxL, dxidxR, jacL, jacR,
                      res1L_dotL, res1L_dotR, res1R_dotL, res1R_dotR)


    # apply Rgk^T, Rgn^T, Dgk^T, Dgn^T
    fill!(res_jacLL, 0); fill!(res_jacRL, 0)
    fill!(res_jacLR, 0); fill!(res_jacRR, 0)
    interiorFaceIntegrate_jac!(mesh.sbpface, iface_red, res1L_dotL, res1R_dotL,
                               res_jacLL, res_jacRL, subtract, subtract, false)
    interiorFaceIntegrate_jac!(mesh.sbpface, iface_red, res1L_dotR, res1R_dotR,
                               res_jacLR, res_jacRR, subtract, subtract, false)

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


"""
  Differentiated version of [`getFaceVariables`](@ref).
"""
function getFaceVariables_diff(params::ParamType,
                          capture::SBPParabolicReduced2SC{Tsol, Tres},
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
                         ) where {Tsol, Tres}

  numDofPerNode, numNodesPerElement = size(wL)
  numNodesPerFace = size(delta_w, 2)
  dim = size(gradwL, 3)

  getFaceVariables(capture, sbpface, iface_red, wL, wR, gradwL, gradwR,
                   nrm_face, delta_w)

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

  return nothing
end


#------------------------------------------------------------------------------
# Dirichlet boundary condition


"""
  Differentiated Dirichlet boundary term
"""
function computeDirichletBoundaryTerm_diff(mesh, sbp, eqn, opts,
                      capture::ParabolicReduceds{Tsol, Tres},
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
                      capture::ParabolicReduceds{Tsol, Tres},
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

  #op = SummationByParts.Subtract()
  subtract = SummationByParts.Subtract()
  op = SummationByParts.Add()

  for i in bc_range
    bndry_i = shockmesh.bndryfaces[i].bndry
    idx_orig = shockmesh.bndryfaces[i].idx_orig
    elnum_orig = shockmesh.elnums_all[bndry_i.element]
    fill!(res_dot, 0)

    w_i = ro_sview(capture.w_el, :, :, bndry_i.element)
    q_i = ro_sview(eqn.q, :, :, elnum_orig)
    coords_i = ro_sview(mesh.coords_bndry, :, :, idx_orig)
    nrm_i = ro_sview(mesh.nrm_bndry, :, :, idx_orig)
    #alpha = capture.alpha_b[i]
    alpha = 1
    dxidxL = ro_sview(mesh.dxidx, :, :, :, elnum_orig)
    jacL = ro_sview(mesh.jac, :, elnum_orig)
    res_i = sview(eqn.res, :, :, elnum_orig)

    computeDirichletDelta_diff(capture, eqn.params, eqn.params_complex,
                          mesh.sbpface, bndry_i,
                          bc_func, entropy_vars, w_i, q_i, coords_i,
                          nrm_i, delta_w, delta_w_dot)

    applyDirichletPenalty_diff(penalty, sbp, eqn.params, mesh.sbpface, diffusion,
                          bndry_i, delta_w, delta_w_dot, q_i, w_i, coords_i,
                          nrm_i, alpha,
                          dxidxL, jacL, res1_dot)

    # apply R^T to T_D * delta_w
    #fill!(res_dot, 0)
    boundaryFaceIntegrate_jac!(mesh.sbpface, bndry_i.face, res1_dot, res_dot,
                               subtract, include_quadrature=false)
#    boundaryFaceIntegrate_jac!(mesh.sbpface, bndry_i.face, delta_w_dot, res_dot,
#                               op, include_quadrature=false)


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
                      nrm_i, dxidxL, jacL, res_dot, res_dotR, subtract)

    assembleElement(assem, mesh, elnum_orig, res_dot)
  end

  return nothing
end


function computeDirichletDelta_diff(capture::ParabolicReduceds{Tsol, Tres},
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




"""
  Jacobian of Dgk * t2L only
"""
function applyDgkTranspose_diff(capture::ParabolicReduceds{Tsol, Tres}, sbp,
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




#------------------------------------------------------------------------------
# Penalty

"""
  Differentiated version of BR2 penalty
"""
function applyReduced2Penalty_diff(penalty::BR2Penalty{Tsol, Tres}, sbp,
                      params::AbstractParamType, sbpface,
                      diffusion::AbstractDiffusion, iface::Interface,
                      delta_w::AbstractMatrix{Tsol},
                      delta_w_dotL::Abstract4DArray, delta_w_dotR::Abstract4DArray,
                      theta::AbstractMatrix{Tres},
                      theta_dotL::Abstract4DArray, theta_dotR::Abstract4DArray,
                      qL_el::AbstractMatrix, qR_el::AbstractMatrix,
                      wL::AbstractMatrix, wR::AbstractMatrix,
                      coordsL::AbstractMatrix, coordsR::AbstractMatrix,
                      nrm_face::AbstractMatrix,
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
      facL = jacL[p]/sbp.w[p]
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
      facR = jacR[p]/sbp.w[p]
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
#  fill!(res1L_dotL, 0); fill!(res1L_dotR, 0)  #TODO: tile this
#  fill!(res1R_dotL, 0); fill!(res1R_dotR, 0)
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

  return false
end


