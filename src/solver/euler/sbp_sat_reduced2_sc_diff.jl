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


