# Modified scheme of Bassi and Rebay (BR2)


"""
  Penalty for BR2
"""
struct BR2Penalty{Tsol, Tres} <: AbstractDiffusionPenalty

  # primal method
  delta_w_n::Matrix{Tsol}
  qL::Array{Tres, 3}
  qR::Array{Tres, 3}
  t1L::Array{Tres, 3}
  t1R::Array{Tres, 3}
  t2L::Array{Tres, 3}
  t2R::Array{Tres, 3}

  # diff method
  t1_dotL::Array{Tres, 5}
  t1_dotR::Array{Tres, 5}

  # the first L/R indicates if this is element kappa or nu
  # the second L/R indices if the derivative is wrt qL or qR
  t2LL::Array{Tres, 5}
  t2LR::Array{Tres, 5}
  t2RL::Array{Tres, 5}
  t2RR::Array{Tres, 5}

  t3LL::Array{Tres, 5}
  t3LR::Array{Tres, 5}
  t3RL::Array{Tres, 5}
  t3RR::Array{Tres, 5}

  # re-use t1_dotL,R
  t4RL::Array{Tres, 5}
  t4RR::Array{Tres, 5}

  function BR2Penalty{Tsol, Tres}(mesh, sbp, opts) where {Tsol, Tres}

    @unpack mesh numDofPerNode numNodesPerFace numNodesPerElement dim

    delta_w_n = zeros(Tsol, numDofPerNode, numNodesPerFace)
    qL = zeros(Tres, numDofPerNode, numNodesPerElement, dim)
    qR = zeros(Tres, numDofPerNode, numNodesPerElement, dim)
    t1L = zeros(Tres, numDofPerNode, numNodesPerElement, dim)
    t1R = zeros(Tres, numDofPerNode, numNodesPerElement, dim)
    t2L = zeros(Tres, numDofPerNode, numNodesPerFace, dim)
    t2R = zeros(Tres, numDofPerNode, numNodesPerFace, dim)

    t1_dotL = zeros(Tres, numDofPerNode, numDofPerNode, dim, numNodesPerFace,
                   numNodesPerElement)
    t1_dotR = zeros(Tres, numDofPerNode, numDofPerNode, dim, numNodesPerFace,
                     numNodesPerElement)

    t2LL = zeros(Tres, numDofPerNode, numDofPerNode, dim, numNodesPerElement,
                      numNodesPerElement)
    t2LR = zeros(Tres, numDofPerNode, numDofPerNode, dim, numNodesPerElement,
                      numNodesPerElement)
   
    t2RL = zeros(Tres, numDofPerNode, numDofPerNode, dim, numNodesPerElement,
                      numNodesPerElement)
    t2RR = zeros(Tres, numDofPerNode, numDofPerNode, dim, numNodesPerElement,
                      numNodesPerElement)
   
    t3LL = zeros(Tres, numDofPerNode, numDofPerNode, dim, numNodesPerElement,
                       numNodesPerElement)
    t3LR = zeros(Tres, numDofPerNode, numDofPerNode, dim, numNodesPerElement,
                       numNodesPerElement)

    t3RL = zeros(Tres, numDofPerNode, numDofPerNode, dim, numNodesPerElement,
                       numNodesPerElement)
    t3RR = zeros(Tres, numDofPerNode, numDofPerNode, dim, numNodesPerElement,
                       numNodesPerElement)


    t4RL = zeros(Tres, numDofPerNode, numDofPerNode, dim, numNodesPerFace,
                     numNodesPerElement)
    t4RR = zeros(Tres, numDofPerNode, numDofPerNode, dim, numNodesPerFace,
                       numNodesPerElement)


    return new(delta_w_n, qL, qR, t1L, t1R, t2L, t2R,
               t1_dotL, t1_dotR, t2LL, t2LR, t2RL, t2RR, t3LL, t3LR, t3RL, t3RR,
               t4RL, t4RR)
  end

end

"""
  Applies the penalty for the SBP-SAT generalization of the modified scheme
  of Bassi and Rebay (BR2).
"""
function applyPenalty(penalty::BR2Penalty{Tsol, Tres}, sbp, sbpface,
                      diffusion::AbstractDiffusion, iface::Interface,
                      delta_w::AbstractMatrix{Tsol}, theta::AbstractMatrix{Tres},
                      wL::AbstractMatrix, wR::AbstractMatrix,
                      nrm_face::AbstractMatrix,
                      alphas::AbstractVector,
                      jacL::AbstractVector, jacR::AbstractVector,
                      res1L::AbstractMatrix, res1R::AbstractMatrix,
                      res2L::AbstractMatrix, res2R::AbstractMatrix
                     ) where {Tsol, Tres}

  fill!(res1L, 0); fill!(res1R, 0)
  fill!(res2L, 0); fill!(res2R, 0)
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
  applyDiffusionTensor(diffusion, wL, iface.elementL, sbpface, iface.faceL,
                       qL, t1L)
  applyDiffusionTensor(diffusion, wR, iface.elementR, sbpface, iface.faceR,
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

  #----------------------------
  # apply T2 and T3
  @simd for j=1:numNodesPerFace
    @simd for k=1:numDofPerNode
      val1 =  0.5*sbpface.wface[j]*theta[k, j]
      val2 = -0.5*sbpface.wface[j]*delta_w[k, j]
      res1L[k, j] += val1
      res1R[k, j] += val1
      res2L[k, j] += val2
      res2R[k, j] -= val2  # delta_w is reversed for elementR
    end
  end

  return nothing
end


"""
  Differentiated version of BR2 penalty
"""
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

  # T2 and T3
  @simd for q=1:numNodesPerElement
    @simd for p=1:numNodesPerFace
      @simd for j=1:numDofPerNOde
        @simd for i=1:numDofPerNode
          res2L_dotL[i, j, p, q] = 0
          res2L_dotR[i, j, p, q] = 0
          res2R_dotL[i, j, p, q] = 0
          res2R_dotR[i, j, p, q] = 0
        end
      end
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
  fill!(res2L, 0); fill!(res2R, 0)
  @simd for j=1:numNodesPerFace
    @simd for k=1:numDofPerNode
      val2 = -0.5*sbpface.wface[j]*delta_w[k, j]
      res2L[k, j] =  val2
      res2R[k, j] = -val2  # delta_w is reversed for elementR
    end
  end


  return false
end


