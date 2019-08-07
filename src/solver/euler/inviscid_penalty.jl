"""
  Penalty for BR2
"""
struct InviscidPenalty{Tsol, Tres} <: AbstractDiffusionPenalty


  function InviscidPenalty{Tsol, Tres}(mesh, sbp, opts) where {Tsol, Tres}

    return new()
  end
end


function applyPenalty(penalty::InviscidPenalty{Tsol, Tres}, sbp, params::AbstractParamType,
                      sbpface,
                      diffusion::AbstractDiffusion, iface::Interface,
                      delta_w::AbstractMatrix{Tsol}, theta::AbstractMatrix{Tres},
                      qL_el::AbstractMatrix, qR_el::AbstractMatrix,
                      wL::AbstractMatrix, wR::AbstractMatrix,
                      coordsL::AbstractMatrix, coordsR::AbstractMatrix,
                      nrm_face::AbstractMatrix,
                      alphas::AbstractVector,
                      dxidxL::Abstract3DArray, dxidxR::Abstract3DArray,
                      jacL::AbstractVector, jacR::AbstractVector,
                      res1L::AbstractMatrix, res1R::AbstractMatrix,
                      res2L::AbstractMatrix, res2R::AbstractMatrix
                     ) where {Tsol, Tres}

  fill!(res1L, 0); fill!(res1R, 0)
  fill!(res2L, 0); fill!(res2R, 0)

  numDofPerNode, numNodesPerFace = size(delta_w)
  numNodesPerElement = length(jacL)
  dim = size(nrm_face, 1)

  for i=1:numNodesPerFace
    fac = calcLength(params, ro_sview(nrm_face, :, i))*sbpface.wface[i]
    for j=1:numDofPerNode
      res1L[j, i] =  fac*delta_w[j]
      res1R[j, i] = -fac*delta_w[j]
    end
  end

  return nothing
end


function applyPenalty_diff(penalty::InviscidPenalty{Tsol, Tres}, sbp,
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
                      alphas::AbstractVector,
                      dxidxL::Abstract3DArray, dxidxR::Abstract3DArray,
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


  for i=1:numNodesPerFace
    fac = calcLength(params, ro_sview(nrm_face, :, i))
    for j=1:numDofPerNode
      #res1L[j, i] =  fac*delta_w[j]
      #res1R[j, i] = -fac*delta_w[j]
    end
  end

  for q=1:numNodesPerElement
    for p=1:numNodesPerFace
      fac = calcLength(params, ro_sview(nrm_face, :, p))*sbpface.wface[p]
      for j=1:numDofPerNode
        for i=1:numDofPerNode
          res1L_dotL[i, j, p, q] =  fac*delta_w_dotL[i, j, p, q]
          res1L_dotR[i, j, p, q] =  fac*delta_w_dotR[i, j, p, q]
          res1R_dotL[i, j, p, q] = -fac*delta_w_dotL[i, j, p, q]
          res1R_dotL[i, j, p, q] = -fac*delta_w_dotR[i, j, p, q]
        end
      end
    end
  end

  fill!(res2L_dotL, 0); fill!(res2L_dotR, 0)
  fill!(res2R_dotL, 0); fill!(res2R_dotR, 0)

  return false
end



