@doc """
  This function calculates the potential flux psi from entropy stability 
  theory.  Methods are available for both 2 and 3D
"""

function getPsi(params::ParamType{2, :conservative}, q_vals::AbstractVector, nrm::AbstractVector)
  return nrm[1]*q_vals[2] + nrm[2]*q_vals[3]
end

function getPsi(params::ParamType{3, :conservative}, q_vals::AbstractVector, nrm::AbstractVector)
  return nrm[1]*q_vals[2] + nrm[2]*q_vals[3] +nrm[3]*q_vals[4]
end

function getPsi{Tsol}(params, qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, nrm::AbstractVector)

  numDofPerNode, numNodesPerElement = size(qL)
  psiL = zeros(Tsol, numNodesPerElement)
  psiR = zeros(psiL)
  dim = length(nrm)
  # calculate psi at volume nodes
  for i=1:numNodesPerElement
    qL_i = qL[:, i]
    qR_i = qR[:, i]

    psi_vecL = psi_vec(params, qL_i)
    psi_vecR = psi_vec(params, qR_i)
    for d=1:dim
      psiL[i] += nrm[d]*psi_vecL[d]
      psiR[i] += nrm[d]*psi_vecR[d]
    end
  end

  return psiL, psiR
end

@doc """
  This function computes the integral of the potential flux over an interface
"""
function computeInterfacePotentialFlux{Tdim, Tsol, Tres}(
                params::ParamType{Tdim, :conservative, Tsol, Tres}, 
                iface::Interface, sbpface, dxidx_face, 
                qL::AbstractMatrix, qR::AbstractMatrix)
# compute the potential flux then compute the reduction with Eface

  rhs = zero(Tres)
  for dim=1:Tdim
    rhs += reduceEface(params, iface, sbpface, dxidx_face, dim, qL, qR)
  end

  return rhs
end

@doc """
  Computes 1^T(E_faceL * psiL + E_faceR * psiR), where E_face performs the 
  face integral over the specified face and psi is the potential flux 
  (from entropy stability theory).  The dir argument indicates whether 
  the normal vector used for the integral is x, y, or z, ie. in 

  integral psi dot n dGamma

  it determines the normal vector n
"""
function reduceEface{Tdim, Tsol, Tres, Tmsh}(params::ParamType{Tdim, :conservative, Tsol, Tres, Tmsh}, iface, sbpface, dxidx_face::Abstract3DArray{Tmsh}, dir::Integer, qL::AbstractMatrix, qR::AbstractMatrix)
  # compute Ex_gamma kappa * psiL + Ex_gamma_nu * psiR, where x is one 
  # of either x or y, as specified by dir

  RHS1 = zero(Tres)
  RHS2 = zero(Tres)

  flux_nrm = params.nrm
  fill!(flux_nrm, 0.0)
  flux_nrm[dir] = 1

  for i=1:sbpface.stencilsize
    for j=1:sbpface.stencilsize
      p_jL = sbpface.perm[j, iface.faceL]
      p_jR = sbpface.perm[j, iface.faceR]
      psiL = getPsi(params, sview(qL, :, p_jL), flux_nrm)
      psiR = getPsi(params, sview(qR, :, p_jR), flux_nrm)
      for k=1:sbpface.numnodes
        nrm_k = zero(Tmsh)
        for d=1:Tdim
          nrm_k += sbpface.normal[d, iface.faceL]*dxidx_face[d, dir, k]
        end
        val = sbpface.interp[i,k]*sbpface.interp[j,k]*sbpface.wface[k]*nrm_k
        RHS1 += val*psiL

        kR = sbpface.nbrperm[k, iface.orient]
        val = sbpface.interp[i, kR]*sbpface.interp[j, kR]*sbpface.wface[k]*nrm_k
        RHS1 -= val*psiR
      end
    end
  end

  return RHS1 + RHS2
end

#TODO: this can be made more efficient once SBP stores E
function computeVolumePotentialFlux{Tdim, Tsol, Tres}(params::ParamType{Tdim, :conservative, Tsol, Tres}, sbp, q_i::AbstractMatrix, dxidx)
  
  numDofPerNode, numNodesPerElement = size(q_i)

  nrm = params.nrm
  # calculate psi vector
  psi = zeros(numNodesPerElement, 2)
  for j=1:numNodesPerElement
    q_j = sview(q_i, :, j)
    for d=1:Tdim
      fill!(nrm, 0.0)
      nrm[d] = 1
      psi[j, d] = getPsi(params, q_j, nrm)
    end
  end

  rhs_reduced = zero(Tres)
  for d=1:Tdim
#    println("dimension = ", d)
    E_d = (sbp.Q[:, :, d] + sbp.Q[:, :, d].')
    psi_nrm = dxidx[d, 1, 1]*psi[:, 1] + dxidx[d, 2, 1]*psi[:, 2]
#    println("psi_nrm = ", psi_nrm)
    val = sum(E_d*psi_nrm)
#    println("rhs_reduced = ", val)

    rhs_reduced += val
  end

  return -rhs_reduced
end

