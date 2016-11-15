# functions that do face integral-like operations, but operate on data from
# the entire element
include("IR_stab.jl")  # stabilization for the IR flux

"""
  Calculate the face integrals in an entropy stable manner for a given
  interface.  Unlike standard face integrals, this requires data from
  the entirety of both elements, not just data interpolated to the face

  resL and resR are updated with the results of the computation for the 
  left and right elements, respectively.

  Note that dxidx_face must contains the scaled metric terms interpolated 
  to the face nodes.

  Aliasing restrictions: none, although its unclear what the meaning of this
                         function would be if resL and resR alias

  Performance note: the version in the tests is the same speed as this one
                    for p=1 Omega elements and about 10% faster for 
                    p=4 elements, but would not be able to take advantage of t
                    he sparsity of R for SBP Gamma elements
"""
                  
function calcESFaceIntegral{Tdim, Tsol, Tres, Tmsh}(
                             params::AbstractParamType{Tdim}, 
                             sbpface::AbstractFace, 
                             iface::Interface,
                             qL::AbstractMatrix{Tsol}, 
                             qR::AbstractMatrix{Tsol}, 
                             aux_vars::AbstractMatrix{Tres}, 
                             dxidx_face::Abstract3DArray{Tmsh},
                             functor::FluxType, 
                             resL::AbstractMatrix{Tres}, 
                             resR::AbstractMatrix{Tres})


  Flux_tmp = params.flux_vals1
  numDofPerNode = length(Flux_tmp)
  nrm = params.nrm
  for dim = 1:Tdim
    fill!(nrm, 0.0)
    nrm[dim] = 1

    # loop over the nodes of "left" element that are in the stencil of interp
    for i = 1:sbpface.stencilsize
      p_i = sbpface.perm[i, iface.faceL]
      qi = sview(qL, :, p_i)
      aux_vars_i = sview(aux_vars, :, p_i)  # !!!! why no aux_vars_j???

      # loop over the nodes of "right" element that are in the stencil of interp
      for j = 1:sbpface.stencilsize
        p_j = sbpface.perm[j, iface.faceR]
        qj = sview(qR, :, p_j)

        # accumulate entry p_i, p_j of E
        Eij = zero(Tres)  # should be Tres
        for k = 1:sbpface.numnodes
          # the computation of nrm_k could be moved outside i,j loops and saved
          # in an array of size [3, sbp.numnodes]
          nrm_k = zero(Tmsh)
          for d = 1:Tdim
            nrm_k += sbpface.normal[d, iface.faceL]*dxidx_face[d, dim, k]
          end
          kR = sbpface.nbrperm[k, iface.orient]
          Eij += sbpface.interp[i,k]*sbpface.interp[j,kR]*sbpface.wface[k]*nrm_k
        end  # end loop k
        
        # compute flux and add contribution to left and right elements
        functor(params, qi, qj, aux_vars_i, nrm, Flux_tmp)
        for p=1:numDofPerNode
          resL[p, p_i] -= Eij*Flux_tmp[p]
          resR[p, p_j] += Eij*Flux_tmp[p]
        end

      end
    end
  end  # end loop Tdim


  return nothing
end
function calcEDissipativeFaceIntegral{Tdim, Tsol, Tres, Tmsh}(
                             params::AbstractParamType{Tdim}, 
                             sbpface::AbstractFace, 
                             iface::Interface,
                             qL::AbstractMatrix{Tsol}, 
                             qR::AbstractMatrix{Tsol}, 
                             aux_vars::AbstractMatrix{Tres}, 
                             dxidx_face::Abstract3DArray{Tmsh},
                             functor::FluxType, 
                             resL::AbstractMatrix{Tres}, 
                             resR::AbstractMatrix{Tres})

  calcESFaceIntegral(params, sbpface, iface, qL, qR, aux_vars, dxidx_face, 
                     functorm, resL, resR)
  calcEntropyPenaltyIntegral(params, sbpface, iface, qL, qR, aux_vars, 
                             dxidx_face, resL, resR)

  return nothing
end

"""

  Aliasing restrictions: params.nrm2, params.A0
"""

function calcEntropyPenaltyIntegral{Tdim, Tsol, Tres, Tmsh}(
             params::ParamType{Tdim, :conservative, Tsol, Tres, Tmsh},
             sbpface::AbstractFace, iface::Interface, 
             qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, 
             aux_vars::AbstractMatrix{Tres}, dxidx_face::Abstract3DArray{Tmsh},
             resL::AbstractMatrix{Tres}, resR::AbstractMatrix{Tres})

#  println("----- entered calcEntropyDissipativeIntegral -----")

  numDofPerNode = size(qL, 1)

  #TODO: add these arrays to params
  # convert qL and qR to entropy variables (only the nodes that will be used)
  wL = params.w_vals_stencil
  wR = params.w_vals2_stencil
#  wL = zeros(Tsol, numDofPerNode, sbpface.stencilsize)
#  wR = zeros(wL)

  for i=1:sbpface.stencilsize
    # apply sbpface.perm here
    p_iL = sbpface.perm[i, iface.faceL]
    p_iR = sbpface.perm[i, iface.faceR]
    qL_i = sview(qL, :, p_iL)
    qR_i = sview(qR, :, p_iR)
    wL_i = sview(wL, :, i)
    wR_i = sview(wR, :, i)
    convertToIR(params, qL_i, wL_i)
    convertToIR(params, qR_i, wR_i)
  end

  # convert to IR entropy variables

  # accumulate wL at the node
  wL_i = zeros(Tsol, numDofPerNode)
  wR_i = zeros(Tsol, numDofPerNode)
  # convert wL at the node back to qL
  qL_i = zeros(Tsol, numDofPerNode)
  qR_i = zeros(Tsol, numDofPerNode)
  dir = params.nrm2
  A0 = params.A0

  for i=1:sbpface.numnodes  # loop over face nodes
    ni = sbpface.nbrperm[i]
    fill!(wL_i, 0.0)
    fill!(wR_i, 0.0)

    # interpolate wL and wR to this node
    for j=1:sbpface.stencilsize
      interpL = sbpface.interp[j, i]
      interpR = sbpface.interp[j, ni]

      for k=1:numDofPerNode
        wL_i[k] += interpL*wL[k, j]
        wR_i[k] += interpR*wR[k, j]
      end
    end

    # get the normal vector (scaled)
    for dim =1:Tdim
      nrm_dim = zero(Tmsh)
      for d = 1:Tdim
        nrm_dim += sbpface.normal[d, iface.faceL]*dxidx_face[d, dim, i]
      end
      dir[dim] = nrm_dim
    end

    convertToConservativeFromIR_(params, wL_i, qL_i)
    convertToConservativeFromIR_(params, wR_i, qR_i)
    # get lambda * IRA0
    lambda_max = getLambdaMax(params, qL_i, qR_i, dir)
    
    # compute average qL
    # also delta w (used later)
    for j=1:numDofPerNode
      qL_i[j] =0.5*( qL_i[j] + qR_i[j])
      wL_i[j] -= wR_i[j]
    end

    getIRA0(params, qL_i, A0)

    # wface[i] * lambda_max * A0 * delta w
    smallmatvec!(A0, wL_i, wR_i)
    scale!(wR_i, sbpface.wface[i]*lambda_max)

    # interpolate back to volume nodes
    for j=1:sbpface.stencilsize
      j_pL = sbpface.perm[j, iface.faceL]
      j_pR = sbpface.perm[j, iface.faceR]

      for p=1:numDofPerNode
        resL[p, j_pL] -= sbpface.interp[j, i]*wR_i[p]
        resR[p, j_pR] += sbpface.interp[j, ni]*wR_i[p]
      end
    end

  end  # end loop i

  return nothing
end


