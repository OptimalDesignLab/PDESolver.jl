include("viscous_penalty.jl")

@doc """

Calculate fluxes at edge cubature points using face-based form.
eqn.flux_face, eqn.xflux, eqn.yflux will be updated.

Input:
  mesh :
  sbp  :
  eqn  :
  opts :
Output:

"""->
function calcViscousFlux_interior(mesh::AbstractDGMesh{Tmsh},
                                  sbp::AbstractOperator,
                                  eqn::NSData{Tsol, Tres, Tdim},
                                  opts) where {Tmsh, Tsol, Tres, Tdim}

  Ma      = eqn.params.euler_params.Ma
  Re      = eqn.params.Re
  gamma_1 = eqn.params.euler_params.gamma_1
  Pr      = 0.72
  coef_nondim = Ma/Re
  interfaces  = sview(mesh.interfaces, :)
  nfaces      = length(mesh.interfaces)
  p    = opts["order"]
  dq   = Array{Tsol}(mesh.numDofPerNode, mesh.numNodesPerFace)
  dqn  = Array{Tsol}(Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)
  GtL  = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, Tdim, Tdim, mesh.numNodesPerFace)
  GtR  = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, Tdim, Tdim, mesh.numNodesPerFace)
  pMat  = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerFace)
  Fv_faceL = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)
  Fv_faceR = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)
  Fv_avg   = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)
  vecfluxL = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)
  vecfluxR = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)

  sbpface = mesh.sbpface
  sat_type = opts["SAT_type"]
  # this one is Harmann's definition
  const_tii = (p + 1.0)*(p + Tdim)/(2.0*Tdim)
  # const_tii = calcTraceInverseInequalityConst(sbp, sbpface)
  area_sum = sview(eqn.area_sum, :)

  params = eqn.params
  # sigma = calcTraceInverseInequalityConst(sbp, sbpface)
  # println("rho_max = ", sigma)

  for f = 1:nfaces    # loop over faces

    flux  = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
    face = interfaces[f]
    elemL = face.elementL
    elemR = face.elementR
    faceL = face.faceL
    faceR = face.faceR
    permL = sview(sbpface.perm, :, faceL)
    permR = sview(sbpface.perm, :, faceR)

    nrm_xy = ro_sview(mesh.nrm_face, :, :, f)

    # We need viscous flux and diffusion tensor on interfaces, and there
    # are different ways to compute them. For viscous flux:
    # 1) Fv = Fv(q, ∇q). Since we already have Q on face nodes, if 𝛻Q is also available 
    # on interface, then we are done. This way is probably more consistent with computation
    # of other terms like Fv(q_b, ∇ϕ)
    # 2) we can compute viscous flux on volume nodes and then interpolate to interface nodes.
    # It's logically simple but computationally expensive.

    # compute viscous flux and diffusion tensor
    # Arrayviews returns a n x 1 x m array (ie. 3D, not 2D), so don't use it
    q_faceL = sview(eqn.q_face, :, 1, :, f)
    q_faceR = sview(eqn.q_face, :, 2, :, f)       # TODO: use sview instead of slice
    q_elemL = sview(eqn.q, :, :, elemL)
    q_elemR = sview(eqn.q, :, :, elemR)
    calcDiffusionTensor(eqn.params, q_faceL, GtL)
    calcDiffusionTensor(eqn.params, q_faceR, GtR)

    # one way to compute Fv_face 
    # calcFvis_interiorFace(mesh, sbp, f, q_elemL, q_elemR, Fv_face)    

    # compute the face derivatives first, i.e., we first compute 
    # the derivatives at element nodes, and then do the interpolation.
    dqdx_face  = Array{Tsol}(Tdim, mesh.numDofPerNode, 2, mesh.numNodesPerFace)
    dqdx_elemL = Array{Tsol}(Tdim, mesh.numDofPerNode, mesh.numNodesPerElement)
    dqdx_elemR = Array{Tsol}(Tdim, mesh.numDofPerNode, mesh.numNodesPerElement)
    calcGradient(mesh, sbp, elemL, q_elemL, dqdx_elemL)
    calcGradient(mesh, sbp, elemR, q_elemR, dqdx_elemR)

    #
    # TODO: we can consider the first 2 dimension as a single dimension,
    # then we will not need slice here any more.
    #
    # for d = 1 : Tdim
      # dqdxL = sview(dqdx_elemL, d, :, :)
      # dqdxR = sview(dqdx_elemR, d, :, :)
      # dqdx_f = sview(dqdx_face, d, :, :, :)
      # interiorfaceinterpolate(sbpface, face, dqdxL, dqdxR, dqdx_f)
    # end

    # # Now both G and dqdx are avaiable at face nodes  
    # dqdx_faceL = sview(dqdx_face, :, :, 1, :)
    # dqdx_faceR = sview(dqdx_face, :, :, 2, :)
    # calcFvis(params, GtL, dqdx_faceL, Fv_faceL)
    # calcFvis(params, GtR, dqdx_faceR, Fv_faceR)

    Fv_face = zeros(Tsol, Tdim, mesh.numDofPerNode, 2, mesh.numNodesPerFace)
    jacL = ro_sview(mesh.jac, :, elemL)
    jacR = ro_sview(mesh.jac, :, elemR)
    dxidxL = ro_sview(mesh.dxidx, :,:,:,elemL)
    dxidxR = ro_sview(mesh.dxidx, :,:,:,elemR)
    calcFaceFvis(params, sbp, sbpface, q_elemL, q_elemR, dxidxL, jacL, dxidxR, jacR, face, Fv_face)
    Fv_faceL = sview(Fv_face, :,:,1,:)
    Fv_faceR = sview(Fv_face, :,:,2,:)
    # diffL = maximum(abs.(real(view(Fv_face, :, :, 1, :) - Fv_faceL)))
    # diffR = maximum(abs.(real(view(Fv_face, :, :, 2, :) - Fv_faceR)))
    # if (diffL > 1.e-8)
        # println(diffL)
    # end
    
    cmptIPMat(mesh, sbp, eqn, opts, f, GtL, GtR, pMat)

    # Start to compute fluxes. We have 3 terms on interfaces:
    # 1) {Fv}⋅[ϕ]
    # 2) {G^T ∇ϕ}⋅[q] 
    # 3) δ{G}[q]:[ϕ]

    # q jump
    for n = 1 : mesh.numNodesPerFace
      for iDof = 1 : mesh.numDofPerNode
        dq[iDof, n] = q_faceL[iDof, n] - q_faceR[iDof, n]
      end
    end

    # average viscous flux on face
    for n = 1 : mesh.numNodesPerFace
      for iDof = 2 : mesh.numDofPerNode
        for d = 1 : Tdim
          Fv_avg[d, iDof, n] = 0.5 * (Fv_faceL[d, iDof, n] + Fv_faceR[d, iDof, n] )
        end
      end
    end

    # Finally, everything is ready, let's compute fluxes, or penalties

    # This part computes the contribution of
    #     ∫ {G^T∇ϕ}:[q] dΓ = ∫ ∇ϕ⋅F dΓ , 
    # where 
    #     [q] = (q+ - q-) ⊗ n = Δq⊗n , 
    # Then we can consider Δq⊗n as ∇q and F as viscous flux.
    fill!(vecfluxL, 0.0)
    fill!(vecfluxR, 0.0)
    for n = 1 : mesh.numNodesPerFace
      for iDof = 1 : mesh.numDofPerNode
        #
        # sum up columns of each row
        #
        for iDim = 1 : Tdim
          # vecfluxL[iDim, iDof, n] = 0.0
          # vecfluxR[iDim, iDof, n] = 0.0
          for jDim = 1 : Tdim
            tmpL = 0.0
            tmpR = 0.0
            for jDof = 1 : mesh.numDofPerNode
              tmpL += GtL[iDof, jDof, iDim, jDim, n]
              tmpR += GtR[iDof, jDof, iDim, jDim, n]
            end
            vecfluxL[iDim, iDof, n] += tmpL * nrm_xy[jDim, n]
            vecfluxR[iDim, iDof, n] += tmpR * nrm_xy[jDim, n]
          end
          vecfluxL[iDim,iDof,n] *=  dq[iDof,n]
          vecfluxR[iDim,iDof,n] *=  dq[iDof,n]
        end
      end
    end

    # δ{G}[q]:n, contributing to  δ{G}[q]:[ϕ]
    for n = 1:mesh.numNodesPerFace
      for iDof = 2 : mesh.numDofPerNode
        for jDof = 1 : mesh.numDofPerNode
          flux[iDof, n] +=  pMat[iDof, jDof, n]*dq[jDof, n]
        end
      end
    end

    # {Fv}⋅n, contributing to {Fv}⋅[ϕ]
    for n = 1:mesh.numNodesPerFace
      for iDof = 2 : mesh.numDofPerNode
        for iDim = 1 : Tdim
          flux[iDof, n] -= Fv_avg[iDim, iDof, n]*nrm_xy[iDim,n] 
        end
      end
    end
    # accumulate fluxes
    for n = 1:mesh.numNodesPerFace
      for iDof = 2 : Tdim+2
        for iDim = 1 : Tdim
          eqn.vecflux_faceL[iDim, iDof, n, f] -=  vecfluxL[iDim, iDof, n]*coef_nondim
          eqn.vecflux_faceR[iDim, iDof, n, f] -=  vecfluxR[iDim, iDof, n]*coef_nondim
        end
        eqn.flux_face[iDof, n, f] += flux[iDof, n]*coef_nondim
      end
    end
  end # end of loop over all interfaces

  return nothing
end # end of function calcViscousFlux_interior



@doc """
Now actually we are integrating 
∫ G∇ϕ:[q] dΓ
where G is the diffusion tensor and q is the solution variable. It is not
immediately in the 2nd form. A alternative (or better) way to do this 
integral is as follows
∫ ∇ϕ⋅(G^t (qn))

Input:
Output:

"""->
function evalFaceIntegrals_vector(mesh::AbstractDGMesh{Tmsh},
                                  sbp::AbstractOperator,
                                  eqn::NSData{Tsol, Tres, Tdim},
                                  opts) where {Tmsh, Tsol, Tres, Tdim}
  # This part computes ∫ ∇ϕ⋅F  dΓ, 
  sbpface = mesh.sbpface
  DxL = Array{Tmsh}(mesh.numNodesPerElement, mesh.numNodesPerElement, Tdim)
  DxR = Array{Tmsh}(mesh.numNodesPerElement, mesh.numNodesPerElement, Tdim)

  GtL = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, Tdim, Tdim, mesh.numNodesPerFace)
  GtR = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, Tdim, Tdim, mesh.numNodesPerFace)

  R = sview(sbpface.interp[:,:])
  w = sview(sbpface.wface, :)
  res = sview(eqn.res, :,:,:)

  numNodes_elem = mesh.numNodesPerElement    # number of Nodes per elemet
  numNodes_face = mesh.numNodesPerFace       # number of nodes on interfaces
  stencilsize = sbpface.stencilsize        # size of stencil for interpolation

  RDxL = Array{Tmsh}(mesh.numNodesPerFace, mesh.numNodesPerElement, Tdim)
  RDxR = Array{Tmsh}(mesh.numNodesPerFace, mesh.numNodesPerElement, Tdim)
  FvL  = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerFace, mesh.numNodesPerElement)
  FvR  = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerFace, mesh.numNodesPerElement)
  GtRDxL = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerFace, mesh.numNodesPerElement)
  GtRDxR = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerFace, mesh.numNodesPerElement)
  nrm    = Array{Tmsh}(Tdim, mesh.numNodesPerFace)
  dq     =  Array{Tsol}(mesh.numDofPerNode, mesh.numNodesPerFace)

  nfaces = length(mesh.interfaces)

  for f = 1 : nfaces
    face = mesh.interfaces[f]

    elemL = face.elementL
    elemR = face.elementR
    faceL = face.faceL
    faceR = face.faceR
    pL = sview(sbpface.perm, :, faceL)
    pR = sview(sbpface.perm, :, faceR)

    # compute RDx
    calcDx(mesh, sbp, elemL, DxL)
    calcDx(mesh, sbp, elemR, DxR)

    for i = 1 : length(RDxL)
      RDxL[i] = 0.0
      RDxR[i] = 0.0
    end

    for d =  1 : Tdim    
      for row = 1 : numNodes_face
        rowR = sbpface.nbrperm[row, face.orient]
        for col = 1 : numNodes_elem
          for s = 1 : stencilsize
            RDxL[row, col, d] +=  R[s, row]  * DxL[pL[s], col, d]     
            RDxR[row, col, d] +=  R[s, rowR] * DxR[pR[s], col, d]     
          end
        end
      end
    end

    vecfluxL = sview(eqn.vecflux_faceL,:,:,:,f)
    vecfluxR = sview(eqn.vecflux_faceR,:,:,:,f)
    for i = 1 : numNodes_elem
      for j = 1 : numNodes_face
        for iDof = 2 : mesh.numDofPerNode
          tmpL = 0.0
          tmpR = 0.0
          for iDim = 1 : Tdim
            tmpL += RDxL[j, i, iDim] * vecfluxL[iDim, iDof, j]
            tmpR += RDxR[j, i, iDim] * vecfluxR[iDim, iDof, j]
          end
          res[iDof, i, elemL] += tmpL * w[j]
          res[iDof, i, elemR] += tmpR * w[j]
        end
      end
    end
  end

  return nothing
end





"""

Integrate ∫ ∇ϕ⋅F dΩ
TODO: consider combine it together with `weakdifferentiate`

Input:
  mesh   : 
  sbp    :
  eqn    :
  res    :
Output:

"""
function weakdifferentiate2!(mesh::AbstractMesh{Tmsh},
                             sbp::AbstractOperator{Tsbp},
                             eqn::NSData{Tsol, Tres, Tdim},
                             res::AbstractArray{Tres,3}) where {Tmsh, Tsbp, Tsol, Tres, Tdim}
  @assert (sbp.numnodes ==  size(res,2))

  dim             = Tdim
  numElems        = mesh.numEl
  numNodesPerElem = mesh.numNodesPerElement
  numDofsPerNode  = mesh.numDofPerNode

  gamma_1 = eqn.params.euler_params.gamma_1
  Pr = 0.72
  Ma = eqn.params.euler_params.Ma
  Re = eqn.params.Re
  coef_nondim = Ma/Re 

  Qx = Array{Tsbp}(numNodesPerElem, numNodesPerElem, dim)
  Fv = zeros(Tres, Tdim, numDofsPerNode, numNodesPerElem)
  w = sview(sbp.w, :)

  for elem = 1 : numElems
    # compute viscous flux
    q      = sview(eqn.q, :, :, elem)
    dxidx = sview(mesh.dxidx, :,:,:,elem)
    jac      = sview(mesh.jac, :, elem)

    calcFvis_elem(eqn.params, sbp, q, dxidx, jac, Fv)

    calcQx(mesh, sbp, elem, Qx)

    for d = 1 : dim
      for i = 1 : sbp.numnodes
        for j = 1 : sbp.numnodes
          for iDof = 2 : numDofsPerNode
            res[iDof, i, elem] -=  coef_nondim * Qx[j,i,d] * Fv[d, iDof, j] / jac[j]
          end
        end
      end
    end
  end
end


"""
  This flux function throws an error. Useful for defaults.
"""
mutable struct ErrorFlux <: FluxType
end

function (obj::ErrorFlux)(params::ParamType,
              uL::AbstractArray{Tsol,1},
              uR::AbstractArray{Tsol,1},
              aux_vars::AbstractVector{Tres},
              nrm::AbstractArray,
              F::AbstractArray{Tres}) where {Tsol, Tres}

  error("ErrorFlux called.")
  return nothing
end



"""
  Calls the [`SIPG`](@ref) (viscous) flux
"""
mutable struct SIPGViscousFlux <: FluxType
end

function (obj::SIPGViscousFlux)(params::ParamType,
              sbp::AbstractOperator,
              sbpface,    # TODO: type
              uL::AbstractArray{Tsol,1},
              uR::AbstractArray{Tsol,1},
              dxidxL,     # TODO: type
              jacL,       # TODO: type
              dxidxR,     # TODO: type
              jacR,       # TODO: type
              face,       # TODO: type
              F::AbstractVector{Tres}) where {Tsol, Tres}

  # calcViscousFlux_SIPG(params, uL, uR, aux_vars, nrm, F)    # this is the inviscid flux signature, needs to be changed
  calcViscousFlux_SIPG(params, sbp, sbpface, uL, uR, dxidxL, jacL, dxidxR, jacR, face, F)
  return nothing
end


"""
### NavierStokesMod.FluxDict

  This dictonary maps the names of the fluxes (Strings) to the
  functor object itself.  All flux functors should be added to the dictionary.

  TODO: document signature of the functors here

"""
global const FluxDict = Dict{String, FluxType}(
"SIPGViscousFlux" => SIPGViscousFlux(),
"ErrorFlux" => ErrorFlux(),
)

"""
### NavierStokesMod.getFluxFunctors

  This function retrieves the flux functors from the dictonary and
  stores them to eqn.viscous_flux_func.

  Inputs:
    mesh: an AbstractDGMesh
    sbp
    eqn
    opts
"""
function getFluxFunctors(mesh::AbstractDGMesh, sbp, eqn, opts)


  name = opts["Viscous_flux_name"]
  eqn.viscous_flux_func = FluxDict[name]

  return nothing
end

