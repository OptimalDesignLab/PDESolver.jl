# export evalElliptic, init, testOperatorGradientInterpolation, dataPrep, interpolateFace

import PDESolver.evalResidual

function init(mesh::AbstractMesh{Tmsh},
              sbp::AbstractSBP,
              eqn::AbstractEllipticData{Tsol, Tres},
              opts,
              pmesh=mesh) where {Tmsh, Tsol, Tres}
  #
  # TODO: initMPIStructures(mesh, opts)
  #
  getBCFunctors(mesh, sbp, eqn, opts)
  getBCFunctors(pmesh, sbp, eqn, opts)
  # getFluxFunctors(mesh, sbp, eqn, opts)
  getSrcFuntors(mesh, sbp, eqn, opts)
  getDiffnFunc(mesh, sbp, eqn, opts)
  if haskey(opts, "Functional")
    getFunctional(mesh, sbp, eqn, opts)
  end

  calcDiffn(mesh, sbp, eqn, eqn.diffusion_func, eqn.lambda)
  dim = size(mesh.coords, 1)

  # interpolate diffusion from elements to faces
  for d2 = 1 : dim
    for d1 = 1 : dim
      lambda_elem = Base.view(eqn.lambda, d1, d2, :, :, :)
      lambda_face = Base.view(eqn.lambda_face, d1, d2, :, :, :, :)
      lambda_bndry = Base.view(eqn.lambda_bndry, d1, d2, :, :, :)
      interpolateFace(mesh, sbp, eqn, opts, lambda_elem, lambda_face)
      interpolateBoundary(mesh, sbp, eqn, opts, lambda_elem, lambda_bndry)


      # lambda_elem[:,:,:] = eqn.lambda[d1, d2, :,:,:]
      # interpolateFace(mesh, sbp, eqn, opts, lambda_elem, lambda_face)
      # interpolateBoundary(mesh, sbp, eqn, opts, lambda_elem, lambda_bndry)
      # eqn.lambda_face[d1, d2, :,:,:,:] = lambda_face[:,:,:,:]
      # eqn.lambda_bndry[d1, d2, :,:,:] = lambda_bndry[:,:,:]
    end
  end
  return nothing
end

function evalResidual(mesh::AbstractMesh,
                      sbp::AbstractSBP,
                      eqn::EllipticData,
                      opts::Dict,
                      t=0.0)
  dataPrep(mesh, sbp, eqn, opts)
  evalVolumeIntegrals(mesh, sbp, eqn, opts)
  evalFaceIntegrals(mesh, sbp, eqn, opts)
  evalBoundaryIntegrals(mesh, sbp, eqn)

  for e = 1:mesh.numEl
    for n = 1:mesh.numNodesPerElement
      for v = 1:mesh.numDofPerNode
        # eqn.res[v,n,e] = -1.0*eqn.res[v,n,e]
      end
    end
  end

  # if haskey(opts, "TimeAdvance") && opts["TimeAdvance"] == "BDF2"
  # BDF2(mesh, sbp, eqn, opts)
  # end
  #
  # TODO: parallel commmunication
  #
  # if mesh.commsize > 1
  #     evalSharedFaceIntegrals(mesh, sbp, eqn, opts)
  # end if
  return nothing
end


function dataPrep(mesh::AbstractMesh{Tmsh},
                         sbp::AbstractSBP,
                         eqn::AbstractEllipticData{Tsol, Tres},
                         opts) where {Tmsh, Tsol, Tres}
  #
  # TODO: Not sure if filtering is needed
  #



  fill!(eqn.res, 0.0)
  # fill!(eqn.res_edge, 0.0)

  # source term is multiplied 1/|J|
  calcSource(mesh, sbp, eqn, eqn.src_func, eqn.src)
  # calcDiffn(mesh, sbp, eqn, eqn.diffusion_func, eqn.lambda)
  calcGradient(mesh, sbp, eqn, eqn.q, eqn.q_grad)

  fill!(eqn.q_face, 0.0)
  fill!(eqn.q_grad_face, 0.0)
  fill!(eqn.q_bndry, 0.0)
  fill!(eqn.q_grad_bndry, 0.0)

  # interpolate diffusion from elements to faces
  # dim = size(mesh.coords, 1)
  # lambda_elem = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
  # lambda_face = Array(Tsol, mesh.numDofPerNode, 2, mesh.numNodesPerFace, mesh.numInterfaces)
  # lambda_bndry = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace, mesh.numBoundaryFaces)

  # for d2 = 1 : dim
  # for d1 = 1 : dim
  # lambda_elem = view(eqn.lambda, d1, d2, :, :, :)
  # lambda_face = view(eqn.lambda_face, d1, d2, :, :, :, :)
  # lambda_bndry = view(eqn.lambda_bndry, d1, d2, :, :, :)
  # interpolateFace(mesh, sbp, eqn, opts, lambda_elem, lambda_face)
  # interpolateBoundary(mesh, sbp, eqn, opts, lambda_elem, lambda_bndry)


  # # lambda_elem[:,:,:] = eqn.lambda[d1, d2, :,:,:]
  # # interpolateFace(mesh, sbp, eqn, opts, lambda_elem, lambda_face)
  # # interpolateBoundary(mesh, sbp, eqn, opts, lambda_elem, lambda_bndry)
  # # eqn.lambda_face[d1, d2, :,:,:,:] = lambda_face[:,:,:,:]
  # # eqn.lambda_bndry[d1, d2, :,:,:] = lambda_bndry[:,:,:]
  # end
  # end

  # interpolate q from volume to interior edge and boundary edge
  interpolateBoundary(mesh, sbp, eqn, opts, eqn.q, eqn.q_bndry)
  interpolateFace(mesh, sbp, eqn, opts, eqn.q, eqn.q_face)
  # interpolate grad(q) from volume to interior edge and boundary edge
  interpolateBoundary(mesh, sbp, eqn, opts, eqn.q_grad, eqn.q_grad_bndry)
  interpolateFace(mesh, sbp, eqn, opts, eqn.q_grad, eqn.q_grad_face)

  fill!(eqn.xflux_face, 0.0)
  fill!(eqn.yflux_face, 0.0)
  fill!(eqn.flux_face,  0.0)
  fill!(eqn.xflux_bndry, 0.0)
  fill!(eqn.yflux_bndry, 0.0)
  fill!(eqn.flux_bndry,  0.0)
  # eqn.flux_func = FluxType["SIPG"]
  calcFaceFlux(mesh, sbp, eqn, opts, mesh.interfaces,
               eqn.xflux_face, eqn.yflux_face, eqn.flux_face)
  getBCFluxes(mesh, sbp, eqn, opts)

  return nothing
end

function evalVolumeIntegrals(mesh::AbstractMesh{Tmsh},
                             sbp::AbstractSBP,
                             eqn::EllipticData{Tsol, Tres, Tdim},
                             opts) where {Tmsh, Tsol, Tres, Tdim}
  #
  # integral of term-1
  #
  weakdifferentiate2!(mesh, sbp, eqn, eqn.q_grad, eqn.res)
  # integral of term-6, i.e., source term
  strongintegrate(mesh, sbp, eqn.src, eqn.res)
end

function strongintegrate(mesh::AbstractMesh{Tmsh},
                         sbp::AbstractSBP{Tsbp},
                         f::AbstractArray{Tflx,3},
                         res::AbstractArray{Tres,3}) where {Tmsh, Tsbp, Tflx, Tres}
  for elem = 1:mesh.numEl
    jac = sview(mesh.jac, :, elem)
    for i = 1:mesh.numNodesPerElement
      res[:, i, elem] -= sbp.w[i]/jac[i]*f[:, i, elem]
    end
  end

end
#
# Volume integration of term looking like:
#     \int ∇v⋅(Λ∇f) dΩ
# where Λ is matrix coefficient
#
function weakdifferentiate2!(mesh::AbstractMesh{Tmsh},
                                    sbp::AbstractSBP{Tsbp},
                                    eqn::EllipticData{Tsol, Tres, Tdim},
                                    q_grad::AbstractArray{Tflx,4},
                                    res::AbstractArray{Tres,3}) where {Tmsh, Tsbp,Tflx,Tsol, Tres, Tdim}
  @assert (sbp.numnodes == size(q_grad, 2))
  @assert (sbp.numnodes == size(res,2))
  @assert (size(q_grad, 1) == size(res,1))    # numDofPerNode
  @assert (size(mesh.dxidx, 1) == size(q_grad, 4))
  dim           = size(q_grad, 4)
  nElems        = size(q_grad, 3)
  nNodesPerElem = size(q_grad, 2)
  nDofsPerNode  = size(q_grad, 1)
  Qx = Array(Tsbp, nNodesPerElem, nNodesPerElem, dim)
  lambda_dqdx = Array(Tres, dim, nDofsPerNode, nNodesPerElem)
  w = sview(sbp.w, :)
  for elem=1:nElems
    # compute Λ∇q
    lambda = sview(eqn.lambda, :,:,:,:, elem)
    lambda_dqdx[:,:,:] = 0.0
    for d1 = 1:dim
      for d2 = 1:dim
        for node = 1:nNodesPerElem
          for ivar = 1:nDofsPerNode
            lambda_dqdx[d1, ivar, node] += lambda[d1, d2, ivar, node]*q_grad[ivar, node, elem, d2]  
          end
        end
      end
    end
    jac = sview(mesh.jac, :, elem)
    calcQx(mesh, sbp, Int(elem), Qx)
    for d = 1:dim
      for i = 1:sbp.numnodes
        for j = 1:sbp.numnodes
          for ivar = 1:nDofsPerNode
            # res[ivar, i, elem] += Qx[j,i,d]*q_grad[ivar, j, elem, d]/jac[j]
            res[ivar, i, elem] += Qx[j,i,d]*lambda_dqdx[d, ivar, j]/jac[j]
          end
        end
      end
    end
  end
end

#
# Integration over interior interfaces. The integration is categorized into 2
# parts.
#
# TODO: put the integral of ∫∇v ⋅ (fx, fy) dΓ into a function,
#       let's say interiorfaceintegrate2
#
function evalFaceIntegrals(
                                  mesh::AbstractDGMesh{Tmsh},
                                  sbp::AbstractSBP,
                                  eqn::EllipticData{Tsol, Tres, Tdim},
                                  opts) where {Tmsh, Tsol, Tres, Tdim}

  #
  # Integrate face_flux
  #
  interiorfaceintegrate!(mesh.sbpface, mesh.interfaces, eqn.flux_face, eqn.res)

  #
  # Integrate face_flux
  #
  sbpface = mesh.sbpface
  DxL = Array(Tmsh, mesh.numNodesPerElement, mesh.numNodesPerElement, Tdim)
  DxR = Array(Tmsh, mesh.numNodesPerElement, mesh.numNodesPerElement, Tdim)

  GxL = Array(Tmsh, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numNodesPerElement, Tdim)
  GxR = Array(Tmsh, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numNodesPerElement, Tdim)
  # interpolation matrix, keep in mind R is transfered for efficiency
  R = sview(sbpface.interp[:,:])
  w = sview(sbpface.wface, :)
  res = sview(eqn.res, :,:,:)

  nv = mesh.numNodesPerElement    # number of Nodes per elemet
  nf = mesh.numNodesPerFace       # number of nodes on interfaces
  ns = sbpface.stencilsize        # size of stencil for interpolation

  # loop over all the boundaries
  nfaces = length(mesh.interfaces)
  for f = 1:nfaces
    face = mesh.interfaces[f]

    xflux = sview(eqn.xflux_face, :,:,f)
    yflux = sview(eqn.yflux_face, :,:,f)
    elemL = face.elementL
    elemR = face.elementR
    faceL = face.faceL
    faceR = face.faceR
    pL = sview(sbpface.perm, :, faceL)
    pR = sview(sbpface.perm, :, faceR)

    calcDx(mesh, sbp, Int(elemL), DxL)
    calcDx(mesh, sbp, Int(elemR), DxR)
    #
    # Compute operator G = Λ∇
    #
    lambdaL = sview(eqn.lambda, :, :, :, :, elemL)
    lambdaR = sview(eqn.lambda, :, :, :, :, elemR)
    GxL[:,:,:,:] = 0.0
    GxR[:,:,:,:] = 0.0
    for d = 1:Tdim
      for i = 1:mesh.numNodesPerElement
        for j  = 1:mesh.numNodesPerElement
          GxL[:, i, j, d] = lambdaL[d, 1, :, i]*DxL[i,j,1] + lambdaL[d, 2, :, i]*DxL[i,j,2] 
          GxR[:, i, j, d] = lambdaR[d, 1, :, i]*DxR[i,j,1] + lambdaR[d, 2, :, i]*DxR[i,j,2] 
        end
      end
    end

    for i = 1 : nv
      for j = 1 : ns
        for k = 1 : nf
          kR = sbpface.nbrperm[k, face.orient]
          for ivar = 1:mesh.numDofPerNode
            res[ivar, i, elemL] += GxL[ivar, pL[j], i, 1]*R[j, k]*w[k]*xflux[ivar, k]
            res[ivar, i, elemL] += GxL[ivar, pL[j], i, 2]*R[j, k]*w[k]*yflux[ivar, k]
            res[ivar, i, elemR] += GxR[ivar, pR[j], i, 1]*R[j,kR]*w[k]*xflux[ivar, k]
            res[ivar, i, elemR] += GxR[ivar, pR[j], i, 2]*R[j,kR]*w[k]*yflux[ivar, k]
          end
        end
      end
    end

  end
  return nothing
end

function evalBoundaryIntegrals(mesh::AbstractMesh{Tmsh},
                               sbp::AbstractSBP,
                               eqn::EllipticData{Tsol, Tres, Tdim}) where {Tmsh, Tsol, Tres, Tdim}
  boundaryintegrate!(mesh.sbpface, mesh.bndryfaces, eqn.flux_bndry, eqn.res)

  sbpface = mesh.sbpface
  Dx = Array(Tmsh, (mesh.numNodesPerElement, mesh.numNodesPerElement, Tdim))
  # interpolation matrix, keep in mind R is transfered for efficiency
  R = sview(sbpface.interp[:,:])
  w = sview(sbpface.wface, :)
  res = sview(eqn.res, :,:,:)

  nv = mesh.numNodesPerElement
  nf = mesh.numNodesPerFace
  ns = sbpface.stencilsize
  Gx = Array(Tmsh, (mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numNodesPerElement, Tdim))
  # loop over all the boundaries
  for bc = 1:mesh.numBC
    indx0 = mesh.bndry_offsets[bc]
    indx1 = mesh.bndry_offsets[bc+1] - 1
    for f = indx0:indx1
      xflux = sview(eqn.xflux_bndry, :,:,f)
      yflux = sview(eqn.yflux_bndry, :,:,f)
      bndry = mesh.bndryfaces[f]
      elem = bndry.element
      face = bndry.face
      p = sview(sbpface.perm, :, face)

      calcDx(mesh, sbp, Int(elem), Dx)

      lambda = sview(eqn.lambda, :, :, :, :, elem)
      Gx[:,:,:] = 0.0
      for d = 1:Tdim
        for i = 1:mesh.numNodesPerElement
          for j  = 1:mesh.numNodesPerElement
            Gx[:, i, j, d] = lambda[d, 1, :, i]*Dx[i,j,1] + lambda[d, 2, :, i]*Dx[i,j,2]  
          end
        end
      end

      for i = 1 : nv      # loop over volume nodes
        for j = 1 : ns      # loop over stencils
          for k = 1 : nf      # loop over face nodes
            for ivar = 1:mesh.numDofPerNode
              res[ivar, i, elem] += Gx[ivar, p[j], i, 1]*R[j ,k]*w[k]*xflux[ivar, k]
              res[ivar, i, elem] += Gx[ivar, p[j], i, 2]*R[j ,k]*w[k]*yflux[ivar, k]
            end
          end
        end
      end
    end
  end

  return nothing
end  # end evalBoundaryIntegrals
#
# Compute Dx
#
function calcDx(mesh::AbstractMesh{Tmsh},
                sbp::AbstractSBP,
                elem::Int,
                Dx::AbstractArray{Tmsh, 3}) where Tmsh
  @assert(size(Dx, 1) == mesh.numNodesPerElement)
  @assert(size(Dx, 2) == mesh.numNodesPerElement)
  @assert(size(Dx, 3) == size(mesh.dxidx, 1))
  dxidx = sview(mesh.dxidx, :,:,:,elem) # (dim, dim, numNodesPerElement)
  jac = mesh.jac[:, elem]
  Dx[:,:,:] = 0.0
  dim = size(Dx, 3)
  for d=1:dim
    for dd=1:dim
      for l=1:mesh.numNodesPerElement
        for m = 1:mesh.numNodesPerElement
          # Since dxidx is scaled by 1/|J|, we need to get it back,
          # that's why jac is here
          Dx[l,m,d] += dxidx[dd,d,l]*jac[l]*sbp.Q[l,m,dd]
        end
      end
    end
    # Until here Dx stores Qx, we need to left multiply H^(-1)
    for l=1:mesh.numNodesPerElement
      Dx[l,:,d] /= sbp.w[l]
    end
  end
end

function calcQx(mesh::AbstractMesh{Tmsh},
                sbp::AbstractSBP,
                elem::Int,
                Qx::AbstractArray{Tmsh, 3}) where Tmsh
  @assert(size(Qx, 1) == mesh.numNodesPerElement)
  @assert(size(Qx, 2) == mesh.numNodesPerElement)
  @assert(size(Qx, 3) == size(mesh.dxidx, 1))
  dxidx = sview(mesh.dxidx, :,:,:,elem) # (dim, dim, numNodesPerElement)
  jac = mesh.jac[:, elem]
  Qx[:,:,:] = 0.0
  dim = size(Qx, 3)
  for d=1:dim
    for dd=1:dim
      for l=1:mesh.numNodesPerElement
        # Since dxidx is scaled by 1/|J|, we need to get it back,
        # that's why jac is here
        for m = 1:mesh.numNodesPerElement
          Qx[l,m,d] += dxidx[dd,d,l]*jac[l]*sbp.Q[l,m,dd]
        end
      end
    end
  end
end

