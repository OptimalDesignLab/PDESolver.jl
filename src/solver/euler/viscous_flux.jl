include("viscous_penalty.jl")

@doc """

Calculate fluxes at edge cubature points using face-based form.
eqn.flux_face, eqn.xflux, eqn.yflux will be updated.

Input:
  mesh :
  sbp  :
  eqn  :
  opts :
  peeridx: 0 for regular interior
             1-n for face on interface 1-n

Output:

# """->
function calcViscousFlux_interior{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractDGMesh{Tmsh},
                                                          sbp::AbstractSBP,
                                                          eqn::EulerData{Tsol, Tres, Tdim},
                                                          opts::Dict,
                                                          peeridx::Int=0)

  # DEBUGAA
  #=
  println(eqn.params.f, "Entered calcViscousFlux_interior")
  println(eqn.params.f, " size(mesh.shared_interfaces):", size(mesh.shared_interfaces))
  println(eqn.params.f, " size(eqn.flux_sharedface): ", size(eqn.flux_sharedface))
  =#

  # note about eqn.shared_data: when this function is called by
  #   calcSharedFaceIntegrals_nopre_element_inner, finishExchangeData will have made the proper
  #   alterations to eqn.shared_data, so it is safe to call directly. See line 192 in Utils/parallel.jl

  Ma      = eqn.params.Ma
  Re      = eqn.params.Re
  gamma_1 = eqn.params.gamma_1
  Pr      = 0.72
  coef_nondim = Ma/Re

  if peeridx == 0      # if this function is being asked to compute viscous fluxes across interior
                        #   interfaces, all local to a partition
    interfaces  = sview(mesh.interfaces, :)
    nfaces      = length(mesh.interfaces)

    flux_location = eqn.flux_face
  else                  # if this function is being asked to compute viscous fluxes across interfaces
                        #   between partitions

    interfaces = mesh.shared_interfaces[peeridx]
    nfaces = length(interfaces)

    flux_location = eqn.flux_sharedface[peeridx]

  end

  p    = opts["order"]
  dq   = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  dqn  = Array(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)
  GtL  = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, Tdim, Tdim, mesh.numNodesPerFace)
  GtR  = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, Tdim, Tdim, mesh.numNodesPerFace)
  pMat  = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerFace)
  Fv_faceL = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)
  Fv_faceR = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)
  Fv_avg   = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)
  vecfluxL = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)
  vecfluxR = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)

  sbpface = mesh.sbpface            # This should be ok in parallel. sbpface does not change
                                    #   across partitions or even across the mesh since we do not
                                    #   have mixed elements. It simply describes the face between
                                    #   two sbp elements.
  sat_type = opts["SAT_type"]
  # this one is Hartmann's definition
  const_tii = (p + 1.0)*(p + Tdim)/(2.0*Tdim)
  # const_tii = calcTraceInverseInequalityConst(sbp, sbpface)
  area_sum = sview(eqn.area_sum, :)

  params = eqn.params
  # sigma = calcTraceInverseInequalityConst(sbp, sbpface)
  # println("rho_max = ", sigma)

  # peeridx explanation:
  #   peeridx contains the index in a range of peers (starting at 1, ending at npeers).
  #   If peeridx is 0, then perform some action on the local part, using serial logic.
  #   If it is anything but 0, then we are working with one of the peers, and then parallel logic
  #     must be used.
  if peeridx != 0         # AAAAA3
    start_elnum = mesh.shared_element_offsets[peeridx]
  end
  # TODO (speed): make start_elnum = -1 in the peeridx == 0 case, so that elemR is face.elementR? see lines 100-104

  for f = 1:nfaces    # loop over faces
    # println(eqn.params.f, "----- in loop f = 1:nfaces, f = ", f, "  nfaces = ", nfaces,"  peeridx = ", peeridx, " -----")

    flux  = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
    face = interfaces[f]          # AA: think this section is OK for parallel, because interfaces is set properly above
    elemL = face.elementL
    if peeridx == 0               # AAAAA3: forgot mesh.shared_element_offsets[peeridx]
      elemR = face.elementR
    else
      elemR = face.elementR - start_elnum + 1
    end
    faceL = face.faceL
    faceR = face.faceR
    permL = sview(sbpface.perm, :, faceL)
    permR = sview(sbpface.perm, :, faceR)

    # AAAAA4: parallelized properly. mesh.nrm_sharedface
    if peeridx == 0
      nrm_xy = ro_sview(mesh.nrm_face, :, :, f)
    else
      nrm_xy = ro_sview(mesh.nrm_sharedface[peeridx], :, :, f)
    end

    # DEBUGAA
    # println(eqn.params.f, "peeridx: ", peeridx, ", elemL: ", elemL, ", elemR: ", elemR, ", f: ", f, ", nrm_xy: ", nrm_xy)

    # We need viscous flux and diffusion tensor on interfaces, and there
    # are different ways to compute them. For viscous flux:
    # 1) Fv = Fv(q, ∇q). Since we already have Q on face nodes, if 𝛻Q is also available 
    # on interface, then we are done. This way is probably more consistent with computation
    # of other terms like Fv(q_b, ∇ϕ)
    # 2) we can compute viscous flux on volume nodes and then interpolate to interface nodes.
    # It's logically simple but computationally expensive.

    # compute viscous flux and diffusion tensor
    q_faceL = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
    q_faceR = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)

    # q_face{L,R} computed below this if statement now

    if peeridx == 0
      q_elemL = sview(eqn.q, :, :, elemL)
      q_elemR = sview(eqn.q, :, :, elemR)

      q_faceL = slice(eqn.q_face, :, 1, :, f)
      q_faceR = slice(eqn.q_face, :, 2, :, f)       # TODO: use sview instead of slice

      # replace the above with this?
      # interiorFaceInterpolate!(mesh.sbpface, face, q_elemL, q_elemR, q_faceL, q_faceR)

    else            # AAAAA parallelized
      data = eqn.shared_data[peeridx]
      q_elemL = ro_sview(eqn.q, :, :, elemL)
      # q_elemR: for comparison, see line 686 & 703 in flux.jl
      q_elemR = ro_sview(data.q_recv, :, :, elemR)
      # q_face*: for comparison, see flux.jl, line 702.

      interiorFaceInterpolate!(mesh.sbpface, face, q_elemL, q_elemR, q_faceL, q_faceR)
      # face: same as iface_j in calcSharedFaceIntegral_nopre_element_inner.
      #   face = interfaces[j] where j = 1:length(mesh.interfaces)
    end
    calcDiffusionTensor(eqn.params, q_faceL, GtL)
    calcDiffusionTensor(eqn.params, q_faceR, GtR)

    # one way to compute Fv_face 
    # calcFvis_interiorFace(mesh, sbp, f, q_elemL, q_elemR, Fv_face)    

    # compute the face derivatives first, i.e., we first compute 
    # the derivatives at element nodes, and then do the interpolation.
    dqdx_face  = Array(Tsol, Tdim, mesh.numDofPerNode, 2, mesh.numNodesPerFace)
    dqdx_elemL = Array(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerElement)
    dqdx_elemR = Array(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerElement)
    calcGradient(mesh, sbp, elemL, q_elemL, dqdx_elemL)
    calcGradient(mesh, sbp, elemR, q_elemR, dqdx_elemR)
    # TODO: elemR? set above as (in peeridx > 0 case) 
    #     elemR = face.elementR - start_elnum + 1

    #
    # TODO: we can consider the first 2 dimension as a single dimension,
    # then we will not need slice here any more.
    #
    # for d = 1 : Tdim
      # dqdxL = slice(dqdx_elemL, d, :, :)
      # dqdxR = slice(dqdx_elemR, d, :, :)
      # dqdx_f = slice(dqdx_face, d, :, :, :)
      # interiorfaceinterpolate(sbpface, face, dqdxL, dqdxR, dqdx_f)
    # end

    # # Now both G and dqdx are available at face nodes
    # dqdx_faceL = slice(dqdx_face, :, :, 1, :)
    # dqdx_faceR = slice(dqdx_face, :, :, 2, :)
    # calcFvis(params, GtL, dqdx_faceL, Fv_faceL)
    # calcFvis(params, GtR, dqdx_faceR, Fv_faceR)

    Fv_face = zeros(Tsol, Tdim, mesh.numDofPerNode, 2, mesh.numNodesPerFace)
    # Here is where we need to decide where dxidxL, dxidxR, jacL, and jacR come from.
    # If this interface is between elements on the same partition, then the former method
    # (using ro_sview of mesh.jac and mesh.dxidx) works. If it is between partitions,
    # we need to obtain the data from RemoteMetrics.
    jacL = ro_sview(mesh.jac, :, elemL)
    dxidxL = ro_sview(mesh.dxidx, :,:,:,elemL)
    # AAAAA parallelized
    if peeridx == 0
      jacR = ro_sview(mesh.jac, :, elemR)
      dxidxR = ro_sview(mesh.dxidx, :,:,:,elemR)
    else
      dxidxR = ro_sview(mesh.remote_metrics[peeridx].dxidx, :, :, :, elemR)
      jacR = ro_sview(mesh.remote_metrics[peeridx].jac, :, elemR)
    end

    # TODO: in peeridx > 0 case, Fv_faceR?
    calcFaceFvis(params, sbp, sbpface, q_elemL, q_elemR, dxidxL, jacL, dxidxR, jacR, face, Fv_face)
    Fv_faceL = sview(Fv_face, :,:,1,:)
    Fv_faceR = sview(Fv_face, :,:,2,:)
    # diffL = maximum(abs(real(slice(Fv_face, :, :, 1, :) - Fv_faceL)))
    # diffR = maximum(abs(real(slice(Fv_face, :, :, 2, :) - Fv_faceR)))
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

    #------------------------------------------------------------
    # Method 2: peeridx if statement outside outermost loop
    # accumulate fluxes
    #=
    if peeridx == 0
      for n = 1:mesh.numNodesPerFace
        for iDof = 2 : Tdim+2
          for iDim = 1 : Tdim
            eqn.vecflux_faceL[iDim, iDof, n, f] -=  vecfluxL[iDim, iDof, n]*coef_nondim
            eqn.vecflux_faceR[iDim, iDof, n, f] -=  vecfluxR[iDim, iDof, n]*coef_nondim
          end   # end 'for iDim = 1 : Tdim'
          flux_location[iDof, n, f] += flux[iDof, n]*coef_nondim
        end   # end 'for iDof = 2 : Tdim+2'
      end   # end 'for n = 1:mesh.numNodesPerFace'
    else    # if peeridx != 0
      for n = 1:mesh.numNodesPerFace
        for iDof = 2 : Tdim+2
          for iDim = 1 : Tdim
            eqn.vecflux_faceL_shared[peeridx][iDim, iDof, n, f] -= vecfluxL[iDim, iDof, n]*coef_nondim
          end   # end 'for iDim = 1 : Tdim'
          flux_location[iDof, n, f] += flux[iDof, n]*coef_nondim
        end   # end 'for iDof = 2 : Tdim+2'
      end   # end 'for n = 1:mesh.numNodesPerFace'
    end   # end 'if peeridx == 0'
    =#
    #------------------------------------------------------------




    #------------------------------------------------------------
    # Method 1: peeridx if statement inside innermost loop
    # accumulate fluxes
    for n = 1:mesh.numNodesPerFace
      for iDof = 2 : Tdim+2
        for iDim = 1 : Tdim
          if peeridx == 0       # TODO (speed) put if statement outside of loop
            # DEBUGAA
            # println(eqn.params.f, "inside cVF_i, acc into vecflux. size(eqn.vecflux_faceL):", size(eqn.vecflux_faceL), " iDim:", 
                    # iDim, " iDof:", iDof, " n:", n, " f:", f, " size(vecfluxL):", size(vecfluxL))

            eqn.vecflux_faceL[iDim, iDof, n, f] -=  vecfluxL[iDim, iDof, n]*coef_nondim
            eqn.vecflux_faceR[iDim, iDof, n, f] -=  vecfluxR[iDim, iDof, n]*coef_nondim
          else
            # AAAAA2: different array for shared case
            eqn.vecflux_faceL_shared[peeridx][iDim, iDof, n, f] -= vecfluxL[iDim, iDof, n]*coef_nondim
          end   # end 'if peeridx == 0'
        end   # end 'for iDim = 1 : Tdim'

        # This line accumulates the scalar part of the numerical viscous flux here.
        #   *eqn.flux_face* which also stores inviscid flux and will be integrated and assembled
        #   into residual. (edited)
        #   remember there are two groups of viscous fluxes, and this is the *scalar* one.
        #   *evalFaceIntegrals_vector* takes care of the *vector* one.

        flux_location[iDof, n, f] += flux[iDof, n]*coef_nondim
      end   # end 'for iDof = 2 : Tdim+2'
    end   # end 'for n = 1:mesh.numNodesPerFace'
    #------------------------------------------------------------

  end     # end of loop over all interfaces

  #=
  if mesh.myrank == 0
    println(" ")
    println(" (in calcViscousFlux_interior)")
    for f_ix = 1:nfaces
      println(" eqn.vecflux_faceL[:, :, :, $f_ix]: ", eqn.vecflux_faceL[:, :, :, f_ix])
    end
    println(" ")
    if mesh.commsize == 1
      for f_ix = 1:nfaces
        println(" eqn.vecflux_faceR[:, :, :, $f_ix]: ", eqn.vecflux_faceR[:, :, :, f_ix])
      end
    elseif mesh.commsize > 1
      if peeridx == 1
      # println(" size(eqn.vecflux_faceL_shared): ", size(eqn.vecflux_faceL_shared))
      # println(" peeridx: ", peeridx)
      # println(" size(eqn.vecflux_faceL_shared[peeridx]): ", size(eqn.vecflux_faceL_shared[peeridx]))
        for f_ix = 1:nfaces
          println(" eqn.vecflux_faceL_shared[peeridx][:, :, :, $f_ix]: ", eqn.vecflux_faceL_shared[peeridx][:, :, :, f_ix])
        end
      end
    end
    println(" ")
  end
  =#

  return nothing
end # end of function calcViscousFlux_interior


function calcViscousFlux_boundary{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh},
                                                          sbp::AbstractSBP,
                                                          eqn::EulerData{Tsol, Tres, Tdim},
                                                          opts::Dict)
  # freestream info
  Ma = eqn.params.Ma
  Re = eqn.params.Re
  gamma_1 = eqn.params.gamma_1
  Pr = 0.72
  coef_nondim = Ma/Re

  p = opts["order"]
  sat_type = opts["SAT_type"]
  const_tii = (p + 1.0)*(p + Tdim)/Tdim
  sbpface = mesh.sbpface
  dq    = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)    
  dqn   = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)    
  q_bnd = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)    
  pMat  = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerFace)
  Gt    = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, Tdim, Tdim, mesh.numNodesPerFace)
  Gt_bnd  = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, Tdim, Tdim, mesh.numNodesPerFace)
  Fv_face = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)
  Fv_bnd  = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)
  vecflux = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)

  nrm1 = Array(Tmsh, Tdim, mesh.numNodesPerFace)
  area = Array(Tmsh, mesh.numNodesPerFace)
  area_sum = sview(eqn.area_sum, :)

  # sigma = calcTraceInverseInequalityConst(sbp, sbpface)
  dqdx_elem = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerElement )
  dqdx_face = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace )
  for iBC = 1 : mesh.numBC
    indx0 = mesh.bndry_offsets[iBC]
    indx1 = mesh.bndry_offsets[iBC+1] - 1

    # specify boundary value function
    # TODO: Put it into a function 
    bnd_functor::AbstractBoundaryValueType
    key_i = string("BC", iBC, "_name")
    val = opts[key_i]
    Gt_functor = calcDiffusionTensor
    if val == "FreeStreamBC"
      bnd_functor = Farfield()
    elseif val == "ExactChannelBC"
      bnd_functor = ExactChannel()
    elseif val == "nonslipBC"
      bnd_functor = AdiabaticWall()
      Gt_functor = calcDiffusionTensorOnAdiabaticWall
    elseif val == "noPenetrationBC"
      continue
    elseif val == "zeroPressGradientBC"
      bnd_functor = Farfield()
    else
      error("iBC = ", iBC, ", Only 'FreeStreamBC' and 'nonslipBC' available")
    end

    for f = indx0 : indx1

      flux  = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
      bndry = mesh.bndryfaces[f]
      elem = bndry.element
      face = bndry.face
      perm = sview(sbpface.perm, :, face)
      xy = sview(mesh.coords_bndry, :, :, f)

      # Compute geometric info on face
      nrm_xy = ro_sview(mesh.nrm_bndry, :, :, f)
      for n = 1 : mesh.numNodesPerFace
        area[n] = norm(ro_sview(nrm_xy, :, n)) 

        for i = 1 : Tdim
          nrm1[i,n] = nrm_xy[i,n] / area[n]
        end
      end

      # We need viscous flux and diffusion tensor on interfaces, and there
      # are different ways to compute them. For viscous flux:
      # 1) since we have Q on face nodes, if 𝛻Q is available on interface, then we are done.
      # 2) we can comoute viscous flux on volume nodes and then interpolate to interface node.
      # It's logically simple but computationally expensive.

      # Compute boundary viscous flux, F(q_b, ∇q) = G(q_b)∇q.
      # so we need viscousity tensor G, and derivative of q.
      q_face = sview(eqn.q_bndry, :, :, f)
      bnd_functor(q_face, xy, nrm1, eqn.params, q_bnd)

      # diffusion matrix used in penalty term should be computed from q_face rather than q_bnd
      if val == "nonslipBC"
        Gt_functor(eqn.params, q_bnd, nrm1, Gt)
      else
        Gt_functor(eqn.params, q_bnd, Gt)
      end
      q_elem = sview(eqn.q, :, :, elem)
      calcGradient(mesh, sbp, elem, q_elem, dqdx_elem)

      #
      # TODO: we can consider the first 2 dimension as a single dimension,
      # then we will not need slice here any more.
      #
      for d = 1 : Tdim
        q_x_node = slice(dqdx_elem, d, :, :)        # TODO: use sview instead of slice, AA
        q_x_face = slice(dqdx_face, d, :, :)
        boundaryinterpolate(sbpface, bndry, q_x_node, q_x_face) 
      end

      calcFvis(eqn.params, Gt, dqdx_face, Fv_face)

      # compute penalty matrix
      cmptBPMat(mesh, sbp, eqn, opts, f, Gt, pMat)

      # Start to compute fluxes.  We have 3 terms on interfaces:
      # 1) -{Fv}⋅[ϕ]
      # 2) -{G^T ∇ϕ}⋅[q] 
      # 3) +δ{G}[q]:[ϕ]
      for n = 1 : mesh.numNodesPerFace
        for iDof = 1 : mesh.numDofPerNode
          dq[iDof, n] = q_face[iDof, n] - q_bnd[iDof, n]
          # dqn[1, iDof, n] = -dq[iDof, n]*nrm_xy[1, n]
          # dqn[2, iDof, n] = -dq[iDof, n]*nrm_xy[2, n]
        end
      end
      # DEBUG BEGIN
      # if maximum(abs(real(dq))) > 1.e-11
        # println(real(dq))
      # end
      # DEBUG END


      # -----------------------------------------------
      # This part computes the contribution of
      # ∫ {G^T∇ϕ}:[q] dΓ = ∫ ∇ϕ⋅F dΓ , 
      # where 
      # [q] = (q+ - q-) ⊗ n, 
      # G = G(q_b) depends on boudanry value.
      # Then we can consider Δq⊗n as ∇q and F as viscous flux.
      # -----------------------------------------------

      # calcFvis(Gt, dqn, vecflux)
      fill!(vecflux, 0.0)
      for n = 1 : mesh.numNodesPerFace
        for iDof = 1 : mesh.numDofPerNode
          for iDim = 1 : Tdim
            for jDim = 1 : Tdim
              tmp = 0.0
              for jDof = 1 : mesh.numDofPerNode
                tmp += Gt[iDof, jDof, iDim, jDim, n]
              end
              vecflux[iDim, iDof, n] += tmp * nrm_xy[jDim, n]
            end
            vecflux[iDim,iDof,n] *=  dq[iDof,n]
          end
        end
      end

      for n = 1 : mesh.numNodesPerFace
        for iDof = 2 : mesh.numDofPerNode
          for iDim = 1 : Tdim
            flux[iDof, n] -= Fv_face[iDim, iDof, n]*nrm_xy[iDim,n] 
          end
        end
      end

      for n = 1 : mesh.numNodesPerFace
        for iDof = 2 : mesh.numDofPerNode
          for jDof = 1: mesh.numDofPerNode
            flux[iDof, n] +=  pMat[iDof, jDof, n]*dq[jDof, n]
          end
        end
      end

      # accumulate fluxes
      for n = 1 : mesh.numNodesPerFace
        for iDof = 2 : Tdim+2
          for iDim = 1 : Tdim
            eqn.vecflux_bndry[iDim, iDof, n, f] -=  vecflux[iDim, iDof, n]*coef_nondim
          end
          eqn.bndryflux[iDof, n, f] += flux[iDof, n]*coef_nondim
        end
      end
    end # loop over faces of one BC
  end # loop over BCs
  return nothing 
end


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
function evalFaceIntegrals_vector{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractDGMesh{Tmsh},
                                                          sbp::AbstractSBP,
                                                          eqn::EulerData{Tsol, Tres, Tdim},
                                                          opts::Dict,
                                                          peeridx::Int=0)

  # DEBUGAA
  # println(eqn.params.f, "Entered evalFaceIntegrals_vector")

  # This part computes ∫ ∇ϕ⋅F  dΓ, 
  sbpface = mesh.sbpface
  DxL = Array(Tmsh, mesh.numNodesPerElement, mesh.numNodesPerElement, Tdim)
  DxR = Array(Tmsh, mesh.numNodesPerElement, mesh.numNodesPerElement, Tdim)

  # TODO (speed): are these needed here?
  #=
  GtL = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, Tdim, Tdim, mesh.numNodesPerFace)
  GtR = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, Tdim, Tdim, mesh.numNodesPerFace)
  =#

  # TODO: why colon????
  R = sview(sbpface.interp[:,:])
  w = sview(sbpface.wface, :)
  res = sview(eqn.res, :,:,:)

  numNodes_elem = mesh.numNodesPerElement    # number of Nodes per elemet
  numNodes_face = mesh.numNodesPerFace       # number of nodes on interfaces
  stencilsize = sbpface.stencilsize        # size of stencil for interpolation

  RDxL = Array(Tmsh, mesh.numNodesPerFace, mesh.numNodesPerElement, Tdim)
  RDxR = Array(Tmsh, mesh.numNodesPerFace, mesh.numNodesPerElement, Tdim)
  # TODO (speed): are these needed?
  #=
  FvL  = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerFace, mesh.numNodesPerElement) # AA: unused?
  FvR  = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerFace, mesh.numNodesPerElement) # AA: unused?
  GtRDxL = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerFace, mesh.numNodesPerElement)   # AA: unused?
  GtRDxR = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerFace, mesh.numNodesPerElement)   # AA: unused?
  nrm    = Array(Tmsh, Tdim, mesh.numNodesPerFace)      # AA: unused?
  dq     =  Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)   # AA: unused?
  =#

  if peeridx == 0
    interfaces = mesh.interfaces
  else      # AAAAA parallelized
    interfaces = mesh.shared_interfaces[peeridx]
  end

  nfaces = length(interfaces)

  if peeridx != 0
    start_elnum = mesh.shared_element_offsets[peeridx]
  end

  for f = 1 : nfaces

    face = interfaces[f]

    elemL = face.elementL
    if peeridx == 0
      elemR = face.elementR
    else
      elemR = face.elementR - start_elnum + 1     # AAAAA4
    end

    faceL = face.faceL
    pL = sview(sbpface.perm, :, faceL)      # WRONG in shared face. correct in local face.
                                            # ser: [1,2]  par: [3,1]
    if peeridx == 0
      faceR = face.faceR
      pR = sview(sbpface.perm, :, faceR)
    end

    # DEBUGAA
    # println(eqn.params.f, "in vec visc flux. elemL:", elemL, " elemR:", elemR, " f:", f)

    # compute RDx
    # Dx: moves SBP's D matrix to the physical domain
    calcDx(mesh, sbp, elemL, DxL)
    if peeridx == 0
      calcDx(mesh, sbp, elemR, DxR)
    end

    for i = 1 : length(RDxL)      # TODO (speed): maybe save on computation by using zeros?
      RDxL[i] = 0.0
      if peeridx == 0
        RDxR[i] = 0.0
      end
    end

    for d =  1 : Tdim    
      for row = 1 : numNodes_face
        if peeridx == 0
          rowR = sbpface.nbrperm[row, face.orient]      # this is an integer. 1 or 2 in the 2D test case
        end
        for col = 1 : numNodes_elem
          for s = 1 : stencilsize
            RDxL[row, col, d] +=  R[s, row]  * DxL[pL[s], col, d]     
            if peeridx == 0
              RDxR[row, col, d] +=  R[s, rowR] * DxR[pR[s], col, d]
            end
          end
        end
      end
    end

    #----------------------------------------------------------------------------
    # DEBUG: bottom center interface comparison
    # DEBUGAA
    #=
    print_now = 0
    if mesh.commsize == 1   # serial
      if f == 4
        print_now = 2
      end
      if f == 5
        print_now = 1
      end
    else    # parallel
      if mesh.myrank == 0 && peeridx == 1     # if we're running parallel and if we're on a shared iface
        if f == 1
          print_now = 2
        end
      end
      if mesh.myrank == 0 && peeridx == 0
        if f == 1
          print_now = 1
        end
      end
    end

    if print_now > 0
      println("----------------------- selected values for this face ------------------------")
      if print_now == 2
        println(" SHARED FACE")
      elseif print_now == 1
        println(" LOCAL FACE")
      end
      println(" mesh.commsize: ", mesh.commsize)
      println(" peeridx: ", peeridx)
      println(" mesh.myrank: ", mesh.myrank)
      println(" f: ", f)
      println(" ")
      println(" w: ", w)
      println(" RDxL: ", RDxL)
      if print_now == 1 || mesh.commsize == 1
        println(" RDxR: ", RDxR)
        # println(" rowR: ", rowR)
      end
      println(" DxL: ", DxL)
      if print_now == 1 || mesh.commsize == 1
        println(" DxR: ", DxR)
      end
      println(" pL: ", pL)
      if print_now == 1 || mesh.commsize == 1
        println(" pR: ", pR)
      end
      println(" stencilsize: ", stencilsize)
      println("------------------------------------------------------------------------------")
    end
    =#
    #----------------------------------------------------------------------------


    # Why we need a shared version of vecflux_faceL:
    #   the local part only needs to care about its local element on the shared interface
    #   This is why there are all these "if peeridx == 0" statements around the R computations
    # TODO (speed): Should we group all R fns together around one if statement, and loop inside that?
    if peeridx == 0
      vecfluxL = sview(eqn.vecflux_faceL,:,:,:,f)     # AAAA2: if statement around this for assigning to the shared vecflux
      vecfluxR = sview(eqn.vecflux_faceR,:,:,:,f)
    else
      vecfluxL = sview(eqn.vecflux_faceL_shared[peeridx],:,:,:,f)     # AAAA2: if statement around this for assigning to the shared vecflux
    end

    # DEBUGAA
    #=
    if mesh.myrank == 0
      println(" ")
      println(" (in viscflux)")
      println(" mesh.commsize: ", mesh.commsize)
      println(" peeridx: ", peeridx)
      println(" mesh.myrank: ", mesh.myrank)
      println(" f: ", f)
      println(" print_now: ", print_now)
      if peeridx == 0
        println(" mesh.interfaces[", f, "].elementL: ", mesh.interfaces[f].elementL)
        println(" mesh.interfaces[", f, "].elementR: ", mesh.interfaces[f].elementR)
        println(" mesh.coords_interface[:, :, ", f, "]: ", mesh.coords_interface[:, :, f])
      else
        println(" mesh.shared_interfaces[1][", f, "].elementL: ", mesh.shared_interfaces[1][f].elementL)
        println(" mesh.shared_interfaces[1][", f, "].elementR: ", mesh.shared_interfaces[1][f].elementR)
        println(" mesh.coords_sharedface[1][:, :, ", f, "]: ", mesh.coords_sharedface[1][:, :, f])
      end
      println(" elemL: ", elemL)
      println(" elemR: ", elemR)
      println(" face.elementR: ", face.elementR)
      println(" ")
      if mesh.commsize > 1
        println(" size(eqn.vecflux_faceL_shared): ", size(eqn.vecflux_faceL_shared))
        println(" size(eqn.vecflux_faceL_shared[1]): ", size(eqn.vecflux_faceL_shared[1]))
        # println(" eqn.vecflux_faceL_shared: ", eqn.vecflux_faceL_shared)
      else
        # println(" eqn.vecflux_faceL: ", eqn.vecflux_faceL)
      end
      println(" ")
      println(" size(eqn.vecflux_faceL): ", size(eqn.vecflux_faceL))
      # println(" eqn.vecflux_faceL: ", eqn.vecflux_faceL)
      println(" size(vecfluxL): ", size(vecfluxL))
      println(" vecfluxL: ", vecfluxL)
      println(" ")
    end
    =#

    for i = 1 : numNodes_elem
      for j = 1 : numNodes_face
        for iDof = 2 : mesh.numDofPerNode
          # DEBUGAA
          # println(eqn.params.f, " vector flux -> res, iDof:", iDof, " i:", i, " elemL:", elemL, " elemR:", elemR)
          tmpL = 0.0
          # AAAAA2: R only needed in local interface case
          if peeridx == 0
            tmpR = 0.0
          end
          for iDim = 1 : Tdim
            tmpL += RDxL[j, i, iDim] * vecfluxL[iDim, iDof, j]
            # AAAAA2: R only needed in local interface case
            if peeridx == 0
              tmpR += RDxR[j, i, iDim] * vecfluxR[iDim, iDof, j]
            end
          end
          res[iDof, i, elemL] += tmpL * w[j]
          # AAAAA2: This will generate a bounds error in the shared case, so need an if statement. 
          if peeridx == 0
            res[iDof, i, elemR] += tmpR * w[j]
          end
        end   # end of loop 'iDof = 2 : mesh.numDofPerNode'
      end   # end of loop 'j = 1 : numNodes_face'
    end   # end of loop 'i = 1 : numNodes_elem'
  end   # end of loop 'f = 1 : nfaces'

  return nothing
end



@doc """
Now actually we are integrating 
  ∫ G∇ϕ:[q] dΓ

Input: 
  mesh
  sbp
  eqn
  opts
Output:

"""->
function evalBoundaryIntegrals_vector{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh},
                                                              sbp::AbstractSBP,
                                                              eqn::EulerData{Tsol, Tres, Tdim},
                                                              opts)

  sbpface = mesh.sbpface
  Dx = Array(Tmsh, (mesh.numNodesPerElement, mesh.numNodesPerElement, Tdim))
  R = sview(sbpface.interp[:,:])
  w = sview(sbpface.wface, :)
  res = sview(eqn.res, :,:,:)

  numNodes_elem = mesh.numNodesPerElement
  numNodes_face = mesh.numNodesPerFace
  stencilsize   = sbpface.stencilsize
  q_bnd = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  dq    = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  Gt = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, Tdim, Tdim, mesh.numNodesPerFace)
  RDx = zeros(Tmsh, mesh.numNodesPerFace, mesh.numNodesPerElement, Tdim)
  GtRDx = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerFace, mesh.numNodesPerElement) # AA: unused?
  nrm     = Array(Tmsh, Tdim, mesh.numNodesPerFace)
  nrm1 = Array(Tmsh, Tdim, mesh.numNodesPerFace)
  area = Array(Tmsh, mesh.numNodesPerFace)

  # loop over all the boundaries
  for bc = 1:mesh.numBC
    indx0 = mesh.bndry_offsets[bc]
    indx1 = mesh.bndry_offsets[bc+1] - 1

    for f = indx0:indx1
      bndry = mesh.bndryfaces[f]
      elem = bndry.element
      face = bndry.face
      p = sview(sbpface.perm, :, face)


      # compute RDx
      calcDx(mesh, sbp, elem, Dx)

      for i = 1 : length(RDx)
        RDx[i] = 0.0
      end

      for d =  1 : Tdim    
        for row = 1 : numNodes_face
          for col = 1 : numNodes_elem
            for s = 1 : stencilsize
              RDx[row, col, d] +=  R[s, row] * Dx[p[s], col, d]     
            end
          end
        end
      end

      vecflux = sview(eqn.vecflux_bndry, :,:,:,f)
      for i = 1 : numNodes_elem
        for j = 1 : numNodes_face
          for iDof = 2 : mesh.numDofPerNode
            # res[iDof, i, elem] +=  ( RDx[j, i, 1] * vecflux[1, iDof, j] 
            # + RDx[j, i, 2] * vecflux[2, iDof, j] ) * w[j]
            tmp = 0.0
            for iDim = 1 : Tdim
              tmp += RDx[j, i, iDim] * vecflux[iDim, iDof, j]
            end
            res[iDof, i, elem] += tmp * w[j]
          end
        end
      end
    end
  end

  return nothing
end  # end evalBoundaryIntegrals_vector



@doc """

Integrate ∫ ∇ϕ⋅F dΩ
TODO: consider combine it together with `weakdifferentiate`

Input:
  mesh   : 
  sbp    :
  eqn    :
  res    :
Output:

# """->
function weakdifferentiate2!{Tmsh, Tsbp, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh},
                                                           sbp::AbstractSBP{Tsbp},
                                                           eqn::EulerData{Tsol, Tres, Tdim},
                                                           res::AbstractArray{Tres,3})
  @assert (sbp.numnodes ==  size(res,2))

  dim             = Tdim
  numElems        = mesh.numEl
  numNodesPerElem = mesh.numNodesPerElement
  numDofsPerNode  = mesh.numDofPerNode

  gamma_1 = eqn.params.gamma_1
  Pr = 0.72
  Ma = eqn.params.Ma
  Re = eqn.params.Re
  coef_nondim = Ma/Re 

  Qx = Array(Tsbp, numNodesPerElem, numNodesPerElem, dim)
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

