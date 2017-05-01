
@doc """

Compute the constant coefficent in inverse trace ineqality, i.e.,
the largest eigenvalue of 
B^{1/2} R H^{-1} R^{T} B^{1/2}

Input:
sbp
Output:
cont_tii
"""->
function calcTraceInverseInequalityConst{Tsbp}(sbp::AbstractSBP{Tsbp},
                                               sbpface::AbstractFace{Tsbp})
  R = sview(sbpface.interp, :,:)
  perm = Array(Tsbp, sbp.numnodes, sbpface.stencilsize)
  Hinv = Array(Tsbp, sbp.numnodes, sbp.numnodes)
  Bsqrt = Array(Tsbp, sbpface.numnodes, sbpface.numnodes)
  HRBRH = Array(Tsbp, sbpface.numnodes, sbp.numnodes)
  BsqrtRHinvRtBsqrt = Array(Tsbp, sbpface.numnodes, sbpface.numnodes)
  for s = 1:sbpface.stencilsize
    perm[sbpface.perm[s, 1], s] = 1.0
  end
  for i = 1:sbp.numnodes
    Hinv[i,i] = 1.0/sbp.w[i]
  end
  for i = 1:sbpface.numnodes
    Bsqrt[i,i] = sqrt(sbpface.wface[i])
  end

  BsqrtRHinvRtBsqrt = Bsqrt*R.'*perm.'*Hinv*perm*R*Bsqrt 
  const_tii = eigmax(BsqrtRHinvRtBsqrt)
  return const_tii
end


@doc """

Calculate fluxes at edge cubature points using face-based form

Input:
mesh :
sbp  :
eqn  :
opts :
Output:

update eqn.flux_face, eqn.xflux, eqn.yflux
# """->
function calcViscousFlux_interior{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractDGMesh{Tmsh},
                                                          sbp::AbstractSBP,
                                                          eqn::EulerData{Tsol, Tres, Tdim},
                                                          opts)

  Ma        = eqn.params.Ma
  Re        = eqn.params.Re
  gamma_1 = eqn.params.gamma_1
  Pr        = 0.72
  coef_nondim = Ma/Re
  interfaces  = sview(mesh.interfaces, :)
  nfaces        = length(mesh.interfaces)
  p     = opts["order"]
  sat_scalar = zeros(Tsol, mesh.numNodesPerFace)
  dq     = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  dqn     = Array(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)
  nrm  = Array(Tmsh, Tdim, mesh.numNodesPerFace)
  nrm0 = Array(Tmsh, Tdim, mesh.numNodesPerFace)
  area = Array(Tmsh, mesh.numNodesPerFace)
  GtL = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, Tdim, Tdim, mesh.numNodesPerFace)
  GtR = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, Tdim, Tdim, mesh.numNodesPerFace)
  penalty = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerFace)
  Fv_faceL = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)
  Fv_faceR = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)
  Fv_avg   = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)
  vecfluxL = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)
  vecfluxR = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)

  sat_type = opts["SAT_type"]
  const_tii = (p + 1.0)*(p + Tdim)/(2.0*Tdim)
  area_sum = sview(eqn.area_sum, :)

  sbpface = mesh.sbpface
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

    # Compute geometric info on face
    for n = 1:mesh.numNodesPerFace
      dxidx = sview(mesh.dxidx_face, :, :, n, f)

      # norm vector in reference element
      nrm_xi = sview(sbpface.normal, :, faceL)
      nrm[1,n] = dxidx[1, 1]*nrm_xi[1] + dxidx[2, 1]*nrm_xi[2]
      nrm[2,n] = dxidx[1, 2]*nrm_xi[1] + dxidx[2, 2]*nrm_xi[2]

      area[n] = sqrt(nrm[1,n]*nrm[1,n] + nrm[2,n]*nrm[2,n])

      # norm vector in physical domain without any scale, ie, |nrm| = 1
      nrm0[1,n] = nrm[1,n]/area[n]
      nrm0[2,n] = nrm[2,n]/area[n]
    end

    # Compute the size of element and face, and then meas(face)/meas(elem)
    elem_volL = 0.0
    elem_volR = 0.0
    face_area = 0.0
    for n = 1:mesh.numNodesPerElement
      elem_volL +=  sbp.w[n]/mesh.jac[n, elemL]
      elem_volR +=  sbp.w[n]/mesh.jac[n, elemR]
    end

    for n = 1:mesh.numNodesPerFace
      face_area +=  sbpface.wface[n]*area[n]
    end

    heL = elem_volL/eqn.area_sum[elemL]
    heR = elem_volR/eqn.area_sum[elemR]
    he = min(heL, heR)
    area_weightL = area_sum[elemL]/face_area
    area_weightR = area_sum[elemR]/face_area

    # We need viscous flux and diffusion tensor on interfaces, and there
    # are different ways to compute them. For viscous flux:
    # 1) Fv = Fv(q, ∇q). Since we already have Q on face nodes, if 𝛻Q is also available 
    # on interface, then we are done. This way is probably more consistent with computation
    # of other terms like Fv(q_b, ∇ϕ)
    # 2) we can compute viscous flux on volume nodes and then interpolate to interface nodes.
    # It's logically simple but computationally expensive.

    # compute viscous flux and diffusion tensor
    q_faceL = slice(eqn.q_face, :, 1, :, f)
    q_faceR = slice(eqn.q_face, :, 2, :, f)
    q_elemL = sview(eqn.q, :, :, elemL)
    q_elemR = sview(eqn.q, :, :, elemR)
    calcDiffusionTensor(q_faceL, GtL)
    calcDiffusionTensor(q_faceR, GtR)

    # one way to compute Fv_face 
    # calcFvis_interiorFace(mesh, sbp, f, q_elemL, q_elemR, Fv_face)    

    # compute the face derivatives first, i.e., we first compute 
    # the derivatives at element nodes, and then do the interpolation.
    dqdx_face  = Array(Tsol, Tdim, mesh.numDofPerNode, 2, mesh.numNodesPerFace)
    dqdx_elemL = Array(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerElement)
    dqdx_elemR = Array(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerElement)
    calcGradient(mesh, sbp, elemL, q_elemL, dqdx_elemL)
    calcGradient(mesh, sbp, elemR, q_elemR, dqdx_elemR)

    for d = 1 : Tdim
      dqdxL = slice(dqdx_elemL, d, :, :)
      dqdxR = slice(dqdx_elemR, d, :, :)
      dqdx_f = slice(dqdx_face, d, :, :, :)
      interiorfaceinterpolate(sbpface, face, dqdxL, dqdxR, dqdx_f)
    end

    # Now both G and dqdx are avaiable at face nodes  
    dqdx_faceL = slice(dqdx_face, :, :, 1, :)
    dqdx_faceR = slice(dqdx_face, :, :, 2, :)
    calcFvis(GtL, dqdx_faceL, Fv_faceL)
    calcFvis(GtR, dqdx_faceR, Fv_faceR)

    # First compute penalty
    if sat_type ==  "SAT-SIPG"
      Cip = opts["Cip"]*p*p
      for n = 1:mesh.numNodesPerFace
        sat_scalar[n] = Cip * area[n] / he
      end
    else
      #    
      # Reference:
      #    Shahbazi, Mavriplis, Multigrid algorithms for high order discontinuous 
      #    Galerkin discretization of the compressible Navier-Stokes equations, 
      #    JCP, volume 228 Issue 21, November, 2009
      #

      #
      # Compute n_i G_{ij} n_j
      #
      for n = 1 : mesh.numNodesPerFace
        # for iDof = 1:mesh.numDofPerNode
        # for jDof = 1:mesh.numDofPerNode
        # penalty[iDof,jDof,n] =  nrm0[1,n]*nrm0[1,n] * (GtL[iDof,jDof,1,1,n] + GtR[iDof,jDof,1,1,n])
        # penalty[iDof,jDof,n] +=  nrm0[1,n]*nrm0[2,n] * (GtL[iDof,jDof,1,2,n] + GtR[iDof,jDof,1,2,n])
        # penalty[iDof,jDof,n] +=  nrm0[2,n]*nrm0[1,n] * (GtL[iDof,jDof,2,1,n] + GtR[iDof,jDof,2,1,n])
        # penalty[iDof,jDof,n] +=  nrm0[2,n]*nrm0[2,n] * (GtL[iDof,jDof,2,2,n] + GtR[iDof,jDof,2,2,n])
        # penalty[iDof,jDof,n] *=  0.5*const_tii * area[n] / he
        # end
        # end

        penalty[2,1,n] = ( nrm0[1,n]*nrm0[1,n] * (GtL[2,1,1,1,n] + GtR[2,1,1,1,n])
                         + nrm0[1,n]*nrm0[2,n] * (GtL[2,1,1,2,n] + GtR[2,1,1,2,n])
                         + nrm0[2,n]*nrm0[1,n] * (GtL[2,1,2,1,n] + GtR[2,1,2,1,n])
                         + nrm0[2,n]*nrm0[2,n] * (GtL[2,1,2,2,n] + GtR[2,1,2,2,n]) )

        penalty[2,2,n] = ( nrm0[1,n]*nrm0[1,n] * (GtL[2,2,1,1,n] + GtR[2,2,1,1,n])
                         + nrm0[2,n]*nrm0[2,n] * (GtL[2,2,2,2,n] + GtR[2,2,2,2,n]) )

        penalty[2,3,n] = ( nrm0[1,n]*nrm0[2,n] * (GtL[2,3,1,2,n] + GtR[2,3,1,2,n])
                         + nrm0[2,n]*nrm0[1,n] * (GtL[2,3,2,1,n] + GtR[2,3,2,1,n]) ) 

        penalty[3,1,n] = ( nrm0[1,n]*nrm0[1,n] * (GtL[3,1,1,1,n] + GtR[3,1,1,1,n])
                         + nrm0[1,n]*nrm0[2,n] * (GtL[3,1,1,2,n] + GtR[3,1,1,2,n])
                         + nrm0[2,n]*nrm0[1,n] * (GtL[3,1,2,1,n] + GtR[3,1,2,1,n])
                         + nrm0[2,n]*nrm0[2,n] * (GtL[3,1,2,2,n] + GtR[3,1,2,2,n]) ) 

        penalty[3,2,n] = ( nrm0[1,n]*nrm0[2,n] * (GtL[3,2,1,2,n] + GtR[3,2,1,2,n])
                         + nrm0[2,n]*nrm0[1,n] * (GtL[3,2,2,1,n] + GtR[3,2,2,1,n]) ) 

        penalty[3,3,n] = ( nrm0[1,n]*nrm0[1,n] * (GtL[3,3,1,1,n] + GtR[3,3,1,1,n])
                         + nrm0[2,n]*nrm0[2,n] * (GtL[3,3,2,2,n] + GtR[3,3,2,2,n]) )

        penalty[4,1,n] = ( nrm0[1,n]*nrm0[1,n] * (GtL[4,1,1,1,n] + GtR[4,1,1,1,n])
                         + nrm0[1,n]*nrm0[2,n] * (GtL[4,1,1,2,n] + GtR[4,1,1,2,n])
                         + nrm0[2,n]*nrm0[1,n] * (GtL[4,1,2,1,n] + GtR[4,1,2,1,n])
                         + nrm0[2,n]*nrm0[2,n] * (GtL[4,1,2,2,n] + GtR[4,1,2,2,n]) )

        penalty[4,2,n] = ( nrm0[1,n]*nrm0[1,n] * (GtL[4,2,1,1,n] + GtR[4,2,1,1,n])
                         + nrm0[1,n]*nrm0[2,n] * (GtL[4,2,1,2,n] + GtR[4,2,1,2,n])
                         + nrm0[2,n]*nrm0[1,n] * (GtL[4,2,2,1,n] + GtR[4,2,2,1,n])
                         + nrm0[2,n]*nrm0[2,n] * (GtL[4,2,2,2,n] + GtR[4,2,2,2,n]) )

        penalty[4,3,n] = ( nrm0[1,n]*nrm0[1,n] * (GtL[4,3,1,1,n] + GtR[4,3,1,1,n])
                         + nrm0[1,n]*nrm0[2,n] * (GtL[4,3,1,2,n] + GtR[4,3,1,2,n])
                         + nrm0[2,n]*nrm0[1,n] * (GtL[4,3,2,1,n] + GtR[4,3,2,1,n])
                         + nrm0[2,n]*nrm0[2,n] * (GtL[4,3,2,2,n] + GtR[4,3,2,2,n]) ) 

        penalty[4,4,n] = ( nrm0[1,n]*nrm0[1,n] * (GtL[4,4,1,1,n] + GtR[4,4,1,1,n])
                         + nrm0[2,n]*nrm0[2,n] * (GtL[4,4,2,2,n] + GtR[4,4,2,2,n]) ) 

        # in order to enhance stability, we try this factor on penalty terms
        factor = opts["Cip"]
        for iDof = 2 : 4
          for jDof = 1 : 4
            penalty[iDof,jDof,n] *=  0.5*const_tii * area[n] / he * factor
          end
        end
      end
    end

    # Start to compute fluxes. We have 3 terms on interfaces:
    # 1) {Fv}⋅[ϕ]
    # 2) {G^T ∇ϕ}⋅[q] 
    # 3) δ{G}[q]:[ϕ]

    # q jump
    for n = 1 : mesh.numNodesPerFace
      for iDof = 1 : mesh.numDofPerNode
        dq[iDof, n] = q_faceL[iDof, n] - q_faceR[iDof, n]
        # factor 0.5 comes from the average operator {G^T ∇ϕ}
        # dqn[1, iDof, n] = -0.5*dq[iDof, n]*nrm[1,n] 
        # dqn[2, iDof, n] = -0.5*dq[iDof, n]*nrm[2,n]
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

    # finally, everything is ready, let's compute fluxes, or penalties

    # This part computes the contribution of
    # ∫ {G^T∇ϕ}:[q] dΓ = ∫ ∇ϕ⋅F dΓ , 
    # where 
    # [q] = (q+ - q-) ⊗ n = Δq⊗n , 
    # Then we can consider Δq⊗n as ∇q and F as viscous flux.
    for n = 1 : mesh.numNodesPerFace
      for iDof = 1 : mesh.numDofPerNode
        #
        # sum up columns of each row
        #
        vecfluxL[1,iDof,n] = ((GtL[iDof,1,1,1,n] + GtL[iDof,2,1,1,n] + GtL[iDof,3,1,1,n]+ GtL[iDof,4,1,1,n]) * nrm[1,n] 
                           + (GtL[iDof,1,1,2,n] + GtL[iDof,2,1,2,n] + GtL[iDof,3,1,2,n]+ GtL[iDof,4,1,2,n]) * nrm[2,n] )
        vecfluxL[2,iDof,n] = ((GtL[iDof,1,2,1,n] + GtL[iDof,2,2,1,n] + GtL[iDof,3,2,1,n]+ GtL[iDof,4,2,1,n]) * nrm[1,n] 
                           + (GtL[iDof,1,2,2,n] + GtL[iDof,2,2,2,n] + GtL[iDof,3,2,2,n]+ GtL[iDof,4,2,2,n]) * nrm[2,n] )
        vecfluxL[1,iDof,n] *=  dq[iDof,n]
        vecfluxL[2,iDof,n] *=  dq[iDof,n]

        vecfluxR[1,iDof,n] = ((GtR[iDof,1,1,1,n] + GtR[iDof,2,1,1,n] + GtR[iDof,3,1,1,n]+ GtR[iDof,4,1,1,n]) * nrm[1,n] 
                           + (GtR[iDof,1,1,2,n] + GtR[iDof,2,1,2,n] + GtR[iDof,3,1,2,n]+ GtR[iDof,4,1,2,n]) * nrm[2,n] )
        vecfluxR[2,iDof,n] = ((GtR[iDof,1,2,1,n] + GtR[iDof,2,2,1,n] + GtR[iDof,3,2,1,n]+ GtR[iDof,4,2,1,n]) * nrm[1,n] 
                           + (GtR[iDof,1,2,2,n] + GtR[iDof,2,2,2,n] + GtR[iDof,3,2,2,n]+ GtR[iDof,4,2,2,n]) * nrm[2,n] )
        vecfluxR[1,iDof,n] *=  dq[iDof,n]
        vecfluxR[2,iDof,n] *=  dq[iDof,n]
      end
    end

    # δ{G}[q]:n, contributing to  δ{G}[q]:[ϕ]
    if sat_type ==  "SAT-SIPG"
      for n = 1 : mesh.numNodesPerFace
        for iDof = 1 : mesh.numDofPerNode
          flux[iDof, n] +=  sat_scalar[n] * dq[iDof, n]
        end
      end
    else
      for n = 1:mesh.numNodesPerFace
        for iDof = 1 : mesh.numDofPerNode
          for jDof = 1:mesh.numDofPerNode
            flux[iDof, n] +=  penalty[iDof, jDof, n]*dq[jDof, n]
          end
        end
      end
    end
    
    # {Fv}⋅n, contributing to {Fv}⋅[ϕ]
    for n = 1:mesh.numNodesPerFace
      for iDof = 2 : mesh.numDofPerNode
        flux[iDof, n]  -=  Fv_avg[1, iDof, n]*nrm[1,n] + Fv_avg[2, iDof, n]*nrm[2,n]
      end
    end
    # accumulate fluxes
    for n = 1:mesh.numNodesPerFace
      for iDof = 1 : Tdim+2
        # eqn.vecflux_face[1, iDof, n, f] -=  dqn[1, iDof, n]*coef_nondim
        # eqn.vecflux_face[2, iDof, n, f] -=  dqn[2, iDof, n]*coef_nondim
        eqn.vecflux_faceL[1, iDof, n, f] -=  vecfluxL[1, iDof, n]*coef_nondim
        eqn.vecflux_faceL[2, iDof, n, f] -=  vecfluxL[2, iDof, n]*coef_nondim
        eqn.vecflux_faceR[1, iDof, n, f] -=  vecfluxR[1, iDof, n]*coef_nondim
        eqn.vecflux_faceR[2, iDof, n, f] -=  vecfluxR[2, iDof, n]*coef_nondim
        eqn.flux_face[iDof, n, f]  +=   flux[iDof, n]*coef_nondim
      end
    end
  end # end of loop over all interfaces

  return nothing
end # end of function calcViscousFlux_interior


function calcViscousFlux_boundary{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh},
                                                          sbp::AbstractSBP,
                                                          eqn::EulerData{Tsol, Tres, Tdim},
                                                          opts)
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
  dq      = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)    
  dqn      = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)    
  q_bnd = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)    
  penalty = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerFace)
  sat_scalar = zeros(Tsol, mesh.numNodesPerFace)

  Gt = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, Tdim, Tdim, mesh.numNodesPerFace)
  Gt_bnd = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, Tdim, Tdim, mesh.numNodesPerFace)
  Fv_face = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)
  Fv_bnd = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)
  vecflux = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)

  nrm     = Array(Tmsh, Tdim, mesh.numNodesPerFace)
  nrm0 = Array(Tmsh, Tdim, mesh.numNodesPerFace)
  area = Array(Tmsh, mesh.numNodesPerFace)
  area_sum = sview(eqn.area_sum, :)

  # sigma = calcTraceInverseInequalityConst(sbp, sbpface)
  dqdx_elem = Array(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerElement )
  dqdx_face = Array(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace )
  for iBC = 1:mesh.numBC
    indx0 = mesh.bndry_offsets[iBC]
    indx1 = mesh.bndry_offsets[iBC+1] - 1

    # specify boundary value function
    # TODO: Put it into a function 
    bnd_functor::AbstractBoundaryValueType
    key_i = string("BC", iBC, "_name")
    val = opts[key_i]
    if "FreeStreamBC" ==  val
      bnd_functor = Farfield()
      calcGt_functor = calcDiffusionTensor
    elseif "nonslipBC" ==  val
      bnd_functor = AdiabaticWall()
      calcGt_functor = calcDiffusionTensor_adiabaticWall
    elseif "noPenetrationBC" ==  val
      continue
    elseif "zeroPressGradientBC" ==  val
      bnd_functor = Farfield()
      calcGt_functor = calcDiffusionTensor
    else
      error("iBC = ", iBC, ", Only 'FreeStreamBC' and 'nonslipBC' available")
    end

    # @bp

    for f = indx0:indx1
      flux  = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)

      bndry = mesh.bndryfaces[f]
      elem = bndry.element
      face = bndry.face
      perm = sview(sbpface.perm, :, face)

      # Compute geometric info on face
      for n = 1 : mesh.numNodesPerFace
        dxidx = sview(mesh.dxidx_bndry, :, :, n, f)
        nrm_xi = sview(sbpface.normal, :, bndry.face)
        nrm[1,n] = dxidx[1, 1]*nrm_xi[1] + dxidx[2, 1]*nrm_xi[2]
        nrm[2,n] = dxidx[1, 2]*nrm_xi[1] + dxidx[2, 2]*nrm_xi[2]
        area[n] = sqrt(nrm[1,n]*nrm[1,n] + nrm[2,n]*nrm[2,n])
        nrm0[1,n] = nrm[1,n]/area[n]
        nrm0[2,n] = nrm[2,n]/area[n]
      end

      # compute element size `he` and face size
      elem_vol = 0.0
      for n = 1:mesh.numNodesPerElement
        elem_vol +=  sbp.w[n]/mesh.jac[n, elem]
      end

      face_area = 0.0
      for n = 1:mesh.numNodesPerFace
        face_area +=  sbpface.wface[n]*area[n]
      end

      # On boundary, the area is weighted twice, that's why we have 0.5
      area_weight = 0.5*area_sum[elem]/face_area

      he = elem_vol/eqn.area_sum[elem]

      # We need viscous flux and diffusion tensor on interfaces, and there
      # are different ways to compute them. For viscous flux:
      # 1) since we have Q on face nodes, if 𝛻Q is available on interface, then we are done.
      # 2) we can comoute viscous flux on volume nodes and then interpolate to interface node.
      # It's logically simple but computationally expensive.

      # Compute boundary viscous flux, F(q_b, ∇q) = G(q_b)∇q.
      # so we need viscousity tensor G, and derivative of q.
      q_face = sview(eqn.q_bndry, :, :, f)
      bnd_functor(q_face, nrm0, eqn.params, q_bnd)

      # diffusion matrix used in penalty term should be computed from q_face rather than q_bnd
      calcGt_functor(q_bnd, Gt)

      q_elem = sview(eqn.q, :, :, elem)
      calcGradient(mesh, sbp, elem, q_elem, dqdx_elem)

      for d = 1 : Tdim
        q_x_node = slice(dqdx_elem, d, :, :)
        q_x_face = slice(dqdx_face, d, :, :)
        boundaryinterpolate(sbpface, bndry, q_x_node, q_x_face) 
      end

      # Li Wang's approach
      calcFvis(Gt, dqdx_face, Fv_face)
      # Hartman's approach
      # calcFvis(q_bnd, dqdx_face, Fv_face)

      # First compute penalty
      if sat_type ==  "SAT-SIPG"
        Cip = 2*opts["Cip"]*p*p
        for n = 2:mesh.numNodesPerFace
          sat_scalar[n] = Cip * area[n] / he
        end
      else
        for n = 1:mesh.numNodesPerFace
          # for iDof = 1:mesh.numDofPerNode
          # for jDof = 1:mesh.numDofPerNode
          # penalty[iDof,jDof,n] =  nrm0[1,n]*nrm0[1,n] * Gt[iDof,jDof,1,1,n]
          # penalty[iDof,jDof,n] +=  nrm0[1,n]*nrm0[2,n] * Gt[iDof,jDof,1,2,n]
          # penalty[iDof,jDof,n] +=  nrm0[2,n]*nrm0[1,n] * Gt[iDof,jDof,2,1,n]
          # penalty[iDof,jDof,n] +=  nrm0[2,n]*nrm0[2,n] * Gt[iDof,jDof,2,2,n]
          # penalty[iDof,jDof,n] *=  const_tii * area[n] / he
          # end
          # end
          penalty[2,1,n]  = ( nrm0[1,n]*nrm0[1,n] * Gt[2,1,1,1,n]
                            + nrm0[1,n]*nrm0[2,n] * Gt[2,1,1,2,n] 
                            + nrm0[2,n]*nrm0[1,n] * Gt[2,1,2,1,n]
                            + nrm0[2,n]*nrm0[2,n] * Gt[2,1,2,2,n] )

          penalty[2,2,n]  = ( nrm0[1,n]*nrm0[1,n] * Gt[2,2,1,1,n]
                            + nrm0[2,n]*nrm0[2,n] * Gt[2,2,2,2,n] )

          penalty[2,3,n]  = ( nrm0[1,n]*nrm0[2,n] * Gt[2,3,1,2,n]
                            + nrm0[2,n]*nrm0[1,n] * Gt[2,3,2,1,n] )

          penalty[3,1,n]  = ( nrm0[1,n]*nrm0[1,n] * Gt[3,1,1,1,n]
                            + nrm0[1,n]*nrm0[2,n] * Gt[3,1,1,2,n]
                            + nrm0[2,n]*nrm0[1,n] * Gt[3,1,2,1,n]
                            + nrm0[2,n]*nrm0[2,n] * Gt[3,1,2,2,n] )

          penalty[3,2,n]  = ( nrm0[1,n]*nrm0[2,n] * Gt[3,2,1,2,n]
                            + nrm0[2,n]*nrm0[1,n] * Gt[3,2,2,1,n] )

          penalty[3,3,n]  = ( nrm0[1,n]*nrm0[1,n] * Gt[3,3,1,1,n]
                            + nrm0[2,n]*nrm0[2,n] * Gt[3,3,2,2,n] )

          penalty[4,1,n]  = ( nrm0[1,n]*nrm0[1,n] * Gt[4,1,1,1,n] 
                            + nrm0[1,n]*nrm0[2,n] * Gt[4,1,1,2,n]
                            + nrm0[2,n]*nrm0[1,n] * Gt[4,1,2,1,n]
                            + nrm0[2,n]*nrm0[2,n] * Gt[4,1,2,2,n] )

          penalty[4,2,n]  = ( nrm0[1,n]*nrm0[1,n] * Gt[4,2,1,1,n]
                            + nrm0[1,n]*nrm0[2,n] * Gt[4,2,1,2,n]
                            + nrm0[2,n]*nrm0[1,n] * Gt[4,2,2,1,n]
                            + nrm0[2,n]*nrm0[2,n] * Gt[4,2,2,2,n] )

          penalty[4,3,n]  = ( nrm0[1,n]*nrm0[1,n] * Gt[4,3,1,1,n]
                            + nrm0[1,n]*nrm0[2,n] * Gt[4,3,1,2,n]
                            + nrm0[2,n]*nrm0[1,n] * Gt[4,3,2,1,n]
                            + nrm0[2,n]*nrm0[2,n] * Gt[4,3,2,2,n] )

          penalty[4,4,n]  = ( nrm0[1,n]*nrm0[1,n] * Gt[4,4,1,1,n] 
                            + nrm0[2,n]*nrm0[2,n] * Gt[4,4,2,2,n] )

          # in order to enhance the stablity, we try this relaxation factor 
          # on penalty terms
          factor = opts["Cip"]
          for iDof = 2 : 4
            for jDof = 1 : 4
              penalty[iDof,jDof,n] *=  const_tii * area[n] / he * factor
            end
          end
        end 
      end

      # Start to compute fluxes.  We have 3 terms on interfaces:
      # 1) -{Fv}⋅[ϕ]
      # 2) -{G^T ∇ϕ}⋅[q] 
      # 3) +δ{G}[q]:[ϕ]
      for n = 1 : mesh.numNodesPerFace
        for iDof = 1 : mesh.numDofPerNode
          dq[iDof, n] = q_face[iDof, n] - q_bnd[iDof, n]
          # dqn[1, iDof, n] = -dq[iDof, n]*nrm[1, n]
          # dqn[2, iDof, n] = -dq[iDof, n]*nrm[2, n]
        end
      end


      # This part computes the contribution of
      # ∫ {G^T∇ϕ}:[q] dΓ = ∫ ∇ϕ⋅F dΓ , 
      # where 
      # [q] = (q+ - q-) ⊗ n, 
      # G = G(q_b) depends on boudanry value.
      # Then we can consider Δq⊗n as ∇q and F as viscous flux.

      # calcFvis(Gt, dqn, vecflux)
      for n = 1 : mesh.numNodesPerFace
        for iDof = 1 : mesh.numDofPerNode
          vecflux[1,iDof,n] = ( (Gt[iDof,1,1,1,n] + Gt[iDof,2,1,1,n] + Gt[iDof,3,1,1,n]+ Gt[iDof,4,1,1,n])*nrm[1,n] 
                            + (Gt[iDof,1,1,2,n] + Gt[iDof,2,1,2,n] + Gt[iDof,3,1,2,n]+ Gt[iDof,4,1,2,n])*nrm[2,n] )
          vecflux[2,iDof,n] = ( (Gt[iDof,1,2,1,n] + Gt[iDof,2,2,1,n] + Gt[iDof,3,2,1,n]+ Gt[iDof,4,2,1,n])*nrm[1,n] 
                            + (Gt[iDof,1,2,2,n] + Gt[iDof,2,2,2,n] + Gt[iDof,3,2,2,n]+ Gt[iDof,4,2,2,n])*nrm[2,n] )
          vecflux[1,iDof,n] *=  dq[iDof,n]
          vecflux[2,iDof,n] *=  dq[iDof,n]
        end
      end

      for n = 1 : mesh.numNodesPerFace
        for iDof = 1 : mesh.numDofPerNode
          flux[iDof, n] -=  ( Fv_face[1, iDof, n]*nrm[1,n] + Fv_face[2, iDof, n]*nrm[2,n] )
        end
      end
      if sat_type ==  "SAT-SIPG"
        for n = 1 : mesh.numNodesPerFace
          for iDof = 1 : mesh.numDofPerNode
            flux[iDof, n] +=  sat_scalar[n] * dq[iDof, n]
          end
        end
      else
        for n = 1 : mesh.numNodesPerFace
          for iDof = 1 : mesh.numDofPerNode
            for jDof = 1: mesh.numDofPerNode
              flux[iDof, n] +=  penalty[iDof, jDof, n]*dq[jDof, n]
            end
          end
        end
      end

      # accumulate fluxes
      for n = 1:mesh.numNodesPerFace
        for iDof = 1 : Tdim+2
          # eqn.vecflux_bndry[1, iDof, n, f] -=  dqn[1, iDof, n]*coef_nondim
          # eqn.vecflux_bndry[2, iDof, n, f] -=  dqn[2, iDof, n]*coef_nondim
          eqn.vecflux_bndry[1, iDof, n, f] -=  vecflux[1, iDof, n]*coef_nondim
          eqn.vecflux_bndry[2, iDof, n, f] -=  vecflux[2, iDof, n]*coef_nondim
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
                                                          opts)
  # This part computes ∫ ∇ϕ⋅F  dΓ, 
  sbpface = mesh.sbpface
  DxL = Array(Tmsh, mesh.numNodesPerElement, mesh.numNodesPerElement, Tdim)
  DxR = Array(Tmsh, mesh.numNodesPerElement, mesh.numNodesPerElement, Tdim)

  GtL = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, Tdim, Tdim, mesh.numNodesPerFace)
  GtR = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, Tdim, Tdim, mesh.numNodesPerFace)

  R = sview(sbpface.interp[:,:])
  w = sview(sbpface.wface, :)
  res = sview(eqn.res, :,:,:)

  numNodes_elem = mesh.numNodesPerElement    # number of Nodes per elemet
  numNodes_face = mesh.numNodesPerFace       # number of nodes on interfaces
  stencilsize = sbpface.stencilsize        # size of stencil for interpolation

  RDxL = Array(Tmsh, mesh.numNodesPerFace, mesh.numNodesPerElement, Tdim)
  RDxR = Array(Tmsh, mesh.numNodesPerFace, mesh.numNodesPerElement, Tdim)
  FvL  = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerFace, mesh.numNodesPerElement)
  FvR  = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerFace, mesh.numNodesPerElement)
  GtRDxL = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerFace, mesh.numNodesPerElement)
  GtRDxR = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerFace, mesh.numNodesPerElement)
  nrm    = Array(Tmsh, Tdim, mesh.numNodesPerFace)
  dq     =  Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)

  nfaces = length(mesh.interfaces)
  # @bp
  for f = 1 : nfaces
    face = mesh.interfaces[f]

    elemL = face.elementL
    elemR = face.elementR
    faceL = face.faceL
    faceR = face.faceR
    pL = sview(sbpface.perm, :, faceL)
    pR = sview(sbpface.perm, :, faceR)


    # xflux = slice(eqn.vecflux_face, 1,:,:,f)
    # yflux = slice(eqn.vecflux_face, 2,:,:,f)

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
    #
    # the old way
    #
    # for n = 1 : mesh.numNodesPerFace
    # dxidx = sview(mesh.dxidx_face, :, :, n, f)
    # nrm_xi = sview(sbp.facenormal, :, faceL)
    # nrm[1,n] = dxidx[1, 1]*nrm_xi[1] + dxidx[2, 1]*nrm_xi[2]
    # nrm[2,n] = dxidx[1, 2]*nrm_xi[1] + dxidx[2, 2]*nrm_xi[2]
    # end

    # q_faceL = slice(eqn.q_face, :, 1, :, f)
    # q_faceR = slice(eqn.q_face, :, 2, :, f)
    # calcDiffusionTensor(q_faceL, GtL)
    # calcDiffusionTensor(q_faceR, GtR)

    # for i = 1 : length(q_faceL)
    # dq[i] = q_faceL[i] - q_faceR[i]
    # end

    # for nf = 1 : mesh.numNodesPerFace
    # for ne = 1 : mesh.numNodesPerElement
    # for iDof = 1 : mesh.numDofPerNode
    # for jDof = 1 : mesh.numDofPerNode
    # GtRDxL[iDof, jDof, nf, ne]  = ( nrm[1,nf]*GtL[iDof, jDof, 1, 1, nf]  
    # + nrm[2,nf]*GtL[iDof, jDof, 1, 2, nf] ) * RDxL[nf, ne, 1] 
    # GtRDxL[iDof, jDof, nf, ne] +=  ( nrm[1,nf]*GtL[iDof, jDof, 2, 1, nf]  
    # + nrm[2,nf]*GtL[iDof, jDof, 2, 2, nf] ) * RDxL[nf, ne, 2] 

    # GtRDxR[iDof, jDof, nf, ne]  = ( nrm[1,nf]*GtR[iDof, jDof, 1, 1, nf]  
    # + nrm[2,nf]*GtR[iDof, jDof, 1, 2, nf] ) * RDxR[nf, ne, 1] 
    # GtRDxR[iDof, jDof, nf, ne] +=  ( nrm[1,nf]*GtR[iDof, jDof, 2, 1, nf]  
    # + nrm[2,nf]*GtR[iDof, jDof, 2, 2, nf] ) * RDxR[nf, ne, 2] 
    # end
    # end
    # end
    # end

    # for ne = 1 : mesh.numNodesPerElement
    # for nf = 1 : mesh.numNodesPerFace
    # for iDof = 1 : mesh.numDofPerNode
    # for jDof = 1 : mesh.numDofPerNode
    # res[iDof, ne, elemL] +=  GtRDxL[iDof, jDof, nf, ne] * dq[iDof, nf] * w[nf]
    # res[iDof, ne, elemR] +=  GtRDxR[iDof, jDof, nf, ne] * dq[iDof, nf] * w[nf]
    # end
    # end
    # end
    # end

    #
    # the new way
    #
    vecfluxL = sview(eqn.vecflux_faceL,:,:,:,f)
    vecfluxR = sview(eqn.vecflux_faceR,:,:,:,f)
    for i = 1 : numNodes_elem
      for j = 1 : numNodes_face
        for iDof = 2 : mesh.numDofPerNode
          res[iDof, i, elemL] +=  ( RDxL[j, i, 1] * vecfluxL[1, iDof, j] 
                                  + RDxL[j, i, 2] * vecfluxL[2, iDof, j] ) * w[j] 
          res[iDof, i, elemR] +=  ( RDxR[j, i, 1] * vecfluxR[1, iDof, j]
                                  + RDxR[j, i, 2] * vecfluxR[2, iDof, j] ) * w[j] 
        end
      end
    end
  end

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
  # Fv = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerFace, mesh.numNodesPerElement)
  GtRDx = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerFace, mesh.numNodesPerElement)
  nrm     = Array(Tmsh, Tdim, mesh.numNodesPerFace)
  nrm0 = Array(Tmsh, Tdim, mesh.numNodesPerFace)
  area = Array(Tmsh, mesh.numNodesPerFace)

  # loop over all the boundaries
  for bc = 1:mesh.numBC
    indx0 = mesh.bndry_offsets[bc]
    indx1 = mesh.bndry_offsets[bc+1] - 1

    # bnd_functor::AbstractBoundaryValueType
    # key_i = string("BC", bc, "_name")
    # val = opts[key_i]
    # if "FreeStreamBC" ==  val
      # bnd_functor = Farfield()
    # elseif "nonslipBC" ==  val
      # bnd_functor = AdiabaticWall()
    # elseif "noPenetrationBC" ==  val
      # continue
    # elseif "zeroPressGradientBC" ==  val
      # bnd_functor = Farfield()
    # else
      # error("iBC = ", bc, ", Only 'FreeStreamBC' and 'nonslipBC' available")
    # end

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

      #
      # old way
      #
      # for n = 1 : mesh.numNodesPerFace
      # dxidx = sview(mesh.dxidx_bndry, :, :, n, f)
      # nrm_xi = sview(sbp.facenormal, :, bndry.face)
      # nrm[1,n] = dxidx[1, 1]*nrm_xi[1] + dxidx[2, 1]*nrm_xi[2]
      # nrm[2,n] = dxidx[1, 2]*nrm_xi[1] + dxidx[2, 2]*nrm_xi[2]
      # area[n] = sqrt(nrm[1,n]*nrm[1,n] + nrm[2,n]*nrm[2,n])
      # nrm0[1,n] = nrm[1,n]/area[n]
      # nrm0[2,n] = nrm[2,n]/area[n]
      # end

      # q_face = sview(eqn.q_bndry, :, :, f)
      # bnd_functor(q_face, nrm0, eqn.params, q_bnd)

      # calcDiffusionTensor(q_bnd, Gt)

      # for i = 1 : length(q_bnd)
      # dq[i] = q_face[i] - q_bnd[i]
      # end

      # for nf = 1 : mesh.numNodesPerFace
      # for ne = 1 : mesh.numNodesPerElement
      # for iDof = 1 : mesh.numDofPerNode
      # for jDof = 1 : mesh.numDofPerNode
      # GtRDx[iDof, jDof, nf, ne]  = ( nrm[1,nf] * Gt[iDof, jDof, 1, 1, nf]  
      # + nrm[2,nf] * Gt[iDof, jDof, 1, 2, nf] ) * RDx[nf, ne, 1] 
      # GtRDx[iDof, jDof, nf, ne] +=  ( nrm[1,nf] * Gt[iDof, jDof, 2, 1, nf]  
      # + nrm[2,nf] * Gt[iDof, jDof, 2, 2, nf] ) * RDx[nf, ne, 2] 
      # end
      # end
      # end
      # end

      # for ne = 1 : mesh.numNodesPerElement
      # for nf = 1 : mesh.numNodesPerFace
      # for iDof = 1 : mesh.numDofPerNode
      # for jDof = 1 : mesh.numDofPerNode
      # res[iDof, ne, elem] +=  GtRDx[iDof, jDof, nf, ne] * dq[iDof, nf] * w[nf] 
      # end
      # end
      # end
      # end

      vecflux = sview(eqn.vecflux_bndry, :,:,:,f)
      for i = 1 : numNodes_elem
        for j = 1 : numNodes_face
          for iDof = 2 : mesh.numDofPerNode
            res[iDof, i, elem] +=  ( RDx[j, i, 1] * vecflux[1, iDof, j] 
                                   + RDx[j, i, 2] * vecflux[2, iDof, j] ) * w[j]
          end
        end
      end
    end
  end

  return nothing
end  # end evalBoundaryIntegrals_vector



@doc """

Integrate ∫ ∇ϕ⋅F dΩ
This function needs to be combined together with `weakdifferentiate`

Input:
mesh    :: 
sbp    ::
eqn    ::
res    ::
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

    calcFvis_elem(sbp, q, dxidx, jac, Fv)

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
