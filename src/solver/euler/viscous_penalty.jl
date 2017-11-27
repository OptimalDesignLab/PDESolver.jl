"""
  compute the interior penalty matrix
  **Input**
   * mesh
   * sbp
   * eqn
   * opts
   * iface: the index of interface
  **Input/Output**
   * pMat: penalty matrix
"""
function cmptIPMat{Tmsh, Tdim, Tsol, Tres}(mesh::AbstractMesh{Tmsh},
                                           sbp::AbstractSBP,
                                           eqn::EulerData{Tsol, Tres, Tdim},
                                           opts,
                                           iface::Int,
                                           GtL::AbstractArray{Tsol, 5},
                                           GtR::AbstractArray{Tsol, 5},
                                           pMat::AbstractArray{Tsol, 3})
  if opts["SAT_type"] == "Hartman"
    cmptIPMat_hartman(mesh, sbp, eqn, opts, iface, GtL, GtR, pMat)
  elseif opts["SAT_type"] == "SAT-SIPG"
    cmptIPMat_SIPG(mesh, sbp, eqn, opts, iface, GtL, GtR, pMat)
  elseif opts["SAT_type"] == "SAT-BR2"
    cmptIPMat_BR2(mesh, sbp, eqn, opts, iface, GtL, GtR, pMat)
  end
end

function cmptIPMat_SIPG{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh},
                                          sbp::AbstractSBP,
                                          eqn::EulerData{Tsol, Tres, 2},
                                          opts,
                                          iface::Int,
                                          GtL::AbstractArray{Tsol, 5},
                                          GtR::AbstractArray{Tsol, 5},
                                          pMat::AbstractArray{Tsol, 3})
  error("SAT-SIPG not available yet")
end

function cmptIPMat_BR2{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh},
                                         sbp::AbstractSBP,
                                         eqn::EulerData{Tsol, Tres, 2},
                                         opts,
                                         iface::Int,
                                         GtL::AbstractArray{Tsol, 5},
                                         GtR::AbstractArray{Tsol, 5},
                                         pMat::AbstractArray{Tsol, 3})
  error("SAT-BR2 not available yet")
end

function cmptIPMat_hartman{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh},
                                             sbp::AbstractSBP,
                                             eqn::EulerData{Tsol, Tres, 2},
                                             opts,
                                             iface::Int,
                                             GtL::AbstractArray{Tsol, 5},
                                             GtR::AbstractArray{Tsol, 5},
                                             pMat::AbstractArray{Tsol, 3})
# Reference:
#    Shahbazi, Mavriplis, Multigrid algorithms for high order discontinuous 
#    Galerkin discretization of the compressible Navier-Stokes equations, 
#    JCP, volume 228 Issue 21, November, 2009
#
  # in order to enhance/test stability, we try this factor on penalty terms
  factor = eqn.params.penalty_relaxation
  const_tii = eqn.params.const_tii
  Tdim = 2
  sbpface = mesh.sbpface

  face = mesh.interfaces[iface]
  elemL = face.elementL
  elemR = face.elementR

  # Compute geometric info on face
  nrm1 = zeros(Tmsh, Tdim, mesh.numNodesPerFace)
  area = zeros(Tmsh, mesh.numNodesPerFace)
  for n = 1 : mesh.numNodesPerFace
    nrm_xy = ro_sview(mesh.nrm_face, :, n, iface)
    area[n] = norm(nrm_xy)

    for i = 1 : Tdim
      nrm1[i,n] = nrm_xy[i] / area[n]
    end
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

  for n = 1 : mesh.numNodesPerFace
    # for iDof = 1:mesh.numDofPerNode
      # for jDof = 1:mesh.numDofPerNode
        # pMat[iDof, jDof, n] = 0.0
        # for jDim = 1 : Tdim
          # for iDim = 1 : Tdim
            # pMat[iDof,jDof,n] +=  nrm1[iDim,n]*nrm1[jDim,n] * (GtL[iDof,jDof,iDim,jDim,n] + GtR[iDof,jDof,iDim,jDim,n])
          # end
        # end
        # pMat[iDof,jDof,n] *=  0.5*const_tii * area[n] / he * factor
      # end

    # end
    coef = 0.5*const_tii * area[n] / he * factor

    pMat[2,1,n] = ( nrm1[1,n]*nrm1[1,n] * (GtL[2,1,1,1,n] + GtR[2,1,1,1,n])
                  + nrm1[1,n]*nrm1[2,n] * (GtL[2,1,1,2,n] + GtR[2,1,1,2,n])
                  + nrm1[2,n]*nrm1[1,n] * (GtL[2,1,2,1,n] + GtR[2,1,2,1,n])
                  + nrm1[2,n]*nrm1[2,n] * (GtL[2,1,2,2,n] + GtR[2,1,2,2,n])
                 ) * coef 

    pMat[2,2,n] = ( nrm1[1,n]*nrm1[1,n] * (GtL[2,2,1,1,n] + GtR[2,2,1,1,n])
                  + nrm1[2,n]*nrm1[2,n] * (GtL[2,2,2,2,n] + GtR[2,2,2,2,n])
                 ) * coef 

    pMat[2,3,n] = ( nrm1[1,n]*nrm1[2,n] * (GtL[2,3,1,2,n] + GtR[2,3,1,2,n])
                  + nrm1[2,n]*nrm1[1,n] * (GtL[2,3,2,1,n] + GtR[2,3,2,1,n]) 
                 ) * coef 

    pMat[3,1,n] = ( nrm1[1,n]*nrm1[1,n] * (GtL[3,1,1,1,n] + GtR[3,1,1,1,n])
                  + nrm1[1,n]*nrm1[2,n] * (GtL[3,1,1,2,n] + GtR[3,1,1,2,n])
                  + nrm1[2,n]*nrm1[1,n] * (GtL[3,1,2,1,n] + GtR[3,1,2,1,n])
                  + nrm1[2,n]*nrm1[2,n] * (GtL[3,1,2,2,n] + GtR[3,1,2,2,n]) 
                 ) * coef 

    pMat[3,2,n] = ( nrm1[1,n]*nrm1[2,n] * (GtL[3,2,1,2,n] + GtR[3,2,1,2,n])
                  + nrm1[2,n]*nrm1[1,n] * (GtL[3,2,2,1,n] + GtR[3,2,2,1,n]) 
                 ) * coef 

    pMat[3,3,n] = ( nrm1[1,n]*nrm1[1,n] * (GtL[3,3,1,1,n] + GtR[3,3,1,1,n])
                  + nrm1[2,n]*nrm1[2,n] * (GtL[3,3,2,2,n] + GtR[3,3,2,2,n]) 
                 ) * coef 

    pMat[4,1,n] = ( nrm1[1,n]*nrm1[1,n] * (GtL[4,1,1,1,n] + GtR[4,1,1,1,n])
                  + nrm1[1,n]*nrm1[2,n] * (GtL[4,1,1,2,n] + GtR[4,1,1,2,n])
                  + nrm1[2,n]*nrm1[1,n] * (GtL[4,1,2,1,n] + GtR[4,1,2,1,n])
                  + nrm1[2,n]*nrm1[2,n] * (GtL[4,1,2,2,n] + GtR[4,1,2,2,n]) 
                 ) * coef 

    pMat[4,2,n] = ( nrm1[1,n]*nrm1[1,n] * (GtL[4,2,1,1,n] + GtR[4,2,1,1,n])
                  + nrm1[1,n]*nrm1[2,n] * (GtL[4,2,1,2,n] + GtR[4,2,1,2,n])
                  + nrm1[2,n]*nrm1[1,n] * (GtL[4,2,2,1,n] + GtR[4,2,2,1,n])
                  + nrm1[2,n]*nrm1[2,n] * (GtL[4,2,2,2,n] + GtR[4,2,2,2,n]) 
                 ) * coef 

    pMat[4,3,n] = ( nrm1[1,n]*nrm1[1,n] * (GtL[4,3,1,1,n] + GtR[4,3,1,1,n])
                  + nrm1[1,n]*nrm1[2,n] * (GtL[4,3,1,2,n] + GtR[4,3,1,2,n])
                  + nrm1[2,n]*nrm1[1,n] * (GtL[4,3,2,1,n] + GtR[4,3,2,1,n])
                  + nrm1[2,n]*nrm1[2,n] * (GtL[4,3,2,2,n] + GtR[4,3,2,2,n]) 
                 ) * coef 

    pMat[4,4,n] = ( nrm1[1,n]*nrm1[1,n] * (GtL[4,4,1,1,n] + GtR[4,4,1,1,n])
                  + nrm1[2,n]*nrm1[2,n] * (GtL[4,4,2,2,n] + GtR[4,4,2,2,n]) 
                 ) * coef 

  end # end of loop over nodes
  return nothing
end

function cmptIPMat_hartman{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh},
                                             sbp::AbstractSBP,
                                             eqn::EulerData{Tsol, Tres, 3},
                                             opts,
                                             iface::Int,
                                             GtL::AbstractArray{Tsol, 5},
                                             GtR::AbstractArray{Tsol, 5},
                                             pMat::AbstractArray{Tsol, 3})
  # in order to enhance/test stability, we try this factor on penalty terms
  factor = eqn.params.penalty_relaxation
  const_tii = eqn.params.const_tii
  Tdim = 3
  sbpface = mesh.sbpface
  face = mesh.interfaces[iface]
  elemL = face.elementL
  elemR = face.elementR

  # Compute geometric info on face
  nrm1 = Array(Tmsh, Tdim, mesh.numNodesPerFace)
  area = Array(Tmsh, mesh.numNodesPerFace)
  for n = 1 : mesh.numNodesPerFace
    nrm_xy = ro_sview(mesh.nrm_face, :, n, iface)
    area[n] = norm(nrm_xy)

    for i = 1 : Tdim
      nrm1[i,n] = nrm_xy[i] / area[n]
    end
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

  for n = 1 : mesh.numNodesPerFace
    # for iDof = 1:mesh.numDofPerNode
      # for jDof = 1:mesh.numDofPerNode
        # pMat[iDof, jDof, n] = 0.0
        # for jDim = 1 : Tdim
          # for iDim = 1 : Tdim
            # pMat[iDof,jDof,n] +=  nrm1[iDim,n]*nrm1[jDim,n] * (GtL[iDof,jDof,iDim,jDim,n] + GtR[iDof,jDof,iDim,jDim,n])
          # end
        # end
        # pMat[iDof,jDof,n] *=  0.5*const_tii * area[n] / he * factor
      # end

    # end
    coef = 0.5*const_tii * area[n] / he * factor
                 

    pMat[2,1,n] = ( nrm1[1,n]*nrm1[1,n] * (GtL[2,1,1,1,n] + GtR[2,1,1,1,n])
                  + nrm1[1,n]*nrm1[2,n] * (GtL[2,1,1,2,n] + GtR[2,1,1,2,n])
                  + nrm1[1,n]*nrm1[3,n] * (GtL[2,1,1,3,n] + GtR[2,1,1,3,n])
                  + nrm1[2,n]*nrm1[1,n] * (GtL[2,1,2,1,n] + GtR[2,1,2,1,n])
                  + nrm1[2,n]*nrm1[2,n] * (GtL[2,1,2,2,n] + GtR[2,1,2,2,n]) 
                  # + nrm1[2,n]*nrm1[3,n] * (GtL[2,1,2,3,n] + GtR[2,1,2,3,n]) 
                  + nrm1[3,n]*nrm1[1,n] * (GtL[2,1,3,1,n] + GtR[2,1,3,1,n]) 
                  # + nrm1[3,n]*nrm1[2,n] * (GtL[2,1,3,2,n] + GtR[2,1,3,2,n]) 
                  + nrm1[3,n]*nrm1[3,n] * (GtL[2,1,3,3,n] + GtR[2,1,3,3,n]) 
                 ) * coef

    
    pMat[2,2,n] = ( nrm1[1,n]*nrm1[1,n] * (GtL[2,2,1,1,n] + GtR[2,2,1,1,n])
                  # + nrm1[1,n]*nrm1[2,n] * (GtL[2,2,1,2,n] + GtR[2,2,1,2,n])
                  # + nrm1[1,n]*nrm1[3,n] * (GtL[2,2,1,3,n] + GtR[2,2,1,3,n])
                  # + nrm1[2,n]*nrm1[1,n] * (GtL[2,2,2,1,n] + GtR[2,2,2,1,n])
                  + nrm1[2,n]*nrm1[2,n] * (GtL[2,2,2,2,n] + GtR[2,2,2,2,n]) 
                  # + nrm1[2,n]*nrm1[3,n] * (GtL[2,2,2,3,n] + GtR[2,2,2,3,n]) 
                  # + nrm1[3,n]*nrm1[1,n] * (GtL[2,2,3,1,n] + GtR[2,2,3,1,n]) 
                  # + nrm1[3,n]*nrm1[2,n] * (GtL[2,2,3,2,n] + GtR[2,2,3,2,n]) 
                  + nrm1[3,n]*nrm1[3,n] * (GtL[2,2,3,3,n] + GtR[2,2,3,3,n]) 
                 ) * coef

    pMat[2,3,n] = ( 
                   # nrm1[1,n]*nrm1[1,n] * (GtL[2,3,1,1,n] + GtR[2,3,1,1,n])
                  + nrm1[1,n]*nrm1[2,n] * (GtL[2,3,1,2,n] + GtR[2,3,1,2,n])
                  # + nrm1[1,n]*nrm1[3,n] * (GtL[2,3,1,3,n] + GtR[2,3,1,3,n])
                  + nrm1[2,n]*nrm1[1,n] * (GtL[2,3,2,1,n] + GtR[2,3,2,1,n])
                  # + nrm1[2,n]*nrm1[2,n] * (GtL[2,3,2,2,n] + GtR[2,3,2,2,n]) 
                  # + nrm1[2,n]*nrm1[3,n] * (GtL[2,3,2,3,n] + GtR[2,3,2,3,n]) 
                  # + nrm1[3,n]*nrm1[1,n] * (GtL[2,3,3,1,n] + GtR[2,3,3,1,n]) 
                  # + nrm1[3,n]*nrm1[2,n] * (GtL[2,3,3,2,n] + GtR[2,3,3,2,n]) 
                  # + nrm1[3,n]*nrm1[3,n] * (GtL[2,3,3,3,n] + GtR[2,3,3,3,n]) 
                 ) * coef

    pMat[2,4,n] = ( 
                    # nrm1[1,n]*nrm1[1,n] * (GtL[2,4,1,1,n] + GtR[2,4,1,1,n])
                  # + nrm1[1,n]*nrm1[2,n] * (GtL[2,4,1,2,n] + GtR[2,4,1,2,n])
                  + nrm1[1,n]*nrm1[3,n] * (GtL[2,4,1,3,n] + GtR[2,4,1,3,n])
                  # + nrm1[2,n]*nrm1[1,n] * (GtL[2,4,2,1,n] + GtR[2,4,2,1,n])
                  # + nrm1[2,n]*nrm1[2,n] * (GtL[2,4,2,2,n] + GtR[2,4,2,2,n]) 
                  # + nrm1[2,n]*nrm1[3,n] * (GtL[2,4,2,3,n] + GtR[2,4,2,3,n]) 
                  + nrm1[3,n]*nrm1[1,n] * (GtL[2,4,3,1,n] + GtR[2,4,3,1,n]) 
                  # + nrm1[3,n]*nrm1[2,n] * (GtL[2,4,3,2,n] + GtR[2,4,3,2,n]) 
                  # + nrm1[3,n]*nrm1[3,n] * (GtL[2,4,3,3,n] + GtR[2,4,3,3,n]) 
                 ) * coef

    # pMat[2,5,n] = ( 
                    # # nrm1[1,n]*nrm1[1,n] * (GtL[2,5,1,1,n] + GtR[2,5,1,1,n])
                  # # + nrm1[1,n]*nrm1[2,n] * (GtL[2,5,1,2,n] + GtR[2,5,1,2,n])
                  # # + nrm1[1,n]*nrm1[3,n] * (GtL[2,5,1,3,n] + GtR[2,5,1,3,n])
                  # # + nrm1[2,n]*nrm1[1,n] * (GtL[2,5,2,1,n] + GtR[2,5,2,1,n])
                  # # + nrm1[2,n]*nrm1[2,n] * (GtL[2,5,2,2,n] + GtR[2,5,2,2,n]) 
                  # # + nrm1[2,n]*nrm1[3,n] * (GtL[2,5,2,3,n] + GtR[2,5,2,3,n]) 
                  # # + nrm1[3,n]*nrm1[1,n] * (GtL[2,5,3,1,n] + GtR[2,5,3,1,n]) 
                  # # + nrm1[3,n]*nrm1[2,n] * (GtL[2,5,3,2,n] + GtR[2,5,3,2,n]) 
                  # # + nrm1[3,n]*nrm1[3,n] * (GtL[2,5,3,3,n] + GtR[2,5,3,3,n]) )

    pMat[3,1,n] = ( nrm1[1,n]*nrm1[1,n] * (GtL[3,1,1,1,n] + GtR[3,1,1,1,n])
                  + nrm1[1,n]*nrm1[2,n] * (GtL[3,1,1,2,n] + GtR[3,1,1,2,n])
                  # + nrm1[1,n]*nrm1[3,n] * (GtL[3,1,1,3,n] + GtR[3,1,1,3,n])
                  + nrm1[2,n]*nrm1[1,n] * (GtL[3,1,2,1,n] + GtR[3,1,2,1,n])
                  + nrm1[2,n]*nrm1[2,n] * (GtL[3,1,2,2,n] + GtR[3,1,2,2,n]) 
                  + nrm1[2,n]*nrm1[3,n] * (GtL[3,1,2,3,n] + GtR[3,1,2,3,n]) 
                  # + nrm1[3,n]*nrm1[1,n] * (GtL[3,1,3,1,n] + GtR[3,1,3,1,n]) 
                  + nrm1[3,n]*nrm1[2,n] * (GtL[3,1,3,2,n] + GtR[3,1,3,2,n]) 
                  + nrm1[3,n]*nrm1[3,n] * (GtL[3,1,3,3,n] + GtR[3,1,3,3,n]) 
                 ) * coef

    pMat[3,2,n] = ( 
                    # nrm1[1,n]*nrm1[1,n] * (GtL[3,2,1,1,n] + GtR[3,2,1,1,n])
                  + nrm1[1,n]*nrm1[2,n] * (GtL[3,2,1,2,n] + GtR[3,2,1,2,n])
                  # + nrm1[1,n]*nrm1[3,n] * (GtL[3,2,1,3,n] + GtR[3,2,1,3,n])
                  + nrm1[2,n]*nrm1[1,n] * (GtL[3,2,2,1,n] + GtR[3,2,2,1,n])
                  # + nrm1[2,n]*nrm1[2,n] * (GtL[3,2,2,2,n] + GtR[3,2,2,2,n]) 
                  # + nrm1[2,n]*nrm1[3,n] * (GtL[3,2,2,3,n] + GtR[3,2,2,3,n]) 
                  # + nrm1[3,n]*nrm1[1,n] * (GtL[3,2,3,1,n] + GtR[3,2,3,1,n]) 
                  # + nrm1[3,n]*nrm1[2,n] * (GtL[3,2,3,2,n] + GtR[3,2,3,2,n]) 
                  # + nrm1[3,n]*nrm1[3,n] * (GtL[3,2,3,3,n] + GtR[3,2,3,3,n]) 
                 ) * coef

    pMat[3,3,n] = ( nrm1[1,n]*nrm1[1,n] * (GtL[3,3,1,1,n] + GtR[3,3,1,1,n])
                  # + nrm1[1,n]*nrm1[2,n] * (GtL[3,3,1,2,n] + GtR[3,3,1,2,n])
                  # + nrm1[1,n]*nrm1[3,n] * (GtL[3,3,1,3,n] + GtR[3,3,1,3,n])
                  # + nrm1[2,n]*nrm1[1,n] * (GtL[3,3,2,1,n] + GtR[3,3,2,1,n])
                  + nrm1[2,n]*nrm1[2,n] * (GtL[3,3,2,2,n] + GtR[3,3,2,2,n]) 
                  # + nrm1[2,n]*nrm1[3,n] * (GtL[3,3,2,3,n] + GtR[3,3,2,3,n]) 
                  # + nrm1[3,n]*nrm1[1,n] * (GtL[3,3,3,1,n] + GtR[3,3,3,1,n]) 
                  # + nrm1[3,n]*nrm1[2,n] * (GtL[3,3,3,2,n] + GtR[3,3,3,2,n]) 
                  + nrm1[3,n]*nrm1[3,n] * (GtL[3,3,3,3,n] + GtR[3,3,3,3,n]) 
                 ) * coef

    pMat[3,4,n] = ( 
                    # nrm1[1,n]*nrm1[1,n] * (GtL[3,4,1,1,n] + GtR[3,4,1,1,n])
                  # + nrm1[1,n]*nrm1[2,n] * (GtL[3,4,1,2,n] + GtR[3,4,1,2,n])
                  # + nrm1[1,n]*nrm1[3,n] * (GtL[3,4,1,3,n] + GtR[3,4,1,3,n])
                  # + nrm1[2,n]*nrm1[1,n] * (GtL[3,4,2,1,n] + GtR[3,4,2,1,n])
                  # + nrm1[2,n]*nrm1[2,n] * (GtL[3,4,2,2,n] + GtR[3,4,2,2,n]) 
                  + nrm1[2,n]*nrm1[3,n] * (GtL[3,4,2,3,n] + GtR[3,4,2,3,n]) 
                  # + nrm1[3,n]*nrm1[1,n] * (GtL[3,4,3,1,n] + GtR[3,4,3,1,n]) 
                  + nrm1[3,n]*nrm1[2,n] * (GtL[3,4,3,2,n] + GtR[3,4,3,2,n]) 
                  # + nrm1[3,n]*nrm1[3,n] * (GtL[3,4,3,3,n] + GtR[3,4,3,3,n]) 
                 ) * coef

    # pMat[3,5,n] = ( 
                    # # nrm1[1,n]*nrm1[1,n] * (GtL[3,5,1,1,n] + GtR[3,5,1,1,n])
                  # # + nrm1[1,n]*nrm1[2,n] * (GtL[3,5,1,2,n] + GtR[3,5,1,2,n])
                  # # + nrm1[1,n]*nrm1[3,n] * (GtL[3,5,1,3,n] + GtR[3,5,1,3,n])
                  # # + nrm1[2,n]*nrm1[1,n] * (GtL[3,5,2,1,n] + GtR[3,5,2,1,n])
                  # # + nrm1[2,n]*nrm1[2,n] * (GtL[3,5,2,2,n] + GtR[3,5,2,2,n]) 
                  # # + nrm1[2,n]*nrm1[3,n] * (GtL[3,5,2,3,n] + GtR[3,5,2,3,n]) 
                  # # + nrm1[3,n]*nrm1[1,n] * (GtL[3,5,3,1,n] + GtR[3,5,3,1,n]) 
                  # # + nrm1[3,n]*nrm1[2,n] * (GtL[3,5,3,2,n] + GtR[3,5,3,2,n]) 
                  # # + nrm1[3,n]*nrm1[3,n] * (GtL[3,5,3,3,n] + GtR[3,5,3,3,n]) 
                 # )

    pMat[4,1,n] = ( nrm1[1,n]*nrm1[1,n] * (GtL[4,1,1,1,n] + GtR[4,1,1,1,n])
                  # + nrm1[1,n]*nrm1[2,n] * (GtL[4,1,1,2,n] + GtR[4,1,1,2,n])
                  + nrm1[1,n]*nrm1[3,n] * (GtL[4,1,1,3,n] + GtR[4,1,1,3,n])
                  # + nrm1[2,n]*nrm1[1,n] * (GtL[4,1,2,1,n] + GtR[4,1,2,1,n])
                  + nrm1[2,n]*nrm1[2,n] * (GtL[4,1,2,2,n] + GtR[4,1,2,2,n]) 
                  + nrm1[2,n]*nrm1[3,n] * (GtL[4,1,2,3,n] + GtR[4,1,2,3,n]) 
                  + nrm1[3,n]*nrm1[1,n] * (GtL[4,1,3,1,n] + GtR[4,1,3,1,n]) 
                  + nrm1[3,n]*nrm1[2,n] * (GtL[4,1,3,2,n] + GtR[4,1,3,2,n]) 
                  + nrm1[3,n]*nrm1[3,n] * (GtL[4,1,3,3,n] + GtR[4,1,3,3,n]) 
                 ) * coef

    pMat[4,2,n] = ( 
                    # nrm1[1,n]*nrm1[1,n] * (GtL[4,2,1,1,n] + GtR[4,2,1,1,n])
                  # + nrm1[1,n]*nrm1[2,n] * (GtL[4,2,1,2,n] + GtR[4,2,1,2,n])
                  + nrm1[1,n]*nrm1[3,n] * (GtL[4,2,1,3,n] + GtR[4,2,1,3,n])
                  # + nrm1[2,n]*nrm1[1,n] * (GtL[4,2,2,1,n] + GtR[4,2,2,1,n])
                  # + nrm1[2,n]*nrm1[2,n] * (GtL[4,2,2,2,n] + GtR[4,2,2,2,n]) 
                  # + nrm1[2,n]*nrm1[3,n] * (GtL[4,2,2,3,n] + GtR[4,2,2,3,n]) 
                  + nrm1[3,n]*nrm1[1,n] * (GtL[4,2,3,1,n] + GtR[4,2,3,1,n]) 
                  # + nrm1[3,n]*nrm1[2,n] * (GtL[4,2,3,2,n] + GtR[4,2,3,2,n]) 
                  # + nrm1[3,n]*nrm1[3,n] * (GtL[4,2,3,3,n] + GtR[4,2,3,3,n]) 
                 ) * coef

    pMat[4,3,n] = ( 
                    # nrm1[1,n]*nrm1[1,n] * (GtL[4,3,1,1,n] + GtR[4,3,1,1,n])
                  # + nrm1[1,n]*nrm1[2,n] * (GtL[4,3,1,2,n] + GtR[4,3,1,2,n])
                  # + nrm1[1,n]*nrm1[3,n] * (GtL[4,3,1,3,n] + GtR[4,3,1,3,n])
                  # + nrm1[2,n]*nrm1[1,n] * (GtL[4,3,2,1,n] + GtR[4,3,2,1,n])
                  # + nrm1[2,n]*nrm1[2,n] * (GtL[4,3,2,2,n] + GtR[4,3,2,2,n]) 
                  + nrm1[2,n]*nrm1[3,n] * (GtL[4,3,2,3,n] + GtR[4,3,2,3,n]) 
                  # + nrm1[3,n]*nrm1[1,n] * (GtL[4,3,3,1,n] + GtR[4,3,3,1,n]) 
                  + nrm1[3,n]*nrm1[2,n] * (GtL[4,3,3,2,n] + GtR[4,3,3,2,n]) 
                  # + nrm1[3,n]*nrm1[3,n] * (GtL[4,3,3,3,n] + GtR[4,3,3,3,n]) 
                 ) * coef

    pMat[4,4,n] = ( nrm1[1,n]*nrm1[1,n] * (GtL[4,4,1,1,n] + GtR[4,4,1,1,n])
                  # + nrm1[1,n]*nrm1[2,n] * (GtL[4,4,1,2,n] + GtR[4,4,1,2,n])
                  # + nrm1[1,n]*nrm1[3,n] * (GtL[4,4,1,3,n] + GtR[4,4,1,3,n])
                  # + nrm1[2,n]*nrm1[1,n] * (GtL[4,4,2,1,n] + GtR[4,4,2,1,n])
                  + nrm1[2,n]*nrm1[2,n] * (GtL[4,4,2,2,n] + GtR[4,4,2,2,n]) 
                  # + nrm1[2,n]*nrm1[3,n] * (GtL[4,4,2,3,n] + GtR[4,4,2,3,n]) 
                  # + nrm1[3,n]*nrm1[1,n] * (GtL[4,4,3,1,n] + GtR[4,4,3,1,n]) 
                  # + nrm1[3,n]*nrm1[2,n] * (GtL[4,4,3,2,n] + GtR[4,4,3,2,n]) 
                  + nrm1[3,n]*nrm1[3,n] * (GtL[4,4,3,3,n] + GtR[4,4,3,3,n]) 
                 ) * coef

    # pMat[4,5,n] = ( 
                    # # nrm1[1,n]*nrm1[1,n] * (GtL[4,5,1,1,n] + GtR[4,5,1,1,n])
                  # # + nrm1[1,n]*nrm1[2,n] * (GtL[4,5,1,2,n] + GtR[4,5,1,2,n])
                  # # + nrm1[1,n]*nrm1[3,n] * (GtL[4,5,1,3,n] + GtR[4,5,1,3,n])
                  # # + nrm1[2,n]*nrm1[1,n] * (GtL[4,5,2,1,n] + GtR[4,5,2,1,n])
                  # # + nrm1[2,n]*nrm1[2,n] * (GtL[4,5,2,2,n] + GtR[4,5,2,2,n]) 
                  # # + nrm1[2,n]*nrm1[3,n] * (GtL[4,5,2,3,n] + GtR[4,5,2,3,n]) 
                  # # + nrm1[3,n]*nrm1[1,n] * (GtL[4,5,3,1,n] + GtR[4,5,3,1,n]) 
                  # # + nrm1[3,n]*nrm1[2,n] * (GtL[4,5,3,2,n] + GtR[4,5,3,2,n]) 
                  # # + nrm1[3,n]*nrm1[3,n] * (GtL[4,5,3,3,n] + GtR[4,5,3,3,n]) 
                 # )

    pMat[5,1,n] = ( nrm1[1,n]*nrm1[1,n] * (GtL[5,1,1,1,n] + GtR[5,1,1,1,n])
                  + nrm1[1,n]*nrm1[2,n] * (GtL[5,1,1,2,n] + GtR[5,1,1,2,n])
                  + nrm1[1,n]*nrm1[3,n] * (GtL[5,1,1,3,n] + GtR[5,1,1,3,n])
                  + nrm1[2,n]*nrm1[1,n] * (GtL[5,1,2,1,n] + GtR[5,1,2,1,n])
                  + nrm1[2,n]*nrm1[2,n] * (GtL[5,1,2,2,n] + GtR[5,1,2,2,n]) 
                  + nrm1[2,n]*nrm1[3,n] * (GtL[5,1,2,3,n] + GtR[5,1,2,3,n]) 
                  + nrm1[3,n]*nrm1[1,n] * (GtL[5,1,3,1,n] + GtR[5,1,3,1,n]) 
                  + nrm1[3,n]*nrm1[2,n] * (GtL[5,1,3,2,n] + GtR[5,1,3,2,n]) 
                  + nrm1[3,n]*nrm1[3,n] * (GtL[5,1,3,3,n] + GtR[5,1,3,3,n]) 
                 ) * coef


    pMat[5,2,n] = ( nrm1[1,n]*nrm1[1,n] * (GtL[5,2,1,1,n] + GtR[5,2,1,1,n])
                  + nrm1[1,n]*nrm1[2,n] * (GtL[5,2,1,2,n] + GtR[5,2,1,2,n])
                  + nrm1[1,n]*nrm1[3,n] * (GtL[5,2,1,3,n] + GtR[5,2,1,3,n])
                  + nrm1[2,n]*nrm1[1,n] * (GtL[5,2,2,1,n] + GtR[5,2,2,1,n])
                  + nrm1[2,n]*nrm1[2,n] * (GtL[5,2,2,2,n] + GtR[5,2,2,2,n]) 
                  # + nrm1[2,n]*nrm1[3,n] * (GtL[5,2,2,3,n] + GtR[5,2,2,3,n]) 
                  + nrm1[3,n]*nrm1[1,n] * (GtL[5,2,3,1,n] + GtR[5,2,3,1,n]) 
                  # + nrm1[3,n]*nrm1[2,n] * (GtL[5,2,3,2,n] + GtR[5,2,3,2,n]) 
                  + nrm1[3,n]*nrm1[3,n] * (GtL[5,2,3,3,n] + GtR[5,2,3,3,n]) 
                 ) * coef

    pMat[5,3,n] = ( nrm1[1,n]*nrm1[1,n] * (GtL[5,3,1,1,n] + GtR[5,3,1,1,n])
                  + nrm1[1,n]*nrm1[2,n] * (GtL[5,3,1,2,n] + GtR[5,3,1,2,n])
                  # + nrm1[1,n]*nrm1[3,n] * (GtL[5,3,1,3,n] + GtR[5,3,1,3,n])
                  + nrm1[2,n]*nrm1[1,n] * (GtL[5,3,2,1,n] + GtR[5,3,2,1,n])
                  + nrm1[2,n]*nrm1[2,n] * (GtL[5,3,2,2,n] + GtR[5,3,2,2,n]) 
                  + nrm1[2,n]*nrm1[3,n] * (GtL[5,3,2,3,n] + GtR[5,3,2,3,n]) 
                  # + nrm1[3,n]*nrm1[1,n] * (GtL[5,3,3,1,n] + GtR[5,3,3,1,n]) 
                  + nrm1[3,n]*nrm1[2,n] * (GtL[5,3,3,2,n] + GtR[5,3,3,2,n]) 
                  + nrm1[3,n]*nrm1[3,n] * (GtL[5,3,3,3,n] + GtR[5,3,3,3,n]) 
                 ) * coef

    pMat[5,4,n] = ( nrm1[1,n]*nrm1[1,n] * (GtL[5,4,1,1,n] + GtR[5,4,1,1,n])
                  # + nrm1[1,n]*nrm1[2,n] * (GtL[5,4,1,2,n] + GtR[5,4,1,2,n])
                  + nrm1[1,n]*nrm1[3,n] * (GtL[5,4,1,3,n] + GtR[5,4,1,3,n])
                  # + nrm1[2,n]*nrm1[1,n] * (GtL[5,4,2,1,n] + GtR[5,4,2,1,n])
                  + nrm1[2,n]*nrm1[2,n] * (GtL[5,4,2,2,n] + GtR[5,4,2,2,n]) 
                  + nrm1[2,n]*nrm1[3,n] * (GtL[5,4,2,3,n] + GtR[5,4,2,3,n]) 
                  + nrm1[3,n]*nrm1[1,n] * (GtL[5,4,3,1,n] + GtR[5,4,3,1,n]) 
                  + nrm1[3,n]*nrm1[2,n] * (GtL[5,4,3,2,n] + GtR[5,4,3,2,n]) 
                  + nrm1[3,n]*nrm1[3,n] * (GtL[5,4,3,3,n] + GtR[5,4,3,3,n]) 
                 ) * coef

    pMat[5,5,n] = ( nrm1[1,n]*nrm1[1,n] * (GtL[5,5,1,1,n] + GtR[5,5,1,1,n])
                  # + nrm1[1,n]*nrm1[2,n] * (GtL[5,5,1,2,n] + GtR[5,5,1,2,n])
                  # + nrm1[1,n]*nrm1[3,n] * (GtL[5,5,1,3,n] + GtR[5,5,1,3,n])
                  # + nrm1[2,n]*nrm1[1,n] * (GtL[5,5,2,1,n] + GtR[5,5,2,1,n])
                  + nrm1[2,n]*nrm1[2,n] * (GtL[5,5,2,2,n] + GtR[5,5,2,2,n]) 
                  # + nrm1[2,n]*nrm1[3,n] * (GtL[5,5,2,3,n] + GtR[5,5,2,3,n]) 
                  # + nrm1[3,n]*nrm1[1,n] * (GtL[5,5,3,1,n] + GtR[5,5,3,1,n]) 
                  # + nrm1[3,n]*nrm1[2,n] * (GtL[5,5,3,2,n] + GtR[5,5,3,2,n]) 
                  + nrm1[3,n]*nrm1[3,n] * (GtL[5,5,3,3,n] + GtR[5,5,3,3,n]) 
                 ) * coef

  end # end of loop over nodes
  return nothing
end

"""
  compute the interior penalty matrix
  **Input**
   * mesh
   * sbp
   * eqn
   * opts
   * iface: the index of interface
  **Input/Output**
   * pMat: penalty matrix
"""
function cmptBPMat{Tmsh, Tdim, Tsol, Tres}(mesh::AbstractMesh{Tmsh},
                                           sbp::AbstractSBP,
                                           eqn::EulerData{Tsol, Tres, Tdim},
                                           opts,
                                           iface::Int,
                                           Gt::AbstractArray{Tsol, 5},
                                           pMat::AbstractArray{Tsol, 3})
  if opts["SAT_type"] == "Hartman"
    cmptBPMat_hartman(mesh, sbp, eqn, opts, iface, Gt, pMat)
  elseif opts["SAT_type"] == "SAT-SIPG"
    cmptBPMat_SIPG(mesh, sbp, eqn, opts, iface, Gt, pMat)
  elseif opts["SAT_type"] == "SAT-BR2"
    cmptBPMat_BR2(mesh, sbp, eqn, opts, iface, Gt, pMat)
  end
end

function cmptBPMat_SIPG{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh},
                                          sbp::AbstractSBP,
                                          eqn::EulerData{Tsol, Tres, 2},
                                          opts,
                                          iface::Int,
                                          Gt::AbstractArray{Tsol, 5},
                                          pMat::AbstractArray{Tsol, 3})
  error("SAT-SIPG not available yet")
end

function cmptBPMat_BR2{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh},
                                         sbp::AbstractSBP,
                                         eqn::EulerData{Tsol, Tres, 2},
                                         opts,
                                         iface::Int,
                                         Gt::AbstractArray{Tsol, 5},
                                         pMat::AbstractArray{Tsol, 3})
  error("SAT-BR2 not available yet")
end

function cmptBPMat_hartman{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh},
                                             sbp::AbstractSBP,
                                             eqn::EulerData{Tsol, Tres, 2},
                                             opts,
                                             iface::Int,
                                             Gt::AbstractArray{Tsol, 5},
                                             pMat::AbstractArray{Tsol, 3})
# Reference:
#    Shahbazi, Mavriplis, Multigrid algorithms for high order discontinuous 
#    Galerkin discretization of the compressible Navier-Stokes equations, 
#    JCP, volume 228 Issue 21, November, 2009
#
  # in order to enhance/test stability, we try this factor on penalty terms
  factor = eqn.params.penalty_relaxation
  const_tii = eqn.params.const_tii
  Tdim = 2
  sbpface = mesh.sbpface

  face = mesh.bndryfaces[iface]
  elem = face.element

  # Compute geometric info on face
  nrm1 = zeros(Tmsh, Tdim, mesh.numNodesPerFace)
  area = zeros(Tmsh, mesh.numNodesPerFace)
  for n = 1 : mesh.numNodesPerFace
    nrm_xy = ro_sview(mesh.nrm_bndry, :, n, iface)
    area[n] = norm(nrm_xy)

    for i = 1 : Tdim
      nrm1[i,n] = nrm_xy[i] / area[n]
    end
  end

  # Compute the size of element and face, and then meas(face)/meas(elem)
  elem_vol = 0.0
  face_area = 0.0
  for n = 1:mesh.numNodesPerElement
    elem_vol +=  sbp.w[n]/mesh.jac[n, elem]
  end

  for n = 1:mesh.numNodesPerFace
    face_area +=  sbpface.wface[n]*area[n]
  end

  he = elem_vol/eqn.area_sum[elem]

  for n = 1 : mesh.numNodesPerFace
    # for iDof = 1:mesh.numDofPerNode
      # for jDof = 1:mesh.numDofPerNode
        # pMat[iDof, jDof, n] = 0.0
        # for jDim = 1 : Tdim
          # for iDim = 1 : Tdim
            # pMat[iDof,jDof,n] +=  nrm1[iDim,n]*nrm1[jDim,n] * (GtL[iDof,jDof,iDim,jDim,n] + GtR[iDof,jDof,iDim,jDim,n])
          # end
        # end
        # pMat[iDof,jDof,n] *=  0.5*const_tii * area[n] / he * factor
      # end

    # end
    coef = const_tii * area[n] / he * factor

    pMat[2,1,n] = ( nrm1[1,n]*nrm1[1,n] * Gt[2,1,1,1,n]
                  + nrm1[1,n]*nrm1[2,n] * Gt[2,1,1,2,n]
                  + nrm1[2,n]*nrm1[1,n] * Gt[2,1,2,1,n]
                  + nrm1[2,n]*nrm1[2,n] * Gt[2,1,2,2,n]
                 ) * coef 

    pMat[2,2,n] = ( nrm1[1,n]*nrm1[1,n] * Gt[2,2,1,1,n]
                  + nrm1[2,n]*nrm1[2,n] * Gt[2,2,2,2,n]
                 ) * coef 

    pMat[2,3,n] = ( nrm1[1,n]*nrm1[2,n] * Gt[2,3,1,2,n]
                  + nrm1[2,n]*nrm1[1,n] * Gt[2,3,2,1,n] 
                 ) * coef 

    pMat[3,1,n] = ( nrm1[1,n]*nrm1[1,n] * Gt[3,1,1,1,n]
                  + nrm1[1,n]*nrm1[2,n] * Gt[3,1,1,2,n]
                  + nrm1[2,n]*nrm1[1,n] * Gt[3,1,2,1,n]
                  + nrm1[2,n]*nrm1[2,n] * Gt[3,1,2,2,n] 
                 ) * coef 

    pMat[3,2,n] = ( nrm1[1,n]*nrm1[2,n] * Gt[3,2,1,2,n]
                  + nrm1[2,n]*nrm1[1,n] * Gt[3,2,2,1,n] 
                 ) * coef 

    pMat[3,3,n] = ( nrm1[1,n]*nrm1[1,n] * Gt[3,3,1,1,n]
                  + nrm1[2,n]*nrm1[2,n] * Gt[3,3,2,2,n] 
                 ) * coef 

    pMat[4,1,n] = ( nrm1[1,n]*nrm1[1,n] * Gt[4,1,1,1,n]
                  + nrm1[1,n]*nrm1[2,n] * Gt[4,1,1,2,n]
                  + nrm1[2,n]*nrm1[1,n] * Gt[4,1,2,1,n]
                  + nrm1[2,n]*nrm1[2,n] * Gt[4,1,2,2,n] 
                 ) * coef 

    pMat[4,2,n] = ( nrm1[1,n]*nrm1[1,n] * Gt[4,2,1,1,n]
                  + nrm1[1,n]*nrm1[2,n] * Gt[4,2,1,2,n]
                  + nrm1[2,n]*nrm1[1,n] * Gt[4,2,2,1,n]
                  + nrm1[2,n]*nrm1[2,n] * Gt[4,2,2,2,n] 
                 ) * coef 

    pMat[4,3,n] = ( nrm1[1,n]*nrm1[1,n] * Gt[4,3,1,1,n]
                  + nrm1[1,n]*nrm1[2,n] * Gt[4,3,1,2,n]
                  + nrm1[2,n]*nrm1[1,n] * Gt[4,3,2,1,n]
                  + nrm1[2,n]*nrm1[2,n] * Gt[4,3,2,2,n] 
                 ) * coef 

    pMat[4,4,n] = ( nrm1[1,n]*nrm1[1,n] * Gt[4,4,1,1,n]
                  + nrm1[2,n]*nrm1[2,n] * Gt[4,4,2,2,n] 
                 ) * coef 

  end # end of loop over nodes
  return nothing
end

function cmptBPMat_hartman{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh},
                                             sbp::AbstractSBP,
                                             eqn::EulerData{Tsol, Tres, 3},
                                             opts,
                                             iface::Int,
                                             Gt::AbstractArray{Tsol, 5},
                                             pMat::AbstractArray{Tsol, 3})
  # in order to enhance/test stability, we try this factor on penalty terms
  factor = eqn.params.penalty_relaxation
  const_tii = eqn.params.const_tii
  Tdim = 3
  sbpface = mesh.sbpface
  face = mesh.bndryfaces[iface]
  elem = face.element

  # Compute geometric info on face
  nrm1 = Array(Tmsh, Tdim, mesh.numNodesPerFace)
  area = Array(Tmsh, mesh.numNodesPerFace)
  for n = 1 : mesh.numNodesPerFace
    nrm_xy = ro_sview(mesh.nrm_bndry, :, n, iface)
    area[n] = norm(nrm_xy)

    for i = 1 : Tdim
      nrm1[i,n] = nrm_xy[i] / area[n]
    end
  end

  # Compute the size of element and face, and then meas(face)/meas(elem)
  elem_vol = 0.0
  face_area = 0.0
  for n = 1:mesh.numNodesPerElement
    elem_vol +=  sbp.w[n]/mesh.jac[n, elem]
  end

  for n = 1:mesh.numNodesPerFace
    face_area +=  sbpface.wface[n]*area[n]
  end

  he = elem_vol/eqn.area_sum[elem]

  for n = 1 : mesh.numNodesPerFace
    # for iDof = 1:mesh.numDofPerNode
      # for jDof = 1:mesh.numDofPerNode
        # pMat[iDof, jDof, n] = 0.0
        # for jDim = 1 : Tdim
          # for iDim = 1 : Tdim
            # pMat[iDof,jDof,n] +=  nrm1[iDim,n]*nrm1[jDim,n] * (GtL[iDof,jDof,iDim,jDim,n] + GtR[iDof,jDof,iDim,jDim,n])
          # end
        # end
        # pMat[iDof,jDof,n] *=  0.5*const_tii * area[n] / he * factor
      # end

    # end
    coef = const_tii * area[n] / he * factor
                 

    pMat[2,1,n] = ( nrm1[1,n]*nrm1[1,n] * Gt[2,1,1,1,n]
                  + nrm1[1,n]*nrm1[2,n] * Gt[2,1,1,2,n]
                  + nrm1[1,n]*nrm1[3,n] * Gt[2,1,1,3,n]
                  + nrm1[2,n]*nrm1[1,n] * Gt[2,1,2,1,n]
                  + nrm1[2,n]*nrm1[2,n] * Gt[2,1,2,2,n] 
                  # + nrm1[2,n]*nrm1[3,n] * Gt[2,1,2,3,n]
                  + nrm1[3,n]*nrm1[1,n] * Gt[2,1,3,1,n]
                  # + nrm1[3,n]*nrm1[2,n] * Gt[2,1,3,2,n]
                  + nrm1[3,n]*nrm1[3,n] * Gt[2,1,3,3,n]
                 ) * coef

    
    pMat[2,2,n] = ( nrm1[1,n]*nrm1[1,n] * Gt[2,2,1,1,n]
                  # + nrm1[1,n]*nrm1[2,n] * Gt[2,2,1,2,n]
                  # + nrm1[1,n]*nrm1[3,n] * Gt[2,2,1,3,n]
                  # + nrm1[2,n]*nrm1[1,n] * Gt[2,2,2,1,n]
                  + nrm1[2,n]*nrm1[2,n] * Gt[2,2,2,2,n] 
                  # + nrm1[2,n]*nrm1[3,n] * Gt[2,2,2,3,n] 
                  # + nrm1[3,n]*nrm1[1,n] * Gt[2,2,3,1,n] 
                  # + nrm1[3,n]*nrm1[2,n] * Gt[2,2,3,2,n] 
                  + nrm1[3,n]*nrm1[3,n] * Gt[2,2,3,3,n] 
                 ) * coef

    pMat[2,3,n] = ( 
                   # nrm1[1,n]*nrm1[1,n] * Gt[2,3,1,1,n]
                  + nrm1[1,n]*nrm1[2,n] * Gt[2,3,1,2,n]
                  # + nrm1[1,n]*nrm1[3,n] * Gt[2,3,1,3,n]
                  + nrm1[2,n]*nrm1[1,n] * Gt[2,3,2,1,n]
                  # + nrm1[2,n]*nrm1[2,n] * Gt[2,3,2,2,n] 
                  # + nrm1[2,n]*nrm1[3,n] * Gt[2,3,2,3,n] 
                  # + nrm1[3,n]*nrm1[1,n] * Gt[2,3,3,1,n] 
                  # + nrm1[3,n]*nrm1[2,n] * Gt[2,3,3,2,n] 
                  # + nrm1[3,n]*nrm1[3,n] * Gt[2,3,3,3,n] 
                 ) * coef

    pMat[2,4,n] = ( 
                    # nrm1[1,n]*nrm1[1,n] * Gt[2,4,1,1,n]
                  # + nrm1[1,n]*nrm1[2,n] * Gt[2,4,1,2,n]
                  + nrm1[1,n]*nrm1[3,n] * Gt[2,4,1,3,n]
                  # + nrm1[2,n]*nrm1[1,n] * Gt[2,4,2,1,n]
                  # + nrm1[2,n]*nrm1[2,n] * Gt[2,4,2,2,n] 
                  # + nrm1[2,n]*nrm1[3,n] * Gt[2,4,2,3,n] 
                  + nrm1[3,n]*nrm1[1,n] * Gt[2,4,3,1,n] 
                  # + nrm1[3,n]*nrm1[2,n] * Gt[2,4,3,2,n] 
                  # + nrm1[3,n]*nrm1[3,n] * Gt[2,4,3,3,n] 
                 ) * coef

    # pMat[2,5,n] = ( 
                    # # nrm1[1,n]*nrm1[1,n] * Gt[2,5,1,1,n]
                  # # + nrm1[1,n]*nrm1[2,n] * Gt[2,5,1,2,n]
                  # # + nrm1[1,n]*nrm1[3,n] * Gt[2,5,1,3,n]
                  # # + nrm1[2,n]*nrm1[1,n] * Gt[2,5,2,1,n]
                  # # + nrm1[2,n]*nrm1[2,n] * Gt[2,5,2,2,n] 
                  # # + nrm1[2,n]*nrm1[3,n] * Gt[2,5,2,3,n] 
                  # # + nrm1[3,n]*nrm1[1,n] * Gt[2,5,3,1,n] 
                  # # + nrm1[3,n]*nrm1[2,n] * Gt[2,5,3,2,n] 
                  # # + nrm1[3,n]*nrm1[3,n] * Gt[2,5,3,3,n] )

    pMat[3,1,n] = ( nrm1[1,n]*nrm1[1,n] * Gt[3,1,1,1,n]
                  + nrm1[1,n]*nrm1[2,n] * Gt[3,1,1,2,n]
                  # + nrm1[1,n]*nrm1[3,n] * Gt[3,1,1,3,n]
                  + nrm1[2,n]*nrm1[1,n] * Gt[3,1,2,1,n]
                  + nrm1[2,n]*nrm1[2,n] * Gt[3,1,2,2,n] 
                  + nrm1[2,n]*nrm1[3,n] * Gt[3,1,2,3,n] 
                  # + nrm1[3,n]*nrm1[1,n] * Gt[3,1,3,1,n] 
                  + nrm1[3,n]*nrm1[2,n] * Gt[3,1,3,2,n] 
                  + nrm1[3,n]*nrm1[3,n] * Gt[3,1,3,3,n] 
                 ) * coef

    pMat[3,2,n] = ( # nrm1[1,n]*nrm1[1,n] * Gt[3,2,1,1,n]
                  + nrm1[1,n]*nrm1[2,n] * Gt[3,2,1,2,n]
                  # + nrm1[1,n]*nrm1[3,n] * Gt[3,2,1,3,n]
                  + nrm1[2,n]*nrm1[1,n] * Gt[3,2,2,1,n]
                  # + nrm1[2,n]*nrm1[2,n] * Gt[3,2,2,2,n] 
                  # + nrm1[2,n]*nrm1[3,n] * Gt[3,2,2,3,n] 
                  # + nrm1[3,n]*nrm1[1,n] * Gt[3,2,3,1,n] 
                  # + nrm1[3,n]*nrm1[2,n] * Gt[3,2,3,2,n] 
                  # + nrm1[3,n]*nrm1[3,n] * Gt[3,2,3,3,n] 
                 ) * coef

    pMat[3,3,n] = ( nrm1[1,n]*nrm1[1,n] * Gt[3,3,1,1,n]
                  # + nrm1[1,n]*nrm1[2,n] * Gt[3,3,1,2,n]
                  # + nrm1[1,n]*nrm1[3,n] * Gt[3,3,1,3,n]
                  # + nrm1[2,n]*nrm1[1,n] * Gt[3,3,2,1,n]
                  + nrm1[2,n]*nrm1[2,n] * Gt[3,3,2,2,n] 
                  # + nrm1[2,n]*nrm1[3,n] * Gt[3,3,2,3,n] 
                  # + nrm1[3,n]*nrm1[1,n] * Gt[3,3,3,1,n] 
                  # + nrm1[3,n]*nrm1[2,n] * Gt[3,3,3,2,n] 
                  + nrm1[3,n]*nrm1[3,n] * Gt[3,3,3,3,n] 
                 ) * coef

    pMat[3,4,n] = ( # nrm1[1,n]*nrm1[1,n] * Gt[3,4,1,1,n]
                  # + nrm1[1,n]*nrm1[2,n] * Gt[3,4,1,2,n]
                  # + nrm1[1,n]*nrm1[3,n] * Gt[3,4,1,3,n]
                  # + nrm1[2,n]*nrm1[1,n] * Gt[3,4,2,1,n]
                  # + nrm1[2,n]*nrm1[2,n] * Gt[3,4,2,2,n] 
                  + nrm1[2,n]*nrm1[3,n] * Gt[3,4,2,3,n] 
                  # + nrm1[3,n]*nrm1[1,n] * Gt[3,4,3,1,n] 
                  + nrm1[3,n]*nrm1[2,n] * Gt[3,4,3,2,n] 
                  # + nrm1[3,n]*nrm1[3,n] * Gt[3,4,3,3,n] 
                 ) * coef

    # pMat[3,5,n] = ( 
                    # # nrm1[1,n]*nrm1[1,n] * Gt[3,5,1,1,n]
                  # # + nrm1[1,n]*nrm1[2,n] * Gt[3,5,1,2,n]
                  # # + nrm1[1,n]*nrm1[3,n] * Gt[3,5,1,3,n]
                  # # + nrm1[2,n]*nrm1[1,n] * Gt[3,5,2,1,n)
                  # # + nrm1[2,n]*nrm1[2,n] * Gt[3,5,2,2,n] 
                  # # + nrm1[2,n]*nrm1[3,n] * Gt[3,5,2,3,n] 
                  # # + nrm1[3,n]*nrm1[1,n] * Gt[3,5,3,1,n] 
                  # # + nrm1[3,n]*nrm1[2,n] * Gt[3,5,3,2,n] 
                  # # + nrm1[3,n]*nrm1[3,n] * Gt[3,5,3,3,n] 
                 # )

    pMat[4,1,n] = ( nrm1[1,n]*nrm1[1,n] * Gt[4,1,1,1,n]
                  # + nrm1[1,n]*nrm1[2,n] * Gt[4,1,1,2,n]
                  + nrm1[1,n]*nrm1[3,n] * Gt[4,1,1,3,n]
                  # + nrm1[2,n]*nrm1[1,n] * Gt[4,1,2,1,n]
                  + nrm1[2,n]*nrm1[2,n] * Gt[4,1,2,2,n] 
                  + nrm1[2,n]*nrm1[3,n] * Gt[4,1,2,3,n] 
                  + nrm1[3,n]*nrm1[1,n] * Gt[4,1,3,1,n] 
                  + nrm1[3,n]*nrm1[2,n] * Gt[4,1,3,2,n] 
                  + nrm1[3,n]*nrm1[3,n] * Gt[4,1,3,3,n] 
                 ) * coef

    pMat[4,2,n] = ( 
                    # nrm1[1,n]*nrm1[1,n] * Gt[4,2,1,1,n]
                  # + nrm1[1,n]*nrm1[2,n] * Gt[4,2,1,2,n]
                  + nrm1[1,n]*nrm1[3,n] * Gt[4,2,1,3,n]
                  # + nrm1[2,n]*nrm1[1,n] * Gt[4,2,2,1,n]
                  # + nrm1[2,n]*nrm1[2,n] * Gt[4,2,2,2,n] 
                  # + nrm1[2,n]*nrm1[3,n] * Gt[4,2,2,3,n] 
                  + nrm1[3,n]*nrm1[1,n] * Gt[4,2,3,1,n] 
                  # + nrm1[3,n]*nrm1[2,n] * Gt[4,2,3,2,n] 
                  # + nrm1[3,n]*nrm1[3,n] * Gt[4,2,3,3,n] 
                 ) * coef

    pMat[4,3,n] = ( 
                    # nrm1[1,n]*nrm1[1,n] * Gt[4,3,1,1,n]
                  # + nrm1[1,n]*nrm1[2,n] * Gt[4,3,1,2,n]
                  # + nrm1[1,n]*nrm1[3,n] * Gt[4,3,1,3,n]
                  # + nrm1[2,n]*nrm1[1,n] * Gt[4,3,2,1,n]
                  # + nrm1[2,n]*nrm1[2,n] * Gt[4,3,2,2,n] 
                  + nrm1[2,n]*nrm1[3,n] * Gt[4,3,2,3,n] 
                  # + nrm1[3,n]*nrm1[1,n] * Gt[4,3,3,1,n] 
                  + nrm1[3,n]*nrm1[2,n] * Gt[4,3,3,2,n] 
                  # + nrm1[3,n]*nrm1[3,n] * Gt[4,3,3,3,n] 
                 ) * coef

    pMat[4,4,n] = ( nrm1[1,n]*nrm1[1,n] * Gt[4,4,1,1,n]
                  # + nrm1[1,n]*nrm1[2,n] * Gt[4,4,1,2,n]
                  # + nrm1[1,n]*nrm1[3,n] * Gt[4,4,1,3,n]
                  # + nrm1[2,n]*nrm1[1,n] * Gt[4,4,2,1,n]
                  + nrm1[2,n]*nrm1[2,n] * Gt[4,4,2,2,n] 
                  # + nrm1[2,n]*nrm1[3,n] * Gt[4,4,2,3,n] 
                  # + nrm1[3,n]*nrm1[1,n] * Gt[4,4,3,1,n] 
                  # + nrm1[3,n]*nrm1[2,n] * Gt[4,4,3,2,n] 
                  + nrm1[3,n]*nrm1[3,n] * Gt[4,4,3,3,n] 
                 ) * coef

    # pMat[4,5,n] = ( 
                    # # nrm1[1,n]*nrm1[1,n] * Gt[4,5,1,1,n]
                  # # + nrm1[1,n]*nrm1[2,n] * Gt[4,5,1,2,n]
                  # # + nrm1[1,n]*nrm1[3,n] * Gt[4,5,1,3,n)
                  # # + nrm1[2,n]*nrm1[1,n] * Gt[4,5,2,1,n]
                  # # + nrm1[2,n]*nrm1[2,n] * Gt[4,5,2,2,n] 
                  # # + nrm1[2,n]*nrm1[3,n] * Gt[4,5,2,3,n] 
                  # # + nrm1[3,n]*nrm1[1,n] * Gt[4,5,3,1,n] 
                  # # + nrm1[3,n]*nrm1[2,n] * Gt[4,5,3,2,n] 
                  # # + nrm1[3,n]*nrm1[3,n] * Gt[4,5,3,3,n] 
                 # )

    pMat[5,1,n] = ( nrm1[1,n]*nrm1[1,n] * Gt[5,1,1,1,n]
                  + nrm1[1,n]*nrm1[2,n] * Gt[5,1,1,2,n]
                  + nrm1[1,n]*nrm1[3,n] * Gt[5,1,1,3,n]
                  + nrm1[2,n]*nrm1[1,n] * Gt[5,1,2,1,n]
                  + nrm1[2,n]*nrm1[2,n] * Gt[5,1,2,2,n] 
                  + nrm1[2,n]*nrm1[3,n] * Gt[5,1,2,3,n] 
                  + nrm1[3,n]*nrm1[1,n] * Gt[5,1,3,1,n] 
                  + nrm1[3,n]*nrm1[2,n] * Gt[5,1,3,2,n] 
                  + nrm1[3,n]*nrm1[3,n] * Gt[5,1,3,3,n] 
                 ) * coef


    pMat[5,2,n] = ( nrm1[1,n]*nrm1[1,n] * Gt[5,2,1,1,n]
                  + nrm1[1,n]*nrm1[2,n] * Gt[5,2,1,2,n]
                  + nrm1[1,n]*nrm1[3,n] * Gt[5,2,1,3,n]
                  + nrm1[2,n]*nrm1[1,n] * Gt[5,2,2,1,n]
                  + nrm1[2,n]*nrm1[2,n] * Gt[5,2,2,2,n] 
                  # + nrm1[2,n]*nrm1[3,n] * Gt[5,2,2,3,n] 
                  + nrm1[3,n]*nrm1[1,n] * Gt[5,2,3,1,n] 
                  # + nrm1[3,n]*nrm1[2,n] * Gt[5,2,3,2,n] 
                  + nrm1[3,n]*nrm1[3,n] * Gt[5,2,3,3,n] 
                 ) * coef

    pMat[5,3,n] = ( nrm1[1,n]*nrm1[1,n] * Gt[5,3,1,1,n]
                  + nrm1[1,n]*nrm1[2,n] * Gt[5,3,1,2,n]
                  # + nrm1[1,n]*nrm1[3,n] * Gt[5,3,1,3,n]
                  + nrm1[2,n]*nrm1[1,n] * Gt[5,3,2,1,n]
                  + nrm1[2,n]*nrm1[2,n] * Gt[5,3,2,2,n] 
                  + nrm1[2,n]*nrm1[3,n] * Gt[5,3,2,3,n] 
                  # + nrm1[3,n]*nrm1[1,n] * Gt[5,3,3,1,n] 
                  + nrm1[3,n]*nrm1[2,n] * Gt[5,3,3,2,n] 
                  + nrm1[3,n]*nrm1[3,n] * Gt[5,3,3,3,n] 
                 ) * coef

    pMat[5,4,n] = ( nrm1[1,n]*nrm1[1,n] * Gt[5,4,1,1,n]
                  # + nrm1[1,n]*nrm1[2,n] * Gt[5,4,1,2,n]
                  + nrm1[1,n]*nrm1[3,n] * Gt[5,4,1,3,n]
                  # + nrm1[2,n]*nrm1[1,n] * Gt[5,4,2,1,n]
                  + nrm1[2,n]*nrm1[2,n] * Gt[5,4,2,2,n] 
                  + nrm1[2,n]*nrm1[3,n] * Gt[5,4,2,3,n] 
                  + nrm1[3,n]*nrm1[1,n] * Gt[5,4,3,1,n] 
                  + nrm1[3,n]*nrm1[2,n] * Gt[5,4,3,2,n] 
                  + nrm1[3,n]*nrm1[3,n] * Gt[5,4,3,3,n] 
                 ) * coef

    pMat[5,5,n] = ( nrm1[1,n]*nrm1[1,n] * Gt[5,5,1,1,n]
                  # + nrm1[1,n]*nrm1[2,n] * Gt[5,5,1,2,n]
                  # + nrm1[1,n]*nrm1[3,n] * Gt[5,5,1,3,n]
                  # + nrm1[2,n]*nrm1[1,n] * Gt[5,5,2,1,n]
                  + nrm1[2,n]*nrm1[2,n] * Gt[5,5,2,2,n] 
                  # + nrm1[2,n]*nrm1[3,n] * Gt[5,5,2,3,n] 
                  # + nrm1[3,n]*nrm1[1,n] * Gt[5,5,3,1,n] 
                  # + nrm1[3,n]*nrm1[2,n] * Gt[5,5,3,2,n] 
                  + nrm1[3,n]*nrm1[3,n] * Gt[5,5,3,3,n] 
                 ) * coef

  end # end of loop over nodes
  return nothing
end
