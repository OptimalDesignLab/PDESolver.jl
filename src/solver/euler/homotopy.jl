# file for homotopy functions
# TODO: add mechanism for selecting which homotopy function to use


"""
  Calculate a first order accurate dissipation to use as a homotopy function


  Inputs:
    mesh: a DG mesh
    sbp: an SBP operator
    eqn: an EulerData object
    opts: options dictionary
  
  Inputs/Outputs:
    res: 3D array to store the homotopy function in

  Note eqn.res is *not* modified by this function.

  Aliasing restrictions: none
"""
function calcHomotopyDiss{Tsol, Tres, Tmsh}(mesh::AbstractDGMesh{Tmsh}, sbp, 
                          eqn::EulerData{Tsol, Tres}, opts, 
                          res::Abstract3DArray{Tres})

  @assert mesh.commsize == 1  # this doesn't work in parallel yet

  #----------------------------------------------------------------------------
  # volume dissipation

  dudxi = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl, mesh.dim)

  for i=1:mesh.dim
    differentiate!(sbp, i, eqn.q, sview(dudxi, :, :, :, i))
  end

  # scale by maximum wave speed
  nrm = zeros(Tmsh, mesh.dim)
  for dim=1:mesh.dim  # loop over parametric directions
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        q_j = sview(eqn.q, :, j, i)

        # get vector in xi direction defind by dim
        for k=1:mesh.dim
          nrm[i] = mesh.dxidx[dim, k, j, i]
        end

        # get maximum eigenvalue
        lambda_max = getLambdaMax(eqn.params, q_j, nrm)

        # scale by maximum eigenvalue
        for k=1:mesh.numDofPerNode
          dudxi[k, j, i, dim] *= lambda_max
        end

      end  # end loop j
    end  # end loop i
  end  # end loop over parametric directions

  # apply D^T * H
  for i=1:mesh.dim
    weakdifferentiate!(sbp, i, sview(dudxi, :, :, :, i), res, trans=true)
  end

  #----------------------------------------------------------------------------
  # interface dissipation

  # interpolate to face rather than using eqn.q_face, in case it hasn't been
  # updated since eqn.q was updated

  q_faceL = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  q_faceR = zeros(q_faceL)
  nrm2 = eqn.params.nrm2
  flux = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  for i=1:mesh.numInterfaces
    iface_i = mesh.interfaces[i]
    qL = sview(eqn.q, :, :, iface_i.elementL)
    qR = sview(eqn.q, :, :, iface_i.elementR)
    resL = sview(res, :, :, iface_i.elementL)
    resR = sview(res, :, :, iface_i.elementR)
    
    interiorFaceInterpolate!(mesh.sbpface, iface_i, qL, qR, q_faceL, q_faceR)

    # calculate numerical flux at each face node
    for j=1:mesh.numNodesPerFace
      qL_j = sview(q_faceL, :, j)
      qR_j = sview(q_faceR, :, j)

      # get the face normal
      dxidx_j = sview(mesh.dxidx_face, :, :, j, i)
      nrm_xi = sview(mesh.sbpface.normal, :, iface_i.faceL)
      calcBCNormal(eqn.params, dxidx_j, nrm_xi, nrm2)

      lambda_max = getLambdaMaxSimple(eqn.params, qL_j, qR_j, nrm2)

      for k=1:mesh.numDofPerNode
        flux[k, j] = lambda_max*(qL_j[k] - qR_j[k])
      end
    end  # end loop j

    # integrate over the face
    interiorFaceIntegrate!(mesh.sbpface, iface_i, flux, resL, resR)
  end  # end loop i

  #----------------------------------------------------------------------------
  # boundary dissipation
  # use q_faceL, nrm2, flux  from interface dissipation
  qg = eqn.params.qg  # boundary state
  func = FreeStreamBC()  # construct the functor
  for i=1:mesh.numBoundaryFaces
    bndry_i = mesh.bndryfaces[i]
    qL = sview(eqn.q, :, :, bndry_i.element)
    resL = sview(res, :, :, bndry_i.element)

    boundaryFaceInterpolate!(mesh.sbpface, bndry_i.face, qL, q_faceL)

    for j=1:mesh.numNodesPerFace
      q_j = sview(q_faceL, :, j)
      dxidx_j = sview(mesh.dxidx_bndry, :, :, j, i)

      # calculate boundary state
      coords = sview(mesh.coords_bndry, :, j, i)
      calcFreeStream(coords, eqn.params, qg)

      # calculate face normal
      nrm_xi = sview(mesh.sbpface.normal, :, bndry_i.face)
      calcBCNormal(eqn.params, dxidx_j, nrm_xi, nrm2)


      # calculate lambda_max
      lambda_max = getLambdaMaxSimple(eqn.params, q_j, qg, nrm2)

      # calculate dissipation
      for k=1:mesh.numNodesPerElement
        flux[k, j] = lambda_max*(q_j[k] - qg[k])
      end
    end  # end loop j

    # integrate over the face
    boundaryFaceIntegrate!(mesh.sbpface, bndry_i.face, flux, resL)
  end  # end loop i


  return nothing
end

"""
  Calculate the maximum wave speed for a given state

  Inputs:
    params: a ParamType
    qL: a vector of conservative variables at a node
    dir: a direction vector, length 2 in 2D and 3 in 3D

  Outputs:
    lambda_max: the maximum wave speed

  Aliasing restrictions: none
"""
function getLambdaMax{Tsol, Tmsh, Tdim}(params::ParamType{Tdim}, 
                      qL::AbstractVector{Tsol}, 
                      dir::AbstractVector{Tmsh})

  Tres = promote_type(Tsol, Tmsh)
  gamma = params.gamma
  Un = zero(Tres)
  dA = zero(Tmsh)
  rhoLinv = 1/qL[1]

  pL = calcPressure(params, qL)
  aL = sqrt(gamma*pL*rhoLinv)  # speed of sound

  for i=1:Tdim
    Un += dir[i]*qL[i+1]*rhoLinv
    dA += dir[i]*dir[i]
  end

  lambda_max = absvalue(Un) + dA*aL

  return lambda_max
end




