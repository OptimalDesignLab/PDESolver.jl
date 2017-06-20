# file for homotopy functions
# TODO: add mechanism for selecting which homotopy function to use

import PDESolver.evalHomotopy

"""
  This function calls the appropriate homotopy function for the Euler module.
"""
function evalHomotopy(mesh::AbstractMesh, sbp::AbstractSBP, eqn::EulerData, opts::Dict, res::Abstract3DArray, t = 0.0)

  calcHomotopyDiss(mesh, sbp, eqn, opts, res)

  return nothing
end


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

#  println("\nentered calcHomotopyDiss")

  # some checks for when parallelism is enabled
  @assert opts["parallel_data"] == "element"
  for i=1:mesh.npeers
    @assert mesh.recv_waited[i]
  end

  fill!(res, 0.0)

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
          nrm[k] = mesh.dxidx[dim, k, j, i]
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
#  nrm2 = eqn.params.nrm2
  flux = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  for i=1:mesh.numInterfaces
    iface_i = mesh.interfaces[i]
    qL = sview(eqn.q, :, :, iface_i.elementL)
    qR = sview(eqn.q, :, :, iface_i.elementR)
    resL = sview(res, :, :, iface_i.elementL)
    resR = sview(res, :, :, iface_i.elementR)
#    fill!(q_faceL, 0.0)
#    fill!(q_faceR, 0.0)

    interiorFaceInterpolate!(mesh.sbpface, iface_i, qL, qR, q_faceL, q_faceR)

    # calculate numerical flux at each face node
    for j=1:mesh.numNodesPerFace
      qL_j = sview(q_faceL, :, j)
      qR_j = sview(q_faceR, :, j)

      # get the face normal
      nrm2 = sview(mesh.nrm_face, :, j, i)

      lambda_max = getLambdaMaxSimple(eqn.params, qL_j, qR_j, nrm2)

      for k=1:mesh.numDofPerNode
        flux[k, j] = 0.5*lambda_max*(qL_j[k] - qR_j[k])
      end
    end  # end loop j

    # integrate over the face
    interiorFaceIntegrate!(mesh.sbpface, iface_i, flux, resL, resR)
  end  # end loop i
#=
# the boundary term makes the predictor-corrector algorithm converge slower
  #----------------------------------------------------------------------------
  # boundary dissipation
  # use q_faceL, nrm2, flux  from interface dissipation
  qg = eqn.params.qg  # boundary state
  for i=1:mesh.numBoundaryFaces
    bndry_i = mesh.bndryfaces[i]
    qL = sview(eqn.q, :, :, bndry_i.element)
    resL = sview(res, :, :, bndry_i.element)
    fill!(q_faceL, 0.0)

    boundaryFaceInterpolate!(mesh.sbpface, bndry_i.face, qL, q_faceL)

#    q_faceL = sview(eqn.q_bndry, :, :, i)
    for j=1:mesh.numNodesPerFace
      q_j = sview(q_faceL, :, j)
#      dxidx_j = sview(mesh.dxidx_bndry, :, :, j, i)

      # calculate boundary state
      coords = sview(mesh.coords_bndry, :, j, i)
      calcFreeStream(coords, eqn.params, qg)

      # calculate face normal
      nrm2 = sview(mesh.nrm_bndry, :, j, i)

      # calculate lambda_max
      lambda_max = getLambdaMaxSimple(eqn.params, q_j, qg, nrm2)

      # calculate dissipation
      for k=1:mesh.numDofPerNode
        flux[k, j] = 0.5*lambda_max*(q_j[k] - qg[k])
      end

    end  # end loop j

    
    # integrate over the face
    boundaryFaceIntegrate!(mesh.sbpface, bndry_i.face, flux, resL)
  end  # end loop i
=#
  
  #----------------------------------------------------------------------------
  
  # shared face integrals
  # use q_faceL, q_faceR, flux from above
  workarr = zeros(q_faceR)
  for peer=1:mesh.npeers
    # get data for this peer
    bndries_local = mesh.bndries_local[peer]
    bndries_remote = mesh.bndries_remote[peer]
    interfaces_peer = mesh.shared_interfaces[peer]

    qR_peer = eqn.q_face_recv[peer]
#    dxidx_peer = mesh.dxidx_sharedface[peer]
    nrm_peer = mesh.nrm_sharedface[peer]
    start_elnum = mesh.shared_element_offsets[peer]

    for i=1:length(bndries_local)
      bndry_i = bndries_local[i]
      bndryR_i = bndries_remote[i]
      iface_i = interfaces_peer[i]
      qL_i = sview(eqn.q, :, :, bndry_i.element)
      qR_i = sview(qR_peer, :, :, iface_i.elementR - start_elnum + 1)
      resL = sview(res, :, :, bndry_i.element)

      # interpolate to face
#      interiorFaceInterpolate!(mesh.sbpface, iface_i, qL_i, qR_i, q_faceL, q_faceR)
      boundaryFaceInterpolate!(mesh.sbpface, bndry_i.face, qL_i, q_faceL)
      boundaryFaceInterpolate!(mesh.sbpface, bndryR_i.face, qR_i, q_faceR)

      # permute elementR
      permvec = sview(mesh.sbpface.nbrperm, :, iface_i.orient)
      SummationByParts.permuteface!(permvec, workarr, q_faceR)

      # compute flux at every face node
      for j=1:mesh.numNodesPerFace
        qL_j = sview(q_faceL, :, j)
        qR_j = sview(q_faceR, :, j)
        nrm2 = sview(nrm_peer, :, j, i)

        # get max wave speed
        lambda_max = getLambdaMaxSimple(eqn.params, qL_j, qR_j, nrm2)

        # calculate flux
        for k=1:mesh.numDofPerNode
          flux[k, j] = 0.5*lambda_max*(qL_j[k] - qR_j[k])
        end
      end  # end loop j

      # integrate over the face
      boundaryFaceIntegrate!(mesh.sbpface, bndry_i.face, flux, resL)
    end  # end loop i
  end  # end loop peer
  


  # negate for consistency with the physics module
  for i=1:length(res)
    res[i] = -res[i]
  end

#  println("homotopy residual norm = ", norm(vec(res)))

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




