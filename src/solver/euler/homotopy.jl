# file for homotopy functions
# TODO: add mechanism for selecting which homotopy function to use
# TODO: better idea: move the homotopy function to NLSolver, because it is
#                    physics agnostic

import PDESolver.evalHomotopy

"""
  This function calls the appropriate homotopy function for the Euler module.
"""
function evalHomotopy(mesh::AbstractMesh, sbp::AbstractOperator, eqn::EulerData, opts::Dict, res::Abstract3DArray, t = 0.0)

  
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
function calcHomotopyDiss(mesh::AbstractDGMesh{Tmsh}, sbp, 
        eqn::EulerData{Tsol, Tres}, opts, 
        res::Abstract3DArray{Tres}) where {Tsol, Tres, Tmsh}

#  println("\nentered calcHomotopyDiss")

  # some checks for when parallelism is enabled
  @assert opts["parallel_data"] == PARALLEL_DATA_ELEMENT
  for i=1:mesh.npeers
    @assert eqn.shared_data[i].recv_waited
  end

  fill!(res, 0.0)

  #----------------------------------------------------------------------------
  # volume dissipation


  t1 = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerElement)
  t2 = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)
  fill!(res, 0.0)
  nrm2 = zeros(Tmsh, mesh.dim)

  # compute the D operator in each direction
  D = zeros(mesh.numNodesPerElement, mesh.numNodesPerElement, mesh.dim)
  for d=1:mesh.dim
    D[:, :, d] = inv(diagm(sbp.w))*sbp.Q[:, :, d]
  end


  for el=1:mesh.numEl
    q_el = sview(eqn.q, :, :, el)
    res_el = sview(res, :, :, el)
    for d1=1:mesh.dim
      fill!(t1, 0.0)
      fill!(t2, 0.0)
      differentiateElement!(sbp, d1, q_el, t1)

      for j=1:mesh.numNodesPerElement
        q_j = sview(eqn.q, :, j, el)
        for k=1:mesh.dim
          nrm2[k] = mesh.dxidx[d1, k, j, el]
        end

        lambda_max = getLambdaMax(eqn.params, q_j, nrm2)
        for k=1:mesh.numDofPerNode
          t2[k, j] = lambda_max*t1[k, j]
        end

      end  # end loop j
     weakDifferentiateElement!(sbp, d1, t2, res_el, SummationByParts.Add(), true)
    end  # end loop dim
  end  # end loop el

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


# the boundary term makes the predictor-corrector algorithm converge slower
  #----------------------------------------------------------------------------
  # boundary dissipation
  # use q_faceL, nrm2, flux  from interface dissipation
  if opts["homotopy_addBoundaryIntegrals"]
    qg = zeros(Tsol, mesh.numDofPerNode)  # boundary state
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
        calcFreeStream(eqn.params, coords, qg)
#        calcInvChannelIC(eqn.params, coords, qg)

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
  end  # end if addBoundaryIntegrals

  
  #----------------------------------------------------------------------------
  
  # shared face integrals
  # use q_faceL, q_faceR, flux from above

  workarr = zeros(q_faceR)
  for peer=1:mesh.npeers
    # get data for this peer
    bndries_local = mesh.bndries_local[peer]
    interfaces_peer = mesh.shared_interfaces[peer]

    qR_peer = eqn.shared_data[peer].q_recv
    nrm_peer = mesh.nrm_sharedface[peer]
    start_elnum = mesh.shared_element_offsets[peer]

    for i=1:length(bndries_local)
      bndry_i = bndries_local[i]
      iface_i = interfaces_peer[i]
      qL_i = sview(eqn.q, :, :, bndry_i.element)
      qR_i = sview(qR_peer, :, :, iface_i.elementR - start_elnum + 1)
      resL = sview(res, :, :, bndry_i.element)

      # interpolate to face
      interiorFaceInterpolate!(mesh.sbpface, iface_i, qL_i, qR_i, q_faceL, q_faceR)
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





