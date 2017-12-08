include("viscous_penalty.jl")

"""
  Similar to euler's calcSharedFaceIntegrals_nopre_element_inner.
  Actually, relatively identical.
"""

function calcSharedFaceIntegrals_viscous{Tmsh, Tsol, Tres}(
                            mesh::AbstractDGMesh{Tmesh},
                            sbp::AbstractSBP,
                            eqn::EulerData{Tsol, Tres},
                            opts,
                            data::SharedFaceData,
                            functor::FluxType)

  q = eqn.q
  params = eqn.params

  q_faceL = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  q_faceR = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)

  idx = data.peeridx                    # index of this peer in mesh.peer_parts
  interfaces = data.interfaces          # vector of Interfaces
  bndries_local = data.bndries_local    # vector of Boundaries describing faces from local side
  bndries_remote = data.bndries_remote  # vector of Boundaries describing faces from remote side

                          # see comment in corresponding inviscid parallel function for why no qL_arr
  qR_arr = data.q_recv
  nrm_arr = mesh.nrm_sharedface[idx]
  aux_vars_arr = eqn.aux_vars_sharedface[idx]
  flux_face = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)   # populated in the inner loop, then
                                                                      # passed into boundaryFaceIntegrate!
  
  start_elnum = mesh.shared_element_offsets[idx]                                                                    

  for j = 1:length(interfaces)
    iface_j = interfaces[j]       # this interface      # TODO: what is the type of iface_j?
    bndryL_j = bndries_local[j]   # this L bndry
    bndryR_j = bndries_remote[j]  # this R bndry

    fL = bndryL_j.face            # the face of this L bndry

    # interpolate to face
    qL = ro_sview(q, :, :, iface_j.elementL)    # q of local element of this interface
    el_r = iface_j.elementR - start_elnum + 1   # remote element of this interface. need to account for
                                                #   offset in mesh.shared_element_offsets
    qR = ro_sview(qR_arr, :, :, el_r)                                            
    interiorFaceInterpolate!(mesh.sbpface, iface_j, qL, qR, q_faceL, q_faceR)

    # calculate flux
    for k = 1:mesh.numNodesPerFace
      qL_k = ro_sview(q_faceL, :, k)      # this node's q data (all dof's) from qL
      qR_k = ro_sview(q_faceR, :, k)      #       "                        from qR
      aux_vars = ro_sview(aux_vars_arr, :, k, j)
      nrm_xy = ro_sview(nrm_arr, :, k, j)

      flux_k = sview(flux_face, :, k)     # this node's flux. will need to be modified by functor below 
                                          #   to populate the full flux_face to pass to
                                          #   boundaryFaceIntegrate!

      parent(aux_vars)[1] = calcPressure(params, qL_k)

      functor(params, qL_k, qR_k, aux_vars, nrm_xy, flux_k)
    end   # end of loop: k = 1:mesh.numNodesPerFace

    # do the integration
    res_j = sview(eqn.res, :, :, bndryL_j.element)     # get the portion of res corresponding to this element
    boundaryFaceIntegrate!(mesh.sbpface, fL, flux_face, res_j, SummationByParts.Subtract())

  end   # end of loop: j = 1:length(interfaces)

  return nothing

end
                      
function evalSharedFaceIntegrals_viscous(mesh::AbstractDGMesh, sbp, eqn, opts)

  # Notes:
  #   1) Viscous terms will always require that full element data is shared, not just faces.
  #      Therefore, opts["parallel_data"] must be "element"
  #   2) Currently, entropy conservative and entropy stable formulations are not implemented 
  #      for viscous terms.

  if face_integral_type == 1    # non-entropy variables
    if opts["parallel_data"] == "face"
      throw(ErrorException("Viscous terms require shared element and face data, but parallel_data = face. Exiting."))
    elseif opts["parallel_data"] == "element"
      finishExchangeData(mesh, sbp, eqn, opts, eqn.shared_data, calcSharedFaceIntegrals_viscous)
    else
      throw(ErrorException("unsupported parallel_data value: $parallel_data"))
    end
  elseif face_integral_type == 2
    throw(ErrorException("face_integral_type $face_integral_type unsupported for viscous terms currently."))
  else
    throw(ErrorException("unsupported face integral type = $face_integral_type"))
  end

  return nothing

end
