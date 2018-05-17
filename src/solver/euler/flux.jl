# functions for calculating the flux for interior face integrals

@doc """
### EulerEquationMod.calcFaceFlux

  This function calculates the DG flux between a specified set of faces,
  using the solution data at the faces stored in eqn.q_face.
  Note that the flux is negated because the face integrals have a
  negative sign in the weak form.

  Conservative variables only!

  **Inputs**:

   * mesh
   * sbp
   * eqn
   * functor: the functor that calculates the flux at a node
   * interfaces: an array of type Interface that specifies which interfaces
                to calculate the flux for

  **Inputs/Outputs**:
   * face_flux: array to store the flux in, numDofPerNode x nnodesPerFace
               x length(interfaces)

  The functor must have the signature:

  `func( uL, qR, aux_vars, dxidx, nrm, flux_j, eqn.params)`

  where uL and uR are the solution values for a node on the left and right
  elements, aux_vars are the auxiliary variables for the node,
  dxidx is the scaled mapping jacobian for elementL, and nrm is the face
  normal in reference space. flux_j is the array of length numDofPerNode to be
  populated with the flux. params is eqn.params.

"""->
function calcFaceFlux{Tmsh, Tsol, Tres, Tdim}( mesh::AbstractDGMesh{Tmsh},
                          sbp::AbstractSBP,
                          eqn::EulerData{Tsol, Tres, Tdim, :conservative},
                          functor::FluxType,
                          interfaces::AbstractArray{Interface,1},
                          face_flux::AbstractArray{Tres, 3})

  nfaces = length(interfaces)
  params = eqn.params
  for i=1:nfaces  # loop over faces
    interface_i = interfaces[i]
    for j = 1:mesh.numNodesPerFace
      eL = interface_i.elementL
      fL = interface_i.faceL
      # get components
      qL = ro_sview(eqn.q_face, :, 1, j, i)
      qR = ro_sview(eqn.q_face, :, 2, j, i)
      aux_vars = ro_sview(eqn.aux_vars_face, :, j, i)
      nrm_xy = ro_sview(mesh.nrm_face, :, j, i)
      flux_j = sview(face_flux, :, j, i)

      functor(params, qL, qR, aux_vars, nrm_xy, flux_j)
    end
  end

  return nothing
end

"""
  Like [`calcFaceFlux`](@ref), but computes the flux for a single element and
  then integrates it immediately, updating eqn.res

  **Inputs**:
   * mesh
   * sbp
   * eqn
   * opts
   * functor: a FluxType that evalutes the flux
   * interfaces: the vector of [`Interface`](@ref)s to compute the integrals for

"""
function calcFaceIntegral_nopre{Tmsh, Tsol, Tres, Tdim}(
                                mesh::AbstractDGMesh{Tmsh},
                                sbp::AbstractSBP,
                                eqn::EulerData{Tsol, Tres, Tdim, :conservative},
                                opts,
                                functor::FluxType,
                                interfaces::AbstractArray{Interface, 1})


  nfaces = length(interfaces)
  params = eqn.params

  q_faceL = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  q_faceR = zeros(q_faceL)

  flux_face = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)

  for i=1:nfaces
    iface_i = interfaces[i]

    qL = ro_sview(eqn.q, :, :, iface_i.elementL)
    qR = ro_sview(eqn.q, :, :, iface_i.elementR)
    interiorFaceInterpolate!(mesh.sbpface, iface_i, qL, qR, q_faceL,
                             q_faceR)

    for j=1:mesh.numNodesPerFace
      qL_j = ro_sview(q_faceL, :, j)
      qR_j = ro_sview(q_faceR, :, j)

      eqn.aux_vars_face[1, j, i] = calcPressure(params, qL_j)
      aux_vars = ro_sview(eqn.aux_vars_face, :, j, i)

      nrm_xy = ro_sview(mesh.nrm_face, :, j, i)
      flux_j = sview(flux_face, :, j)

      functor(params, qL_j, qR_j, aux_vars, nrm_xy, flux_j)
    end  # end loop j

    resL = sview(eqn.res, :, :, iface_i.elementL)
    resR = sview(eqn.res, :, :, iface_i.elementR)
    interiorFaceIntegrate!(mesh.sbpface, iface_i, flux_face, resL, resR,
                           SummationByParts.Subtract())
  end  # end loop i

  return nothing
end

"""
  This function loops over interfaces and computes a face integral that
  uses data from all volume nodes. See [`FaceElementIntegralType`](@ref)
  for details on the integral performed.
"""
function getFaceElementIntegral{Tmsh, Tsol, Tres, Tdim}(
                           mesh::AbstractDGMesh{Tmsh},
                           sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim},
                           face_integral_functor::FaceElementIntegralType,
                           flux_functor::FluxType,
                           sbpface::AbstractFace,
                           interfaces::AbstractArray{Interface, 1})

  params = eqn.params
#  sbpface = mesh.sbpface
  nfaces = length(interfaces)
#  resL2 = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)
#  resR2 = zeros(resL2)
  for i=1:nfaces
    iface = interfaces[i]
    elL = iface.elementL
    elR = iface.elementR
    qL = ro_sview(eqn.q, :, :, elL)
    qR = ro_sview(eqn.q, :, :, elR)
    aux_vars = ro_sview(eqn.aux_vars, :, :, elL)
    nrm_face = ro_sview(mesh.nrm_face, :, :, i)
    resL = sview(eqn.res, :, :, elL)
    resR = sview(eqn.res, :, :, elR)

    face_integral_functor(params, sbpface, iface, qL, qR, aux_vars,
                       nrm_face, flux_functor, resL, resR)
  end

  return nothing
end

"""
  This function loops over interfaces and computes a face integral that
  uses data from all volume nodes. See FaceElementIntegralType for details on
  the integral performed. This function works on the staggered grid.

  Inputs:
    mesh_s: mesh object on solution mesh
    mesh_f: mesh object on flux mesh
    sbp_s: SBP operator on solution mesh
    sbp_f: SBP operator on flux mesh
    face_integral_functor: a [`FaceElementIntegralType`](@doc) functor
    flux_functor: a [`FluxType`] functor that is passed to face_integral_functor
                  to use as the numerical flux function
    sbpface: an AbstractFace for the flux grid
    interfaces: vector of [`Interfaces`](@doc) that the integral will be
                computed for.

  Inputs/Outputs:
    eqn: equation object (implicitly lives on solution grid).  eqn.res is
          updated with the results.
"""
function getFaceElementIntegral{Tmsh, Tsol, Tres, Tdim}(
                           mesh_s::AbstractDGMesh{Tmsh},
                           mesh_f::AbstractDGMesh{Tmsh},
                           sbp_s::AbstractSBP, sbp_f::AbstractSBP,
                           eqn::EulerData{Tsol, Tres, Tdim},
                           face_integral_functor::FaceElementIntegralType,
                           flux_functor::FluxType,
                           sbpface::AbstractFace,
                           interfaces::AbstractArray{Interface, 1})
  params = eqn.params
  qf = eqn.q_flux
  res = eqn.res
  nfaces = length(interfaces)
  resL_s = params.ress_el1
  resR_s = params.ress_el2
  resL_f = params.res_el1
  resR_f = params.res_el2

  aux_vars = zeros(Tres, 1, mesh_f.numNodesPerElement)

  for i=1:nfaces
    iface = interfaces[i]
    elL = iface.elementL
    elR = iface.elementR
    nrm_face = ro_sview(mesh_f.nrm_face, :, :, i)
    
    qf_L = ro_sview(qf, :, :, elL)
    qf_R = ro_sview(qf, :, :, elR)

    for j=1:mesh_f.numNodesPerElement
      aux_vars[1, j] = calcPressure(params, ro_sview(qf_L, :, j))
    end

    fill!(resL_f, 0.0)
    fill!(resR_f, 0.0)

    face_integral_functor(params, sbpface, iface, qf_L, qf_R, aux_vars,
                       nrm_face, flux_functor, resL_f, resR_f)

    # interpolate back
    smallmatmat!(resL_f, mesh_s.I_S2F, resL_s)
    smallmatmat!(resR_f, mesh_s.I_S2F, resR_s)

    # accumulate into res
    @simd for j=1:mesh_s.numNodesPerElement
      @simd for k=1:mesh_s.numDofPerNode
        res[k, j, elL] += resL_s[k, j]
        res[k, j, elR] += resR_s[k, j]
      end
    end

  end  # end loop i

  return nothing
end

                         
"""
  This function is a thin wrapper around
  [`calcSharedFaceElementIntegrals_element_inner`](@ref),
  presenting the interface needed by [`finishExchangeData`](@ref).
  See that function for the interface details.
"""
function calcSharedFaceElementIntegrals_element{Tmsh, Tsol, Tres}(
                            mesh::AbstractDGMesh{Tmsh},
                            sbp::AbstractSBP, eqn::EulerData{Tsol, Tres},
                            opts, data::SharedFaceData)

  if opts["use_staggered_grid"]
    calcSharedFaceElementIntegralsStaggered_element_inner(mesh, mesh.mesh2, 
                            sbp, mesh.sbp2, eqn, opts, data,
                            eqn.face_element_integral_func, eqn.flux_func)

  else
    calcSharedFaceElementIntegrals_element_inner(mesh, sbp, eqn, opts, data,
                            eqn.face_element_integral_func,  eqn.flux_func)
  end

  return nothing
end

"""
  This function loops over given set of shared faces and computes a face
  integral that
  uses data from all volume nodes.  See [`FaceElementIntegralType`](@ref)
  for details on the integral performed.

  **Inputs**:

   * mesh
   * sbp
   * eqn
   * opts
   * data: a SharedFaceData specifying which shared faces to compute
   * face_integral_functor
   * flux_functor

"""
function calcSharedFaceElementIntegrals_element_inner{Tmsh, Tsol, Tres}(
                            mesh::AbstractDGMesh{Tmsh},
                            sbp::AbstractSBP, eqn::EulerData{Tsol, Tres},
                            opts, data::SharedFaceData,
                            face_integral_functor::FaceElementIntegralType,
                            flux_functor::FluxType)

  if opts["parallel_data"] != "element"
    throw(ErrorException("cannot use calcSharedFaceIntegrals_element without parallel element data"))
  end

  q = eqn.q
  params = eqn.params
  # we don't care about elementR here, so use this throwaway array
  resR = Array(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)

  # get the data for the parallel interface
  idx = data.peeridx
  interfaces = data.interfaces
  bndries_local = data.bndries_local
  bndries_remote = data.bndries_remote
#    qL_arr = eqn.q_face_send[i]
  qR_arr = data.q_recv
  nrm_face_arr = mesh.nrm_sharedface[idx]

  start_elnum = mesh.shared_element_offsets[idx]

  # compute the integrals
  for j=1:length(interfaces)
    iface_j = interfaces[j]
    elL = iface_j.elementL
    elR = iface_j.elementR - start_elnum + 1  # is this always equal to j?
    qL = ro_sview(q, :, :, elL)
    qR = ro_sview(qR_arr, :, :, elR)
    aux_vars = ro_sview(eqn.aux_vars, :, :, elL)
    nrm_face = ro_sview(nrm_face_arr, :, :, j)
    resL = sview(eqn.res, :, :, elL)

    face_integral_functor(eqn.params, mesh.sbpface, iface_j, qL, qR, aux_vars,
                          nrm_face, flux_functor, resL, resR)
  end  # end loop j

  return nothing
end

"""
  Like [`calcSharedFaceElementIntegrals_element_inner`](@ref), but for staggered grid.

  data.q_recv is the solution grid data.  This function interpolates it to the
  flux grid on the fly.

  **Inputs**:

   * mesh_s: solution grid mesh
   * mesh_f: flux grid mesh
   * sbp_s: SBP operator for solution grid
   * sbp_f: SBP operator for the flux grid
   * opts: options dictionary
   * data: [`SharedFaceData`](@ref)
   * face_integral_functor: [`FaceElementIntegralType`](@ref)
   * flux_functor: [`FluxType`](@ref) passed to face_integral functor

  **Inputs/Outputs**:

    * eqn: equation object (lives on solution grid)

  Aliasing restrictions: none
"""
function calcSharedFaceElementIntegralsStaggered_element_inner{Tmsh, Tsol, Tres}(
                            mesh_s::AbstractDGMesh{Tmsh},
                            mesh_f::AbstractDGMesh{Tmsh},
                            sbp_s::AbstractSBP, sbp_f::AbstractSBP,
                            eqn::EulerData{Tsol, Tres},
                            opts, data::SharedFaceData,
                            face_integral_functor::FaceElementIntegralType,
                            flux_functor::FluxType)

  if opts["parallel_data"] != "element"
    throw(ErrorException("cannot use calcSharedFaceIntegrals_element without parallel element data"))
  end

  q = eqn.q_flux
  res = eqn.res
  params = eqn.params

  # temporary arrays needed for interpolating to the flux grid
  qvars_f = params.q_el2
  resL_f = params.res_el1
  resR_f = params.res_el2
  resL_s = params.ress_el1
  # we don't care about resR here
  aux_vars = Array(Tres, 1, mesh_f.numNodesPerElement)

  # get the data for the parallel interface
  idx = data.peeridx
  interfaces = data.interfaces
  qR_arr = data.q_recv
  nrm_face_arr = mesh_f.nrm_sharedface[idx]
  start_elnum = mesh_f.shared_element_offsets[idx]

  # compute the integrals
  for j=1:length(interfaces)
    iface_j = interfaces[j]
    elL = iface_j.elementL
    elR = iface_j.elementR - start_elnum + 1
    qL = ro_sview(q, :, :, elL)
    qR = ro_sview(qR_arr, :, :, elR)
#    aux_vars = ro_sview(eqn.aux_vars, :, :, elL)
    nrm_face = ro_sview(nrm_face_arr, :, :, j)

    # interpolate to flux grid
    interpolateElementStaggered(params, mesh_s, qR, aux_vars, qvars_f)

    fill!(resL_f, 0.0)
    # we dont carea bout resR, so dont waste time zeroing it out
    face_integral_functor(eqn.params, mesh_f.sbpface, iface_j, qL, qvars_f,
                          aux_vars, nrm_face, flux_functor, resL_f, resR_f)


    # interpolate residual back to solution grid
    smallmatmat!(resL_f, mesh_s.I_S2F, resL_s)

    @simd for k=1:mesh_s.numNodesPerElement
      @simd for p=1:mesh_s.numDofPerNode
        res[p, k, elL] += resL_s[p, k]
      end
    end

  end  # end loop j


  return nothing
end

"""
  This function is a thin wrapper around calcSharedFaceIntegrals_inner.
  It present the interface needed by [`finishExchangeData`](@ref).
"""
function calcSharedFaceIntegrals{Tmsh, Tsol}( mesh::AbstractDGMesh{Tmsh},
                            sbp::AbstractSBP, eqn::EulerData{Tsol},
                            opts, data::SharedFaceData)

  if opts["precompute_face_flux"]
    calcSharedFaceIntegrals_inner(mesh, sbp, eqn, opts, data, eqn.flux_func)
  else
    calcSharedFaceIntegrals_nopre_inner(mesh, sbp, eqn, opts, data, eqn.flux_func)
  end

  return nothing
end


@doc """
### EulerEquationMod.calcSharedFaceIntegrals

  This function calculates the shared face integrals over a given set of
  faces.

  **Inputs**:

   * mesh
   * sbp
   * eqn
   * opts
   * data: the SharedFaceData specifying which faces to compute
   * functor: the FluxType to use for the face flux

"""->
function calcSharedFaceIntegrals_inner{Tmsh, Tsol}( mesh::AbstractDGMesh{Tmsh},
                            sbp::AbstractSBP, eqn::EulerData{Tsol},
                            opts, data::SharedFaceData, functor::FluxType)
# calculate the face flux and do the integration for the shared interfaces

  if opts["parallel_data"] != "face"
    throw(ErrorException("cannot use calcSharedFaceIntegrals without parallel face data"))
  end


  params = eqn.params

  # calculate the flux
  idx = data.peeridx
  interfaces = data.interfaces
  bndries_local = data.bndries_local
  qL_arr = data.q_send
  qR_arr = data.q_recv
  aux_vars_arr = eqn.aux_vars_sharedface[idx]
  nrm_arr = mesh.nrm_sharedface[idx]
  flux_arr = eqn.flux_sharedface[idx]

  for j=1:length(interfaces)
    interface_i = interfaces[j]
    for k=1:mesh.numNodesPerFace
      eL = interface_i.elementL
      fL = interface_i.faceL

      qL = ro_sview(qL_arr, :, k, j)
      qR = ro_sview(qR_arr, :, k, j)
      aux_vars = ro_sview(aux_vars_arr, :, k, j)
      nrm_xy = ro_sview(nrm_arr, :, k, j)
      flux_j = sview(flux_arr, :, k, j)
      functor(params, qL, qR, aux_vars, nrm_xy, flux_j)
    end
  end
  # end flux calculation

  # do the integration
  boundaryintegrate!(mesh.sbpface, bndries_local, flux_arr, eqn.res, SummationByParts.Subtract())

  @debug2 sharedFaceLogging(mesh, sbp, eqn, opts, data, qL_arr, qR_arr)

  return nothing
end

"""
  Like [`calcSharedFaceIntegrals_inner`](@ref), but performs the integration and
  updates eqn.res rather than computing the flux and storing it in
  eqn.flux_sharedface
"""
function calcSharedFaceIntegrals_nopre_inner{Tmsh, Tsol, Tres}(
                            mesh::AbstractDGMesh{Tmsh},
                            sbp::AbstractSBP, eqn::EulerData{Tsol, Tres},
                            opts, data::SharedFaceData, functor::FluxType)

  if opts["parallel_data"] != "face"
    throw(ErrorException("cannot use calcSharedFaceIntegrals without parallel face data"))
  end


  params = eqn.params

  # calculate the flux
  idx = data.peeridx
  interfaces = data.interfaces
  bndries_local = data.bndries_local
  qL_arr = data.q_send
  qR_arr = data.q_recv
  aux_vars_arr = eqn.aux_vars_sharedface[idx]
  nrm_arr = mesh.nrm_sharedface[idx]
  flux_j = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)

  for j=1:length(interfaces)
    interface_i = interfaces[j]
    eL = interface_i.elementL
    fL = interface_i.faceL

    # compute the flux
    for k=1:mesh.numNodesPerFace
      qL = ro_sview(qL_arr, :, k, j)
      qR = ro_sview(qR_arr, :, k, j)
      aux_vars = ro_sview(aux_vars_arr, :, k, j)
      parent(aux_vars)[1] = calcPressure(params, qL)
      nrm_xy = ro_sview(nrm_arr, :, k, j)
      flux_k = sview(flux_j, :, k)
      functor(params, qL, qR, aux_vars, nrm_xy, flux_k)
    end
    
    # do the integration
    res_j = sview(eqn.res, :, :, eL)
    boundaryFaceIntegrate!(mesh.sbpface, fL, flux_j, res_j, SummationByParts.Subtract())

  end

  return nothing
end
 
"""
  This function is a thin wrapper around
  [`calcSharedFaceIntegrals_element_inner`](@ref).
  It presents the interface required by [`finishExchangeData`](@ref)
"""
function calcSharedFaceIntegrals_element{Tmsh, Tsol}(mesh::AbstractDGMesh{Tmsh},
                            sbp::AbstractSBP, eqn::EulerData{Tsol},
                            opts, data::SharedFaceData)

  if opts["precompute_face_flux"]
    calcSharedFaceIntegrals_element_inner(mesh, sbp, eqn, opts, data, eqn.flux_func)
  else
    calcSharedFaceIntegrals_nopre_element_inner(mesh, sbp, eqn, opts, data, eqn.flux_func)
  end

  return nothing
end



# element parallel version
"""
  This function calculates the shared face integrals for a given set of faces.
  It uses the MPI send and receive buffers that contain the solution for the
  elements on the boundary (rather than the data on the faces).  This
  enables calculating a sparse jacobian with minimal parallel communication.

  **Inputs**:

   * mesh
   * sbp
   * eqn
   * opts
   * data: a SharedFaceData specifying the faces to calculate
   * functor: the flux functor

"""
function calcSharedFaceIntegrals_element_inner{Tmsh, Tsol}(
                            mesh::AbstractDGMesh{Tmsh},
                            sbp::AbstractSBP, eqn::EulerData{Tsol},
                            opts, data::SharedFaceData, functor::FluxType)

  q = eqn.q
  params = eqn.params

  @debug2 begin
    qL_face_arr[i] = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace,
                                     mesh.peer_face_counts[i])
    qR_face_arr[i] = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace,
                                     mesh.peer_face_counts[i])
  end 

  # TODO: make these fields of params
  q_faceL = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  q_faceR = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  workarr = zeros(q_faceR)


  idx = data.peeridx
  interfaces = data.interfaces
  bndries_local = data.bndries_local
  bndries_remote = data.bndries_remote
#    qL_arr = eqn.q_face_send[i]
  qR_arr = data.q_recv
  nrm_arr = mesh.nrm_sharedface[idx]
  aux_vars_arr = eqn.aux_vars_sharedface[idx]
  flux_arr = eqn.flux_sharedface[idx]

  start_elnum = mesh.shared_element_offsets[idx]

  @debug2 flush(params.f)
  for j=1:length(interfaces)
    iface_j = interfaces[j]
    bndryL_j = bndries_local[j]
    bndryR_j = bndries_remote[j]
    fL = bndryL_j.face

    # interpolate to face
    qL = ro_sview(q, :, :, iface_j.elementL)
    el_r = iface_j.elementR - start_elnum + 1
    qR = ro_sview(qR_arr, :, :, el_r)
    interiorFaceInterpolate!(mesh.sbpface, iface_j, qL, qR, q_faceL, q_faceR)

    @debug2 qL_face_arr[:, :, j] = q_faceL
    @debug2 qR_face_arr[:, :, j] = q_faceR

    # calculate flux
    for k=1:mesh.numNodesPerFace
      qL_k = ro_sview(q_faceL, :, k)
      qR_k = ro_sview(q_faceR, :, k)
      aux_vars = ro_sview(aux_vars_arr, :, k, j)
      nrm_xy = ro_sview(nrm_arr, :, k, j)
      flux_k = sview(flux_arr, :, k, j)

      parent(aux_vars)[1] = calcPressure(params, qL_k)

      functor(params, qL_k, qR_k, aux_vars, nrm_xy, flux_k)
     end
   end  # end loop over interfaces

  # evaluate integral
  boundaryintegrate!(mesh.sbpface, bndries_local, flux_arr, eqn.res, SummationByParts.Subtract())

  @debug2 sharedFaceLogging(mesh, sbp, eqn, opts, data, qL_face_arr, qR_face_arr)

  return nothing
end

"""
  Like [`calcSharedFaceIntegrals_element_inner`](@ref), but performs the integration and
  updates eqn.res rather than computing the flux only and storing it in
  eqn.flux_sharedface
"""
function calcSharedFaceIntegrals_nopre_element_inner{Tmsh, Tsol, Tres}(
                            mesh::AbstractDGMesh{Tmsh},
                            sbp::AbstractSBP, eqn::EulerData{Tsol, Tres},
                            opts, data::SharedFaceData, functor::FluxType)

  q = eqn.q
  params = eqn.params

  @debug2 begin
    qL_face_arr[i] = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace,
                                     mesh.peer_face_counts[i])
    qR_face_arr[i] = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace,
                                     mesh.peer_face_counts[i])
  end 

  # TODO: make these fields of params
  q_faceL = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  q_faceR = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  workarr = zeros(q_faceR)


  idx = data.peeridx
  interfaces = data.interfaces
  bndries_local = data.bndries_local
  bndries_remote = data.bndries_remote
#    qL_arr = eqn.q_face_send[i]
  qR_arr = data.q_recv
  nrm_arr = mesh.nrm_sharedface[idx]
  aux_vars_arr = eqn.aux_vars_sharedface[idx]
  flux_face = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)

  start_elnum = mesh.shared_element_offsets[idx]

  @debug2 flush(params.f)
  for j=1:length(interfaces)
    iface_j = interfaces[j]
    bndryL_j = bndries_local[j]
    bndryR_j = bndries_remote[j]
    fL = bndryL_j.face

    # interpolate to face
    qL = ro_sview(q, :, :, iface_j.elementL)
    el_r = iface_j.elementR - start_elnum + 1
    qR = ro_sview(qR_arr, :, :, el_r)
    interiorFaceInterpolate!(mesh.sbpface, iface_j, qL, qR, q_faceL, q_faceR)

    @debug2 qL_face_arr[:, :, j] = q_faceL
    @debug2 qR_face_arr[:, :, j] = q_faceR

    # calculate flux
    for k=1:mesh.numNodesPerFace
      qL_k = ro_sview(q_faceL, :, k)
      qR_k = ro_sview(q_faceR, :, k)
      aux_vars = ro_sview(aux_vars_arr, :, k, j)
      nrm_xy = ro_sview(nrm_arr, :, k, j)
      flux_k = sview(flux_face, :, k)

      parent(aux_vars)[1] = calcPressure(params, qL_k)

      functor(params, qL_k, qR_k, aux_vars, nrm_xy, flux_k)
     end

     # do the integration
     res_j = sview(eqn.res, :, :, bndryL_j.element)
     boundaryFaceIntegrate!(mesh.sbpface, fL, flux_face, res_j, SummationByParts.Subtract())
   end  # end loop over interfaces

  @debug2 sharedFaceLogging(mesh, sbp, eqn, opts, data, qL_face_arr, qR_face_arr)

  return nothing
end




@doc """
### EulerEquationMod.interpolateFace

  This function interpolates the solution values from the internal nodes
  to the face flux points of the elements

  **Inputs**:

   * mesh: an AbstractDGMesh
   * sbp
   * eqn
   * opts
   * q: a 3D array of solution values at the nodes, numDofPerNode x
       numNodesPerElement x numEl

  **Inputs/Outputs**:

   * q_face: a 4D array of solution values at each interface,
            numDofPerNode x 2 x numfacenodes x numInterface
            q_face[:, 1, j, i] stores the q values for elementL of interface
            i node j and q_face[:, 2, j, i] stores the values for elementR

  eqn.aux_vars_face is also populated
"""->
function interpolateFace{Tsol}(mesh::AbstractDGMesh, sbp, eqn, opts,
                         q::Abstract3DArray, q_face::AbstractArray{Tsol, 4})

  # interpolate solution
  interiorfaceinterpolate!(mesh.sbpface, mesh.interfaces, q, q_face)

  # recalculte aux_vars
  for i=1:mesh.numInterfaces
    for j=1:mesh.numNodesPerFace
      q_vals = ro_sview(q_face, :, 1, j, i) # always use elementL
      eqn.aux_vars_face[1, j, i] = calcPressure(eqn.params, q_vals)
    end
  end

  if opts["parallel_data"] == 1
    for peer=1:mesh.npeers
      q_vals_p = eqn.q_face_send[peer]
      aux_vars = eqn.aux_vars_sharedface[peer]
      for i=1:mesh.peer_face_counts[peer]
        for j=1:mesh.numNodesPerFace
          q_vals = ro_sview(q_vals_p, :, j, i)
          aux_vars[1, j, i] = calcPressure(eqn.params, q_vals)
        end
      end
    end
  end  # else populate aux_vars_sharedface in evalSharedFaceIntegrals

  return nothing
end




#------------------------------------------------------------------------------
# Different fluxes
#------------------------------------------------------------------------------

"""
  This flux function sets F = q.  Useful for testing
"""
type IdentityFlux <: FluxType
end

function call{Tsol, Tres, Tmsh}(obj::IdentityFlux, params::ParamType,
              uL::AbstractArray{Tsol,1},
              uR::AbstractArray{Tsol,1},
              aux_vars::AbstractVector{Tres},
              nrm::AbstractVector{Tmsh},
              F::AbstractVector{Tres})

  for i=1:length(F)
    F[i] = q[i]
  end

  return nothing
end

function call{Tsol, Tres, Tmsh}(obj::IdentityFlux, params::ParamType,
              uL::AbstractArray{Tsol,1},
              uR::AbstractArray{Tsol,1},
              aux_vars::AbstractVector{Tres},
              nrm::AbstractMatrix{Tmsh},
              F::AbstractMatrix{Tres})

  for j=1:size(F, 2)
    for i=1:length(F)
      F[i, j] = q[i]
    end
  end

  return nothing
end



"""
  Calls the [`RoeSolver`](@ref)
"""

type RoeFlux <: FluxType
end

function call{Tsol, Tres, Tmsh}(obj::RoeFlux, params::ParamType,
              uL::AbstractArray{Tsol,1},
              uR::AbstractArray{Tsol,1},
              aux_vars::AbstractVector{Tres},
              nrm::AbstractVector{Tmsh},
              F::AbstractVector{Tres})

  RoeSolver(params, uL, uR, aux_vars, nrm, F)
  return nothing
end

"""
  Calls [`calcEulerFlux_standard`](@ref)
"""
type LFFlux <: FluxType
end

function call{Tsol, Tres, Tmsh}(obj::LFFlux, params::ParamType,
              uL::AbstractArray{Tsol,1},
              uR::AbstractArray{Tsol,1},
              aux_vars::AbstractVector{Tres},
              nrm::AbstractVector{Tmsh},
              F::AbstractVector{Tres})

  calcLFFlux(params, uL, uR, aux_vars, nrm, F)
  return nothing
end


type StandardFlux <: FluxType
end

function call{Tsol, Tres}(obj::StandardFlux, params::ParamType,
              uL::AbstractArray{Tsol,1},
              uR::AbstractArray{Tsol,1},
              aux_vars::AbstractVector{Tres},
              nrm::AbstractArray, 
              F::AbstractArray{Tres})

  calcEulerFlux_standard(params, uL, uR, aux_vars, nrm, F)
  return nothing
end

"""
  Calls [`calcEulerFlux_Ducros`](@ref)
"""
type DucrosFlux <: FluxType
end

function call{Tsol, Tres}(obj::DucrosFlux, params::ParamType,
              uL::AbstractArray{Tsol,1},
              uR::AbstractArray{Tsol,1},
              aux_vars::AbstractVector{Tres},
              nrm::AbstractVector,
              F::AbstractVector{Tres})

  calcEulerFlux_Ducros(params, uL, uR, aux_vars, nrm, F)
  return nothing
end

"""
  Calls [`calcEulerFlux_IR`](@ref)
"""
type IRFlux <: FluxType
end

function call{Tsol, Tres}(obj::IRFlux, params::ParamType,
              uL::AbstractArray{Tsol,1},
              uR::AbstractArray{Tsol,1},
              aux_vars::AbstractVector{Tres},
              nrm::AbstractArray, 
              F::AbstractArray{Tres})

  # this will dispatch to either the sinlge director or multi-dimension method
  calcEulerFlux_IR(params, uL, uR, aux_vars, nrm, F)
  return nothing
end

"""
  Calls [`calcEulerFlux_IRSLF`](@ref)
"""
type IRSLFFlux <: FluxType
end

function call{Tsol, Tres}(obj::IRSLFFlux, params::ParamType,
              uL::AbstractArray{Tsol,1},
              uR::AbstractArray{Tsol,1},
              aux_vars::AbstractVector{Tres},
              nrm::AbstractVector,
              F::AbstractVector{Tres})

  calcEulerFlux_IRSLF(params, uL, uR, aux_vars, nrm, F)
  return nothing
end




@doc """
### EulerEquationMod.FluxDict

  This dictonary maps the names of the fluxes (ASCIIStrings) to the
  functor object itself.  All flux functors should be added to the dictionary.

  All fluxes have one method that calculates the flux in a particular direction
  at a node.  Some fluxes have an additional method that computes the flux
  in several directions at a node in a single function call, which can be
  more efficient.  See calcEulerFlux_standard for an example.

  In general, these functors call similarly-named function in bc_solvers.jl.
  It is recommened to look at the documentation for those functions.

  TODO: document signature of the functors here

"""->
global const FluxDict = Dict{ASCIIString, FluxType}(
"IdentityFlux" => IdentityFlux(),
"RoeFlux" => RoeFlux(),
"LFFlux" => LFFlux(),
"StandardFlux" => StandardFlux(),
"DucrosFlux" => DucrosFlux(),
"IRFlux" => IRFlux(),
"IRSLFFlux" => IRSLFFlux()
)

@doc """
### EulerEquationMod.getFluxFunctors

  This function retrieves the flux functors from the dictonary and
  stores them to eqn.flux_func.

  Inputs:
    mesh: an AbstractDGMesh
    sbp
    eqn
    opts
"""->
function getFluxFunctors(mesh::AbstractDGMesh, sbp, eqn, opts)

  name = opts["Flux_name"]
  eqn.flux_func = FluxDict[name]
  name = opts["Volume_flux_name"]
  eqn.volume_flux_func = FluxDict[name]
  return nothing
end

type RoeFlux_revm <: FluxType_revm
end

function call{Tsol, Tres, Tmsh}(obj::RoeFlux_revm, params::ParamType,
              uL::AbstractArray{Tsol,1}, uR::AbstractArray{Tsol, 1}, aux_vars,
              nrm::AbstractVector{Tmsh},
              flux_bar::AbstractVector{Tres}, nrm_bar::AbstractArray{Tmsh, 1})

  RoeSolver_revm(params, uL, uR, aux_vars, nrm, flux_bar, nrm_bar)

  return nothing
end

global const FluxDict_revm = Dict{ASCIIString, FluxType_revm}(
"RoeFlux" => RoeFlux_revm(),
)

function getFluxFunctors_revm(mesh::AbstractDGMesh, sbp, eqn, opts)

  name = opts["Flux_name"]
  eqn.flux_func_bar = FluxDict_revm[name]

  return nothing
end # End function getFluxFunctors_revm
