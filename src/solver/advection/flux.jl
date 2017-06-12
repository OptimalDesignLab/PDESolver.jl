# this file contains the defitions of all the fluxes used for DG face integrals
@doc """
### AdvectionEquationMod.calcFaceFlux

  This function calculates the DG flux between a specified set of faces,
  using the solution data at the faces stored in eqn.q_face.
  Note that the flux is negated because the face integrals have a 
  negative sign in the weak form.

  Inputs:
    mesh
    sbp
    eqn
    functor: the functor that calculates the flux at a node
    interfaces: an array of type Interface that specifies which interfaces
                to calculate the flux for

  Inputs/Outputs:
    face_flux: array to store the flux in, numDofPerNode x nnodesPerFace
               x length(interfaces)

  The functor must have the signature:
  func( uL, qR, alpha_x, alpha_y, dxidx, nrm, params)

  where uL and uR are the solution values for a node on the left and right
  elements, alpha_x and alpha_y are the x and y advection velocities,
  dxidx is the scaled mapping jacobian for elementL, and nrm is the face
  normal in reference space.  params is eqn.params

"""->
function calcFaceFlux{Tmsh,  Tsol, Tres}( mesh::AbstractDGMesh{Tmsh}, 
                          sbp::AbstractSBP, eqn::AdvectionData{Tsol}, 
                          functor::FluxType, 
                          interfaces::AbstractArray{Interface,1}, 
                          face_flux::AbstractArray{Tres, 3})
 
  nfaces = length(interfaces)
  for i=1:nfaces  # loop over faces
    interface_i = interfaces[i]
    for j = 1:mesh.numNodesPerFace
      eL = interface_i.elementL
      fL = interface_i.faceL

      # get components
      qL = eqn.q_face[1, 1, j, i]
      qR = eqn.q_face[1, 2, j, i]
      nrm_scaled = sview(mesh.nrm_face, :, j, i)
#      dxidx = sview(mesh.dxidx_face, :, :, j, i)
#      nrm = sview(sbp.facenormal, :, fL)

      face_flux[1, j, i] = -functor(qL, qR, nrm_scaled, eqn.params)
    end
  end

  return nothing
end

"""
  Compute the face integrals without using eqn.q_face or eqn.flux_face.
  The integral is computed directly and res is updated

  Inputs:
    mesh
    sbp
    eqn
    opts
    flux_func: the flux functor that computes the face flux at a node
"""
function calcFaceIntegrals_nopre{Tsol, Tres, Tmsh, Tdim}(
                                 mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP,
                                 eqn::AdvectionData{Tsol, Tres, Tdim}, opts,
                                 flux_func::FluxType)

  q_faceL = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  q_faceR = zeros(q_faceL)
  flux_face = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)

  for i=1:mesh.numInterfaces
    iface_i = mesh.interfaces[i]
    qL = sview(eqn.q, :, :, iface_i.elementL)
    qR = sview(eqn.q, :, :, iface_i.elementR)

    # interpolate to face
    interiorFaceInterpolate!(mesh.sbpface, iface_i, qL, qR, q_faceL, q_faceR)

    for j=1:mesh.numNodesPerFace
      dxidx_j = sview(mesh.dxidx_face, :, :, j, i)
      nrm = sview(sbp.facenormal, :, iface_i.faceL)

      flux_face[1, j] = -flux_func(q_faceL[j], q_faceR[j], dxidx_j, nrm, eqn.params)
    end

    resL = sview(eqn.res, :, :, iface_i.elementL)
    resR = sview(eqn.res, :, :, iface_i.elementR)
    interiorFaceIntegrate!(mesh.sbpface, iface_i, flux_face, resL, resR)
  end

  return nothing
end


"""
  Thin wrapper around calcSharedFaceIntegrals_inner.  This function is passed
  to finishDataExchange, and internally calls calcSharedFaceIntegrals_inner.
  See finishDataExchange for details on the interface and 
  calcSharedFaceIntegrals_inner for the integral that is computed.
"""
function calcSharedFaceIntegrals{Tmsh, Tsol}( mesh::AbstractDGMesh{Tmsh},
                            sbp::AbstractSBP, eqn::AdvectionData{Tsol},
                            opts, data::SharedFaceData)

  if opts["precompute_face_flux"]
    calcSharedFaceIntegrals_inner(mesh, sbp, eqn, opts, data, eqn.flux_func)
  else
    calcSharedFaceIntegrals_inner_nopre(mesh, sbp, eqn, opts, data,
                                        eqn.flux_func)
  end

  return nothing
end

@doc """
### AdvectionEquationMod.calcSharedFaceIntegrals

  This function calculates the shared face integrals for the faces shared
  with a single peer process.  This function is for
  opts["parallel_type"] == "face" and regular face integrals (ie. not the
  entropy-stable face integrals) only.

  Inputs:
    mesh
    sbp
    eqn
    opts:
    data: a SharedFaceData specifying which faces to compute
    functor: the FluxType to use for the face flux

"""->
function calcSharedFaceIntegrals_inner{Tmsh, Tsol}( mesh::AbstractDGMesh{Tmsh},
                            sbp::AbstractSBP, eqn::AdvectionData{Tsol},
                            opts, data::SharedFaceData, functor::FluxType)
# calculate the face flux and do the integration for the shared interfaces

  if opts["parallel_data"] != "face"
    throw(ErrorException("cannot use calcSharedFaceIntegrals without parallel face data"))
  end

  params = eqn.params
  idx = data.peeridx

    # calculate the flux
  interfaces = data.interfaces
  bndries_local = data.bndries_local
  qL_arr = data.q_send
  qR_arr = data.q_recv
  nrm_arr = mesh.nrm_sharedface[idx]
#  dxidx_arr = mesh.dxidx_sharedface[idx]
  flux_arr = eqn.flux_sharedface[idx]

  for j=1:length(interfaces)
    interface_i = interfaces[j]
    for k=1:mesh.numNodesPerFace
      eL = interface_i.elementL
      fL = interface_i.faceL

      qL = qL_arr[1, k, j]
      qR = qR_arr[1, k, j]
#      dxidx = sview(dxidx_arr, :, :, k, j)
#      nrm = sview(sbp.facenormal, :, fL)
      nrm_scaled = sview(nrm_arr, :, k, j)
      flux_arr[1,k,j] = -functor(qL, qR, nrm_scaled, eqn.params)
    end
  end
  # end flux calculation

  # do the integration
  boundaryintegrate!(mesh.sbpface, bndries_local, flux_arr, eqn.res)

  #TODO: fix this
  @debug1 sharedFaceLogging(mesh, sbp, eqn, opts, data,  qL_arr, qr_arr)

  return nothing
end

"""
  Like calcSharedFaceIntegrals_inner_nopre, but it computes the integral one
  face at a time rather than computing all the flux, storing it in
  eqn.flux_sharedface and then doing the integral

  See calcSharedFaceIntegrals_inner() for a description of the arguments
"""
function calcSharedFaceIntegrals_inner_nopre{Tmsh, Tsol}(
                            mesh::AbstractDGMesh{Tmsh},
                            sbp::AbstractSBP, eqn::AdvectionData{Tsol},
                            opts, data::SharedFaceData, functor::FluxType)
# calculate the face flux and do the integration for the shared interfaces

  if opts["parallel_data"] != "face"
    throw(ErrorException("cannot use calcSharedFaceIntegrals without parallel face data"))
  end

  params = eqn.params
  idx = data.peeridx

    # calculate the flux
  interfaces = data.interfaces
  bndries_local = data.bndries_local
  qL_arr = data.q_send
  qR_arr = data.q_recv
  dxidx_arr = mesh.dxidx_sharedface[idx]
#  flux_arr = eqn.flux_sharedface[idx]
  flux_face = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)

  for j=1:length(interfaces)
    interface_i = interfaces[j]
    for k=1:mesh.numNodesPerFace
      eL = interface_i.elementL
      fL = interface_i.faceL

      qL = qL_arr[1, k, j]
      qR = qR_arr[1, k, j]
      dxidx = sview(dxidx_arr, :, :, k, j)
      nrm = sview(sbp.facenormal, :, fL)
      flux_face[1, k] = -functor(qL, qR, dxidx, nrm, eqn.params)
    end

    res_i = sview(eqn.res, :, :, interface_i.elementL)
#    println("size(flux_face) = ", size(flux_face))
#    println("size(mesh.sbpface.interp) = ", size(mesh.sbpface.interp))
    boundaryFaceIntegrate!(mesh.sbpface, interface_i.faceL, flux_face, res_i)
  end
  
  @debug1 sharedFaceLogging(mesh, sbp, eqn, opts, data,  qL_arr, qr_arr)

  return nothing
end


"""
  Thin wrapper around calcSharedFaceIntegrals_inner.  This function is passed
  to finishDataExchange, and internally calls calcSharedFaceIntegrals_inner.
  See finishDataExchange for details on the interface and 
  calcSharedFaceIntegrals_inner for the integral that is computed.
"""
function calcSharedFaceIntegrals_element{Tmsh, Tsol}(mesh::AbstractDGMesh{Tmsh},
                            sbp::AbstractSBP, eqn::AdvectionData{Tsol},
                            opts, data::SharedFaceData)

  if opts["precompute_face_flux"]
    calcSharedFaceIntegrals_element_inner(mesh, sbp, eqn, opts, data, eqn.flux_func)
  else
    calcSharedFaceIntegrals_element_inner_nopre(mesh, sbp, eqn, opts, data, eqn.flux_func)
  end

  return nothing
end


# element parallel version
function calcSharedFaceIntegrals_element_inner{Tmsh, Tsol}( 
                            mesh::AbstractDGMesh{Tmsh},
                            sbp::AbstractSBP, eqn::AdvectionData{Tsol},
                            opts, data::SharedFaceData, functor::FluxType)

  @debug1 println(eqn.params.f, "entered calcSharedFaceIntegrals_element")
  @debug1 flush(eqn.params.f)

  if opts["parallel_data"] != "element"
    throw(ErrorException("cannot use getSharedFaceIntegrals_elemenet without parallel element data"))
  end



  # TODO: make these fields of params
  q_faceL = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  q_faceR = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  workarr = zeros(q_faceR)
  q = eqn.q
  params = eqn.params


  idx = data.peeridx
  interfaces = data.interfaces
  bndries_local = data.bndries_local
  bndries_remote = data.bndries_remote
#    qL_arr = data.q_send
  qR_arr = data.q_recv
#  dxidx_arr = mesh.dxidx_sharedface[idx]
  nrm_arr = mesh.nrm_sharedface[idx]
  flux_arr = eqn.flux_sharedface[idx]

  start_elnum = mesh.shared_element_offsets[idx]

  @debug1 begin
    qL_face_arr[i] = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace, 
                                 mesh.peer_face_counts[i])
    qR_face_arr[i] = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace, 
                                 mesh.peer_face_counts[i])
    flush(params.f)
  end

  for j=1:length(interfaces)
    iface_j = interfaces[j]
    bndryL_j = bndries_local[j]
    bndryR_j = bndries_remote[j]
    fL = bndryL_j.face

    # interpolate to face
    qL = sview(q, :, :, iface_j.elementL)
    el_r = iface_j.elementR - start_elnum + 1
    qR = sview(qR_arr, :, :, el_r)

    interiorFaceInterpolate!(mesh.sbpface, iface_j, qL, qR, q_faceL, q_faceR)
    
    @debug1 qL_face_arr[:, :, j] = q_faceL
    @debug1 qR_face_arr[:, :, j] = q_faceR

    # calculate flux
    for k=1:mesh.numNodesPerFace
      qL_k = q_faceL[k]
      qR_k = q_faceR[k]
      nrm_scaled = sview(nrm_arr, :, k, j)
#      dxidx = sview(dxidx_arr, :, :, k, j)
#      nrm = sview(sbp.facenormal, :, fL)

      flux_tmp = -functor(qL_k, qR_k, nrm_scaled, eqn.params)
      flux_arr[1,k,j] = flux_tmp
    end
  end  # end loop over interfaces

  # evaluate integral
  boundaryintegrate!(mesh.sbpface, bndries_local, flux_arr, eqn.res)

  @debug1 sharedFaceLogging(mesh, sbp, eqn, opts, data, qL_face_arr, qR_face_arr)

  return nothing
end

# element parallel version
"""
  Like calcSharedFaceIntegarl_element_inner, but computes the integral one
  face at a time instead of computing the entire flux and then integrating.

  See calcSharedFaceIntegrals_element_inner for the meaning of the arguments
"""
function calcSharedFaceIntegrals_element_inner_nopre{Tmsh, Tsol, Tres}( 
                            mesh::AbstractDGMesh{Tmsh},
                            sbp::AbstractSBP, eqn::AdvectionData{Tsol, Tres},
                            opts, data::SharedFaceData, functor::FluxType)

  @debug1 println(eqn.params.f, "entered calcSharedFaceIntegrals_element")
  @debug1 flush(eqn.params.f)

  if opts["parallel_data"] != "element"
    throw(ErrorException("cannot use getSharedFaceIntegrals_elemenet without parallel element data"))
  end

  # TODO: make these fields of params
  q_faceL = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  q_faceR = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  workarr = zeros(q_faceR)
  q = eqn.q
  params = eqn.params


  idx = data.peeridx
  interfaces = data.interfaces
  bndries_local = data.bndries_local
  bndries_remote = data.bndries_remote
#    qL_arr = data.q_send
  qR_arr = data.q_recv
  dxidx_arr = mesh.dxidx_sharedface[idx]
#  flux_arr = eqn.flux_sharedface[idx]
  flux_face = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)

  start_elnum = mesh.shared_element_offsets[idx]

  @debug1 begin
    qL_face_arr[i] = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace, 
                                 mesh.peer_face_counts[i])
    qR_face_arr[i] = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace, 
                                 mesh.peer_face_counts[i])
    flush(params.f)
  end

  for j=1:length(interfaces)
    iface_j = interfaces[j]
    bndryL_j = bndries_local[j]
    bndryR_j = bndries_remote[j]
    fL = bndryL_j.face

    # interpolate to face
    qL = sview(q, :, :, iface_j.elementL)
    el_r = iface_j.elementR - start_elnum + 1
    qR = sview(qR_arr, :, :, el_r)


    interiorFaceInterpolate!(mesh.sbpface, iface_j, qL, qR, q_faceL, q_faceR)

    @debug1 qL_face_arr[:, :, j] = q_faceL
    @debug1 qR_face_arr[:, :, j] = q_faceR

    # calculate flux
    for k=1:mesh.numNodesPerFace
      qL_k = q_faceL[k]
      qR_k = q_faceR[k]
      dxidx = sview(dxidx_arr, :, :, k, j)
      nrm = sview(sbp.facenormal, :, fL)

      flux_face[1, k] = -functor(qL_k, qR_k, dxidx, nrm, eqn.params)
    end

    res_j = sview(eqn.res, :, :, bndryL_j.element)
    boundaryFaceIntegrate!(mesh.sbpface, bndryL_j.face, flux_face, res_j)
  end  # end loop over interfaces

  @debug1 sharedFaceLogging(mesh, sbp, eqn, opts, data, qL_face_arr, qR_face_arr)

  return nothing
end



@doc """
### AdvectionEquationMod.avgFlux

  This flux function averages the two states to calculate the flux

  Inputs:
    uL: the left state
    uR: the right state
    nrm: the scaled normal vector for elementL in x-y space

  Outputs:
    the flux
"""->
type avgFlux <: FluxType
end

function call{Tmsh, Tsol}(obj::avgFlux, uL::Tsol, uR::Tsol,
              nrm::AbstractArray{Tmsh,1}, params::ParamType2)

  alpha_x = params.alpha_x
  alpha_y = params.alpha_y
  #=
  alpha_xi = dxidx[1,1]*alpha_x + dxidx[1,2]*alpha_y
  alpha_eta = dxidx[2,1]*alpha_x + dxidx[2,2]*alpha_y
  alpha_n  = alpha_xi*nrm[1] + alpha_eta*nrm[2]
  =#
  alpha_n = alpha_x*nrm[1] + alpha_y*nrm[2]
  u = alpha_n*(uL + uR)*0.5
  return u
end

function call{Tmsh, Tsol}(obj::avgFlux, uL::Tsol, uR::Tsol,
              nrm::AbstractArray{Tmsh,1}, params::ParamType3)

  alpha_x = params.alpha_x
  alpha_y = params.alpha_y
  alpha_z = params.alpha_z
  #=
  alpha_xi = dxidx[1,1]*alpha_x + dxidx[1,2]*alpha_y + dxidx[1,3]*alpha_z
  alpha_eta = dxidx[2,1]*alpha_x + dxidx[2,2]*alpha_y + dxidx[2,3]*alpha_z
  alpha_psi = dxidx[3,1]*alpha_x + dxidx[3,2]*alpha_y + dxidx[3,3]*alpha_z
  alpha_n  = alpha_xi*nrm[1] + alpha_eta*nrm[2] + alpha_psi*nrm[3]
  =#
  alpha_n = alpha_x*nrm[1] + alpha_y*nrm[2] + alpha_z*nrm[3]
  u = alpha_n*(uL + uR)*0.5
  return u
end



@doc """
### AdvectionEquationMod.LFFlux

  The Lax-Friedrich flux, using a parameter alpha to control upwinding.
  alpha = 0 -> centered flux
  alpha = 1 -> completely upwinded
  The alpha parameter 

  Inputs:
    uL: the left state
    uR: the right state
    nrm_scaled: the scaled normal vector for elementL in x-y space

  Outputs:
    the flux
"""->
type LFFlux <: FluxType
end

function call{Tmsh, Tsol}(obj::LFFlux, uL::Tsol, uR::Tsol,
              nrm_scaled::AbstractArray{Tmsh,1}, params::ParamType2)

  alpha_x = params.alpha_x
  alpha_y = params.alpha_y
  #=
  alpha_xi = dxidx[1,1]*alpha_x + dxidx[1,2]*alpha_y
  alpha_eta = dxidx[2,1]*alpha_x + dxidx[2,2]*alpha_y
  alpha_n  = alpha_xi*nrm[1] + alpha_eta*nrm[2]
  =#
  alpha_n = alpha_x*nrm_scaled[1] + alpha_y*nrm_scaled[2]
  alpha_LF = params.LFalpha
  u = alpha_n*(uL + uR)*0.5 + absvalue(alpha_n)*(1 - alpha_LF)*0.5*(uL - uR)
  return u
end

function call{Tmsh, Tsol}(obj::LFFlux, uL::Tsol, uR::Tsol,
              nrm::AbstractArray{Tmsh,1}, params::ParamType3)

  alpha_x = params.alpha_x
  alpha_y = params.alpha_y
  alpha_z = params.alpha_z
  #=
  alpha_xi = dxidx[1,1]*alpha_x + dxidx[1,2]*alpha_y + dxidx[1,3]*alpha_z
  alpha_eta = dxidx[2,1]*alpha_x + dxidx[2,2]*alpha_y + dxidx[2,3]*alpha_z
  alpha_psi = dxidx[3,1]*alpha_x + dxidx[3,2]*alpha_y + dxidx[3,3]*alpha_z
  alpha_n  = alpha_xi*nrm[1] + alpha_eta*nrm[2] + alpha_psi*nrm[3]
  =#
  alpha_n = alpha_x*nrm[1] + alpha_y*nrm[2] + alpha_z*nrm[3]
  alpha_LF = params.LFalpha
  u = alpha_n*(uL + uR)*0.5 + absvalue(alpha_n)*(1 - alpha_LF)*0.5*(uL - uR)
  return u
end



global const FluxDict = Dict{ASCIIString, FluxType}(
"avgFlux" => avgFlux(),
"LFFlux" => LFFlux(),
)

function getFluxFunctors(mesh::AbstractDGMesh, sbp, eqn, opts)

  name = opts["Flux_name"]
  eqn.flux_func = FluxDict[name]
  return nothing
end
