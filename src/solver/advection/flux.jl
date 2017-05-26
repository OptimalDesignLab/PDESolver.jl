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
 
  @debug1 println(eqn.params.f, "entered calcFaceFlux")
  @debug1 flush(eqn.params.f)

  nfaces = length(interfaces)
  for i=1:nfaces  # loop over faces
    interface_i = interfaces[i]
    @debug1 println(eqn.params.f, "i: $i  interface_i.elementL: ", interface_i.elementL, "  elementR: ", interface_i.elementR)
    @debug1 flush(eqn.params.f)
    for j = 1:mesh.numNodesPerFace
      eL = interface_i.elementL
      fL = interface_i.faceL

      # get components
      qL = eqn.q_face[1, 1, j, i]
      qR = eqn.q_face[1, 2, j, i]
      dxidx = sview(mesh.dxidx_face, :, :, j, i)
      nrm = sview(sbp.facenormal, :, fL)

      face_flux[1, j, i] = -functor(qL, qR, dxidx, nrm, eqn.params)
      @debug1 println(eqn.params.f, "  j: $j  face_flux[1, j, i]: ", face_flux[1, j, i], "  qL: $qL  qR: $qR")
      @debug1 flush(eqn.params.f)
    end
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

  calcSharedFaceIntegrals_inner(mesh, sbp, eqn, opts, data, eqn.flux_func)

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
  dxidx_arr = mesh.dxidx_sharedface[idx]
  flux_arr = eqn.flux_sharedface[idx]

  # permute the received nodes to be in the elementR orientation
  # TODO: move this into finishExchangeData, to avoid possible double
  # permutation of the data
  permuteinterface!(mesh.sbpface, interfaces, qR_arr)

  for j=1:length(interfaces)
    interface_i = interfaces[j]
    for k=1:mesh.numNodesPerFace
      eL = interface_i.elementL
      fL = interface_i.faceL

      qL = qL_arr[1, k, j]
      qR = qR_arr[1, k, j]
      dxidx = sview(dxidx_arr, :, :, k, j)
      nrm = sview(sbp.facenormal, :, fL)
      flux_arr[1,k,j] = -functor(qL, qR, dxidx, nrm, 
                                  eqn.params)
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
  Thin wrapper around calcSharedFaceIntegrals_inner.  This function is passed
  to finishDataExchange, and internally calls calcSharedFaceIntegrals_inner.
  See finishDataExchange for details on the interface and 
  calcSharedFaceIntegrals_inner for the integral that is computed.
"""
function calcSharedFaceIntegrals_element{Tmsh, Tsol}(mesh::AbstractDGMesh{Tmsh},
                            sbp::AbstractSBP, eqn::AdvectionData{Tsol},
                            opts, data::SharedFaceData)

  calcSharedFaceIntegrals_element_inner(mesh, sbp, eqn, opts, data, eqn.flux_func)

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
  dxidx_arr = mesh.dxidx_sharedface[idx]
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

    boundaryFaceInterpolate!(mesh.sbpface, bndryL_j.face, qL, q_faceL)
    boundaryFaceInterpolate!(mesh.sbpface, bndryR_j.face, qR, q_faceR)

    # permute elementR
    permvec = sview(mesh.sbpface.nbrperm, :, iface_j.orient)
    SummationByParts.permuteface!(permvec, workarr, q_faceR)

    @debug1 qL_face_arr[:, :, j] = q_faceL
    @debug1 qR_face_arr[:, :, j] = q_faceR

    # calculate flux
    for k=1:mesh.numNodesPerFace
      qL_k = q_faceL[k]
      qR_k = q_faceR[k]
      dxidx = sview(dxidx_arr, :, :, k, j)
      nrm = sview(sbp.facenormal, :, fL)

      flux_tmp = -functor(qL_k, qR_k, dxidx, nrm, 
                                   eqn.params)
      flux_arr[1,k,j] = flux_tmp
    end
  end  # end loop over interfaces

  # evaluate integral
  boundaryintegrate!(mesh.sbpface, bndries_local, flux_arr, eqn.res)

  @debug1 sharedFaceLogging(mesh, sbp, eqn, opts, data, qL_face_arr, qR_face_arr)

  return nothing
end



@doc """
### AdvectionEquationMod.avgFlux

  This flux function averages the two states to calculate the flux

  Inputs:
    uL: the left state
    uR: the right state
    alpha_x: the advection velocity in the x direction
    alpha_y: the advection velocity in the y direction
    dxidx: the scaled mapping jacobian for elementL
    nrm: the face normal vector for elementL in parametric space

  Outputs:
    the flux
"""->
type avgFlux <: FluxType
end

function call{Tmsh, Tsol}(obj::avgFlux, uL::Tsol, uR::Tsol,
              dxidx::AbstractArray{Tmsh,2}, 
              nrm::AbstractArray{Tmsh,1}, params::ParamType2)

  alpha_x = params.alpha_x
  alpha_y = params.alpha_y
  alpha_xi = dxidx[1,1]*alpha_x + dxidx[1,2]*alpha_y
  alpha_eta = dxidx[2,1]*alpha_x + dxidx[2,2]*alpha_y
  alpha_n  = alpha_xi*nrm[1] + alpha_eta*nrm[2]

  u = alpha_n*(uL + uR)*0.5
  return u
end

function call{Tmsh, Tsol}(obj::avgFlux, uL::Tsol, uR::Tsol,
              dxidx::AbstractArray{Tmsh,2}, 
              nrm::AbstractArray{Tmsh,1}, params::ParamType3)

  alpha_x = params.alpha_x
  alpha_y = params.alpha_y
  alpha_z = params.alpha_z
  alpha_xi = dxidx[1,1]*alpha_x + dxidx[1,2]*alpha_y + dxidx[1,3]*alpha_z
  alpha_eta = dxidx[2,1]*alpha_x + dxidx[2,2]*alpha_y + dxidx[2,3]*alpha_z
  alpha_psi = dxidx[3,1]*alpha_x + dxidx[3,2]*alpha_y + dxidx[3,3]*alpha_z
  alpha_n  = alpha_xi*nrm[1] + alpha_eta*nrm[2] + alpha_psi*nrm[3]

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
    alpha_x: the advection velocity in the x direction
    alpha_y: the advection velocity in the y direction
    dxidx: the scaled mapping jacobian for elementL
    nrm: the face normal vector for elementL in parametric space

  Outputs:
    the flux
"""->
type LFFlux <: FluxType
end

function call{Tmsh, Tsol}(obj::LFFlux, uL::Tsol, uR::Tsol,
              dxidx::AbstractArray{Tmsh,2}, 
              nrm::AbstractArray{Tmsh,1}, params::ParamType2)

  alpha_x = params.alpha_x
  alpha_y = params.alpha_y
  alpha_xi = dxidx[1,1]*alpha_x + dxidx[1,2]*alpha_y
  alpha_eta = dxidx[2,1]*alpha_x + dxidx[2,2]*alpha_y
  alpha_n  = alpha_xi*nrm[1] + alpha_eta*nrm[2]
  alpha_LF = params.LFalpha
  u = alpha_n*(uL + uR)*0.5 + absvalue(alpha_n)*(1 - alpha_LF)*0.5*(uL - uR)
  return u
end

function call{Tmsh, Tsol}(obj::LFFlux, uL::Tsol, uR::Tsol,
              dxidx::AbstractArray{Tmsh,2}, 
              nrm::AbstractArray{Tmsh,1}, params::ParamType3)

  alpha_x = params.alpha_x
  alpha_y = params.alpha_y
  alpha_z = params.alpha_z
  alpha_xi = dxidx[1,1]*alpha_x + dxidx[1,2]*alpha_y + dxidx[1,3]*alpha_z
  alpha_eta = dxidx[2,1]*alpha_x + dxidx[2,2]*alpha_y + dxidx[2,3]*alpha_z
  alpha_psi = dxidx[3,1]*alpha_x + dxidx[3,2]*alpha_y + dxidx[3,3]*alpha_z

  alpha_n  = alpha_xi*nrm[1] + alpha_eta*nrm[2] + alpha_psi*nrm[3]
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
