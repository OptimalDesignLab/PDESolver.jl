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
      alpha_x = eqn.alpha_x
      alpha_y = eqn.alpha_y
      dxidx = sview(mesh.dxidx_face, :, :, j, i)
      nrm = sview(sbp.facenormal, :, fL)

      face_flux[1, j, i] = -functor(qL, qR, alpha_x, alpha_y, dxidx, nrm, eqn.params)
    end
  end

  return nothing
end

@doc """
### AdvectionEquationMod.calcSharedFaceIntegrals

  This function waits for the MPI receives to complete and then calculates
  the integrals over the shared interfaces.

  Inputs:
    mesh
    sbp
    eqn
    functor: the FluxType to use for the face flux

"""->
function calcSharedFaceIntegrals{Tmsh, Tsol}( mesh::AbstractDGMesh{Tmsh},
                            sbp::AbstractSBP, eqn::AdvectionData{Tsol},
                            opts, functor::FluxType)
# calculate the face flux and do the integration for the shared interfaces

  alpha_x = eqn.alpha_x
  alpha_y = eqn.alpha_y
  params = eqn.params

  for i=1:mesh.npeers
    params.t_wait += @elapsed idx, stat = MPI.Waitany!(mesh.recv_reqs)
    mesh.recv_stats[idx] = stat
    mesh.recv_reqs[idx] = MPI.REQUEST_NULL  # make sure this request is not used
    mesh.recv_waited[idx] = true
    # calculate the flux
    interfaces = mesh.shared_interfaces[idx]
    qL_arr = eqn.q_face_send[idx]
    qR_arr = eqn.q_face_recv[idx]
    dxidx_arr = mesh.dxidx_sharedface[idx]
    flux_arr = eqn.flux_sharedface[idx]

    # permute the received nodes to be in the elementR orientation
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
        flux_arr[1,k,j] = -functor(qL, qR, alpha_x, alpha_y, dxidx, nrm, 
                                    eqn.params)
      end
    end
    # end flux calculation

    # do the integration
    boundaryintegrate!(mesh.sbpface, mesh.bndries_local[idx], flux_arr, eqn.res)
  end  # end loop over npeers

  @debug1 sharedFaceLogging(mesh, sbp, eqn, opts, eqn.q_face_send, eqn.q_face_recv)

  return nothing
end


# element parallel version
function calcSharedFaceIntegrals_element{Tmsh, Tsol}( mesh::AbstractDGMesh{Tmsh},
                            sbp::AbstractSBP, eqn::AdvectionData{Tsol},
                            opts, functor::FluxType)

  q = eqn.q
  alpha_x = eqn.alpha_x
  alpha_y = eqn.alpha_y
  params = eqn.params

  @debug1 begin
    qL_face_arr = Array(Array{Tsol, 3}, mesh.npeers)
    qR_face_arr = Array(Array{Tsol, 3}, mesh.npeers)
    for i=1:mesh.npeers
      qL_face_arr[i] = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace, 
                                   mesh.peer_face_counts[i])
      qR_face_arr[i] = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace, 
                                   mesh.peer_face_counts[i])
    end
  end

  npeers = mesh.npeers
  val = sum(mesh.recv_waited)
  if val !=  mesh.npeers && val != 0
    throw(ErrorException("Receive waits in inconsistent state: $val / $npeers already waited on"))
  end

  # TODO: make these fields of params
  q_faceL = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  q_faceR = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  workarr = zeros(q_faceR)
  for i=1:mesh.npeers
    if val == 0
      params.t_wait += @elapsed idx, stat = MPI.Waitany!(mesh.recv_reqs)
      mesh.recv_stats[idx] = stat
      mesh.recv_reqs[idx] = MPI.REQUEST_NULL  # make sure this request is not used
      mesh.recv_waited[idx] = true
    else
      idx = i
    end

    interfaces = mesh.shared_interfaces[idx]
    bndries_local = mesh.bndries_local[idx]
    bndries_remote = mesh.bndries_remote[idx]
#    qL_arr = eqn.q_face_send[i]
    qR_arr = eqn.q_face_recv[idx]
    dxidx_arr = mesh.dxidx_sharedface[idx]
    flux_arr = eqn.flux_sharedface[idx]

    start_elnum = mesh.shared_element_offsets[idx]

    @debug1 qL_face_arr_i = qL_face_arr[i]
    @debug1 qR_face_arr_i = qR_face_arr[i]

    flush(params.f)
    for j=1:length(interfaces)
      iface_j = interfaces[j]
      bndryL_j = bndries_local[j]
      bndryR_j = bndries_remote[j]
      fL = bndryL_j.face

      # interpolate to face
      qL = sview(q, :, :, iface_j.elementL)
      el_r = iface_j.elementR - start_elnum + 1
      qR = sview(qR_arr, :, :, el_r)

      boundaryinterpolate!(mesh.sbpface, bndryL_j.face, qL, q_faceL)
      boundaryinterpolate!(mesh.sbpface, bndryR_j.face, qR, q_faceR)

      # permute elementR
      permvec = sview(mesh.sbpface.nbrperm, :, iface_j.orient)
      SummationByParts.permuteface!(permvec, workarr, q_faceR)

      @debug1 qL_face_arr_i[:, :, j] = q_faceL
      @debug1 qR_face_arr_i[:, :, j] = q_faceR

      # calculate flux
      for k=1:mesh.numNodesPerFace
        qL_k = q_faceL[k]
        qR_k = q_faceR[k]
        dxidx = sview(dxidx_arr, :, :, k, j)
        nrm = sview(sbp.facenormal, :, fL)

         flux_tmp = -functor(qL_k, qR_k, alpha_x, alpha_y, dxidx, nrm, 
                                      eqn.params)
         flux_arr[1,k,j] = flux_tmp
       end
     end  # end loop over interfaces

    # evaluate integral
    boundaryintegrate!(mesh.sbpface, bndries_local, flux_arr, eqn.res)
  end  # end loop over peers

  @debug1 sharedFaceLogging(mesh, sbp, eqn, opts, qL_face_arr, qR_face_arr)

  return nothing
end

function sharedFaceLogging{Tsol}(mesh, sbp, eqn::AdvectionData{Tsol}, opts, qL_arr, qR_arr)

  if opts["writeqface"]
    myrank = mesh.myrank
    for i=1:mesh.npeers
      tmp_arr = zeros(Tsol, mesh.numDofPerNode, 2, mesh.numNodesPerFace, mesh.peer_face_counts[i])
      qL_arr_i = qL_arr[i]
      qR_arr_i = qR_arr[i]
      for j = 1:mesh.peer_face_counts[i]
        for k=1:mesh.numNodesPerFace
          tmp_arr[:, 1, k, j] = qL_arr_i[:, k, j]
          tmp_arr[:, 2, k, j] = qR_arr_i[:, k, j]
        end
      end
      println(eqn.params.f, "q_sharedface $i = \n", tmp_arr)
      fname = string("qsharedface_", i, "_", myrank, ".dat")
      writedlm(fname, tmp_arr)
    end  # end loop over peers

  end  # end if

  if opts["write_fluxface"]
    for i=1:mesh.npeers
      fname = string("fluxsharedface_", i, "_", myrank, ".dat")
      writedlm(fname, eqn.flux_sharedface[i])
    end
  end

  flush(params.f)
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
              alpha_x, alpha_y, dxidx::AbstractArray{Tmsh,2}, 
              nrm::AbstractArray{Tmsh,1}, params::ParamType)

  alpha_xi = dxidx[1,1]*alpha_x + dxidx[1,2]*alpha_y
  alpha_eta = dxidx[2,1]*alpha_x + dxidx[2,2]*alpha_y
  alpha_n  = alpha_xi*nrm[1] + alpha_eta*nrm[2]

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
              alpha_x, alpha_y, dxidx::AbstractArray{Tmsh,2}, 
              nrm::AbstractArray{Tmsh,1}, params::ParamType)

  alpha_xi = dxidx[1,1]*alpha_x + dxidx[1,2]*alpha_y
  alpha_eta = dxidx[2,1]*alpha_x + dxidx[2,2]*alpha_y
  alpha_n  = alpha_xi*nrm[1] + alpha_eta*nrm[2]
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
