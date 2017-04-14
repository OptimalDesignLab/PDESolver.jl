# functions for calculating the flux for interior face integrals

@doc """
### EulerEquationMod.calcFaceFlux

  This function calculates the DG flux between a specified set of faces,
  using the solution data at the faces stored in eqn.q_face.
  Note that the flux is negated because the face integrals have a
  negative sign in the weak form.

  Conservative variables only!

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
  func( uL, qR, aux_vars, dxidx, nrm, flux_j, eqn.params

  where uL and uR are the solution values for a node on the left and right
  elements, aux_vars are the auxiliary variables for the node,
  dxidx is the scaled mapping jacobian for elementL, and nrm is the face
  normal in reference space. flux_j is the array of length numDofPerNode to be
  populated with the flux. params is eqn.params.

"""->
function calcFaceFlux{Tmsh,  Tsol, Tres, Tdim}( mesh::AbstractDGMesh{Tmsh},
                          sbp::AbstractSBP,
                          eqn::EulerData{Tsol, Tres, Tdim, :conservative},
                          functor::FluxType,
                          interfaces::AbstractArray{Interface,1},
                          face_flux::AbstractArray{Tres, 3})

  nfaces = length(interfaces)
  nrm = zeros(Tmsh, size(sbp.facenormal,1))
  for i=1:nfaces  # loop over faces
    interface_i = interfaces[i]
    for j = 1:mesh.numNodesPerFace
      eL = interface_i.elementL
      fL = interface_i.faceL

      # get components
      qL = sview(eqn.q_face, :, 1, j, i)
      qR = sview(eqn.q_face, :, 2, j, i)
      dxidx = sview(mesh.dxidx_face, :, :, j, i)
      aux_vars = sview(eqn.aux_vars_face, :, j, i)
      # nrm = sview(sbp.facenormal, :, fL)
      nrm[:] = sbp.facenormal[:,fL]

      flux_j = sview(face_flux, :, j, i)
      functor(eqn.params, qL, qR, aux_vars, dxidx, nrm, flux_j)
    end
  end

  return nothing
end

"""
  This function loops over interfaces and computes a face integral that
  uses data from all volume nodes. See FaceElementIntegralType for details on
  the integral performed.
"""
function getFaceElementIntegral{Tmsh, Tsol, Tres, Tdim}(
                           mesh::AbstractDGMesh{Tmsh},
                           sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim},
                           face_integral_functor::FaceElementIntegralType,
                           flux_functor::FluxType,
                           interfaces::AbstractArray{Interface, 1})

#  println("----- entered getECFaceIntegral -----")
  nfaces = length(interfaces)
  for i=1:nfaces
    iface = interfaces[i]
    elL = iface.elementL
    elR = iface.elementR
    qL = sview(eqn.q, :, :, elL)
    qR = sview(eqn.q, :, :, elR)
    aux_vars = sview(eqn.aux_vars, :, :, elL)
    dxidx_face = sview(mesh.dxidx_face, :, :, :, i)
    resL = sview(eqn.res, :, :, elL)
    resR = sview(eqn.res, :, :, elR)

    face_integral_functor(eqn.params, mesh.sbpface, iface, qL, qR, aux_vars,
                       dxidx_face, flux_functor, resL, resR)
  end

  return nothing
end

"""
  This function loops over shared interfaces and computes a face integral that
  uses data from all volume nodes.  See FaceElementIntegralType for details on
  the integral performed.
"""
function getSharedFaceElementIntegrals_element{Tmsh, Tsol, Tres}(
                            mesh::AbstractDGMesh{Tmsh},
                            sbp::AbstractSBP, eqn::EulerData{Tsol, Tres},
                            opts,
                            face_integral_functor::FaceElementIntegralType,
                            flux_functor::FluxType)

  if opts["parallel_data"] != "element"
    throw(ErrorException("cannot use getESSharedFaceIntegrals_element without parallel element data"))
  end

  q = eqn.q
  params = eqn.params
  # we don't care about elementR here, so use this throwaway array
  resR = Array(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)

  npeers = mesh.npeers
  val = sum(mesh.recv_waited)
  if val !=  mesh.npeers && val != 0
    throw(ErrorException("Receive waits in inconsistent state: $val / $npeers already waited on"))
  end

  for i=1:mesh.npeers
    if val == 0
      params.time.t_wait += @elapsed idx, stat = MPI.Waitany!(mesh.recv_reqs)
      mesh.recv_stats[idx] = stat
      mesh.recv_reqs[idx] = MPI.REQUEST_NULL  # make sure this request is not used
      mesh.recv_waited[idx] = true
    else
      idx = i
    end

    # get the data for the parallel interface
    interfaces = mesh.shared_interfaces[idx]
    bndries_local = mesh.bndries_local[idx]
    bndries_remote = mesh.bndries_remote[idx]
#    qL_arr = eqn.q_face_send[i]
    qR_arr = eqn.q_face_recv[idx]
    dxidx_face_arr = mesh.dxidx_sharedface[idx]
#    aux_vars_arr = eqn.aux_vars_sharedface[idx]
#    flux_arr = eqn.flux_sharedface[idx]

    start_elnum = mesh.shared_element_offsets[idx]

    # compute the integrals
    for j=1:length(interfaces)
      iface_j = interfaces[j]
      elL = iface_j.elementL
      elR = iface_j.elementR - start_elnum + 1  # is this always equal to j?
      qL = sview(q, :, :, elL)
      qR = sview(qR_arr, :, :, elR)
      aux_vars = sview(eqn.aux_vars, :, :, elL)
      dxidx_face = sview(dxidx_face_arr, :, :, :, j)
      resL = sview(eqn.res, :, :, elL)

      face_integral_functor(eqn.params, mesh.sbpface, iface_j, qL, qR, aux_vars,
                         dxidx_face, flux_functor, resL, resR)
    end  # end loop j

  end  # end loop over peers

  return nothing
end


@doc """
### EulerEquationMod.calcSharedFaceIntegrals

  This function waits for the MPI receives to complete and then calculates
  the integrals over the shared interfaces.

  Inputs:
    mesh
    sbp
    eqn
    functor: the FluxType to use for the face flux

"""->
function calcSharedFaceIntegrals{Tmsh, Tsol}( mesh::AbstractDGMesh{Tmsh},
                            sbp::AbstractSBP, eqn::EulerData{Tsol},
                            opts, functor::FluxType)
# calculate the face flux and do the integration for the shared interfaces

  if opts["parallel_data"] != "face"
    throw(ErrorException("cannot use calcSharedFaceIntegrals without parallel face data"))
  end


  params = eqn.params

  npeers = mesh.npeers
  val = sum(mesh.recv_waited)
  if val !=  mesh.npeers && val != 0
    throw(ErrorException("Receive waits in inconsistent state: $val / $npeers already waited on"))
  end


  for i=1:mesh.npeers
    if val == 0
      params.time.t_wait += @elapsed idx, stat = MPI.Waitany!(mesh.recv_reqs)
      mesh.recv_stats[idx] = stat
      mesh.recv_reqs[idx] = MPI.REQUEST_NULL  # don't use this request again
      mesh.recv_waited[idx] = true
    else
      idx = i
    end

    # calculate the flux
    interfaces = mesh.shared_interfaces[idx]
    qL_arr = eqn.q_face_send[idx]
    qR_arr = eqn.q_face_recv[idx]
    aux_vars_arr = eqn.aux_vars_sharedface[idx]
    dxidx_arr = mesh.dxidx_sharedface[idx]
    flux_arr = eqn.flux_sharedface[idx]

    # permute the received nodes to be in the elementR orientation
    permuteinterface!(mesh.sbpface, interfaces, qR_arr)
    for j=1:length(interfaces)
      interface_i = interfaces[j]
      for k=1:mesh.numNodesPerFace
        eL = interface_i.elementL
        fL = interface_i.faceL

        qL = sview(qL_arr, :, k, j)
        qR = sview(qR_arr, :, k, j)
        dxidx = sview(dxidx_arr, :, :, k, j)
        aux_vars = sview(aux_vars_arr, :, k, j)
        nrm = sview(sbp.facenormal, :, fL)
        flux_j = sview(flux_arr, :, k, j)
        functor(params, qL, qR, aux_vars, dxidx, nrm, flux_j)
      end
    end
    # end flux calculation

    # do the integration
    boundaryintegrate!(mesh.sbpface, mesh.bndries_local[idx], flux_arr, eqn.res, SummationByParts.Subtract())
  end  # end loop over npeers

  @debug1 sharedFaceLogging(mesh, sbp, eqn, opts, eqn.q_face_send, eqn.q_face_recv)

  return nothing
end


# element parallel version
function calcSharedFaceIntegrals_element{Tmsh, Tsol}( mesh::AbstractDGMesh{Tmsh},
                            sbp::AbstractSBP, eqn::EulerData{Tsol},
                            opts, functor::FluxType)

  if opts["parallel_data"] != "element"
    throw(ErrorException("cannot use calcSharedFaceIntegrals_element without parallel element data"))
  end


  q = eqn.q
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
      params.time.t_wait += @elapsed idx, stat = MPI.Waitany!(mesh.recv_reqs)
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
    aux_vars_arr = eqn.aux_vars_sharedface[idx]
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
      boundaryFaceInterpolate!(mesh.sbpface, bndryL_j.face, qL, q_faceL)
      boundaryFaceInterpolate!(mesh.sbpface, bndryR_j.face, qR, q_faceR)

      # permute elementR
      permvec = sview(mesh.sbpface.nbrperm, :, iface_j.orient)
      SummationByParts.permuteface!(permvec, workarr, q_faceR)

      @debug1 qL_face_arr_i[:, :, j] = q_faceL
      @debug1 qR_face_arr_i[:, :, j] = q_faceR

      # calculate flux
      for k=1:mesh.numNodesPerFace
        qL_k = sview(q_faceL, :, k)
        qR_k = sview(q_faceR, :, k)
        aux_vars = sview(aux_vars_arr, :, k, j)
        dxidx = sview(dxidx_arr, :, :, k, j)
        nrm = sview(sbp.facenormal, :, fL)
        flux_k = sview(flux_arr, :, k, j)

        aux_vars[1] = calcPressure(params, qL_k)

        functor(params, qL_k, qR_k, aux_vars, dxidx, nrm, flux_k)
       end
     end  # end loop over interfaces

    # evaluate integral
    boundaryintegrate!(mesh.sbpface, bndries_local, flux_arr, eqn.res, SummationByParts.Subtract())
  end  # end loop over peers

  @debug1 sharedFaceLogging(mesh, sbp, eqn, opts, qL_face_arr, qR_face_arr)

  return nothing
end


@doc """
### EulerEquationMod.interpolateFace

  This function interpolates the solution values from the internal nodes
  to the face flux points of the elements

  Inputs:
    mesh: an AbstractDGMesh
    sbp
    eqn
    opts
    q: a 3D array of solution values at the nodes, numDofPerNode x
       numNodesPerElement x numEl

  Inputs/Outputs:
    q_face: a 4D array of solution values at each interface,
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
      q_vals = sview(q_face, :, 1, j, i) # always use elementL
      eqn.aux_vars_face[1, j, i] = calcPressure(eqn.params, q_vals)
    end
  end

  if opts["parallel_data"] == 1
    for peer=1:mesh.npeers
      q_vals_p = eqn.q_face_send[peer]
      aux_vars = eqn.aux_vars_sharedface[peer]
      for i=1:mesh.peer_face_counts[peer]
        for j=1:mesh.numNodesPerFace
          q_vals = sview(q_vals_p, :, j, i)
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

type RoeFlux <: FluxType
end

function call{Tsol, Tres, Tmsh}(obj::RoeFlux, params::ParamType,
              uL::AbstractArray{Tsol,1},
              uR::AbstractArray{Tsol,1},
              aux_vars, dxidx::AbstractArray{Tmsh, 2}, nrm::AbstractVector,
              F::AbstractVector{Tres})

  RoeSolver(params, uL, uR, aux_vars, dxidx, nrm, F)
end

function call{Tsol, Tres}(obj::RoeFlux, params::ParamType,
              uL::AbstractArray{Tsol,1},
              uR::AbstractArray{Tsol,1},
              aux_vars::AbstractVector{Tres},
              nrm::AbstractVector,
              F::AbstractVector{Tres})

  RoeSolver(params, uL, uR, aux_vars, nrm, F)
  return nothing
end

type StandardFlux <: FluxType
end

function call{Tsol, Tres, Tmsh}(obj::StandardFlux, params::ParamType,
              uL::AbstractArray{Tsol,1},
              uR::AbstractArray{Tsol,1},
              aux_vars, dxidx::AbstractArray{Tmsh, 2}, nrm::AbstractVector,
              F::AbstractVector{Tres})

  calcEulerFlux_standard(params, uL, uR, aux_vars, dxidx, nrm, F)
end

function call{Tsol, Tres}(obj::StandardFlux, params::ParamType,
              uL::AbstractArray{Tsol,1},
              uR::AbstractArray{Tsol,1},
              aux_vars::AbstractVector{Tres},
              nrm::AbstractVector,
              F::AbstractVector{Tres})

  calcEulerFlux_standard(params, uL, uR, aux_vars, nrm, F)
  return nothing
end


type DucrosFlux <: FluxType
end

function call{Tsol, Tres, Tmsh}(obj::DucrosFlux, params::ParamType,
              uL::AbstractArray{Tsol,1},
              uR::AbstractArray{Tsol,1},
              aux_vars, dxidx::AbstractArray{Tmsh, 2}, nrm::AbstractVector,
              F::AbstractVector{Tres})

  calcEulerFlux_Ducros(params, uL, uR, aux_vars, dxidx, nrm, F)
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

type IRFlux <: FluxType
end

function call{Tsol, Tres, Tmsh}(obj::IRFlux, params::ParamType,
              uL::AbstractArray{Tsol,1},
              uR::AbstractArray{Tsol,1},
              aux_vars, dxidx::AbstractArray{Tmsh, 2}, nrm::AbstractVector,
              F::AbstractVector{Tres})

  calcEulerFlux_IR(params, uL, uR, aux_vars, dxidx, nrm, F)
end

function call{Tsol, Tres}(obj::IRFlux, params::ParamType,
              uL::AbstractArray{Tsol,1},
              uR::AbstractArray{Tsol,1},
              aux_vars::AbstractVector{Tres},
              nrm::AbstractVector,
              F::AbstractVector{Tres})

  calcEulerFlux_IR(params, uL, uR, aux_vars, nrm, F)
  return nothing
end

type IRSLFFlux <: FluxType
end

function call{Tsol, Tres, Tmsh}(obj::IRSLFFlux, params::ParamType,
              uL::AbstractArray{Tsol,1},
              uR::AbstractArray{Tsol,1},
              aux_vars, dxidx::AbstractArray{Tmsh, 2}, nrm::AbstractVector,
              F::AbstractVector{Tres})

  calcEulerFlux_IRSLF(params, uL, uR, aux_vars, dxidx, nrm, F)
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
"""->
global const FluxDict = Dict{ASCIIString, FluxType}(
"RoeFlux" => RoeFlux(),
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
              dxidx::AbstractArray{Tmsh, 2}, nrm::AbstractVector,
              flux_bar::AbstractVector{Tres}, dxidx_bar::AbstractArray{Tmsh, 2})

  RoeSolver_revm(params, uL, uR, aux_vars, dxidx, nrm, flux_bar, dxidx_bar)

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
