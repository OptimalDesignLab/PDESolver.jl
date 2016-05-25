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
      nrm = sview(sbp.facenormal, :, fL)

      flux_j = sview(face_flux, :, j, i)
      functor(qL, qR, aux_vars, dxidx, nrm, flux_j, eqn.params)
#      # add the negative sign
#      for k=1:mesh.numDofPerNode
#        flux_j[k] = -flux_j[k]
#      end
    end
  end

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
        aux_vars = sview(eqn.aux_vars_sharedface, :, k, j)
        nrm = sview(sbp.facenormal, :, fL)
        flux_j = sview(flux_arr, :, k, j)
#        flux_arr[1,k,j] = -functor(qL, qR, alpha_x, alpha_y, dxidx, nrm, 
#                                    eqn.params)
        functor(qL, qR, aux_vars, dxidx, nrm, flux_j, eqn.params)
      end
    end
    # end flux calculation

    # do the integration
    boundaryintegrate!(mesh.sbpface, mesh.bndries_local[idx], flux_arr, eqn.res)
  end  # end loop over npeers

  @debug1 sharedFaceLogging(mesh, sbp, eqn, opts, eqn.q_face_send, eqn.q_face_recv)

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
function interpolateFace{Tsol}(mesh::AbstractDGMesh, sbp, eqn, opts, q::Abstract3DArray, q_face::AbstractArray{Tsol, 4})

  # interpolate solution
  interiorfaceinterpolate!(mesh.sbpface, mesh.interfaces, q, q_face)

  # recalculte aux_vars
  for i=1:mesh.numInterfaces
    for j=1:mesh.numNodesPerFace
      q_vals = sview(q_face, :, 1, j, i) # always use elementL
      eqn.aux_vars_face[1, j, i] = calcPressure(eqn.params, q_vals)
    end
  end

  if opts["paralle_type"] == 1  
    for peer=1:mesh.npeers
      q_vals_p = eqn.q_face_send[peer]
      aux_vars = eqn.aux_vars_sharedface[peer]
      for i=1:mesh.peer_face_counts[peer]
        for j=1:mesh.numNodesPerFace
          q_vals = sview(q_vals_p, :, j, i)
          eqn.aux_vars_sharedface[1, j, i] = calcPressure(eqn.params, q_vals)
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

function call{Tsol, Tres, Tmsh}(obj::RoeFlux, uL::AbstractArray{Tsol,1}, 
              uR::AbstractArray{Tsol,1}, 
              aux_vars, dxidx::AbstractArray{Tmsh, 2}, nrm::AbstractVector, 
              F::AbstractVector{Tres}, params::ParamType)

  RoeSolver(uL, uR, aux_vars, dxidx, nrm, F, params)
end

type AvgFlux <: FluxType
end

function call{Tsol, Tres, Tmsh}(obj::AvgFlux, uL::AbstractArray{Tsol,1}, 
              uR::AbstractArray{Tsol,1}, 
              aux_vars, dxidx::AbstractArray{Tmsh, 2}, nrm::AbstractVector, 
              F::AbstractVector{Tres}, params::ParamType)

  AvgSolver(uL, uR, aux_vars, dxidx, nrm, F, params)
end


@doc """
### EulerEquationMod.FluxDict

  This dictonary maps the names of the fluxes (ASCIIStrings) to the
  functor object itself.  All flux functors should be added to the dictionary.
"""->
global const FluxDict = Dict{ASCIIString, FluxType}(
"RoeFlux" => RoeFlux(),
"AvgFlux" => AvgFlux(),
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
  return nothing
end
