# main file for evalResidual_revm

import PDESolver.evalResidual_revm

@doc """
### EulerEquationMod.evalrevm_transposeproduct

Reverse mode of evalResidual with respect to the mesh metrics ∂ξ/∂x

**Arguments**

* mesh  : a mesh object
* sbp   : SBP operator object
* eqn   : an EulerData object
* opts  : options dictionary

"""->

function evalResidual_revm(mesh::AbstractMesh, sbp::AbstractSBP, eqn::EulerData,
                     opts::Dict, t::Number=0.0)

  time = eqn.params.time
  eqn.params.t = t  # record t to params
  myrank = mesh.myrank

  # Forward sweep
  time.t_dataprep += @elapsed dataPrep_for_revm(mesh, sbp, eqn, opts)


  time.t_volume += @elapsed if opts["addVolumeIntegrals"]
    evalVolumeIntegrals_revm(mesh, sbp, eqn, opts)
  end

  if opts["use_GLS"]
    error("GLS not supported for revm product")
  end

  time.t_bndry += @elapsed if opts["addBoundaryIntegrals"]
    evalBoundaryIntegrals_revm(mesh, sbp, eqn, opts)
  end

  # time.t_stab += @elapsed if opts["addStabilization"]
  #   addStabilization(mesh, sbp, eqn, opts)
  # end

  time.t_face += @elapsed if mesh.isDG && opts["addFaceIntegrals"]
    evalFaceIntegrals_revm(mesh, sbp, eqn, opts)
  end

  time.t_sharedface += @elapsed if mesh.commsize > 1
    evalSharedFaceIntegrals_revm(mesh, sbp, eqn, opts)
  end

  time.t_source += @elapsed evalSourceTerm_revm(mesh, sbp, eqn, opts)

  # apply inverse mass matrix to eqn.res, necessary for CN
  if opts["use_Minv"]
    error("use_Minv not supported for revm product")
  end


  time.t_dataprep += @elapsed dataPrep_revm(mesh, sbp, eqn, opts)

  return nothing
end  # end evalResidual


"""
  Dataprep-type function that should be called at the beginning of a reverse
  mode product, rather than the regular `dataPrep` function.
"""
function dataPrep_for_revm(mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP,
                   eqn::AbstractEulerData{Tsol, Tres}, opts) where {Tmsh, Tsol, Tres}

  # apply filtering to input
  if eqn.params.use_filter
    applyFilter(mesh, sbp, eqn, eqn.q, opts)
  end

  getAuxVars(mesh, eqn)

  if opts["use_staggered_grid"]
    error("staggered grid not supported for revm product")
    aux_vars = zeros(Tres, 1, mesh.mesh2.numNodesPerElement)
    for i=1:mesh.numEl
      qs = ro_sview(eqn.q, :, :, i)
      qf = sview(eqn.q_flux, :, :, i)
      interpolateElementStaggered(eqn.params, mesh, qs, aux_vars, qf)
    end
  end


  # only eqn.q_bndry is required, all other methods are no-precompute
  if mesh.isDG
    interpolateBoundary(mesh, sbp, eqn, opts, eqn.q, eqn.q_bndry, eqn.aux_vars_bndry)
  end

  if eqn.params.use_edgestab
    stabscale(mesh, sbp, eqn)
  end

  return nothing
end



@doc """
### EulerEquationMod.dataPrep_revm

Reverse mode of dataPrep w.r.t mesh metrics

* mesh  : a mesh object
* sbp   : SBP operator object
* eqn   : an EulerData object
* opts  : options dictionary

"""->
function dataPrep_revm(mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP,
                   eqn::AbstractEulerData{Tsol, Tres}, opts) where {Tmsh, Tsol, Tres}

  # nothing to do here

  return nothing
end

@doc """
### EulerEquationMod.evalVolumeIntegrals_revm

Reverse mode of evalVolumeIntegrals with respect to the mesh metrics ∂ξ/∂x

**Arguments**

* mesh  : a mesh object
* sbp   : SBP operator object
* eqn   : an EulerData object
* opts  : options dictionary

"""->
function evalVolumeIntegrals_revm(mesh::AbstractMesh{Tmsh},
                             sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim}, opts) where {Tmsh,  Tsol, Tres, Tdim}

  integral_type = opts["volume_integral_type"]
  if integral_type == 1
    if opts["Q_transpose"] == true
      calcVolumeIntegrals_nopre_revm(mesh, sbp, eqn, opts)
    else
      error("Q_transpose = false not supported for revm product")
    end  # end if
  elseif integral_type == 2
    calcVolumeIntegralsSplitForm_revm(mesh, sbp, eqn, opts, eqn.volume_flux_func,
                                 eqn.volume_flux_func_revm)
  else
    throw(ErrorException("Unsupported volume integral type = $integral_type"))
  end

  return nothing
end  # end evalVolumeIntegrals

@doc """
### EulerEquationMod.evalBoundaryIntegrals_revm

Reverse mode of evalBoundaryIntegrals with respect to the mesh metrics ∂ξ/∂x

**Arguments**

* mesh  : a mesh object
* sbp   : SBP operator object
* eqn   : an EulerData object
* opts  : options dictionary

"""->

function evalBoundaryIntegrals_revm(mesh::AbstractMesh{Tmsh},
                               sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim}, opts) where {Tmsh, Tsol, Tres, Tdim}

  getBCFluxes_revm(mesh, sbp, eqn, opts)

  return nothing

end  # end evalBoundaryIntegrals

@doc """
### EulerEquationMod.evalFaceIntegrals_revm

Reverse mode of evalFaceIntegrals with respect to the mesh metrics ∂ξ/∂x

**Arguments**

* mesh  : a mesh object
* sbp   : SBP operator object
* eqn   : an EulerData object
* opts  : options dictionary

"""->

function evalFaceIntegrals_revm(mesh::AbstractDGMesh{Tmsh},
                           sbp::AbstractSBP, eqn::EulerData{Tsol}, opts) where {Tmsh, Tsol}

  face_integral_type = opts["face_integral_type"]
  if face_integral_type == 1
    calcFaceIntegral_nopre_revm(mesh, sbp, eqn, opts, eqn.flux_func_revm,
                         mesh.interfaces)

  elseif face_integral_type == 2

    getFaceElementIntegral_revm(mesh, sbp, eqn, eqn.face_element_integral_func,
                           eqn.flux_func, mesh.sbpface, mesh.interfaces)

  else
    throw(ErrorException("Unsupported face integral type = $face_integral_type"))
  end

  # do some output here?
  return nothing
end

@doc """

Reverse mode evalSharedFaceIntegrals with respect to the mesh metrics

"""

function evalSharedFaceIntegrals_revm(mesh::AbstractDGMesh, sbp, eqn, opts)

  if getParallelData(eqn.shared_data) != "element"
    error("""parallel data setting must be "element" """)
  end


  face_integral_type = opts["face_integral_type"]
  if face_integral_type == 1

    finishExchangeData(mesh, sbp, eqn, opts, eqn.shared_data, calcSharedFaceIntegrals_element_revm)

  elseif face_integral_type == 2
    
    finishExchangeData(mesh, sbp, eqn, opts, eqn.shared_data, calcSharedFaceElementIntegrals_element_revm)
#    getSharedFaceElementIntegrals_element(mesh, sbp, eqn, opts, eqn.face_element_integral_func,  eqn.flux_func)
  else
    throw(ErrorException("unsupported face integral type = $face_integral_type"))
  end

  return nothing
end


#=
#TODO: what is this and why is it here?
"""
### EulerEquationMod.calcSharedFaceIntegrals_revm

Reverse mode of calcSharedFaceIntegrals w.r.t mesh metrics, ∂ξ/∂x.

"""

function calcSharedFaceIntegrals_revm( mesh::AbstractDGMesh{Tmsh},
                            sbp::AbstractSBP, eqn::EulerData{Tsol},
                            opts, functor_revm::FluxType) where {Tmsh, Tsol}
  # calculate the face flux and do the integration for the shared interfaces

  #TODO: update this to use the new parallel communication primatives
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
#    dxidx_arr = mesh.dxidx_sharedface[idx]
    nrm_arr = mesh.nrm_sharedface[idx]
    flux_arr_bar = eqn.flux_sharedface_bar[idx]

    # permute the received nodes to be in the elementR orientation
    permuteinterface!(mesh.sbpface, interfaces, qR_arr)
    for j=1:length(interfaces)
      interface_i = interfaces[j]
      for k=1:mesh.numNodesPerFace
        eL = interface_i.elementL
        fL = interface_i.faceL

        qL = sview(qL_arr, :, k, j)
        qR = sview(qR_arr, :, k, j)
        aux_vars = sview(aux_vars_arr, :, k, j)
        nrm_xy = sview(nrm_arr, :, k, j)
        flux_j = sview(flux_arr_bar, :, k, j)
        functor_revm(params, qL, qR, aux_vars, nrm_xy, flux_j)
      end
    end
    # end flux calculation

    # do the integration
    boundaryintegrate_rev!(mesh.sbpface, mesh.bndries_local[idx], flux_arr_bar, eqn.res_bar, SummationByParts.Subtract())
  end  # end loop over npeers

  @debug1 sharedFaceLogging(mesh, sbp, eqn, opts, eqn.q_face_send, eqn.q_face_recv)

  return nothing
end
=#


"""
### EulerEquationMod.evalSourceTerm_revm

Reverse mode of evalSourceTerm w.r.t mesh metrics.

*This function has not been written properly. The developer must update this
documentation as and when this code is developed*

"""
function evalSourceTerm_revm(mesh::AbstractMesh{Tmsh},
                     sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim},
                     opts) where {Tmsh, Tsol, Tres, Tdim}


  # placeholder for multiple source term functionality (similar to how
  # boundary conditions are done)
  if opts["use_src_term"]
    error("source terms not supported for revm product")
  end

  return nothing
end  # end function

#=
@doc """
###EulerEquationMod.updatePumiMesh

Given the volume nodes in the MeshMovement `volNodes` data structure, update the
Pumi mesh object

* mesh  : a mesh object
* sbp   : SBP operator object
* volNodes : A 2D array containing coordinates of all the mesh vertices/nodes
             owned by a perticular rank. There are no repeated vertices, i.e.
             all the coordinates are unique. size(volNodes) = (3, numVert)

"""->

function updatePumiMesh(mesh::AbstractDGMesh{Tmsh}, sbp::AbstractSBP,
                        volNodes::AbstractArray{Tmsh, 2}) where Tmsh

  for i = 1:mesh.numEl
    for j = 1:size(mesh.vert_coords,2)
      # Get the vertex numbering on the portion of mesh owned by the processor
      local_vertnum = mesh.element_vertnums[j,i]
      for k = 1:mesh.dim
        mesh.vert_coords[k,j,i] = volNodes[k,local_vertnum]
      end
    end
  end

  for i = 1:mesh.numEl
    update_coords(mesh, i, mesh.vert_coords[:,:,i])
  end
  commit_coords(mesh, sbp)

  return nothing
end
=#
