@doc """
### EulerEquationMod.evalrevm_transposeproduct

Reverse mode of evalResidual with respect to the mesh metrics ∂ξ/∂x

**Arguments**

* mesh  : a mesh object
* sbp   : SBP operator object
* eqn   : an EulerData object
* opts  : options dictionary

"""->

function evalrevm_transposeproduct{Tsol}(mesh::AbstractMesh, sbp::AbstractSBP, eqn::EulerData,
                     opts::Dict, input_array::AbstractArray{Tsol, 1}, t=0.0)

  array1DTo3D(mesh, sbp, eqn, opts, eqn.res_bar, input_array)

  time = eqn.params.time
  eqn.params.t = t  # record t to params
  myrank = mesh.myrank

  # TODO: Is this needed??
  # time.t_send += @elapsed if opts["parallel_type"] == 1
  #   println(eqn.params.f, "starting data exchange")
  #   #TODO: update this to use new parallel primatives
  #   startDataExchange(mesh, opts, eqn.res_bar,  eqn.q_face_send, eqn.q_face_recv, eqn.params.f)
  # end


  # !!!! MAKE SURE TO DO DATA EXCHANGE BEFORE !!!!

  # Forward sweep
  time.t_dataprep += @elapsed dataPrep(mesh, sbp, eqn, opts)


  time.t_volume += @elapsed if opts["addVolumeIntegrals"]
    evalVolumeIntegrals_revm(mesh, sbp, eqn, opts)
  end

  if opts["use_GLS"]
    println("adding boundary integrals")
    GLS(mesh,sbp,eqn)
  end

  time.t_bndry += @elapsed if opts["addBoundaryIntegrals"]
   evalBoundaryIntegrals_revm(mesh, sbp, eqn)
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

  # time.t_source += @elapsed evalSourceTerm_revm(mesh, sbp, eqn, opts)


  # # apply inverse mass matrix to eqn.res, necessary for CN
  # if opts["use_Minv"]
  #   applyMassMatrixInverse3D(mesh, sbp, eqn, opts, eqn.res)
  # end


  time.t_dataprep += @elapsed dataPrep_revm(mesh, sbp, eqn, opts)


  return nothing
end  # end evalResidual

@doc """
### EulerEquationMod.dataPrep_revm

Reverse mode of dataPrep w.r.t mesh metrics

* mesh  : a mesh object
* sbp   : SBP operator object
* eqn   : an EulerData object
* opts  : options dictionary

"""->
function dataPrep_revm{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP,
                                     eqn::AbstractEulerData{Tsol, Tres}, opts)

  getBCFluxes_revm(mesh, sbp, eqn, opts)

  if mesh.isDG

    if opts["face_integral_type"] == 1
      calcFaceFlux_revm(mesh, sbp, eqn, eqn.flux_func_bar, mesh.interfaces, eqn.flux_face_bar)
      # interpolateFace_revm(mesh, sbp, eqn, opts, eqn.q, eqn.q_face)
    end

    # fill!(eqn.q_bndry, 0.0)
    # fill!(eqn.q_face, 0.0)
    # fill!(eqn.flux_face, 0.0)
    # interpolateBoundary_revm(mesh, sbp, eqn, opts, eqn.q, eqn.q_bndry)

  end

  getEulerFlux_revm(mesh, sbp,  eqn, opts)

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

function evalVolumeIntegrals_revm{Tmsh,  Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh},
                             sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim}, opts)

  integral_type = opts["volume_integral_type"]
  if integral_type == 1
    if opts["Q_transpose"] == true
      fill!(eqn.flux_parametric_bar, 0.0) # zero out for first use
      for i=1:Tdim
        # Input: eqn.res_bar
        # Output: flux_parametric_bar
        weakdifferentiate_rev!(sbp, i, sview(eqn.flux_parametric_bar, :, :, :, i),
                               eqn.res_bar, trans=true)
      end
    else
      fill!(eqn.flux_parametric_bar, 0.0)
      for i=1:Tdim
        weakdifferentiate_rev!(sbp, i, sview(eqn.flux_parametric_bar, :, :, :, i),
                               eqn.res_bar, SummationByParts.Subtract(), trans=false)
      end
    end  # end if
  elseif integral_type == 2
    error("integral_type == 2 not supported")
    calcVolumeIntegralsSplitForm(mesh, sbp, eqn, opts, eqn.volume_flux_func)
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

function evalBoundaryIntegrals_revm{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh},
                               sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim})

  #TODO: remove conditional
  fill!(eqn.bndryflux_bar, 0.0)
  if mesh.isDG
    boundaryintegrate_rev!(mesh.sbpface, mesh.bndryfaces, eqn.bndryflux_bar,
                           eqn.res_bar, SummationByParts.Subtract())
  else
    boundaryintegrate_rev!(mesh.sbpface, mesh.bndryfaces, eqn.bndryflux_bar,
                           eqn.res_bar, SummationByParts.Subtract())
  end


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

function evalFaceIntegrals_revm{Tmsh, Tsol}(mesh::AbstractDGMesh{Tmsh},
                           sbp::AbstractSBP, eqn::EulerData{Tsol}, opts)

  face_integral_type = opts["face_integral_type"]
  fill!(eqn.flux_face_bar, 0.0)
  if face_integral_type == 1
    # Output to the call below = eqn.flux_face_bar
    interiorfaceintegrate_rev!(mesh.sbpface, mesh.interfaces, eqn.flux_face_bar,
                               eqn.res_bar, SummationByParts.Subtract())

  elseif face_integral_type == 2

    error("integral_type == 2 not supported")
    getFaceElementIntegral_rev(mesh, sbp, eqn, eqn.face_element_integral_func,
                           eqn.flux_func, mesh.interfaces)

  else
    throw(ErrorException("Unsupported face integral type = $face_integral_type"))
  end

  # do some output here?
  return nothing
end

@doc """
### EulerEquationMod.evalSharedFaceIntegrals_revm

Reverse mode evalSharedFaceIntegrals with respect to the mesh metrics ∂ξ/∂x

"""

function evalSharedFaceIntegrals_revm(mesh::AbstractDGMesh, sbp, eqn, opts)

  face_integral_type = opts["face_integral_type"]
  if face_integral_type == 1

    if opts["parallel_data"] == "face"
      calcSharedFaceIntegrals_revm(mesh, sbp, eqn, opts, eqn.flux_func_bar)
    elseif opts["parallel_data"] == "element"
      calcSharedFaceIntegrals_element_revm(mesh, sbp, eqn, opts, eqn.flux_func_bar)
    else
      throw(ErrorException("unsupported parallel data type"))
    end

  elseif face_integral_type == 2

      error("integral_type == 2 not supported")
    getSharedFaceElementIntegrals_element(mesh, sbp, eqn, opts,
                                eqn.face_element_integral_func,  eqn.flux_func)
  else
    throw(ErrorException("unsupported face integral type = $face_integral_type"))
  end

  return nothing
end

"""
### EulerEquationMod.calcSharedFaceIntegrals_revm

Reverse mode of calcSharedFaceIntegrals w.r.t mesh metrics, ∂ξ/∂x.

"""

function calcSharedFaceIntegrals_revm{Tmsh, Tsol}( mesh::AbstractDGMesh{Tmsh},
                            sbp::AbstractSBP, eqn::EulerData{Tsol},
                            opts, functor_revm::FluxType)
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

"""
### EulerEquationMod.evalSourceTerm_revm

Reverse mode of evalSourceTerm w.r.t mesh metrics.

*This function has not been written properly. The developer must update this
documentation as and when this code is developed*

"""

function evalSourceTerm_revm{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh},
                     sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim},
                     opts)


  # placeholder for multiple source term functionality (similar to how
  # boundary conditions are done)
  if opts["use_src_term"]
    applySourceTerm(mesh, sbp, eqn, opts, eqn.src_func)
  end

  return nothing
end  # end function


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

function updatePumiMesh{Tmsh}(mesh::AbstractDGMesh{Tmsh}, sbp::AbstractSBP,
                              volNodes::AbstractArray{Tmsh, 2})

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
