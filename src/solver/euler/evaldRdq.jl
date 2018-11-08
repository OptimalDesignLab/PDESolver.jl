# main file for computing dRdq transposed products

import PDESolver.evalResidual_revq

function evalResidual_revq(mesh::AbstractMesh, sbp::AbstractSBP, eqn::EulerData,
                     opts::Dict, t::Number=0.0)

  @assert mesh.commsize == 1

  time = eqn.params.time
  eqn.params.t = t  # record t to params
  myrank = mesh.myrank

  # TODO: Is this needed??
  # time.t_send += @elapsed if opts["parallel_type"] == 1
  #   println(eqn.params.f, "starting data exchange")
  #   #TODO: update this to use new parallel primatives
  #   startDataExchange(mesh, opts, eqn.res_bar,  eqn.q_face_send, eqn.q_face_recv, eqn.params.f)
  # end

  println("initially:")
  println("  max imag(q) = ", maximum(abs.(imag(eqn.q))))
  println("  max imag(res_bar) = ", maximum(abs.(imag(eqn.res_bar))))
  println("  max imag(q_bar) = ", maximum(abs.(imag(eqn.q_bar))))

  # !!!! MAKE SURE TO DO DATA EXCHANGE BEFORE !!!!

  # Forward sweep
  time.t_dataprep += @elapsed dataPrep_for_revq(mesh, sbp, eqn, opts)


  time.t_volume += @elapsed if opts["addVolumeIntegrals"]
    evalVolumeIntegrals_revq(mesh, sbp, eqn, opts)
  end
  println("after volume integrals:")
  println("  max imag(q) = ", maximum(abs.(imag(eqn.q))))
  println("  max imag(res_bar) = ", maximum(abs.(imag(eqn.res_bar))))
  println("  max imag(q_bar) = ", maximum(abs.(imag(eqn.q_bar))))

 
  if opts["use_GLS"]
    error("GLS not supported for revq product")
  end

  time.t_bndry += @elapsed if opts["addBoundaryIntegrals"]
    evalBoundaryIntegrals_revq(mesh, sbp, eqn, opts)
  end

  println("after boundary integrals:")
  println("  max imag(q) = ", maximum(abs.(imag(eqn.q))))
  println("  max imag(res_bar) = ", maximum(abs.(imag(eqn.res_bar))))
  println("  max imag(q_bar) = ", maximum(abs.(imag(eqn.q_bar))))


  # time.t_stab += @elapsed if opts["addStabilization"]
  #   addStabilization(mesh, sbp, eqn, opts)
  # end

  time.t_face += @elapsed if mesh.isDG && opts["addFaceIntegrals"]
    evalFaceIntegrals_revq(mesh, sbp, eqn, opts)
  end

  println("after face integrals:")
  println("  max imag(q) = ", maximum(abs.(imag(eqn.q))))
  println("  max imag(res_bar) = ", maximum(abs.(imag(eqn.res_bar))))
  println("  max imag(q_bar) = ", maximum(abs.(imag(eqn.q_bar))))

 
  time.t_sharedface += @elapsed if mesh.commsize > 1
    evalSharedFaceIntegrals_revq(mesh, sbp, eqn, opts)
  end

  time.t_source += @elapsed evalSourceTerm_revq(mesh, sbp, eqn, opts)


  # apply inverse mass matrix to eqn.res, necessary for CN
  if opts["use_Minv"]
    error("use_Minv not supported for revq product")
    applyMassMatrixInverse3D(mesh, sbp, eqn, opts, eqn.res)
  end


  time.t_dataprep += @elapsed dataPrep_revq(mesh, sbp, eqn, opts)

  return nothing
end  # end evalResidual


"""
  Dataprep-type function that should be called at the beginning of a reverse
  mode product, rather than the regular `dataPrep` function.
"""
function dataPrep_for_revq(mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP,
                   eqn::AbstractEulerData{Tsol, Tres}, opts) where {Tmsh, Tsol, Tres}

  # apply filtering to input
  if eqn.params.use_filter
    applyFilter(mesh, sbp, eqn, eqn.q, opts)
  end

  getAuxVars(mesh, eqn)

  if opts["use_staggered_grid"]
    error("staggered grid not supported for revq product")
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
### EulerEquationMod.dataPrep_revq

Reverse mode of dataPrep w.r.t mesh metrics

* mesh  : a mesh object
* sbp   : SBP operator object
* eqn   : an EulerData object
* opts  : options dictionary

"""->
function dataPrep_revq(mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP,
                   eqn::AbstractEulerData{Tsol, Tres}, opts) where {Tmsh, Tsol, Tres}

  # nothing to do here

  return nothing
end

@doc """
### EulerEquationMod.evalVolumeIntegrals_revq

Reverse mode of evalVolumeIntegrals with respect to the mesh metrics ∂ξ/∂x

**Arguments**

* mesh  : a mesh object
* sbp   : SBP operator object
* eqn   : an EulerData object
* opts  : options dictionary

"""->
function evalVolumeIntegrals_revq(mesh::AbstractMesh{Tmsh},
                             sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim}, opts) where {Tmsh,  Tsol, Tres, Tdim}

  integral_type = opts["volume_integral_type"]
  if integral_type == 1
    if opts["Q_transpose"] == true
      calcVolumeIntegrals_nopre_revq(mesh, sbp, eqn, opts)
    else
      error("Q_transpose = false not supported for revq product")
    end  # end if
  elseif integral_type == 2
    calcVolumeIntegralsSplitForm_revq(mesh, sbp, eqn, opts, eqn.volume_flux_func,
                                 eqn.volume_flux_func_revq)
  else
    throw(ErrorException("Unsupported volume integral type = $integral_type"))
  end

  return nothing
end  # end evalVolumeIntegrals

@doc """
### EulerEquationMod.evalBoundaryIntegrals_revq

Reverse mode of evalBoundaryIntegrals with respect to the mesh metrics ∂ξ/∂x

**Arguments**

* mesh  : a mesh object
* sbp   : SBP operator object
* eqn   : an EulerData object
* opts  : options dictionary

"""->

function evalBoundaryIntegrals_revq(mesh::AbstractMesh{Tmsh},
                               sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim}, opts) where {Tmsh, Tsol, Tres, Tdim}

  getBCFluxes_revq(mesh, sbp, eqn, opts)

  return nothing

end  # end evalBoundaryIntegrals

@doc """
### EulerEquationMod.evalFaceIntegrals_revq

Reverse mode of evalFaceIntegrals with respect to the mesh metrics ∂ξ/∂x

**Arguments**

* mesh  : a mesh object
* sbp   : SBP operator object
* eqn   : an EulerData object
* opts  : options dictionary

"""->
function evalFaceIntegrals_revq(mesh::AbstractDGMesh{Tmsh},
                           sbp::AbstractSBP, eqn::EulerData{Tsol}, opts) where {Tmsh, Tsol}

  face_integral_type = opts["face_integral_type"]
  if face_integral_type == 1
    calcFaceIntegral_nopre_revq(mesh, sbp, eqn, opts, eqn.flux_func_revq,
                                mesh.interfaces)

  elseif face_integral_type == 2

    error("not supported yet")
    getFaceElementIntegral_revq(mesh, sbp, eqn, eqn.face_element_integral_func,
                           eqn.flux_func, mesh.sbpface, mesh.interfaces)

  else
    throw(ErrorException("Unsupported face integral type = $face_integral_type"))
  end

  # do some output here?
  return nothing
end


"""
### EulerEquationMod.evalSourceTerm_revq

Reverse mode of evalSourceTerm w.r.t mesh metrics.

*This function has not been written properly. The developer must update this
documentation as and when this code is developed*

"""
function evalSourceTerm_revq(mesh::AbstractMesh{Tmsh},
                     sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim},
                     opts) where {Tmsh, Tsol, Tres, Tdim}


  # placeholder for multiple source term functionality (similar to how
  # boundary conditions are done)
  if opts["use_src_term"]
    error("source terms not supported for revq product")
  end

  return nothing
end  # end function


