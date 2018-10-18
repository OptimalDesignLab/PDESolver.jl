import PDESolver: evalJacobian, evalJacobianStrong

"""
  Euler implementation of `evalJacobian`.  Currently only supports the
  Roe scheme.
"""
function evalJacobian(mesh::AbstractMesh, sbp::AbstractSBP, eqn::EulerData, 
                      opts::Dict, assembler::AssembleElementData, t=0.0;
                      start_comm=false)
  time = eqn.params.time
  eqn.params.t = t  # record t to params
  myrank = mesh.myrank

  @assert opts["parallel_type"] == 2
#  println("entered evalResidual")
  time.t_send += @elapsed if start_comm
    startSolutionExchange(mesh, sbp, eqn, opts)
  end


  time.t_dataprep_diff += @elapsed dataPrep_diff(mesh, sbp, eqn, opts)
#  println("dataPrep @time printed above")

  time.t_volume_diff += @elapsed if opts["addVolumeIntegrals"]
    evalVolumeIntegrals_diff(mesh, sbp, eqn, opts, assembler)
#    println("volume integral @time printed above")
  end

  if opts["use_GLS"]
    error("use_GLS not supported by evalJacobian")
  end

  time.t_bndry_diff += @elapsed if opts["addBoundaryIntegrals"]
    evalBoundaryIntegrals_diff(mesh, sbp, eqn, opts, assembler)
#   println("boundary integral @time printed above")
  end


  time.t_stab_diff += @elapsed if opts["addStabilization"]
    addStabilization_diff(mesh, sbp, eqn, opts, assembler)
#    println("stabilizing @time printed above")
  end



  time.t_face_diff += @elapsed if mesh.isDG && opts["addFaceIntegrals"]
    evalFaceIntegrals_diff(mesh, sbp, eqn, opts, assembler)
#    println("face integral @time printed above")
  end



  time.t_sharedface_diff += @elapsed if mesh.commsize > 1
    evalSharedFaceIntegrals_diff(mesh, sbp, eqn, opts, assembler)
#    println("evalSharedFaceIntegrals @time printed above")
  end



  time.t_source_diff += @elapsed evalSourceTerm_diff(mesh, sbp, eqn, opts, assembler)

#  println("source integral @time printed above")

  # apply inverse mass matrix to eqn.res, necessary for CN
  #TODO: this will have to be done at the element level
  if opts["use_Minv"]
#    error("use_Minv must be applied at the element level")
#    applyMassMatrixInverse3D(mesh, sbp, eqn, opts, eqn.res)
  end

  return nothing
end



function evalJacobianStrong(mesh::AbstractMesh, sbp::AbstractSBP,
                      eqn::EulerData, 
                      opts::Dict, assembler::AssembleElementData, t=0.0;
                      start_comm=false)
# currently this function neglects the SAT terms (including boundary conditions)

  time = eqn.params.time
  eqn.params.t = t  # record t to params
  myrank = mesh.myrank


  time.t_volume_diff += @elapsed if opts["addVolumeIntegrals"]
    calcVolumeIntegralsStrong_nopre_diff(mesh, sbp, eqn, opts, assembler)
#    println("volume integral @time printed above")
  end

  return nothing
end



"""
  Differentiated version of [`dataPrep`](@ref).  The Jacobian computation
  does use any of the arrays of precomputed values (options keys 
  `precompute_*`).  Those keys can be specified by the user and will affect
  [`evalResidual`](@ref), but not [`evalJacobian`](@ref).

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
"""
function dataPrep_diff(mesh::AbstractMesh{Tmsh},
     sbp::AbstractSBP,
     eqn::AbstractEulerData{Tsol, Tres}, opts) where {Tmsh, Tsol, Tres}
# gather up all the data needed to do vectorized operatinos on the mesh
# calculates all mesh wide quantities in eqn
# this is almost the exact list of everything we *shouldn't* be storing, but
# rather recalculating on the fly

#println("Entered dataPrep()")

#  println("typeof(eqn) = ", typeof(eqn))
#  println("typeof(eqn.params) = ", typeof(eqn.params))

  # apply filtering to input
  if eqn.params.use_filter
    applyFilter(mesh, sbp, eqn, eqn.q, opts)
  end

  # zero out res
#  fill!(eqn.res, 0.0)
#  fill!(eqn.res_edge, 0.0)

  getAuxVars(mesh, eqn)
#  println("  getAuxVars @time printed above")

  if opts["check_density"]
    checkDensity(eqn, mesh)
#    println("  checkDensity @time printed above")
  end

  if opts["check_pressure"]
#    throw(ErrorException("I'm done"))
    checkPressure(eqn, mesh)
#    println("  checkPressure @time printed above")
  end

  # calculate fluxes

  if opts["use_staggered_grid"]
    error("staggered grid not supported by calcJacobian()")
  end

  # the jacobian computation doesn't precompute anything

  if eqn.params.use_edgestab
    error("Use_edgestab not supported by calcJacobian()")
    stabscale(mesh, sbp, eqn)
  end
#  println("  stabscale @time printed above")

  return nothing
end # end function dataPrep


"""
  Differentiated version of [`evalVolumeIntegrals`](@ref).  Throws an error
  for unsupported schemes.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * assembler
"""
function evalVolumeIntegrals_diff(mesh::AbstractMesh{Tmsh},
                             sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim},
                             opts, assembler::AssembleElementData) where {Tmsh,  Tsol, Tres, Tdim}

  integral_type = opts["volume_integral_type"]

  if opts["Q_transpose"] == false
    error("Q_transpose == false not supported by evalJacobian")
  end  # end if Q_transpose


  if integral_type == 1  # regular volume integrals
    calcVolumeIntegrals_nopre_diff(mesh, sbp, eqn, opts, assembler)
  elseif integral_type == 2  # entropy stable formulation
    calcVolumeIntegralsSplitForm_diff(mesh, sbp, eqn, opts,
                                      eqn.volume_flux_func_diff, assembler)
  else
    throw(ErrorException("Unsupported volume integral type = $integral_type"))
  end

end  # end evalVolumeIntegrals



"""
  Differentiated version of [`evalBoundaryIntegrals`](@ref).

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * assembler.
"""
function evalBoundaryIntegrals_diff(mesh::AbstractMesh{Tmsh},
                               sbp::AbstractSBP,
                               eqn::EulerData{Tsol, Tres, Tdim}, opts,
                               assembler::AssembleElementData) where {Tmsh, Tsol, Tres, Tdim}

  if mesh.isDG
    # when precompute_boundary_flux == false, this fuunction does the
    # integration too, updating res
    getBCFluxes_diff(mesh, sbp, eqn, opts, assembler)
  else
    error("CG meshes not supported by evalJacobian()")
  end

  return nothing

end  # end evalBoundaryIntegrals

"""
  Differentiated version of [`addStabilization`](@ref).

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * assembler
"""
function addStabilization_diff(mesh::AbstractMesh{Tmsh},
                          sbp::AbstractSBP, eqn::EulerData{Tsol}, opts,
                          assembler::AssembleElementData) where {Tmsh,  Tsol}

#  println("==== start of addStabilization ====")

  # ----- Edge Stabilization -----#
  if eqn.params.use_edgestab  #TODO: use options dictionary instead
#    println("applying edge stabilization")
    error("edge stabilization not supported by evalJacobian()")
  end

  # ----- Filtering -----
  if eqn.params.use_res_filter
    error("use_res_filter not supported by evalJacobian()")
  end

  # ----- Artificial Dissipation -----
  if eqn.params.use_dissipation
    error("use_dissipation not supported by evalJacobian()")
  end

  if opts["use_GLS2"]
    error("use_GLS2 not supported by evalJacobian()")
     applyGLS3(mesh, sbp, eqn, opts)
#    test_GLS(mesh, sbp, eqn, opts)
  end

#  println("==== end of addStabilization ====")



  return nothing
end

"""
  Differentiated version of [`evalFaceIntegrals_diff`](@ref).  Throws an
  error for unsupported schemes.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * assembler
"""
function evalFaceIntegrals_diff(mesh::AbstractDGMesh{Tmsh},
                    sbp::AbstractSBP,
                    eqn::EulerData{Tsol}, opts,
                    assembler::AssembleElementData) where {Tmsh, Tsol}

  face_integral_type = opts["face_integral_type"]
  if face_integral_type == 1
      calcFaceIntegral_nopre_diff(mesh, sbp, eqn, opts, eqn.flux_func_diff, mesh.interfaces, assembler)

  elseif face_integral_type == 2
    error("face_integral_type == 2 not supported by evalJacobian")
  else
    throw(ErrorException("Unsupported face integral type = $face_integral_type"))
  end

  return nothing
end

"""
  Differentiated version of [`evalSharedFaceIntegrals`](@ref).
  Currently only supports `opts[parallel_data] == element` because it is
  used for Jacobian calculation.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * assembler
"""
function evalSharedFaceIntegrals_diff(mesh::AbstractDGMesh, sbp, eqn, opts,
                                      assembler::AssembleElementData)

#  println(eqn.params.f, "evaluating shared face integrals")
  face_integral_type = opts["face_integral_type"]
  eqn.assembler = assembler  # stash this here to get retrieved later
  if face_integral_type == 1

    finishExchangeData(mesh, sbp, eqn, opts, eqn.shared_data, calcSharedFaceIntegrals_element_diff)

  elseif face_integral_type == 2
    error("face_integral_type == 2 not supported")
  else
    throw(ErrorException("unsupported face integral type = $face_integral_type"))
  end

  # we don't want to keep a reference to the matrix inside the assembler
  # that would prevent it from getting freed
  eqn.assembler = NullAssembleElementData
  return nothing
end

"""
  Differentiated version of [`evalSourceTerm`](@ref).  Currently, source
  terms are not supported.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * assembler
"""
function evalSourceTerm_diff(mesh::AbstractMesh{Tmsh},
                     sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim},
                     opts, assembler::AssembleElementData) where {Tmsh, Tsol, Tres, Tdim}


  # source terms are not a function of q, so they don't affect the Jacobian

  #=
  # placeholder for multiple source term functionality (similar to how
  # boundary conditions are done)
  if opts["use_src_term"]
    error("source terms not supported by evalJacobian()")
  end
  =#
  return nothing
end  # end function


