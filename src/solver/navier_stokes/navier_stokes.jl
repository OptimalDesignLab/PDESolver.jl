# main file for evaluating the residual

import PDESolver: updateMetricDependents, evalResidual


"""
### EulerEquationMod.init

  This function performs any operations that need to be performed before
  the first residual evaluation.
  Any operations that need to be performed at the beginning of *each*
  residual evaluation should go in dataPrep

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * pmesh (equal to `mesh` by default)

  **Keyword Arguments**

   * init_mesh: set the boundary condition functors, default true,
                should false if `eqn`is nested inside another physics
                module's [`AbstractSolutionData`](@ref) object.
"""
# high level functions
function init(mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP,
              eqn::NSData{Tsol, Tres}, opts, pmesh=mesh;
              init_mesh=true) where {Tmsh, Tsol, Tres}

  # initialize the Euler equation object.  Tell it not to use a source term,
  # because the Navier Stokes module will handle the it
  src_orig = opts["SRCname"]
  opts["SRCname"] = "SRC0"
  EulerEquationMod.init(mesh, sbp, eqn.euler_eqn, opts, init_mesh=false)
  opts["SRCName"] = src_orig

  if init_mesh
    getBCFunctors(mesh, sbp, eqn, opts)
    getBCFunctors(pmesh, sbp, eqn, opts)
  end

  getSRCFunctors(mesh, sbp, eqn, opts)
  if mesh.isDG
    getFluxFunctors(mesh, sbp, eqn, opts)
  end

  return nothing
end



function evalResidual(mesh::AbstractMesh, sbp::AbstractSBP, eqn::NSData,
                     opts::Dict, t=0.0)

  @assert eqn.commsize == 1

  params = eqn.params
  params.t_dataprep += @elapsed dataPrep(mesh, sbp, eqn, opts)

  time.t_volume += @elapsed if opts["addVolumeIntegrals"]
    evalVolumeIntegrals(mesh, sbp, eqn, opts)
#    println("volume integral @time printed above")
  end

  time.t_bndry += @elapsed if opts["addBoundaryIntegrals"]
    # do in inviscid-type boundary integral
    # Because we only support non-precompute, all this function does is
    # the integration, so we can re-use the inviscid one
    EulerEquationMod.evalBoundaryIntegrals(mesh, sbp, eqn, opts)
#   println("boundary integral @time printed above")
  end

  time.t_face += @elapsed if mesh.isDG && opts["addFaceIntegrals"]
    # invisicd face integrals
    EulerEquationMod.evalFaceIntegrals(mesh, sbp, eqn, opts)
#    println("face integral @time printed above")
  end


  time.t_source += @elapsed evalSourceTerm(mesh, sbp, eqn, opts)
#  println("source integral @time printed above")

 
  if eqn.params.isViscous == true
    params.t_face += @elapsed evalFaceIntegrals_vector(mesh, sbp, eqn, opts)
    # do the non-inviscid-type boundary integral
    params.t_bndry += @elapsed evalBoundaryIntegrals_vector(mesh, sbp, eqn, opts)
  end

  # apply inverse mass matrix to eqn.res, necessary for CN
  if opts["use_Minv"]
    EulerEquationMod.applyMassMatrixInverse3D(mesh, sbp, eqn.euler_eqn, opts, eqn.res)
  end



  return nothing
end


"""
### NavierStokesMod.dataPrep

  This function calculates all the quantities that are stored over the entire
  mesh, and performs a few sanity check such as ensuring density and pressure
  are positive.

  This function also calls `dataPrep` for Euler.

  This is a high level function.
"""
function dataPrep(mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP,
                   eqn::NSData{Tsol, Tres}, opts) where {Tmsh, Tsol, Tres}

  # do the Euler part
  # don't do the inviscid boundary conditions
  precompute_bc_orig = opts["precompute_boundary_flux"]
  EulerEquationMod.dataPrep(mesh, sbp, eqn.euler_eqn, opts)
  opts["precompute_boundary_flux"] = precompute_bc_orig

  # do the viscous part
  if eqn.params.isViscous == true
    # fill!(eqn.vecflux_face, 0.0)
    fill!(eqn.vecflux_faceL, 0.0)
    fill!(eqn.vecflux_faceR, 0.0)

    calcViscousFlux_interior(mesh, sbp, eqn, opts)

    fill!(eqn.vecflux_bndry, 0.0)
    calcViscousFlux_boundary(mesh, sbp, eqn, opts)
  end

  return nothing
end


"""
  Does the volume integrals
"""
function evalVolumeIntegrals(mesh::AbstractMesh{Tmsh},
                             sbp::AbstractSBP, eqn::NSData{Tsol, Tres, Tdim},
                             opts) where {Tmsh,  Tsol, Tres, Tdim}

  EulerEquationMod.evalVolumeIntegrals(mesh, sbp, eqn.euler_eqn, opts)

  # now differentiate the volume flux agains for the viscous part
  if eqn.params.isViscous == true
    weakdifferentiate2!(mesh, sbp, eqn, eqn.res)
  end

  return nothing
end

"""
  Updates the residual with the source term contribution.

  **Inputs**:

   * mesh : Abstract mesh type
   * sbp  : Summation-by-parts operator
   * eqn  : Euler equation object
   * opts : options dictonary

  Outputs: none

  Aliasing restrictions: none

"""
function evalSourceTerm(mesh::AbstractMesh{Tmsh},
                     sbp::AbstractSBP, eqn::NSData{Tsol, Tres, Tdim},
                     opts) where {Tmsh, Tsol, Tres, Tdim}


  if opts["use_src_term"]
    # the source term signature are compatible, so use the Euler function
    # to evaluate it
    EulerEquationMod.applySourceTerm(mesh, sbp, eqn.euler_eqn, opts, eqn.src_func)
  end

  return nothing
end  # end function


function updateMetricDependents(mesh::AbstractMesh, sbp::AbstractSBP,
                                eqn::NSData, opts)

  # all the arrays are aliased, so use the Euler version
  EulerEquationMod.updateMetricDependents(mesh, sbp, eqn.euler_eqn, opts)

  return nothing
end


function majorIterationCallback(itr::Integer,
                                mesh::AbstractMesh,
                                sbp::AbstractSBP,
                                eqn::NSData, opts, f::IO)

  # the Euler callback does useful things and doesn't have anything incorrect
  # for viscous, so use it.  Eventually Navier Stokes should have its own
  EulerEquationMod.majorIterationCallback(itr, mesh, sbp, eqn.euler_eqn, opts, f)

  return nothing
end
