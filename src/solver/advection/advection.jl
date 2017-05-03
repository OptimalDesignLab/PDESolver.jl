# advection.jl
import PDESolver.evalResidual

@doc """
### AdvectionEquationMod.evalResidual

This function evaluates the Advection equation.

**Inputs**

*  `mesh` : Abstract mesh object
*  `sbp`  : Summation-by-parts operator
*  `eqn`  : Advection equation object
*  `opts` : Options dictionary
*  `t`    :

Effectively updates eqn.res -- not eqn.res_vec. To make them consistent, use assembleSolution on eqn.res and eqn.res_vec

**Outputs**

*  None

"""->
function evalResidual{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, 
                       sbp::AbstractSBP, eqn::AdvectionData{Tsol, Tres, Tdim},
                       opts::Dict, t=0.0)

  myrank = mesh.myrank
  params = eqn.params
  @debug1 println(params.f, "-----entered evalResidual -----")
  @debug1 printbacktrace(params.f)
  #f = open("pfout_$myrank.dat", "a+")
  #println(f, "----- entered evalResidual -----")
  #close(f)

  eqn.t = t
#  params.time.t_barriers[1] += @elapsed MPI.Barrier(mesh.comm) 
  eqn.res = fill!(eqn.res, 0.0)  # Zero eqn.res for next function evaluation

  @debug1 println(params.f, "== parallel_type: ", opts["parallel_type"])

  # start communication right away
  params.time.t_send += @elapsed if opts["parallel_type"] == 1

    startDataExchange(mesh, opts, eqn.q, eqn.q_face_send, eqn.q_face_recv, 
                      params.f, wait=true)
    @debug1 println(params.f, "-----entered if statement around startDataExchange -----")
    #  println("send parallel data @time printed above")
  end

  params.time.t_volume += @elapsed evalVolumeIntegrals(mesh, sbp, eqn)
#  println("evalVolumeIntegrals @time printed above")

#  params.time.t_barriers[2] += @elapsed MPI.Barrier(mesh.comm) 
  params.time.t_face += @elapsed if mesh.isDG
    evalFaceIntegrals(mesh, sbp, eqn, opts)
  end
#  println("evalFaceIntegrals @time printed above")

#  params.time.t_barriers[3] += @elapsed MPI.Barrier(mesh.comm) 
  params.time.t_source += @elapsed evalSRCTerm(mesh, sbp, eqn, opts)
#  println("evalSRCTerm @time printed above")

#  params.time.t_barriers[4] += @elapsed MPI.Barrier(mesh.comm) 
  params.time.t_bndry += @elapsed evalBoundaryIntegrals(mesh, sbp, eqn)
#  println("evalBoundaryIntegrals @time printed above")

#  params.time.t_barriers[5] += @elapsed MPI.Barrier(mesh.comm) 

  if opts["use_GLS2"]
    applyGLS2(mesh, sbp, eqn, opts, eqn.src_func)
  end
#  println("applyGLS2 @time printed above")


#  params.time.t_barriers[6] += @elapsed MPI.Barrier(mesh.comm) 
  # do parallel computation last
  params.time.t_sharedface += @elapsed if mesh.commsize > 1
    evalSharedFaceIntegrals(mesh, sbp, eqn, opts)
  end

#  params.time.t_barriers[7] += @elapsed MPI.Barrier(mesh.comm) 
#=
  f = open("pfout_$myrank.dat", "a+")
  println(f, "----- finished evalResidual -----")
  close(f)
=#

  if opts["use_Minv"]
    applyMassMatrixInverse3D(mesh, sbp, eqn, opts, eqn.res)
  end

  @debug1 flush(params.f)
#  println(params.f, "----- finished evalResidual -----")
  return nothing
end

@doc """
### AdvectionEquationMod.evalVolumeIntegrals

Evaluate the residual using summation by parts (not including boundary 
integrals) this only works for triangular meshes, where are elements are same

**Inputs**

*  `mesh` : mesh type
*  `sbp`  : Summation-by-parts operator
*  `eqn`  : Advection equation object

**Outputs**

*  None

"""->
function evalVolumeIntegrals{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh},
                                           sbp::AbstractSBP,
                                           eqn::AdvectionData{Tsol, Tres, Tdim})

  # storing flux_parametric in eqn, rather than re-allocating it every time
  flux_parametric = eqn.flux_parametric

  alphas_xy = zeros(Float64, Tdim)      # advection coefficients in the xy directions
  alphas_param = zeros(Tmsh, Tdim)      # advection coefficients in the parametric directions
  dxidx = mesh.dxidx                    # metric term
  q = eqn.q
  alphas_xy[1] = eqn.params.alpha_x
  alphas_xy[2] = eqn.params.alpha_y
  if Tdim == 3
    alphas_xy[3] = eqn.params.alpha_z
  end
  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      q_val = q[1, j, i]
      for k=1:Tdim  # loop over parametric dimensions
        alpha_k = zero(Tmsh)
        for p=1:Tdim  # sum up alpha in the current parametric dimension
          alpha_k += dxidx[k, p, j, i]*alphas_xy[p]
        end
        # the first index is the 'equation number'; advection has only one equation so it's always 1
        # for a vector PDE, there would be more
        flux_parametric[1,j,i,k] = alpha_k*q_val
      end
    end
  end

  # for each dimension, grabbing everything in the mesh and applying weakdifferentiate!
  for i=1:Tdim
    # multiplies flux_parametric by the SBP Q matrix (transposed), stores result in res
    # i is parametric direction
    weakdifferentiate!(sbp, i, sview(flux_parametric, :, :, :, i), eqn.res, trans=true)
  end

  return nothing
end


@doc """
### AdvectionEquationMod.evalBoundaryIntegrals

Evaluate boundary integrals for advection equation

**Inputs**

*  `mesh` : Abstract mesh type
*  `sbp`  : Summation-by-parts operator
*  `eqn`  : Advection equation object

**Outputs**

*  None

"""->
function evalBoundaryIntegrals{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, 
                   sbp::AbstractSBP, eqn::AdvectionData{Tsol, Tres, Tdim})

#  println("----- Entered evalBoundaryIntegrals -----")

  if mesh.isDG
    boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
  end
#  println("    boundaryinterpolate @time printed above")

  for i=1:mesh.numBC
    functor_i = mesh.bndry_funcs[i]
    start_index = mesh.bndry_offsets[i]
    end_index = mesh.bndry_offsets[i+1]
    idx_range_i = start_index:(end_index-1)
    bndry_facenums_i = sview(mesh.bndryfaces, idx_range_i)
    bndryflux_i = sview(eqn.bndryflux, :, :, idx_range_i)
 
    # call the function that calculates the flux for this boundary condition
    # passing the functor into another function avoid type instability
   calcBoundaryFlux(mesh, sbp, eqn, functor_i, idx_range_i, bndry_facenums_i, bndryflux_i)
#   println("    calcBoundaryflux @time printed above")
  end

  if mesh.isDG
    boundaryintegrate!(mesh.sbpface, mesh.bndryfaces, eqn.bndryflux, eqn.res)
  else
    boundaryintegrate!(mesh.sbpface, mesh.bndryfaces, eqn.bndryflux, eqn.res)
  end
#  println("    boundaryintegrate! @time printed above")

#  println("----- Finished evalBoundaryIntegrals -----")
  return nothing
end # end function evalBoundaryIntegrals

@doc """
### AdvectionEquationMod.evalFaceIntegrals

  This function evaluates the interior face integrals for DG methods, using
  the flux function from eqn.flux_func.  The solution variables are interpolated
  to the faces, the flux computed, and then interpolated back to the
  solution points.

  Inputs:
    mesh:  an AbstractDGMesh
    sbp
    eqn
    opts

"""->
function evalFaceIntegrals(mesh::AbstractDGMesh, sbp::AbstractSBP, eqn::AdvectionData,
                      opts)

#  println("----- Entered evalFaceIntegrals -----")
  # interpolate solution to faces
  interiorfaceinterpolate!(mesh.sbpface, mesh.interfaces, eqn.q, eqn.q_face)
#  println("    interiorface interpolate @time printed above")

  myrank = mesh.myrank
  if opts["writeqface"]
    writedlm("qface_$myrank.dat", eqn.q_face)
  end

  # calculate face fluxes
  calcFaceFlux(mesh, sbp, eqn, eqn.flux_func, mesh.interfaces, eqn.flux_face)
#  println("    calcFaceFlux @time printed above")

  if opts["write_fluxface"]
    writedlm("fluxface_$myrank.dat", eqn.flux_face)
  end

  # integrate and interpolate back to solution points
  if mesh.isDG
    interiorfaceintegrate!(mesh.sbpface, mesh.interfaces, eqn.flux_face, eqn.res)
  else
    error("cannot evalFaceIntegrals for non DG mesh")
  end
#  println("    interiorfaceintegrate @time printed above")

#  println("----- Finished evalFaceIntegrals -----")
  return nothing
end

@doc """
### AdvectionEquationMod.evalSRCTerm

  This function performs all the actions necessary to update eqn.res
  with the source term.  The source term is stored in eqn.src_func.  It is
  an abstract field, so it cannot be accessed (performantly) direction, so
  it is passed to an inner function.

  Inputs:
    mesh : Abstract mesh type
    sbp  : Summation-by-parts operator
    eqn  : Advection equation object
    opts : options dictonary

  Outputs: none

  Aliasing restrictions: none

"""->
function evalSRCTerm{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh},
                     sbp::AbstractSBP, eqn::AdvectionData{Tsol, Tres, Tdim}, 
                     opts)


  # placeholder for multiple source term functionality (similar to how
  # boundary conditions are done)
  if opts["use_src_term"]
    applySRCTerm(mesh, sbp, eqn, opts, eqn.src_func)
  end

  return nothing
end  # end function

@doc """
### AdvectionEquationMod.applySRCTerm

  This function updates eqn.res with the source term.  

  Inputs: 
    mesh
    sbp
    eqn
    opts
    src_func:  the functor that returns the value of the source term at a node
               This functor must have the signature:
               src_func(coords, params, t)
               where coords is a vector of length 2 containing the x and y 
               coordinates of the node, alpha_x and alpha_y are the advection
               coefficients, and t is the current time.

  Outputs: none

  Aliasing restrictions: none

"""->
function applySRCTerm(mesh,sbp, eqn, opts, src_func)

  weights = sbp.w
  alpha_x = eqn.params.alpha_x
  alpha_y = eqn.params.alpha_y

  t = eqn.t
  for i=1:mesh.numEl
    jac_i = sview(mesh.jac, :, i)
    res_i = sview(eqn.res, :, :, i)
    for j=1:mesh.numNodesPerElement
      coords_j = sview(mesh.coords, :, j, i)
      src_val = src_func(coords_j, eqn.params, t)
      res_i[j] += weights[j]*src_val/jac_i[j]
    end
  end

  return nothing
end


@doc """
### AdvectionEquationMod.evalSharedFaceIntegrals

  This function does the computation that needs the parallel
  communication to have finished already, namely the face integrals
  for the shared faces

  Inputs:
    mesh
    sbp
    eqn
    opts
"""->
function evalSharedFaceIntegrals(mesh::AbstractDGMesh, sbp, eqn, opts)

  if opts["face_integral_type"] != 1
    throw(ErrorException("unsupported face integral type"))
  end

  if opts["parallel_data"] == "face"
#    println(eqn.params.f, "doing face integrals using face data")
    calcSharedFaceIntegrals(mesh, sbp, eqn, opts, eqn.flux_func)
  elseif opts["parallel_data"] == "element"
#    println(eqn.params.f, "doing face integrals using element data")
    calcSharedFaceIntegrals_element(mesh, sbp, eqn, opts, eqn.flux_func)
  else
    throw(ErrorException("unsupported parallel_data setting"))
  end

  return nothing
end

@doc """
### AdvectionEquationMod.init

  This function initializes the mesh and equation objects with any module
  specific data, such as boundary condition and source term functors.

  Inputs:
    mesh
    sbp
    eqn
    opts

  Outputs: none

  Aliasing restrictions: none
"""->
function init{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP, 
              eqn::AbstractAdvectionData{Tsol, Tres}, opts)

  println("Entering Advection Module")
  getBCFunctors(mesh, sbp, eqn, opts)
  getSRCFunctors(mesh, sbp, eqn, opts)
  if mesh.isDG
    getFluxFunctors(mesh, sbp, eqn, opts)
  end
  initMPIStructures(mesh, opts)
  return nothing
end

@doc """
### AdvectionEquationMod.majorIterationCallback

  This function gets called by the NonlinearSolvers during every major 
  iteration.  This provides the opportunity to do output or monitoring.

  Inputs:
    itr: major iteration number
    mesh
    sbp
    eqn
    opts

    Outputs: none

    Aliasing restrictions: none

"""->
function majorIterationCallback(itr::Integer, mesh::AbstractMesh, 
                                sbp::AbstractSBP, eqn::AbstractAdvectionData, opts, f::IO)
#=
  if mesh.myrank == 0
    println("Performing major Iteration Callback for iteration ", itr)
  end
=#
  if opts["write_vis"] && (((itr % opts["output_freq"])) == 0 || itr == 1)
    vals = real(eqn.q_vec)  # remove unneded imaginary part
    saveSolutionToMesh(mesh, vals)
#    cd("./SolutionFiles")
    fname = string("solution_", itr)
    writeVisFiles(mesh, fname)
    println(f, "finished writing vis file"); flush(f)

    #=
    # DEBUGGING: write error to file
    q_exact = zeros(eqn.q_vec)
    ex_func = ICDict[opts["IC_name"]]
    ex_func(mesh, sbp, eqn, opts, q_exact)
    q_err = real(eqn.q_vec) - q_exact
    saveSolutionToMesh(mesh, q_err)
    fname = string("error_", itr)
    writeVisFiles(mesh, fname)
    =#

#    cd("../")
  end
 
  return nothing
end

@doc """
### AdvectionEquationMod.assembleArray

  This function performs an assignment reduction of a 3D array to a vector.
  Note that because this is an assignment reduction, the order in which 
  3D array is read matters, because only the last value assigned to a location 
  in a vector remains.

  In most cases, what you really want is assembleSolution().

  Inputs:
    mesh
    sbp
    eqn
    opts
    arr: the 3D array to be reduced into a vector

  Inputs/Outputs:
    res_vec: the vector to reduce the array into

  Keywords:
    zeros_resvec: whether or not to zero res_vec before performing the
                  reduction, default true

   Aliasing restrictions: none

"""->
#TODO: replace this with arrToVecAssign?
function assembleArray{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh}, 
                         sbp::AbstractSBP, eqn::AbstractAdvectionData{Tsol}, opts, 
                         arr::Abstract3DArray, res_vec::AbstractArray{Tres,1}, 
                         zero_resvec=true)
# arr is the array to be assembled into res_vec

#  println("in assembleSolution")

  if zero_resvec
    fill!(res_vec, 0.0)
  end


  for i=1:mesh.numEl  # loop over elements
    for j=1:mesh.numNodesPerElement
      for k=size(arr, 1)  # loop over dofs on the node

        dofnum_k = mesh.dofs[k, j, i]

        res_vec[dofnum_k] = arr[k,j,i]
      end
    end
  end
  
  return nothing
end

@doc """
### AdvectionEquationMod.applyMassMatrixInverse3D

  This function applies the 3D inverse mass matrix to an array. 
    The array passed in should always be eqn.res

  Inputs:
    mesh: mesh object, needed for numEl and numDofPerNode fields
    sbp: sbp object, needed for numnodes field
    eqn: equation object, needed for Minv3D field
    opts
    arr: the 3D array to have the 3D mass matrix inverse applied to it

"""->
function applyMassMatrixInverse3D(mesh, sbp, eqn, opts, arr)

  for i = 1:mesh.numEl
    for j = 1:sbp.numnodes
      for k = 1:mesh.numDofPerNode
        arr[k, j, i] = eqn.Minv3D[k, j, i] * arr[k, j, i]
      end
    end
  end

  return arr
end


# functions needed to make it compatible with the NonLinearSolvers module
function matVecA0inv{Tmsh, Tsol, Tdim, Tres}(mesh::AbstractMesh{Tmsh}, 
                     sbp::AbstractSBP, eqn::AdvectionData{Tsol, Tres, Tdim},
                     opts, res_arr::AbstractArray{Tsol, 3})

  return nothing
end

function matVecA0{Tmsh, Tsol, Tdim, Tres}(mesh::AbstractMesh{Tmsh},
                  sbp::AbstractSBP, eqn::AdvectionData{Tsol, Tres, Tdim}, opts,
                  res_arr::AbstractArray{Tsol, 3})

  return nothing
end


