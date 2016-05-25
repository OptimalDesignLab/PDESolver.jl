# advectionFunctions.jl

@doc """
### AdvectionEquationMod.evalAdvection

This function evaluates the Advection equation and preps it for the RK4 solver.
Pass this function as an input argument to the RK4 solver just like evalAdvection.

**Inputs**

*  `mesh` : Abstract mesh object
*  `sbp`  : Summation-by-parts operator
*  `eqn`  : Advection equation object
*  `opts` : Options dictionary
*  `t`    :

**Outputs**

*  None

"""->

function evalAdvection{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, 
                       sbp::AbstractSBP, eqn::AdvectionData{Tsol, Tres, Tdim},
                       opts, t = 0.0)

  myrank = mesh.myrank
  params = eqn.params
#  println(params.f, "-----entered evalAdvection -----")
  #f = open("pfout_$myrank.dat", "a+")
  #println(f, "----- entered evalAdvection -----")
  #close(f)

  eqn.t = t
#  params.time.t_barriers[1] += @elapsed MPI.Barrier(mesh.comm) 
  eqn.res = fill!(eqn.res, 0.0)  # Zero eqn.res for next function evaluation

  # start communication right away
  if opts["parallel_type"] == 1
    params.time.t_send += @elapsed if mesh.commsize > 1
      sendParallelData(mesh, sbp, eqn, opts)
    end
    #  println("send parallel data @time printed above")
  end

  params.time.t_volume += @elapsed evalSCResidual(mesh, sbp, eqn)
#  println("evalSCResidual @time printed above")

#  params.time.t_barriers[2] += @elapsed MPI.Barrier(mesh.comm) 
  params.time.t_face += @elapsed if mesh.isDG
    evalFaceTerm(mesh, sbp, eqn, opts)
  end
#  println("evalFaceTerm @time printed above")

#  params.time.t_barriers[3] += @elapsed MPI.Barrier(mesh.comm) 
  params.time.t_source += @elapsed evalSRCTerm(mesh, sbp, eqn, opts)
#  println("evalSRCTerm @time printed above")

#  params.time.t_barriers[4] += @elapsed MPI.Barrier(mesh.comm) 
  params.time.t_bndry += @elapsed evalBndry(mesh, sbp, eqn)
#  println("evalBndry @time printed above")

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
  println(f, "----- finished evalAdvection -----")
  close(f)
=#
  @debug1 flush(params.f)
#  println(params.f, "----- finished evalAdvection -----")
  return nothing
end

@doc """
### AdvectionEquationMod.evalSCResidual

Evaluate the residual using summation by parts (not including boundary 
integrals) this only works for triangular meshes, where are elements are same

**Inputs**

*  `mesh` : mesh type
*  `sbp`  : Summation-by-parts operator
*  `eqn`  : Advection equation object
*  `alpha_x` and `alpha_y` : advection velocities in x & y directions

**Outputs**

*  None

"""->
function evalSCResidual{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, 
                                    sbp::AbstractSBP, 
                                    eqn::AdvectionData{Tsol, Tres, Tdim}) 


#    println("----- Entered evalSCResidual -----")
  Adq_dxi = eqn.flux_parametric
  alpha_x = eqn.alpha_x
  alpha_y = eqn.alpha_y

#  Adq_dxi = zeros(Tsol, 1, mesh.numNodesPerElement, mesh.numEl, 2)
  for i=1:mesh.numEl  # loop over elements
    for j=1:mesh.numNodesPerElement
      alpha_xi = mesh.dxidx[1, 1, j, i]*alpha_x + 
                 mesh.dxidx[1, 2, j, i]*alpha_y
      alpha_eta = mesh.dxidx[2, 1, j, i]*alpha_x + 
                  mesh.dxidx[2, 2, j, i]*alpha_y
      Adq_dxi[1,j,i,1] = alpha_xi*eqn.q[1,j,i]
      Adq_dxi[1,j,i,2] = alpha_eta*eqn.q[1,j,i]
    end
  end  # end loop over elements

  for i = 1:Tdim
    weakdifferentiate!(sbp,i,sview(Adq_dxi,:,:,:,i), eqn.res, trans = true)
  end

  #  println("----- Finished evalSCResidual -----")
  return nothing
end  # end function

@doc """
### AdvectionEquationMod.evalBndry

Evaluate boundary integrals for advection equation

**Inputs**

*  `mesh` : Abstract mesh type
*  `sbp`  : Summation-by-parts operator
*  `eqn`  : Advection equation object

**Outputs**

*  None

"""->
#TODO: get rid of alpha_x and alpha_y arguments: they are not used
function evalBndry{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, 
                   sbp::AbstractSBP, eqn::AdvectionData{Tsol, Tres, Tdim})

#  println("----- Entered evalBndry -----")

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
    boundaryintegrate!(sbp, mesh.bndryfaces, eqn.bndryflux, eqn.res)
  end
#  println("    boundaryintegrate! @time printed above")

#  println("----- Finished evalBndry -----")
  return nothing
end # end function evalBndry

@doc """
### AdvectionEquationMod.evalFaceTerm

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
function evalFaceTerm(mesh::AbstractDGMesh, sbp::AbstractSBP, eqn::AdvectionData,
                      opts)

#  println("----- Entered evalFaceTerm -----")
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
    error("cannot evalFaceTerm for non DG mesh")
  end
#  println("    interiorfaceintegrate @time printed above")

#  println("----- Finished evalFaceTerm -----")
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
               src_func(coords, alpha_x, alpha_y, t)
               where coords is a vector of length 2 containing the x and y 
               coordinates of the node, alpha_x and alpha_y are the advection
               coefficients, and t is the current time.

  Outputs: none

  Aliasing restrictions: none

"""->
function applySRCTerm(mesh,sbp, eqn, opts, src_func)

  weights = sbp.w
  alpha_x = eqn.alpha_x
  alpha_y = eqn.alpha_y

  t = eqn.t
  for i=1:mesh.numEl
    jac_i = sview(mesh.jac, :, i)
    res_i = sview(eqn.res, :, :, i)
    for j=1:mesh.numNodesPerElement
      coords_j = sview(mesh.coords, :, j, i)
      src_val = src_func(coords_j, alpha_x, alpha_y, t)
      res_i[j] += weights[j]*src_val/jac_i[j]
    end
  end

  return nothing
end


@doc """
### AdvectionEquationMod.sendParallelData

  This function interpolates the data into the send buffer and post
  the Isends and Irecvs.  It does not wait for them to finish

  Inputs:
    mesh
    sbp
    eqn
    opts
"""->
function sendParallelData(mesh::AbstractDGMesh, sbp, eqn, opts)

  for i=1:mesh.npeers
    # interpolate
    mesh.send_waited[i] = getSendData(mesh, opts, eqn.q, mesh.bndries_local[i], eqn.q_face_send[i], mesh.send_reqs[i], mesh.send_waited[i])
  end

  exchangeFaceData(mesh, opts, eqn.q_face_send, eqn.q_face_recv)

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

  if opts["parallel_type"] == 1
#    println(eqn.params.f, "doing face integrals using face data")
    calcSharedFaceIntegrals(mesh, sbp, eqn, opts, eqn.flux_func)
  else
#    println(eqn.params.f, "doing face integrals using element data")
    calcSharedFaceIntegrals_element(mesh, sbp, eqn, opts, eqn.flux_func)
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
  eqn.alpha_x = 1.0
  eqn.alpha_y = 1.0

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

  In most cases, what you really wnat is assembleSolution().

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
