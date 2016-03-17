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

#  println("----- entered evalAdvection -----")
  eqn.t = t
 
  eqn.res = fill!(eqn.res, 0.0)  # Zero eqn.res for next function evaluation
  
  evalSCResidual(mesh, sbp, eqn)
  println("after volume integrals, res = \n", eqn.res)

  # Does not work, should remove
#  if opts["use_GLS"]
#    GLS(mesh, sbp, eqn)
#  end
  evalSRCTerm(mesh, sbp, eqn, opts)
  println("\nafter source term, res = \n", eqn.res)

  evalBndry(mesh, sbp, eqn)
  println("\nafter boundary integrals, res = \n", eqn.res)

  if mesh.isDG
    evalFaceTerm(mesh, sbp, eqn, opts)
    println("\nafter face integrals, res = \n", eqn.res)
  end



  if opts["use_GLS2"]
    applyGLS2(mesh, sbp, eqn, opts, eqn.src_func)
  end

#  println("----- finished evalAdvection -----")
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
  Adq_dxi = zeros(Tsol, 1, mesh.numNodesPerElement, mesh.numEl, 2)
  for i=1:mesh.numEl  # loop over elements
    for j=1:mesh.numNodesPerElement
      alpha_x = eqn.alpha_x
      alpha_y = eqn.alpha_y
      alpha_xi = mesh.dxidx[1, 1, j, i]*alpha_x + 
                 mesh.dxidx[1, 2, j, i]*alpha_y
      alpha_eta = mesh.dxidx[2, 1, j, i]*alpha_x + 
                  mesh.dxidx[2, 2, j, i]*alpha_y
      Adq_dxi[1,j,i,1] = alpha_xi*eqn.q[1,j,i]
      Adq_dxi[1,j,i,2] = alpha_eta*eqn.q[1,j,i]
    end
  end  # end loop over elements

  for i = 1:Tdim
    weakdifferentiate!(sbp,i,view(Adq_dxi,:,:,:,i), eqn.res, trans = true)
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

  for i=1:mesh.numBC
    functor_i = mesh.bndry_funcs[i]
    start_index = mesh.bndry_offsets[i]
    end_index = mesh.bndry_offsets[i+1]
    idx_range_i = start_index:(end_index-1)
    bndry_facenums_i = view(mesh.bndryfaces, idx_range_i)
    bndryflux_i = view(eqn.bndryflux, :, :, idx_range_i)
 
    # call the function that calculates the flux for this boundary condition
    # passing the functor into another function avoid type instability
  calcBoundaryFlux(mesh, sbp, eqn, functor_i, idx_range_i, bndry_facenums_i, bndryflux_i)
  end

  if mesh.isDG
    boundaryintegrate!(mesh.sbpface, mesh.bndryfaces, eqn.bndryflux, eqn.res)
  else
    boundaryintegrate!(sbp, mesh.bndryfaces, eqn.bndryflux, eqn.res)
  end

#  println("----- Finished evalBndry -----")
  return nothing
end # end function evalBndry

function evalFaceTerm(mesh::AbstractMesh, sbp::AbstractSBP, eqn::AdvectionData,
                      opts)

#  println("----- Entered evalFaceTerm -----")
  # interpolate solution to faces
  interiorfaceinterpolate!(mesh.sbpface, mesh.interfaces, eqn.q, eqn.q_face)

  if opts["writeqface"]
    writedlm("qface.dat", eqn.q_face)
  end

  # calculate face fluxes
  calcFaceFlux(mesh, sbp, eqn, eqn.flux_func, mesh.interfaces, eqn.flux_face)

  if opts["write_fluxface"]
    writedlm("fluxface.dat", eqn.flux_face)
  end

  # integrate and interpolate back to solution points
  if mesh.isDG
    #TODO: undo the reshaping once SBP is fixed
    flux_face_reshape = reshape(eqn.flux_face, mesh.sbpface.numnodes, mesh.numInterfaces)
    res_reshape = reshape(eqn.res, mesh.numNodesPerElement, mesh.numEl)
    interiorfaceintegrate!(mesh.sbpface, mesh.interfaces, flux_face_reshape, res_reshape)
  else
    error("cannot evalFaceTerm for non DG mesh")
  end

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
  t = eqn.t
  for i=1:mesh.numEl
    jac_i = view(mesh.jac, :, i)
    res_i = view(eqn.res, :, :, i)
    for j=1:mesh.numNodesPerElement
      coords_j = view(mesh.coords, :, j, i)
      alpha_x = eqn.alpha_x
      alpha_y = eqn.alpha_y

      src_val = src_func(coords_j, alpha_x, alpha_y, t)
      res_i[j] += weights[j]*src_val/jac_i[j]
    end
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
                                sbp::AbstractSBP, eqn::AbstractAdvectionData, opts)

  println("Performing major Iteration Callback for iteration ", itr)

  if opts["write_vis"] && ((itr % opts["output_freq"])) == 0 || itr == 1
    println("writing vtk file")
    vals = real(eqn.q_vec)  # remove unneded imaginary part
    saveSolutionToMesh(mesh, vals)
#    cd("./SolutionFiles")
    fname = string("solution_", itr)
    writeVisFiles(mesh, fname)
#    cd("../")
  end
 
  return nothing
end


