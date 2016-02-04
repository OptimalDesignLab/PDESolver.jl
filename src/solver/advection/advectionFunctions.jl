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
                       sbp::SBPOperator, eqn::AdvectionData{Tsol, Tres, Tdim},
                       opts, t = 0.0)

#  println("entered evalAdvection")
  eqn.t = t
 
  eqn.res = fill!(eqn.res, 0.0)  # Zero eqn.res for next function evaluation
  
  evalSCResidual(mesh, sbp, eqn, eqn.alpha_x, eqn.alpha_y)

  evalBndry(mesh, sbp, eqn, eqn.alpha_x, eqn.alpha_y)
  # println("evalBndry complete")

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

function evalSCResidual{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, 
                                    eqn::AdvectionData{Tsol, Tres, Tdim}, 
                                    alpha_x::AbstractArray{Tsol, 3}, 
                                    alpha_y::AbstractArray{Tsol, 3})


#  println("----- Entered evalSCResidual -----")
  Adq_dxi = zeros(Tsol, 1, mesh.numNodesPerElement, mesh.numEl, 2)
  for i=1:mesh.numEl  # loop over elements
    for j=1:mesh.numNodesPerElement
      alpha_x = eqn.alpha_x[1, j, i]
      alpha_y = eqn.alpha_y[1, j, i]
      alpha_xi = mesh.dxidx[1, 1, j, i]*alpha_x + mesh.dxidx[1, 2, j, i]*alpha_y
      alpha_eta = mesh.dxidx[2, 1, j, i]*alpha_x + mesh.dxidx[2, 2, j, i]*alpha_y
      Adq_dxi[1,j,i,1] = alpha_xi*eqn.q[1,j,i]
      Adq_dxi[1,j,i,2] = alpha_eta*eqn.q[1,j,i]
    end
  end  # end loop over elements


  res2 = zeros(Tres, mesh.numNodesPerElement, mesh.numEl)
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
*  `alpha_x` & `alpha_y` : advection veocities in x & y directions

**Outputs**

*  None

"""->
#TODO: get rid of alpha_x and alpha_y arguments: they are not used
function evalBndry{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, 
                   sbp::SBPOperator, eqn::AdvectionData{Tsol, Tres, Tdim},
                   alpha_x::AbstractArray{Tsol, 3}, alpha_y::AbstractArray{Tsol, 3})

#  println("----- Entered evalBndry -----")
  for i=1:mesh.numBC
  #  println("computing flux for boundary condition ", i)
    functor_i = mesh.bndry_funcs[i]
    start_index = mesh.bndry_offsets[i]
    end_index = mesh.bndry_offsets[i+1]
    bndry_facenums_i = view(mesh.bndryfaces, start_index:(end_index - 1))
    bndryflux_i = view(eqn.bndryflux, :, :, start_index:(end_index - 1))
 
    # call the function that calculates the flux for this boundary condition
    # passing the functor into another function avoid type instability
    calcBoundaryFlux(mesh, sbp, eqn, functor_i, bndry_facenums_i, bndryflux_i)
  end

  boundaryintegrate!(sbp, mesh.bndryfaces, eqn.bndryflux, eqn.res)

#  println("----- Finished evalBndry -----")
  return nothing
end # end function evalBndry


function evalSRCTerm{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh},
                     sbp::SBPOperator, eqn::AdvectionData{Tsol, Tres, Tdim}, 
                     opts)


  # placeholder for multiple source term functionality (similar to how
  # boundary conditions are done)
  if opts["use_src_term"]
    applySRCTerm(mesh, sbp, eqn, opts, eqn.src_func)
  end

  return nothing
end  # end function

function applySRCTerm(mesh,sbp, eqn, opts, src_func)

  weights = sbp.w
  t = eqn.t
  for i=1:mesh.numEl
    jac_i = view(mesh.jac, :, i)
    res_i = view(eqn.res, :, :, i)
    for j=1:mesh.numNodesPerElement
      coords_j = view(mesh.coords, :, j, i)
      alpha_x = eqn.alpha_x[1, j, i]
      alpha_y, = eqn.alpha_y[1, j, i]

      src_val = src_func(coords_j, alpha_x, alpha_y, t)
      res_i[j] += weights[j]*src_val/jac_i[j]
    end
  end

  return nothing
end



@doc """
### AdvectionEquationMod.init
"""->
function init{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, 
              eqn::AbstractAdvectionData{Tsol, Tres}, opts)

  println("Entering Advection Module")
  getBCFunctors(mesh, sbp, eqn, opts)
  getSRCFunctors(mesh, sbp, eqn, opts)
  fill!(eqn.alpha_x, 1.0) # advection velocity in x direction
  fill!(eqn.alpha_y, 0.0) # advection velocity in y direction
  
  return nothing
end


function majorIterationCallback(itr::Integer, mesh::AbstractMesh, 
                                sbp::SBPOperator, eqn::AbstractAdvectionData, opts)

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


function assembleArray{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh}, 
                         sbp::SBPOperator, eqn::AbstractAdvectionData{Tsol}, opts, 
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



