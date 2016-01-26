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
  
  # const eqn.alpha_x = 1.0 
  # const eqn.alpha_y = 1.0 

  eqn.alpha_x = fill!(eqn.alpha_x, 1.0) # advection velocity in x direction
  eqn.alpha_y = fill!(eqn.alpha_y, 1.0) # advection velocity in y direction
  
  eqn.res = fill!(eqn.res, 0.0)  # Zero eqn.res for next function evaluation
  # disassembleSolution(mesh, sbp, eqn, opts, eqn.u_vec)
  # println(eqn.u)
  
  evalSCResidual(mesh, sbp, eqn, eqn.alpha_x, eqn.alpha_y)
  # assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
  # println(eqn.res_vec)
  # println("evalSCResidual complete")
  evalBndry(mesh, sbp, eqn, eqn.alpha_x, eqn.alpha_y)
  # println("evalBndry complete")
  # assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
  # println(eqn.res_vec)
  
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

	                      
  ndof = mesh.numDof  # Total number of dofs
  numEl = mesh.numEl  # Total number of elements
  # nnodes = mesh.numNodesPerElement # count the number of nodes per element 
  ub = zeros(mesh.numNodesPerElement) # holds original solution at for an element
  fluxes = zeros(mesh.numNodesPerElement, 2)  # jacobian term times advection velocity divided
                             # by jac
  dxi_dxq = zeros(Tsol, 1, mesh.numNodesPerElement, numEl, 2) 
  for i=1:numEl  # loop over element
    for j=1:mesh.numNodesPerElement
      dxi_dxq[1,j,i,1] = (mesh.dxidx[1,1,j,i]*alpha_x[1,j,i] + 
                        mesh.dxidx[1,2,j,i]*alpha_y[1,j,i])*eqn.q[1,j,i]
      dxi_dxq[1,j,i,2] = (mesh.dxidx[2,1,j,i]*alpha_x[1,j,i] + 
                        mesh.dxidx[2,2,j,i]*alpha_y[1,j,i])*eqn.q[1,j,i]
    end
  end  # end loop over elements
  
  for i = 1:Tdim
    weakdifferentiate!(sbp,i,view(dxi_dxq,:,:,:,i), eqn.res, trans = true)
  end

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

function evalBndry{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, 
                   sbp::SBPOperator, eqn::AdvectionData{Tsol, Tres, Tdim},
                   alpha_x::AbstractArray{Tsol, 3}, alpha_y::AbstractArray{Tsol, 3})

  # get arguments needed for sbp boundaryintegrate!

  #=bndry_edges = mesh.bndryfaces

  if length(mesh.bndryfaces) != mesh.numBoundaryEdges
    println("Error with Boundary!!!!")
  end
  
  for i = 1:mesh.numBoundaryEdges
    bndry_i = mesh.bndryfaces[i]
    for j = 1:sbp.numfacenodes
      k = sbp.facenodes[j, bndry_i.face]
      u = view(eqn.u, :, k, bndry_i.element)
      x = view(mesh.coords, :, k, bndry_i.element)
      dxidx = view(mesh.dxidx, :, :, k, bndry_i.element)
      nrm = view(sbp.facenormal, :, bndry_i.face)
      bndryflux_i = view(eqn.bndryflux, :, j, i)
      flux1(u, dxidx, nrm, bndryflux_i, alpha_x, alpha_y) # calculate the boundary flux
    end # for j = 1:sbp.numfacenodes
  end # end for i = 1:mesh.numBoundaryEdges =#

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
  
  return nothing
end # end function evalBndry


@doc """
### AdvectionEquationMod.init
"""->
function init{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, 
              eqn::AbstractAdvectionData{Tsol, Tres}, opts)

  println("Entering Advection Module")
  getBCFunctors(mesh, sbp, eqn, opts)
  
  return nothing
end


function majorIterationCallback(itr::Integer, mesh::AbstractMesh, 
                                sbp::SBPOperator, eqn::AbstractAdvectionData, opts)

#  println("Performing major Iteration Callback")

  if opts["write_vis"] && ((itr % opts["output_freq"])) == 0 || itr == 1
    vals = abs(real(eqn.u_vec))  # remove unneded imaginary part
    saveSolutionToMesh(mesh, vals)
    fname = string("solution_", itr)
    writeVisFiles(mesh, fname)
  end
 
  return nothing
end