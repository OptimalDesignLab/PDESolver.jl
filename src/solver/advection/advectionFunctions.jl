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

  println("t = ", t)
  println("entered evalAdvection")
  println("eqn.q_vec = \n", eqn.q_vec)
  println("eqn.q = \n", eqn.q)
  println("centerline diff = ", eqn.q_vec[1] - eqn.q_vec[4])
#=
  for i = 1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      dofnum_j = mesh.dofs[1, j, i]
      x_j = mesh.coords[1, j, i]
      eqn.res_vec[dofnum_j] = cos(-x_j + t)
    end
  end

  # now multiply by M  - actually don't, because the residuals should only be
  # equal after multiplication by Minv
  for i=1:mesh.numDof
    eqn.res_vec[i] = eqn.res_vec[i]
  end

  println("eqn.res_vec = \n", eqn.res_vec)
=#

  # const eqn.alpha_x = 1.0 
  # const eqn.alpha_y = 1.0 

  eqn.t = t
  println("t = ", eqn.t)
  eqn.alpha_x = fill!(eqn.alpha_x, 1.0) # advection velocity in x direction
  eqn.alpha_y = fill!(eqn.alpha_y, 0.0) # advection velocity in y direction
  
  eqn.res = fill!(eqn.res, 0.0)  # Zero eqn.res for next function evaluation
  # disassembleSolution(mesh, sbp, eqn, opts, eqn.u_vec)
  # println(eqn.u)
  
  evalSCResidual(mesh, sbp, eqn, eqn.alpha_x, eqn.alpha_y)
   assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
   println("after evalSCResidual, eqn.res_vec = \n", eqn.res_vec)
   println("eqn.res = \n", eqn.res)
   println("centerline diff = ", eqn.res_vec[1] - eqn.res_vec[4])
  # println("evalSCResidual complete")
#=
  evalBndry(mesh, sbp, eqn, eqn.alpha_x, eqn.alpha_y)
  # println("evalBndry complete")

   assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
   res_scaled = eqn.Minv.*eqn.res_vec
   println("after evalBndry, scaled eqn.res_vec = \n", res_scaled)
   println("centerline diff = ", res_scaled[1] - res_scaled[4])
  # assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
  # println(eqn.res_vec)
=#

  coords = [-1.0, 0]
  u_bc = calc_sinwave(coords, t)
  assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
  for i=1:mesh.numDof
    eqn.res_vec[i] *= eqn.Minv[i]
  end

  alpha_x = 1.0
  eqn.res_vec[8] = -alpha_x*u_bc
  eqn.res_vec[7] = -alpha_x*u_bc
  eqn.res_vec[9] = -alpha_x*u_bc

  coords = [1.0, 0]
  u_bc = calc_sinwave(coords, t)
  eqn.res_vec[5] = -alpha_x*u_bc
  eqn.res_vec[3] = -alpha_x*u_bc
  eqn.res_vec[2] = -alpha_x*u_bc
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
  alpha_x = 1.0
  alpha_y = 0.0
  Adq_dxi = zeros(Tsol, 1, mesh.numNodesPerElement, mesh.numEl, 2)
  for i=1:mesh.numEl  # loop over element
    for j=1:mesh.numNodesPerElement
      alpha_xi = mesh.dxidx[1, 1, j, i]*alpha_x + mesh.dxidx[1, 2, j, i]*alpha_y
      alpha_eta = mesh.dxidx[2, 1, j, i]*alpha_x + mesh.dxidx[2, 2, j, i]*alpha_y
      Adq_dxi[1,j,i,1] = alpha_xi*eqn.q[1,j,i]
      Adq_dxi[1,j,i,2] = alpha_eta*eqn.q[1,j,i]
    end
  end  # end loop over elements

  println("Adq_dxi = \n", Adq_dxi[:, :, :, 1])

  println("Adq_deta = \n", Adq_dxi[:, :, :, 2])

  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      flux_vals = zeros(Tsol, 2)
      flux_vals[1] = Adq_dxi[1, j, i, 1]
      flux_vals[2] = Adq_dxi[1, j, i, 2]
      flux_xy = inv(mesh.dxidx[:, :, j, i])*flux_vals

      println("for element $i, node $j dxidx = \n", mesh.dxidx[:, :, j, i])
      println("for element $i, node $j  q = \n", eqn.q[:, j, i])

      println("flux for element $i node $j direction x = \n", flux_xy[1])
      println("flux for element $i node $j direction y = \n", flux_xy[2])

    end
  end

  for i = 1:Tdim
    weakdifferentiate!(sbp,i,view(Adq_dxi,:,:,:,i), eqn.res, trans = true)
  end

#  println("residual = \n", eqn.res)
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

function evalBndry{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, 
                   sbp::SBPOperator, eqn::AdvectionData{Tsol, Tres, Tdim},
                   alpha_x::AbstractArray{Tsol, 3}, alpha_y::AbstractArray{Tsol, 3})

#  println("----- Entered evalBndry -----")
#  println("bndryfaces = \n", mesh.bndryfaces)
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

#  println("bndryflux = \n", eqn.bndryflux)
  boundaryintegrate!(sbp, mesh.bndryfaces, eqn.bndryflux, eqn.res)
#  println("after adding boundary component, res = \n", eqn.res)

#  println("----- Finished evalBndry -----")
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
    println("writing vtk file")
    vals = real(eqn.q_vec)  # remove unneded imaginary part
    saveSolutionToMesh(mesh, vals)
    cd("./SolutionFiles")
    fname = string("solution_", itr)
    writeVisFiles(mesh, fname)
    cd("../")
  end
 
  return nothing
end
