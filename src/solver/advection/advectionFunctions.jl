# advectionFunctions.jl

@doc """
### AdvectionEquationMod.evalAdvection

This function evaluates the Advection equation and preps it for the RK4 solver.
Pass this function as an input argument to the RK4 solver just like evalEuler.

**Inputs**

*  `mesh` : Abstract mesh object
*  `sbp`  : Summation-by-parts operator
*  `eqn`  : Advection equation object
*  `opts` : Options dictionary
*  `t`    :

**Outputs**

*  None

"""->

function evalAdvection(mesh::AbstractMesh, sbp::SBPOperator,
                       eqn::AdvectionData, opts, t)

  const alpha_x = 1.0 # advection velocity in x direction
  const alpha_y = 1.0 # advection velocity in y direction
  
  eqn.res = fill!(eqn.res, 0.0)  # Zero eqn.res for next function evaluation
  # disassembleSolution(mesh, sbp, eqn, opts, eqn.q_vec)
  # println(eqn.u)
  
  evalSCResidual(mesh, sbp, eqn, alpha_x, alpha_y)
  # assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
  # println(eqn.res_vec)
  # println("evalSCResidual complete")
  evalBndry(mesh, sbp, eqn, alpha_x, alpha_y)
  # println("evalBndry complete")
  # assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
  # println(eqn.res_vec)
  
  eqn.res_vec = eqn.M\eqn.res_vec
  
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

function evalSCResidual{Tsol, Tdim}(mesh::AbstractMesh, sbp::SBPOperator, 
                                    eqn::AdvectionData{Tsol, Tdim}, 
                                    alpha_x::FloatingPoint, 
                                    alpha_y::FloatingPoint)

	                      
  ndof = mesh.numDof  # Total number of dofs
  numEl = mesh.numEl  # Total number of elements
  # nnodes = mesh.numNodesPerElement # count the number of nodes per element 
  ub = zeros(nnodes) # holds original solution at for an element
  fluxes = zeros(nnodes, 2)  # jacobian term times advection velocity divided
                             # by jac
  dxi_dxu = zeros(Tsol, mesh.numNodesPerElement, numEl, Tdim)
  for i=1:numEl  # loop over element
    for j=1:nnodes
      dxi_dxu[j,i,1] = (mesh.dxidx[1,1,j,i]*alpha_x + mesh.dxidx[1,2,j,i]*alpha_y)*eqn.u[1,j,i]
      dxi_dxu[j,i,2] = (mesh.dxidx[2,1,j,i]*alpha_x + mesh.dxidx[2,2,j,i]*alpha_y)*eqn.u[1,j,i]
    end
  end  # end loop over elements
  
  for i = 1:Tdim
    weakdifferentiate!(sbp,i,view(dxi_dxu,:,:,i), eqn.res, trans = true)
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

function evalBndry{Tsol, Tdim}(mesh::PumiMesh2, sbp::SBPOperator, eqn::AdvectionData{Tsol, Tdim},
                   alpha_x::FloatingPoint, alpha_y::FloatingPoint)

  # get arguments needed for sbp boundaryintegrate!

  bndry_edges = mesh.bndryfaces

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
  end # end for i = 1:mesh.numBoundaryEdges

  boundaryintegrate!(sbp, mesh.bndryfaces, eqn.bndryflux, eqn.res)
  
  return nothing
end # end function evalBndry


@doc """
### AdvectionEquationMod.init
"""->
function init{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, 
              eqn::AbstractEulerData{Tsol, Tres}, opts)

  println("Entering Advection Module")
  getBCFunctors(mesh, sbp, eqn, opts)
  
  return nothing
end