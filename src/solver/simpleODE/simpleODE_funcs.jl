# functions to make simpleODE work

@doc """
### SimpleODEMod.evalSimpleODE

This function evaluates the simple ODE equation.

** Inputs **

*  `mesh` : Abstract mesh object
*  `sbp`  : Summation-by-parts operator
*  `eqn`  : Advection equation object
*  `opts` : Options dictionary
*  `t`    :

**Outputs**

*  None

"""->
function evalSimpleODE{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, 
                       sbp::AbstractSBP, eqn::SimpleODEData{Tsol, Tres, Tdim},
                       opts, t = 0.0)

  myrank = mesh.myrank
  params = eqn.params

  # TODO: difference between eqn.params.time and eqn.t?
  eqn.t = t

  eqn.res = fill!(eqn.res, 0.0)         # zero out the residual for the next function evaluation

  # start communication right away
  if opts["parallel_type"] == 1
    params.time.t_send += @elapsed if mesh.commsize > 1
      sendParallelData(mesh, sbp, eqn, opts)
    end
  end

  params.time.t_volume += @elapsed evalSCResidual(mesh, sbp, eqn, t)
#   params.time.t_face += @elapsed if mesh.isDG
#     evalFaceTerm(mesh, sbp, eqn, opts)
#   end
#   params.time.t_source += @elapsed evalSRCTerm(mesh, sbp, eqn, opts)
#   params.time.t_bndry += @elapsed evalBndry(mesh, sbp, eqn)
#   params.time.t_sharedface += @elapsed if mesh.commsize > 1
#     evalSharedFaceIntegrals(mesh, sbp, eqn, opts)
#   end


end   # end of function evalSimpleODE

function evalSCResidual{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh},
                                           sbp::AbstractSBP,
                                           eqn::SimpleODEData{Tsol, Tres, Tdim}, t)

  dxidx = mesh.dxidx
  q = eqn.q

  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      for k=1:Tdim  # loop over parametric dimensions

        # TODO How to get x from mesh -> dxidx?

        # TODO
#         eqn.q[1, j, i] = x^2 + t^2
        eqn.q[1, j, i] = calc_x2_t2(mesh.coords[:, j, i], eqn.params, t)

      end

    end
  end

  return nothing

end

@doc """
### SimpleODEMod.init

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
              eqn::SimpleODEData{Tsol, Tres}, opts)

  println("Entering Simple ODE Module")
  return nothing

end

@doc """
### SimpleODEMod.majorIterationCallback

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
function majorIterationCallback(itr::Integer, mesh::AbstractMesh, sbp::AbstractSBP,
                                eqn::AbstractSimpleODEData, opts, f::IO)

  return nothing

end

