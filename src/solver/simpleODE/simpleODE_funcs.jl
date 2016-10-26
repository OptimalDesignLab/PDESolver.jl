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

  println("===== (in physics) ===== t = $t")

  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      for k=1:Tdim  # loop over parametric dimensions

        # TODO How to get x from mesh -> dxidx?

        # TODO
#         eqn.q[1, j, i] = x^2 + t^2

#         println("===== setting , mesh.coords[:, j, i]

        eqn.res[1, j, i] = calc_4t3(mesh.coords[:, j, i], eqn.params, t)
#         eqn.res[1, j, i] = calc_3t2(mesh.coords[:, j, i], eqn.params, t)
#         eqn.res[1, j, i] = calc_2t(mesh.coords[:, j, i], eqn.params, t)

        # this is for u = x^2 + t OR u = t
#         eqn.res[1, j, i] = 1.0

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

@doc """
### NonlinearSolvers.pde_pre_func

  The pre-function for solving partial differential equations with a physics
  module.  The only operation it performs is disassembling eqn.q_vec into
  eqn.q

  Inputs:
    mesh
    sbp
    eqn
    opts
"""->
function ode_pre_func(mesh, sbp, eqn, opts)
  
  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
end

# TODO: update these comments for pde -> ode

@doc """
### NonlinearSolvers.pde_post_func

  The post-function for solving partial differential equations with a physics
  module.  This function multiplies by A0inv, assembles eqn.res into
  eqn.res_vec, multiplies by the inverse mass matrix, and calculates
  the SBP approximation to the integral L2 norm

  Inputs:
    mesh
    sbp
    eqn
    opts

"""->
function ode_post_func(mesh, sbp, eqn, opts; calc_norm=true)
#   eqn.multiplyA0inv(mesh, sbp, eqn, opts, eqn.res)
  eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
#   for j=1:length(eqn.res_vec) eqn.res_vec[j] = eqn.Minv[j]*eqn.res_vec[j] end
  if calc_norm
#     local_norm = calcNorm(eqn, eqn.res_vec)

    # since simpleODE doesn't use SBP, need to use standard 2-norm instead of SBP norm
    local_norm = norm(eqn.res_vec, 2)/mesh.numDof
    eqn.params.time.t_allreduce += @elapsed global_norm = MPI.Allreduce(local_norm*local_norm, MPI.SUM, mesh.comm)
    return sqrt(global_norm)
  end

  return nothing
end



