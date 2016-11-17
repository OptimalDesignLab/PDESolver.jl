# residual_evaluation.jl: some functions used by NonlinearSolvers relating to 
# residual evaluation

# TODO: calcResidual shouldn't need exporting
export calcResidual, physicsRhs
import Utils.disassembleSolution
@doc """
### NonlinearSolvers.calcResidual

  This function takes the eqn object with the solution varaibles stored in
    q_vec, scatters them into q, evaluates the residual, and then gathers the 
    residual values into res_vec.

    Effectively, this is a wrapper around the physics module residual evaluation
      function (which performs eqn.q -> eqn.res) that performs eqn.q_vec ->
      eqn.res_vec.

    The norm of the strong residual (using the SBP norm) is calculated and returned
    (even though the weak form residual is stores in eq.res_vec).

    Inputs:
      mesh:   an AbstractMesh object
      sbp:    an SBP operator
      eqn:    an AbstractSolutionData object
      opts:   options dictonary
      func_rhs: residual evaluation function
      rhs_vec:  rhs vector
      ctx_residual:  context tuple for func_rhs
      t:      simulation time

    Outputs:
      rhs_norm_global:  norm of residual
      

    To solve steady problems:
      func_rhs should be physicsRhs
      ctx_residual should be (func, ...) where func is the physics function such as evalEuler

    To solve unsteady problems:
      func_rhs should be some function that typically uses physicsRhs as a building block for the time-marching residual
      ctx_residual should be whatever is required by func_rhs

    Note about eqn.q & eqn.q_vec consistency:
      Within this function, when func_rhs is called, eqn.q and eqn.q_vec are consistent.
      At exit from func_rhs, rhs_vec will have the residual in it.

    Aliasing restrictions: none
"""->
function calcResidual(mesh, sbp, eqn, opts, func_rhs, rhs_vec, ctx_residual, t=0.0)
# calculate the residual and its norm

  disassembleSolution(mesh, sbp, eqn, opts, eqn.q_vec)
  time = eqn.params.time
  time.t_send += @elapsed if opts["parallel_type"] == 2 && mesh.npeers > 0
    exchangeElementData(mesh, opts, eqn.q, eqn.q_face_send, eqn.q_face_recv, eqn.params.f)
  end

  func_rhs(mesh, sbp, eqn, opts, rhs_vec, ctx_residual, t)

  rhs_0_norm = calcNorm(eqn, rhs_vec, strongres=true)
  time.t_allreduce += @elapsed rhs_norm_global = MPI.Allreduce(rhs_0_norm*rhs_0_norm, MPI.SUM, mesh.comm)

  return sqrt(rhs_norm_global)

end

@doc """
### NonlinearSolvers.calcResidual

  This method of the function is used for calculating the residual of the physics.
  It is done using the state stored in eqn.q_vec .

"""->
function calcResidual(mesh, sbp, eqn, opts, func_physics, t=0.0)

  # TODO: this equals... deepcopy?
  rhs_vec = eqn.res_vec
  ctx_residual = (func_physics, )

  res_0_norm = calcResidual(mesh, sbp, eqn, opts, physicsRhs, rhs_vec, ctx_residual, t)

  return res_0_norm

end

@doc """
###NonlinearSolver.physicsRhs
  
  RHS (of the physics) calculation, separate from the Newton function

  ctx_residual: func must be the first element. func is the residual evaluation function, i.e. evalEuler
  t:            simulation time

"""->
function physicsRhs(mesh, sbp, eqn, opts, rhs_vec, ctx_residual, t=0.0)

  func = ctx_residual[1]

  func(mesh, sbp, eqn, opts, t)

  assembleResidual(mesh, sbp, eqn, opts, rhs_vec, assemble_edgeres=false)

end


@doc """
### NonlinearSolvers.assembleResidual

  This function takes the residual in eqn.res and performs an additive reduction
  into res_vec.
  This function wraps assembleSolution.

  Inputs:
    mesh: AbstractMesh
    sbp:  SBP operator
    eqn:  AbstractEquation object
    opts: options dictionary
    res_vec: residual vector to put the residual into

  Outputs:
    none

  Aliasing restrictions: none

"""->
function assembleResidual{T}(mesh, sbp, eqn, opts, res_vec::AbstractArray{T, 1}; 
                             assemble_edgeres=true, zero_resvec=true)
# assembles all of the residuals into res_vec
# no aliaising concerns

  eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, res_vec)

  if assemble_edgeres

    for i=1:size(eqn.res_edge, 4)
      eqn.assembleSolution(mesh, sbp, eqn, opts, sview(eqn.res_edge, :, :, :, i),
                           res_vec, zero_resvec=zero_resvec)
    end
  end

  return nothing
end


@doc """
### NonlinearSolvers.disassembleSolution

  This function performs the scatter q_vec -> eqn.q

  Inputs:
    mesh: AbstractMesh
    sbp:  SBP operator
    eqn:  AbstractEquation object
    opts: options dictionary
    q_vec: vector containing solution variables 

  Outputs:
    none

  Aliasing Restrictions: none

"""->
function disassembleSolution{T}(mesh, sbp, eqn, opts, q_vec::AbstractArray{T, 1})
# scatters the q_vec to the 3d array eqn.q
# no aliasing concerns here
  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, q_vec)

  return nothing
end
