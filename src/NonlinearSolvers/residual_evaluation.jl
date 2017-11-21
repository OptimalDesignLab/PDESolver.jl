# residual_evaluation.jl: some functions used by NonlinearSolvers relating to 
# residual evaluation

# TODO: calcResidual shouldn't need exporting
export calcResidual, physicsRhs
import Utils.disassembleSolution
@doc """
### NonlinearSolvers.calcResidual

    This function takes the eqn object with the solution variables stored in
    q_vec, scatters them into q, evaluates the residual, and then gathers the 
    residual values into res_vec.

    This function has 2 methods to support 2 different use cases.
    The first use case is to wrap the evalResidual function defined by each
    physics module.  evalResidual performs eqn.q -> eqn.res.  calcResidual
    performs eqn.q_vec -> eqn.res_vec.  This is accomplished by the method:

    calcResidual(mesh, sbp, eqn, opts, func_physics, t=0.0)

    with Inputs: 
      mesh: an AbstractMesh object
      sbp: an SBP operator
      eqn: an AbstractSolutionData object
      opts: options dictonary
      func_physics: the function that evaluates the residual of the physics,
                    typically evalResidual


    Recall that calcResidual evaluates func_physics as the state eqn.q_vec
    (not eqn.q), and stores the residual in eqn.res_vec

    The second method of calcResidual is a bit more involved.  Its purpose
    is calculate the residual vector for functions which are not necessarily
    a physics function, but internally call a physics function.  This is useful,
    for example, for implicit time marching methods.  This method has the
    signature:

    calcResidual(mesh, sbp, eqn, opts, func_rhs, rhs_vec, tx_residual, t=0.0)

    with Inputs:
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
      

    For this method, func_rhs is the function that computes the residuali vector
    from eqn.q (not eqn.q_vec).  It has the signature

    func_rhs(mesh, sbp, eqn, opts, rhs_vec, ctx_residual, t)
    
    rhs_vec is used to store the residual.  rhs_vec can be the same as 
    eqn.res_vec if func_rhs supports this.
    ctx_residual is a way to pass any additional arguments needed by func_rhs.
    It is a black box as far as calcResidual is concerned.

    The norm of the strong residual (using the SBP norm) is calculated and returned
    (even though the weak form residual is stores in eq.res_vec).


    As a brief example, the first method of calcResidual is actually implemented using
    the second method.  In this case:
      func_rhs should be physicsRhs
      ctx_residual should be (func,) where func is the physics function such as evalEuler

    Note about eqn.q & eqn.q_vec consistency:
      Within this function, when func_rhs is called, eqn.q and eqn.q_vec are consistent.
      At exit from func_rhs, rhs_vec will have the residual in it.  The contents of
      eqn.res and eqn.res_vec are undefind.

    Aliasing restrictions: none
"""->
function calcResidual(mesh, sbp, eqn, opts, func_physics, t=0.0)

  # TODO: this equals... deepcopy?
  rhs_vec = eqn.res_vec
  ctx_residual = (func_physics, )

  res_0_norm = calcResidual(mesh, sbp, eqn, opts, physicsRhs, rhs_vec, ctx_residual, t)

  return res_0_norm

end

"""
  The second method of calcResidual.  See documentation of the first method.
"""
function calcResidual(mesh, sbp, eqn, opts, func_rhs, rhs_vec, ctx_residual, t=0.0)
# calculate the residual and its norm

  disassembleSolution(mesh, sbp, eqn, opts, eqn.q_vec)
  time = eqn.params.time
  time.t_send += @elapsed if opts["parallel_type"] == 2 && mesh.npeers > 0
    startSolutionExchange(mesh, sbp, eqn, opts)
  end

  func_rhs(mesh, sbp, eqn, opts, rhs_vec, ctx_residual, t)

  # TODO: have func_rhs compute res (the 3D array), and call assembleResidual here?
  #       It is much more symmetric to have func_rhs compute eqn.q -> eqn.res
  #       rather than eqn.q -> eqn.res_vec
  rhs_0_norm = calcNorm(eqn, rhs_vec, strongres=true)

  return rhs_0_norm

end



@doc """
###NonlinearSolver.physicsRhs

  This is a stub function used by the first method of calcResidual to make
  evaluating the physics residual look like evaluating some arbitrary rhs_func.

  ctx_residual: func must be the first element. func is the residual evaluation
                function, i.e. evalEuler
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
# no aliasing concerns

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
