# residual_evaluation.jl: some functions used by NonlinearSolvers relating to 
# residual evaluation

# TODO: calcResidual shouldn't need exporting
export calcResidual, physicsRhs
import Utils.disassembleSolution

"""
  This function computes the vector form of of the residual from the vector
  form of the solution, ie. q_vec -> rhs_vec, for a given physics.
  This is one of the two functions required by [`newtonInner`](@ref).

  Users of newtonInner that are not physics modules (for example, implicit
  time marching schemes) will need to implement their own version of this
  function.  See the Implimentation Notes section below.

  **Inputs**
  
   * mesh: an AbstractMesh
   * sbp: an AbstractSBP
   * eqn: an AbstractSolutionData (may be modified during this function)
   * opts: the options dictionary
   * ctx_residual: a tuple of values.  ctx_residual[1] must be a function
                   that computes (q -> res).  Typically this is evalResidual.
                   The other entries of the tuple (if any) are not used.

                   The purpose of this argument is to make the signature of
                   the function generic enough so that implict time marching
                   methods can use it.  (If you are confused about this
                   programming pattern, google how callback are implemented
                   in C, this ctx is like a void* in C).

   * t: the time at which to evalute the residual

  **Inputs/Outputs**

   * rhs_vec: vector to put the residual in

  **Output**

   * norm of the residual vector

  **Implementation Notes:**

  This function is really a wrapper around an evalResidual-like function.
  It has to do 5 things:

    1. scatter eqn.q_vec -> eqn.q
    2. start parallel communication if needed
    3. call the evalResidual-like function to compute eqn.q -> eqn.res
    4. assemble the residual into the output vector, ie. eqn.res -> rhs_vec
    5. compute the norm of the rhs_vec

  Any implementation of this function (used with newtonInner) must have the
  following properties:

    1. on exit, eqn.q and eqn.q_vec must be consistent
    2. this function start parallel communication if needed
    3. allow for the possibility that rhs_vec and eqn.res_vec alias
   
  Is is recommended to use [`calcNorm`](@ref Utils.calcNorm) to compute the norm.

  Other implementations of this function are encouraged to use this function
  to help construct their rhs_vec.
"""
function physicsRhs(mesh, sbp, eqn, opts, rhs_vec, ctx_residual, t=0.0)

  

  # q_vec -> q
  disassembleSolution(mesh, sbp, eqn, opts, eqn.q_vec)

  # start parallel communication if needed
  time = eqn.params.time
  time.t_send += @elapsed if opts["parallel_type"] == 2 && mesh.npeers > 0
    startSolutionExchange(mesh, sbp, eqn, opts)
  end

  # evaluate (q -> res)
  func = ctx_residual[1]
  func(mesh, sbp, eqn, opts, t)

  # res -> rhs_vec
  assembleResidual(mesh, sbp, eqn, opts, rhs_vec, assemble_edgeres=false)

  # calculate the norm
  rhs_0_norm = calcNorm(eqn, rhs_vec, strongres=true)

  return rhs_0_norm
end


"""
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

"""
function assembleResidual{T}(mesh, sbp, eqn, opts, res_vec::AbstractArray{T, 1}; 
                             assemble_edgeres=true, zero_resvec=true)
# assembles all of the residuals into res_vec
# no aliasing concerns
#TODO: clarify zero_resvec (why not used for both calls?)

  assembleSolution(mesh, sbp, eqn, opts, eqn.res, res_vec)

  if assemble_edgeres

    for i=1:size(eqn.res_edge, 4)
      assembleSolution(mesh, sbp, eqn, opts, sview(eqn.res_edge, :, :, :, i),
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
  disassembleSolution(mesh, sbp, eqn, opts, eqn.q, q_vec)

  return nothing
end
