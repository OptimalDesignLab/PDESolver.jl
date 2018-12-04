# These objects are used to control the recalculation of the jacobian
# and preconditioner

"""
  This is the supertype of all types that tell a method when to recalculate
  the Jacobian.  Each subtype defines a different policy for when to
  recalculate.
"""
abstract type RecalculationPolicy end

"""
  This function returns the enum specifying whether the PC and/or LO should
  be recalculated.

  Every implementation of [`RecalculationPolicy`](@ref) should extend this
  function with a new method.

  **Inputs**

   * policy: a [`RecalculationPolicy`](@ref)
   * itr: the current iteration

  **Outputs**
  
   * the num: see [`RECALC_BOTH`](@ref)
"""
function decideRecalculation(policy::RecalculationPolicy, itr::Integer)

  error("abstract decideRecalculation() reached.  Did you forget to extend this function for your RecalculationPolicy type?")
end


"""
  This function uses [`decideRecalculation`](@ref) and the [`LinearSolver`](@ref)
  interface to recalculate the PC and/or LO as specified by the policy.

  This is the main interface the different methods should use for
  recalculating the PC and/or LO.  The only case when
  [`decideRecalculation`](@ref) is better is when the method is not using
  the LinearSolvers module to define the PC and LO.

  **Inputs**

   * policy: a [`RecalculationPolicy`](@ref)
   * itr: current iteration number

  The following arguments exactly match the signature of [`calcPCandLO`](@ref)

   * ls
   * mesh
   * sbp
   * eqn
   * opts
   * ctx_residual
   * t

  When this function exits, the PC and/or LO will have been updated, if
  specified by the RecalculationPolicy.
"""
function doRecalculation(policy::RecalculationPolicy, itr::Integer,
                         ls::LinearSolver,
                         mesh::AbstractMesh, sbp::AbstractOperator,
                         eqn::AbstractSolutionData, opts,
                         ctx_residual, t)


  recalc_type = decideRecalculation(policy, itr)
  if recalc_type == RECALC_BOTH
    calcPCandLO(ls, mesh, sbp, eqn, opts, ctx_residual, t)
  elseif recalc_type == RECALC_PC
    calcPC(ls, mesh, sbp, eqn, opts, ctx_residual, t)
  elseif recalc_type == RECALC_LO
    calcLinearOperator(ls, mesh, sbp, eqn, opts, ctx_residual, t)
  end

  return nothing
end




"""
  This function should reset the RecalculationPolicy to its initial state.
  This function is called at the beginning of Newton's method (making it
  safe to reuse the RecalculationPolicy for repeated calls to Newton)

  Every [`RecalculationPolicy`](@ref) should extend function with a new
  method.

  **Inputs**

   * policy: a [`RecalculationPolicy`](@ref)

  **Outputs**

   * none
"""
function resetRecalculationPolicy(policy::RecalculationPolicy)

  error("abstract resetRecalculationPolicy() reached.  Did you forget to extend this function for your RecalculationPolicy type?")

end
"""
  Helper function to [`decideRecalculation`](@ref)

  **Inputs**

   * recalc_pc: bool specifying whether the PC should be recalculated
   * recalc_lo: bool specifying whether the LO should be recalculated

  **Outputs**

   * recalculation enum (see [`RECALC_BOTH`](@ref)

"""
function getRecalculationEnum(recalc_pc::Bool, recalc_lo::Bool)

  # TODO: move this to own function for convenience
  if recalc_pc && recalc_lo
    return RECALC_BOTH
  elseif recalc_pc
    return RECALC_PC
  elseif recalc_lo
    return RECALC_LO
  else
    return RECALC_NONE
  end

  return nothing
end


"""
  Enums returned by [`decideRecalculation`](@ref) telling the caller what to
  recalculate.

  **Enum Names**

   * RECALC_BOTH: recalculate PC and LO
   * RECALC_PC: recalculate PC only
   * RECALC_LO: recalculate LO only
   * RECALC_NONE: don't recalculate anything
"""
global const RECALC_BOTH = 1
global const RECALC_PC = 2
global const RECALC_LO = 3
global const RECALC_NONE = 4

#------------------------------------------------------------------------------
# RecalculateFixedIntervals

"""
  This [`RecalculationPolicy`](@ref) recalculates the LO and PC every `x` 
  number if iterations.

  **Options Keys**

   * "\$prefix_prec_recalc_freq": frequency of PC recalculation
   * "\$prefix_jac_recalc_freq": frequency of LO recalculation
   * "\$prefix_recalc_first": if true, recalculate the PC and LO on the first
                              iteration, even if other criteria not met.

  where "\$prefix" is the prefix passed into [`getRecalculationPolicy`](@ref)
"""
mutable struct RecalculateFixedIntervals <: RecalculationPolicy
  prec_freq::Int  # recalculate this preconditioner every this many iterations
  jac_freq::Int   # recalculate the jacobian every this many iterations
  last_prec::Int  # last iteration the prec was recalculated
  last_jac::Int    # last iteration the jac was recalculated
  recalc_first::Bool  # true if everything should be recalculated on the
                      # first iteration, even if other criteria not met
end

function RecalculateFixedIntervals(opts, prefix="")

  prec_freq = opts[string(prefix, "_prec_recalc_freq")]
  jac_freq = opts[string(prefix, "_jac_recalc_freq")]
  recalc_first = opts[string(prefix, "_recalc_first")]

  @assert prec_freq > 0
  @assert prec_freq > 0

  last_prec = 0
  last_jac = 0

  return RecalculateFixedIntervals(prec_freq, jac_freq, last_prec, last_jac, recalc_first)
end

function decideRecalculation(policy::RecalculateFixedIntervals, itr::Integer)

  recalc_pc = false
  recalc_lo = false

  if policy.recalc_first && itr == 1
    recalc_pc = true
    recalc_lo = true
  else
    recalc_pc = ( itr - policy.last_prec == policy.prec_freq )
    recalc_lo = ( itr - policy.last_jac == policy.jac_freq )
  end

  if recalc_pc
    policy.last_prec = itr
  end
  if recalc_lo
    policy.last_jac = itr
  end

  return getRecalculationEnum(recalc_pc, recalc_lo)
end

function resetRecalculationPolicy(policy::RecalculateFixedIntervals)

  policy.last_prec = 0
  policy.last_jac = 0

  return nothing
end



#------------------------------------------------------------------------------
# Recalculate Never

"""
 This [`RecalculationPolicy`](@ref) never recalculates the PC or LO.
 This is useful when Newton's method is used inside other methods to solve
 easy problems.
"""
mutable struct RecalculateNever <: RecalculationPolicy
end

function RecalculateNever(opts, prefix="")

  return RecalculateNever()
end


function decideRecalculation(policy::RecalculateNever, itr::Integer)

  recalc_pc = false
  recalc_lo = false

  return getRecalculationEnum(recalc_pc, recalc_lo)
end

function resetRecalculationPolicy(policy::RecalculateNever)

  # nothing to do here

  return nothing
end

#------------------------------------------------------------------------------
# Constructors

"""
  Maps names to [`RecalculationPolicy`](@ref) constructors.  All new
  recalculation policies must be added to this list.
  
  The constructors must have the signature:

    MyPolicyName(opts::Dict, prefix::String)

  where `opts` is the options dictionary and `prefix` is the prefix passed into
  [`getRecalculationPolicy`](@ref).
"""
global const RecalculationPolicyDict = Dict{String, Any}(
"RecalculateFixedIntervals" => RecalculateFixedIntervals,
"RecalculateNever"          => RecalculateNever,
)
# Note: the second type parameter can't be Any, because Julia seems to confuse
#       types and their outer constructors


"""
  This function constructs and returns the recalculation policy specified by
  the options dictionary.

  Methods that want to use a RecalculationPolicy should this this function
  to get it (do not call the constructors directly).

  **Inputs**

   * opts: options dictionary
   * prefix: prefix of options to use.  This is usually the name of the method
             the returned object will be used in, eg. "newton"

  **Outputs**

   * policy: a RecalculationPolicy of some kind

  **Options Keys**

   * "\$prefix_recalculation_policy": name of recalculation policy to get


  **Current Policies**

  $(keys(RecalculationPolicyDict))

"""
function getRecalculationPolicy(opts, prefix="")

  # right now there is only one, but there may be more in the future

  keyname = string(prefix, "_recalculation_policy")
  name = opts[keyname]
  return RecalculationPolicyDict[name](opts, prefix)
end


