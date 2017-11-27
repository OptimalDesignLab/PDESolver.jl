# auxiliary types and functions for Newton's method
# linear operator, newton data etc.

@doc """
  This type holds all the data the might be needed for Newton's method,
  including globalization.  This simplifies the data handling and 
  the C callback used by Petsc
"""->
type NewtonData{Tsol, Tres, Tsolver <: LinearSolver}

  # MPI info
  myrank::Int
  commsize::Int

  #TODO: add newton tolerances here, remove them as keyword args to newtonInner

  # working variables
  res_norm_i::Float64  # current step residual norm
  res_norm_i_1::Float64  # previous step residual norm
  step_norm_i::Float64
  step_norm_i_1::Float64
  res_norm_rel::Float64  # norm of the residual used for the relative residual
                         # tolerance
  step_fac::Float64

  # tolerances (newton)
  res_reltol::Float64
  res_abstol::Float64
  step_tol::Float64
  itermax::Int

  # inexact Newton-Krylov parameters
  krylov_gamma::Float64  # update parameter for krylov tolerance

  ls::Tsolver

end

#TODO: see if the static parameters are still needed
function NewtonData{Tsol, Tres}(mesh, sbp,  
                    eqn::AbstractSolutionData{Tsol, Tres}, opts,
                    ls::LinearSolver)

  println("entered NewtonData constructor")

  myrank = mesh.myrank
  commsize = mesh.commsize



  res_norm_i = 0.0
  res_norm_i_1 = 0.0
  step_norm_i = 0.0
  step_norm_i_1 = 0.0
  res_norm_rel = opts["res_reltol0"]
  step_fac = 1.0

  res_reltol = opts["res_reltol"]
  res_abstol = opts["res_abstol"]
  step_tol = opts["step_tol"]
  itermax = opts["itermax"]

  krylov_gamma = opts["krylov_gamma"]

  return NewtonData{Tsol, Tres, typeof(ls)}(myrank, commsize, 
                    res_norm_i, res_norm_i_1, step_norm_i, step_norm_i_1,
                    res_norm_rel, step_fac,
                    res_reltol, res_abstol, step_tol, itermax,
                    krylov_gamma, ls)
end

include("residual_evaluation.jl")  # some functions for residual evaluation
include("jacobian.jl")

@doc """
### NonlinearSolvers.setupNewton
  Performs setup work for [`newtonInner`](@ref), including creating a 
  [`NewtonData`](@ref) object.

  This function also reset the implicit Euler globalization.

  alloc_rhs: keyword arg to allocate a new object or not for rhs_vec
                true (default) allocates a new vector
                false will use eqn.res_vec

  rhs_func: only used for Petsc in matrix-free mode to do Jac-vec products
            should be the rhs_func passed into [`newtonInner`](@ref)
  ctx_residual: ctx_residual passed into [`newtonInner`](@ref)

  Allocates Jac & RHS

  See [`cleanupNewton`](@ref) to the cleanup function

"""->
function setupNewton{Tsol, Tres}(mesh, pmesh, sbp,
                     eqn::AbstractSolutionData{Tsol, Tres}, opts,
                     ls::LinearSolver; alloc_rhs=true)

  newton_data = NewtonData(mesh, sbp, eqn, opts, ls)

  clearEulerConstants()
  # For simple cases, especially for Newton's method as a steady solver,
  #   having rhs_vec and eqn.res_vec pointing to the same memory
  #   saves us from having to copy back and forth
  if alloc_rhs 
    rhs_vec = zeros(Tsol, size(eqn.res_vec))
  else
    rhs_vec = eqn.res_vec
  end

  # should be all zeros if alloc_rhs is true
#   writedlm("rhs_initial.dat", rhs_vec)

  return newton_data, rhs_vec

end   # end of setupNewton

"""
  Cleans up after running Newton's method.

  **Inputs**

   * newton_data: the NewtonData object

"""
function free(newton_data::NewtonData)

  free(newton_data.ls)

  return nothing
end

"""
  Records the most recent nonlinear residual norm in the NewtonData object.
  Also updates the implicit Euler globalization

  **Inputs**

   * newton_data: the NewtonData
   * res_norm: the residual norm
"""
function recordResNorm(newton_data::NewtonData, res_norm::Number)

  newton_data.res_norm_i_1 = newton_data.res_norm_i
  newton_data.res_norm_i = res_norm
  
  # update implicit Euler globalization
  recordEulerResidual(res_norm)

  return nothing
end

"""
  Records norm of the most recent newton step (ie. the norm of delta q)
  in the NewtonData object

  **Inputs**

   * newton_data: the NewtonData object
   * step_norm: the norm of the step
"""
function recordStepNorm(newton_data::NewtonData, step_norm::Number)

  newton_data.step_norm_i_1 = newton_data.step_norm_i
  newton_data.step_norm_i = step_norm

  return nothing
end

#------------------------------------------------------------------------------
# getter for PC and LO

"""
  Returns the Newton precondtioner and linear operator specified by the options
  dictionary

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * rhs_func: rhs_func required by [`newtonInner`](@ref)
"""
function getNewtonPCandLO(mesh, sbp, eqn, opts)

  # get PC
  if opts["jac_type"] <= 2
    pc = PCNone(mesh, sbp, eqn, opts)
  else
    if opts["use_volume_preconditioner"]
      pc = NewtonVolumePC(mesh, sbp, eqn, opts)
    else
      pc = NewtonMatPC(mesh, sbp, eqn, opts)
    end
  end 

  jactype = opts["jac_type"]
  if jactype == 1
    lo = NewtonDenseLO(pc, mesh, sbp, eqn, opts)
  elseif jactype == 2
    lo = NewtonSparseDirectLO(pc, mesh, sbp, eqn, opts)
  elseif jactype == 3
    lo = NewtonPetscMatLO(pc, mesh, sbp, eqn, opts)
  elseif jactype == 4
    lo = NewtonPetscMatFreeLO(pc, mesh, sbp, eqn, opts)
  end

  return pc, lo
end

#------------------------------------------------------------------------------
# preconditioner


"""
  Adds fields required by all Newton PCs and LOs
"""
macro newtonfields()
  # compute something c = f(a, b)
  return quote
    res_norm_i::Float64  # current step residual norm
    res_norm_i_1::Float64  # previous step residual norm
    # Pseudo-transient continuation Euler
    tau_l::Float64  # current pseudo-timestep
    tau_vec::Array{Float64, 1}  # array of solution at previous pseudo-timestep
  end
end

type NewtonMatPC <: AbstractPetscMatPC
  pc_inner::PetscMatPC
  @newtonfields

end

function NewtonMatPC(mesh::AbstractMesh, sbp::AbstractSBP,
                    eqn::AbstractSolutionData, opts::Dict)


  pc_inner = PetscMatPC(mesh, sbp, eqn, opts)
    res_norm_i = 0.0
    res_norm_i_1 = 0.0
    if opts["newton_globalize_euler"]
      tau_l, tau_vec = initEuler(mesh, sbp, eqn, opts)
    else
      tau_l = opts["euler_tau"]
      tau_vec = []
    end


  return NewtonMatPC(pc_inner, res_norm_i, res_norm_i_1, tau_l, tau_vec)
end

function calcPC(pc::NewtonMatPC, mesh::AbstractMesh, sbp::AbstractSBP,
                eqn::AbstractSolutionData, opts::Dict, ctx_residual, t)

  calcPC(pc.pc_inner, mesh, sbp, eqn, opts, ctx_residual, t)
  physicsJac(mesh, sbp, eqn, opts. pc.A, ctx_residual, t)

  if opts["newton_globalize_euler"]
    # TODO: updating the Euler parameter here is potentially wrong if we
    #       are not updating the Jacobian at every newton step
    updateEuler(pc)
    applyEuler(mesh, sbp, eqn, opts, pc)
  end


  return nothing
end

#------------------------------------------------------------------------------
# linear operator

# because Julia lack multiple inheritance, we have to define 4 of these
# make sure they share the same fields whenever needed

type NewtonDenseLO <: AbstractDenseLO
  lo_inner::DenseLO
  @newtonfields 
end

function NewtonDenseLO(pc::PCNone, mesh::AbstractMesh,
                    sbp::AbstractSBP, eqn::AbstractSolutionData, opts::Dict)

  lo_inner = DenseLO(pc, mesh, sbp, eqn, opts)
  res_norm_i = 0.0
  res_norm_i_1 = 0.0
  if opts["newton_globalize_euler"]
    tau_l, tau_vec = initEuler(mesh, sbp, eqn, opts)
  else
    tau_l = opts["euler_tau"]
    tau_vec = []
  end


  return NewtonDenseLO(lo_inner, res_norm_i, res_norm_i_1, tau_l, tau_vec)
end

type NewtonSparseDirectLO <: AbstractSparseDirectLO
  lo_inner::SparseDirectLO
  @newtonfields
end

function NewtonSparseDirectLO(pc::PCNone, mesh::AbstractMesh,
                    sbp::AbstractSBP, eqn::AbstractSolutionData, opts::Dict)

  lo_inner = SparseDirectLO(pc, mesh, sbp, eqn, opts)

  res_norm_i = 0.0
  res_norm_i_1 = 0.0
  if opts["newton_globalize_euler"]
    tau_l, tau_vec = initEuler(mesh, sbp, eqn, opts)
  else
    tau_l = opts["euler_tau"]
    tau_vec = []
  end


  return NewtonSparseDirectLO(lo_inner, res_norm_i, res_norm_i_1, tau_l, tau_vec)
end

type NewtonPetscMatLO <: AbstractPetscMatLO
  lo_inner::PetscMatLO
  @newtonfields
end

function NewtonPetscMatLO(pc::AbstractPetscPC, mesh::AbstractMesh,
                    sbp::AbstractSBP, eqn::AbstractSolutionData, opts::Dict)

  lo_inner = PetscMatLO(pc, mesh, sbp, eqn, opts)

  res_norm_i = 0.0
  res_norm_i_1 = 0.0
  if opts["newton_globalize_euler"]
    tau_l, tau_vec = initEuler(mesh, sbp, eqn, opts)
  else
    tau_l = opts["euler_tau"]
    tau_vec = []
  end


  return NewtonPetscMatLO(lo_inner, res_norm_i, res_norm_i_1, tau_l, tau_vec)
end

type NewtonPetscMatFreeLO <: AbstractPetscMatFreeLO
  lo_inner::PetscMatFreeLO
  @newtonfields
end

"""
  Newton mat-free linear operator constructor

  **Inputs**

   * pc
   * mesh
   * sbp
   * eqn
   * opts
   * rhs_func: rhs_func from [`newtonInner`](@ref)
"""
function NewtonPetscMatFreeLO(pc::AbstractPetscPC, mesh::AbstractMesh,
                    sbp::AbstractSBP, eqn::AbstractSolutionData, opts::Dict)

  lo_inner = PetscMatFreeLO(pc, mesh, sbp, eqn, opts)
  res_norm_i = 0.0
  res_norm_i_1 = 0.0
  if opts["newton_globalize_euler"]
    tau_l, tau_vec = initEuler(mesh, sbp, eqn, opts)
  else
    tau_l = opts["euler_tau"]
    tau_vec = []
  end


  return NewtonPetscMatFreeLO(lo_inner, rhs_func, res_norm_i, res_norm_i_1, tau_l, tau_vec)
end

"""
  All Newton linear operators
"""
typealias NewtonLO Union{NewtonDenseLO, NewtonSparseDirectLO, NewtonPetscMatLO, NewtonPetscMatFreeLO}

"""
  Newton matrix-explicit linear operators
"""
typealias NewtonMatLO Union{NewtonDenseLO, NewtonSparseDirectLO, NewtonPetscMatLO}

"""
  Any PC or LO that has a matrix in the field `A`
"""
typealias NewtonHasMat Union{NewtonMatPC, NewtonDenseLO, NewtonSparseDirectLO, NewtonPetscMatLO}

"""
  Any Newton PC or LO.
"""
typealias NewtonLinearObject Union{NewtonDenseLO, NewtonSparseDirectLO, NewtonPetscMatLO, NewtonPetscMatFreeLO, NewtonMatPC}


function calcLinearOperator(lo::NewtonMatLO, mesh::AbstractMesh,
                            sbp::AbstractSBP, eqn::AbstractSolutionData,
                            opts::Dict, ctx_residual, t)

   
  calcLinearOperator(lo.lo_inner, mesh, sbp, eqn, opts, ctx_residual, t)

  lo2 = getBaseLO(lo)
  physicsJac(mesh, sbp, eqn, opts, lo2.A, ctx_residual, t)

  if opts["newton_globalize_euler"]
    # TODO: updating the Euler parameter here is potentially wrong if we
    #       are not updating the Jacobian at every newton step
    updateEuler(lo)
    applyEuler(mesh, sbp, eqn, opts, lo)
  end

  return nothing
end

function calcLinearOperator(lo::NewtonPetscMatFreeLO, mesh::AbstractMesh,
                            sbp::AbstractSBP, eqn::AbstractSolutionData,
                            opts::Dict, ctx_residual, t)

  if opts["newton_globalize_euler"]
    updateEuler(lo)
  end

  setLOCtx(lo, mesh, sbp, eqn, opts, ctx_residual, t)

  return nothing
end


function applyLinearOperator{Tsol}(lo::NewtonPetscMatFreeLO, mesh::AbstractMesh,
                             sbp::AbstractSBP, eqn::AbstractSolutionData{Tsol},
                             opts::Dict, ctx_residual, t, x::AbstractVector, 
                             b::AbstractVector)

  @assert !(Tsol <: AbstractFloat)  # complex step only!

  epsilon =  opts["epsilon"]::Float64
  pert = Tsol(0, epsilon)

  # apply perturbation
  for i=1:mesh.numDof
    eqn.q_vec[i] += pert*vec[i]
  end

  physicsRhs(mesh, sbp, eqn, opts, eqn.res_vec, ctx_residual, t)
  
  # calculate derivatives, store into b
  calcJacCol(b, eqn.res_vec, epsilon)

  if opts["newton_globalize_euler"]
    applyEuler(mesh, sbp, eqn, opts, x, lo, b)
  end

  # undo perturbation
  for i=1:mesh.numDof
    eqn.q_vec[i] -= pert*vec[i]
  end

  return nothing
end

function applyLinearOperatorTranspose{Tsol}(lo::NewtonPetscMatFreeLO, 
                             mesh::AbstractMesh,
                             sbp::AbstractSBP, eqn::AbstractSolutionData{Tsol},
                             opts::Dict, ctx_residual, t, x::AbstractVector, 
                             b::AbstractVector)

  error("applyLinearOperatorTranspose() not supported by NewtonPetscMatFreeLO")

end
