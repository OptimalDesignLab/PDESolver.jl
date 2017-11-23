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

  # inexact Newton-Krylov parameters
  reltol::Float64
  abstol::Float64
  dtol::Float64
  itermax::Int
  krylov_gamma::Float64  # update parameter for krylov tolerance

  res_norm_i::Float64  # current step residual norm
  res_norm_i_1::Float64  # previous step residual norm

  ls::Tsolver  # non-concrete type

  #TODO; make this an outer constructor
end

function NewtonData(mesh, sbp, eqn, opts, ls::LinearSolver)

  println("entered NewtonData constructor")
  println("typeof(eqn) = ", typeof(eqn))

  myrank = mesh.myrank
  commsize = mesh.commsize

  reltol = opts["krylov_reltol"]
  abstol = opts["krylov_abstol"]
  dtol = opts["krylov_dtol"]
  itermax = opts["krylov_itermax"]
  krylov_gamma = opts["krylov_gamma"]

  res_norm_i = 0.0
  res_norm_i_1 = 0.0

  return NewtonData{Tsol, Tres, typeof(ls)}(myrank, commsize, reltol, abstol,
                    dtol, itermax, krylov_gamma, res_norm_i, res_norm_i_1, ls)
end


include("residual_evaluation.jl")  # some functions for residual evaluation
include("jacobian.jl")
include("petsc_funcs.jl")  # Petsc related functions

@doc """
### NonlinearSolvers.setupNewton
  
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

  newton_data = NewtonData{Tsol, Tres}(mesh, sbp, eqn, opts, ls)

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

#------------------------------------------------------------------------------
# getter for PC and LO

"""
  Returns the Newton precondtioner and linear operator specified by the options
  dictionary
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

  calcPC(pc.pc_inner)
  #TODO; zero the jacobian here?
  physicsJac(mesh, sbp, eqn, opts. pc.A, ctx_residual, t)

  #TODO: globalization here

  return nothing
end


type NewtonMatFreePC <: AbstractPetscMatFreePC
  pc_inner::PetscMatFreePC
  @newtonfields
end

function NewtonMatFreePC(mesh::AbstractMesh, sbp::AbstractSBP,
                    eqn::AbstractSolutionData, opts::Dict)


  pc_inner = PetscMatFreePC(mesh, sbp, eqn, opts)
  res_norm_i = 0.0
  res_norm_i_1 = 0.0
  if opts["newton_globalize_euler"]
    tau_l, tau_vec = initEuler(mesh, sbp, eqn, opts)
  else
    tau_l = opts["euler_tau"]
    tau_vec = []
  end


  return NewtonMatFreePC(pc_inner, res_norm_i, res_norm_i_1, tau_l, tau_vec)
end

function calcPC(pc::NewtonMatFreePC, mesh::AbstractMesh, sbp::AbstractSBP,
                eqn::AbstractSolutionData, opts::Dict, ctx_residual, t)

  calcPC(pc.pc_inner)
  #TODO; zero the jacobian here?

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

  lo_inner = DenseLO(mesh, sbp, eqn, opts)
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

  lo_inner = SparseDirectLO(mesh, sbp, eqn, opts)

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

  lo_inner = PetscMatLO(mesh, sbp, eqn, opts)

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

function NewtonPetscMatFreeLO(pc::AbstractPetscPC, mesh::AbstractMesh,
                    sbp::AbstractSBP, eqn::AbstractSolutionData, opts::Dict)

  lo_inner = PetscMatFreeLO(mesh, sbp, eqn, opts)
  res_norm_i = 0.0
  res_norm_i_1 = 0.0
  if opts["newton_globalize_euler"]
    tau_l, tau_vec = initEuler(mesh, sbp, eqn, opts)
  else
    tau_l = opts["euler_tau"]
    tau_vec = []
  end


  return NewtonPetscMatFreeLO(lo_inner, res_norm_i, res_norm_i_1, tau_l, tau_vec)
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
  Any PC or LO that does not have an explicit matrix
"""
typealias NewtonHasNoMat Union{NewtonMatFreePC, NewtonPetscMatFreeLO}

"""
  Any Newton PC or LO.
"""
typealias NewtonLinearObject Union{NewtonDenseLO, NewtonSparseDirectLO, NewtonPetscMatLO, NewtonPetscMatFreeLO, NewtonMatPC, NewtonMatFreePC}


function calcLinearOperator(lo::NewtonMatLO, mesh::AbstractMesh,
                            sbp::AbstractSBP, eqn::AbstractSolutionData,
                            opts::Dict, ctx_residual, t)

   
   calcLinearOperator(pc.pc_inner)

  #TODO; zero the jacobian here?
  physicsJac(mesh, sbp, eqn, opts. lo.A, ctx_residual, t)

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

  rhs_func(mesh, sbp, eqn, opts, eqn.res_vec, ctx_residual, t)
  
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
