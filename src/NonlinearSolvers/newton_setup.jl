# auxiliary types and functions for Newton's method
# linear operator, newton data etc.

@doc """
  This type holds all the data the might be needed for Newton's method,
  including globalization.  This simplifies the data handling and 
  the C callback used by Petsc
"""->
type NewtonData{Tsol, Tres}

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
  # Pseudo-transient continuation Euler
  tau_l::Float64  # current pseudo-timestep
  tau_vec::Array{Float64, 1}  # array of solution at previous pseudo-timestep

  lo::StandardLinearSolver  # non-concrete type

  # volume preconditioning
  vol_prec::VolumePreconditioner

  # tuple of values needed to cleanup Petsc data structures
  ctx_newton
  ksp::KSP
  pc::PC

  # ctx needed by MatShell
  # TODO: get list of contents from jac-vec prod wrapper
  ctx_petsc
  ctx_petsc_pc  # for the matrix-free preconditioner

  #TODO; make this an outer constructor
  function NewtonData(mesh, sbp, eqn, opts)

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
    if opts["newton_globalize_euler"]
      tau_l, tau_vec = initEuler(mesh, sbp, eqn, opts)
    else
      tau_l = opts["euler_tau"]
      tau_vec = []
    end

    if opts["use_volume_preconditioner"]
      vol_prec = VolumePreconditioner(mesh, sbp, eqn, opts)
    else
      vol_prec = VolumePreconditioner()
    end

    # NOTE: we are leaving ctx_newton uninitialized because 
    #   createPetscData needs it, but ctx_newton depends on its returned values

    # initialize these to something, will be replaced in setupNewton
    ctx_newton = ()
    # these get replaced later, if needed
    ksp = PETSc.KSP_NULL
    pc = PETSc.PC_NULL
    ctx_petsc = ()
    ctx_petsc_pc = ()

    return new(myrank, commsize, reltol, abstol, dtol, 
                      itermax, krylov_gamma, 
                      res_norm_i, res_norm_i_1, tau_l, tau_vec, 
                      vol_prec, ctx_newton, ksp, pc, ctx_petsc, ctx_petsc_pc)
  end

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
                     eqn::AbstractSolutionData{Tsol, Tres}, opts;
                     alloc_rhs=true)

  jac_type = opts["jac_type"]
  Tjac = typeof(real(eqn.res_vec[1]))  # type of jacobian, residual
  m = mesh.numDof

  # ctx_newton is not defined yet
  newton_data = NewtonData{Tsol, Tres}(mesh, sbp, eqn, opts)

  # Allocation of Jacobian, depending on type of matrix
  eqn.params.time.t_alloc += @elapsed if jac_type == 1  # dense
    jac = zeros(Tjac, m, m)  # storage of the jacobian matrix
    ctx_newton = ()
  elseif jac_type == 2  # sparse
    if typeof(mesh) <: AbstractCGMesh
      println("creating CG SparseMatrix")
      jac = SparseMatrixCSC(mesh.sparsity_bnds, Tjac)
    else
      println("Creating DG sparse matrix")
      jac = SparseMatrixCSC(mesh, Tjac)
    end
    ctx_newton = ()
  elseif jac_type == 3 || jac_type == 4 # petsc
    jac, jacp, x, b, ksp = createPetscData(mesh, pmesh, sbp, eqn, opts, newton_data)
    ctx_newton = (jac, jacp, x, b, ksp)
  end

  # now put ctx_newton into newton_data
  newton_data.ctx_newton = ctx_newton

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

  return newton_data, jac, rhs_vec

end   # end of setupNewton


#------------------------------------------------------------------------------
# preconditioner

type NewtonMatPC <: AbstractPetscMatPC
  pc_inner::PetscMatPC
  #TODO: globalization fields
end

function NewtonMatPC(mesh::AbstractMesh, sbp::AbstractSBP,
                    eqn::AbstractSolutionData, opts::Dict)


  pc_inner = PetscMatPC(mesh, sbp, eqn, opts)

  return NewtonMatPC(pc_inner)
end

function calcPC(pc::NewtonMatPC, mesh::AbstractMesh, sbp::AbstractSBP,
                eqn::AbstractSolutionData, opts::Dict, ctx_residual, t)

  calcPC(pc.pc_inner)
  #TODO; zero the jacobian here?
  physicsJac(mesh, sbp, eqn, opts. pc.Ap, ctx_residual, t)

  #TODO: globalization here

  return nothing
end


type NewtonMatFreePC <: AbstractPetscMatFreePC
  pc_inner::PetscMatFreePC
  #TODO: globalization fields
end

function NewtonMatFreePC(mesh::AbstractMesh, sbp::AbstractSBP,
                    eqn::AbstractSolutionData, opts::Dict)


  pc_inner = PetscMatFreePC(mesh, sbp, eqn, opts)

  return NewtonMatFreePC(pc_inner)
end

function calcPC(pc::NewtonMatFreePC, mesh::AbstractMesh, sbp::AbstractSBP,
                eqn::AbstractSolutionData, opts::Dict, ctx_residual, t)

  calcPC(pc.pc_inner)
  #TODO; zero the jacobian here?

  #TODO: globalization here

  return nothing
end


#------------------------------------------------------------------------------
# linear operator

# because Julia lack multiple inheritance, we have to define 4 of these
# make sure they share the same fields whenever needed

type NewtonDenseLO <: AbstractDenseLO
  lo_inner::DenseLO
end

function NewtonDenseLO(pc::PCNone, mesh::AbstractMesh,
                    sbp::AbstractSBP, eqn::AbstractSolutionData, opts::Dict)

  lo_inner = DenseLO(mesh, sbp, eqn, opts)

  return NewtonDenseLO(lo_inner)
end

type NewtonSparseDirectLO <: AbstractSparseDirectLO
  lo_inner::SparseDirectLO
end

function NewtonSparseDirectLO(pc::PCNone, mesh::AbstractMesh,
                    sbp::AbstractSBP, eqn::AbstractSolutionData, opts::Dict)

  lo_inner = SparseDirectLO(mesh, sbp, eqn, opts)

  return NewtonSparseDirectLO(lo_inner)
end

type NewtonPetscMatLO <: AbstractPetscMatLO
  lo_inner::PetscMatLO
end

function NewtonPetscMatLO(pc::PCNone, mesh::AbstractMesh,
                    sbp::AbstractSBP, eqn::AbstractSolutionData, opts::Dict)

  lo_inner = PetscMatLO(mesh, sbp, eqn, opts)

  return NewtonPetscMatLO(lo_inner)
end

type NewtonPetscMatFreeLO <: AbstractPetscMatFreeLO
  lo_inner::PetscMatFreeLO
end

function NewtonPetscMatFreeLO(pc::PCNone, mesh::AbstractMesh,
                    sbp::AbstractSBP, eqn::AbstractSolutionData, opts::Dict)

  lo_inner = PetscMatFreeLO(mesh, sbp, eqn, opts)

  return NewtonPetscMatFreeLO(lo_inner)
end

"""
  All Newton linear operators
"""
typealias NewtonLO Union{NewtonDenseLO, NewtonSparseDirectLO, NewtonPetscMatLO, newtonPetscMatFreeLO}

"""
  Newton matrix-explicit linear operators
"""
typealias NewtonMatLO Union{NewtonDenseLO, NewtonSparseDirectLO, NewtonPetscMatLO}


function calcLinearOperator(lo::NewtonMatLO, mesh::AbstractMesh,
                            sbp::AbstractSBP, eqn::AbstractSolutionData,
                            opts::Dict, ctx_residual, t)

   calcLinearOperator(pc.pc_inner)

  #TODO; zero the jacobian here?
  physicsJac(mesh, sbp, eqn, opts. lo.A, ctx_residual, t)

  #TODO: globalization here

  return nothing
end

function calcLinearOperator(lo::NewtonMatLO, mesh::AbstractMesh,
                            sbp::AbstractSBP, eqn::AbstractSolutionData,
                            opts::Dict, ctx_residual, t)

  # nothing to do here

  return nothing
end


function applyLinearOperator{Tsol}(lo::NewtonMatFreeLO, mesh::AbstractMesh,
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

  if globalize_euler
    applyEuler(mesh, sbp, eqn, opts, vec, newton_data, b)
  end

  # undo perturbation
  for i=1:mesh.numDof
    eqn.q_vec[i] -= pert*vec[i]
  end

  return nothing
end


