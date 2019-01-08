# This file provide abstraction around performing linear solves (Ax = b)
# and concrete implementations using Dense matrices, SparseMatrixCSC, and
# Petsc (both matrix-explicit and matrix-free)

module LinearSolvers

export LinearSolver, StandardLinearSolver,  # Linear Solver types
       calcPC, calcLinearOperator, calcPCandLO, # Linear Solver interface
       applyPC, applyPCTranspose, 
       linearSolve, linearSolveTranspose, 
       isLOMatFree, isPCMatFree, setTolerances, free,  # utility functions
       applyLinearOperator, applyLinearOperatorTranspose,  # LO functions
       getBaseLO, getBasePC, getBaseObject,
       PCNone, PetscMatPC, PetscMatFreePC,  # PC types
       DenseLO, SparseDirectLO, PetscMatLO, PetscMatFreeLO,  # LO types
       AbstractPC, AbstractPetscMatPC, AbstractPetscMatFreePC,# Abstract PC types
       AbstractPetscPC, getInnerPC,
       AbstractLO, AbstractDenseLO, AbstractSparseDirectLO, # Abstrac LO types
       AbstractPetscMatLO, AbstractPetscMatFreeLO, getInnerLO,
       setPCCtx, setLOCtx,  # matrix-free specific functions
       needParallelData




using ODLCommonTools
using Utils
using SummationByParts
using Base.LinAlg.BLAS
using MPI
using PETSc2

import Utils.free
# import SuiteSparse stuff
import Base.SparseArrays.UMFPACK: UmfpackLU, umfpack_free_numeric,
                                  umfpack_free_symbolic, umfpack_symbolic!,
                                  umfpack_numeric!


"""
  Abstract supertype of all linear solvers.  The [`StandardLinearSolver`](@ref)
  implementation should be general enough for everything we do.

  The purpose of this type is to provide a unified interface for managing
  a linear solve, including preconditioning if needed.  Many of the
  operations on this type are delegated to the pc and lo objects.

  The [`AbstractPC`](@ref) and [`AbstractLO`](@ref) types
  defined in this module are suitible for solving Ax = b when A is the
  Jacobian of the physics.  Other methods (for example, unsteady time marching
  or homotopy methods) should build their own AbstractPC and
  AbstractLO objects and use them with
  [`StandardLinearSolver`](@ref).

  **Required Fields**

   * pc: an [`AbstractPC`](@ref) object
   * lo: an [`AbstractLO`](@ref) object
   * shared_mat:  Bool, true if pc and lo share the same matrix object false
                  otherwise (either different matrix objects or matrix-free)

  **Static Parameters**

   * T1: type of pc
   * T2: type of lo
"""
abstract type LinearSolver{T1, T2} end

"""
  The most commonly used implementation of [`LinearSolver`](@ref).
"""
mutable struct StandardLinearSolver{T1, T2} <: LinearSolver{T1, T2}
  pc::T1
  lo::T2
  shared_mat::Bool
  comm::MPI.Comm
  myrank::Int
  commsize::Int
  ksp::KSP  # used only for Petsc matrices
  is_finalized::Bool

  # tolerances for iterative solve
  reltol::PetscReal
  abstol::PetscReal
  dtol::PetscReal
  itermax::PetscInt
end


"""
  Constructor for StandardLinearSolver.

  **Inputs**

   * pc: an [`AbstractPC`](@ref), fully initialized
   * lo: an [`AbstractLO`](@ref), fully initialized
   * comm: the MPI communicator the pc and lo are defined on
   * opts: options dictionary

  This function throws exceptions if incompatible pc and lo types are used.

  **Options Keys**

  This function uses: "krylov_reltol", "krylov_abstol", "krylov_dtol", and
  "krylov_itermax"

"""
function StandardLinearSolver(pc::T1, lo::T2, comm::MPI.Comm, opts) where {T1, T2}

  if typeof(lo) <: DirectLO
    @assert typeof(pc) <: PCNone
  end

  if typeof(lo) <: PetscLO
    @assert typeof(pc) <: Union{AbstractPetscMatPC, AbstractPetscMatFreePC}
  end

  pc2 = getBasePC(pc)
  lo2 = getBaseLO(lo)
  if typeof(pc2) <: PetscMatPC && typeof(lo2) <: PetscMatLO
    if pc2.A.pobj == lo2.A.pobj
      shared_mat = true
    else
      shared_mat = false
    end
  else
    shared_mat = false
  end

  myrank = MPI.Comm_rank(comm)
  commsize = MPI.Comm_size(comm)

  if typeof(lo) <: PetscLO
    ksp = createKSP(pc2, lo2, comm)
  else
    ksp = KSP_NULL
  end

  is_finalized = false

  reltol = opts["krylov_reltol"]
  abstol = opts["krylov_abstol"]
  dtol = opts["krylov_dtol"]
  itermax = opts["krylov_itermax"]

  ls = StandardLinearSolver{T1, T2}(pc, lo, shared_mat, comm, myrank,
                                    commsize, ksp, is_finalized,
                                    reltol, abstol, dtol, itermax)

  finalizer(ls, free)

  return ls
end


"""
  Abstract supertype of all preconditioner types.  Preconditioners can be used
  with iterative methods only.  When using a direct method, [`PCNone`](@ref)
  should be used.

  The purpose of this type is to provide a consistent interface for different
  types of preconditioners.  In particular, it provides control over when the
  preconditioner is recomputed.

  Users should implement new preconditioners using composition with
  one of the preconditioners defined here, namely [`PetscMatPC`](@ref) or
  [`PetscMatFreePC`](@ref), ie. they should define a new subtype of
  [`AbstractPC`](@ref) that has either `PetscMatPC` or `PetscMatFrePC` as
  a field.  This allows calling the existing functions for these types
  (which compute a preconditioner for the Jacobian of the physics) and
  modifying the preconditioner as needed.
  User defined preconditioners should subtype either [`AbstractPetscMatPC`](@ref)
  of [`AbstractPetscMatFreePC`](@ref).

  **Fields**

   * pc_inner: another [`AbstractPC`](@ref)

  Note that arbitrarily deep nesting of preconditioners is allowed.
  The `pc_inner` field can be one of [`PCNone`](@ref), [`PetscMatPC`](@ref),
  or [`PetscMatFreePC`](@ref) for a non-nested preconditioner, or some
  other [`AbstractPC`](@ref) for a nested preconditioner.
"""
abstract type AbstractPC end

"""
  Abstract supertype of all Petsc matrix-explicit preconditioners
"""
abstract type AbstractPetscMatPC <: AbstractPC end

"""
  Abstract supertype of all Petsc matrix-free preconditioners.
"""
abstract type AbstractPetscMatFreePC <: AbstractPC end

"""
  Alias for any kind of Petsc PC (matrix-explicit or matrix-free)
"""
const AbstractPetscPC = Union{AbstractPetscMatPC, AbstractPetscMatFreePC}

# PC interface
"""
  This function calculates the preconditioner.  Every implementation of
  [`AbstractPC`](@ref) should extend this function with a new method

  When creating new preconditioners, this function should generally be called
  first on pc_inner, and modifications to the Jacobian should be made 
  subsequently.

  **Inputs**

   * pc: the AbstractPC implementation
   * mesh
   * sbp
   * eqn
   * opts
   * ctx_residual: the ctx required by [`physicsRhs`](@ref)
   * t: current time

  Implementation Notes:
    For matrix-explicit preconditioners, this might not actually calculate the
    preconditioner.  Rather, it calculates the matrix the preconditioner is
    based on, and the solve function calcultes the preconditioner from it.
    Nevertheless, it supports the proper semantics for when the PC is updated
    even when the PC matrix and LinearOperator matrix are the same (as long as
    the solve function is called regularly).
"""
function calcPC(pc::AbstractPC, mesh::AbstractMesh, sbp::AbstractOperator,
                eqn::AbstractSolutionData, opts::Dict, ctx_residual, t)

  error("reached AbstractPC calcPC(), did you forget to define calcPC() for your AbstractPC implementation?")
end

"""
  Applies the preconditioner, ie. x = inv(Ap)*b, where Ap is the approximation
  to the matrix A.  Note that Ap itself may not be available for some
  preconditioners, hence there is only an API for applying inv(Ap),
  not Ap itself.

    Matrix-free preconditioners need to extend this function with a new method,
    matrix-explicit preconditioners do not.

  **Inputs**

   * pc: the [`AbstractPC`](@ref) implementation.
   * mesh
   * sbp
   * eqn: this argument should generally not be used because all the solution
          related data should be stored in the pc object by [`calcPC`](@ref)
   * opts
   * t: current time
   * b: a AbstractVector representing the local part of the solution (ie
        eqn.q_vec)

  **Inputs/Outputs**

   * x: AbstractVector overwritten with result (same size as b).
        Some preconditioners may use the value of `x` on entry
        as the initial guess.
"""
function applyPC(pc::AbstractPC, mesh::AbstractMesh, sbp::AbstractOperator,
                 eqn::AbstractSolutionData, opts::Dict, t, b::AbstractVector, 
                 x::AbstractVector)


  error("reached AbstractPC applyPC(), did you forget to define applyPC() for your AbstractPC implementation?")

end

"""
  Applies the transpose of the preconditioner, ie. x = inv(Ap).'*b.
  Similar to [`applyPC`](@ref), see that function for details.

  Note that not every preconditioning method supports this.
"""
function applyPCTranspose(pc::AbstractPC, mesh::AbstractMesh, sbp::AbstractOperator,
                 eqn::AbstractSolutionData, opts::Dict, t, b::AbstractVector, 
                 x::AbstractVector)

  
  error("reached AbstractPC applyPCTranspose(), did you forget to define applyPCTranspose() for your AbstractPC implementation?")

end

"""
  This function returns the underlying preconditioner object, ie.
  [`PCNone`](@ref), [`PetscMatPC`](@ref), or [`PetscMatFreePC`](@ref).

  Note that arbitrarily deep nesting of preconditioners is allowed.
  Users do not have to implement as long as the nested preconditioner is
  stored in a field called `pc_inner`.

  For matrix-explicit preconditioners, this function is useful for getting
  the [`PetscMatPC`](@ref) object, which contains the preconditioning
  Jacobian matrix.

  **Inputs**

   * pc: the users [`AbstractPC`](@ref)
"""
function getBasePC(pc::AbstractPC)

  # this will recurse all the way down to the underyling pc, which returns
  # itself
  return getBasePC(pc.pc_inner)
end

"""
  This function allows extracting a specific type of pc_inner.

  Note: this is not type-stable in julia 0.4, but can be rewritten to be so
        in later version of Julia

  **Inputs**

   * pc: a PC
   * T2: a (possibly abstract) type that is contained as an inner PC somewhere
         inside pc.

  **Outputs**

   * pc2: a PC that is a subtype of PC2
"""

function getInnerPC(pc::T1, ::Type{T2}) where {T1 <: AbstractPC, T2}
# can't have T2 <: AbstractPC because that would preclude Unions like
# NewtonLinearObject
#
  if T1 <: T2
    return pc
  else
    return getInnerPC(pc.pc_inner, T2)
  end
  error("unreachable reached")
end


#=
# This is broken on Julia 0.6
# T1 <: T2 implies T1 == T2 when T1 and T2 are concrete
function getInnerPC(pc::T1, ::Type{T2}) where {T2 <: AbstractPC, T1 <: T2}
  return pc
end

function getInnerPC(pc::T1, ::Type{T2}) where {T1 <:AbstractPC, T2 <: AbstractPC}
  return getInnerPC(pc.pc_inner, T2)
end
=#



"""
  Returns `true` if `calcPC` needs parallel communication started before
  it is called, false otherwise.  Defaults to true.  Users should extend this
  function with a new method if that particular PC does not require parallel
  communication.

  **Inputs**

   * pc:  an `AbstractPC

  **Outputs**

   * Bool
"""
function needParallelData(pc::AbstractPC)
  return true
end

"""
  This function frees any memory belonging to external libraries.  Users must
  call this function when they are finished with an AbstractPC
  object.

  Users do not have to define this function for their
  [`AbstractPC`](@ref) types.

  **Inputs**

   * pc: the AbstractPC object
"""
function free(pc::AbstractPC)

  free(getBasePC(pc))
end


# LinearOperator interface
"""
  Abstract supertype of all linear operators used for A when solving Ax = b.

  The purpose of this type is to provide a consistent interface for different
  types of linear operators.  This type really combines two notions: what the
  type of the linear operator is, and how it should be solved.

  Any implementation of this type should subtype the appropriate catagory
  of: [`AbstractDenseLO`](@ref), [`AbstractSparseDirectLO`](@ref),
  [`AbstractPetscMatLO`](@ref), [`AbstractPetscMatFreeLO`](@ref)

  Note that matrix-explicit implementations can often write a single function
  for all these cases if using an matrix interface functions that are defined
  for all the matrix types.  See [`MatExplicitLO`](@ref)

  **Required Fields**

   * lo_inner: another `AbstractPC`.  Can be one of [`DenseLO`](@ref),
               [`SparseDirectLO`](@ref), [`PetscMatLO`](@ref), or
               [`PetscMatFreeLO`](@ref), or any other user defined linear
               operator.
"""
abstract type AbstractLO end

#TODO: doc these

"""
  Linear operator type for Dense matrices.  This is generally used only for
  debugging.
"""
abstract type AbstractDenseLO <: AbstractLO end

"""
  Linear operator type for `SparseMatrixCSC` matrices, which use a direct
  solver.
"""
abstract type AbstractSparseDirectLO <: AbstractLO end

"""
  Linear operator type for Petsc matrix-explicit.
"""
abstract type AbstractPetscMatLO <: AbstractLO end

"""
  Linear operator type for Petsc matrix-free.
"""
abstract type AbstractPetscMatFreeLO <: AbstractLO end

"""
  Useful union for all the matrix-explicit linear operator types.
  Because matrices have a small set of common interface functions, it is
  often possible to write a single function that works on all the different
  types of matrices.
"""
const MatExplicitLO = Union{AbstractDenseLO, AbstractSparseDirectLO, AbstractPetscMatLO}

"""
  Union of Petsc linear operator types
"""
const PetscLO = Union{AbstractPetscMatLO, AbstractPetscMatFreeLO}

"""
  Union of linear operators that do direct solves
"""
const DirectLO = Union{AbstractDenseLO, AbstractSparseDirectLO}


"""
  This function calculates the linear operator.  Every implementation of
  [`AbstractLO`](@ref) should extend this function with a new
  method.  For matrix-free operators, this function must exist but need
  not perform any actions.

  For matrix-explicit implementations, this function should be called on
  `lo_inner` first and modifications to the Jacobian made subsequently.

  **Inputs**

   * lo: the AbstractLO implementation (fields may be updated)
   * mesh
   * sbp
   * eqn
   * opts
   * ctx_residual: the ctx required by [`physicsRhs`](@ref)
   * t: current time


  Implementation Notes:
    For matrix-free operation, this function sets the Petsc ctx for the
    PetscMat, which contains a reference to the mesh, sbp, eqn, opts arguments.
    This could lead to unexpected behavior if those arguments are modified and
    this function is not called again before the next solve.

"""
function calcLinearOperator(lo::AbstractLO, mesh::AbstractMesh,
                            sbp::AbstractOperator, eqn::AbstractSolutionData,
                            opts::Dict, ctx_residual, t)

  error("reached AbstractLO calcLinearOperator(), did you forget to define calcLinearOperator() for your AbstractLO implementation?")
end

"""
  Applies the linear operator, ie. , Ax = b
  
  Matrix-explicit implementations [`AbstractLO`](@ref) do not have
  to implement this function, though matrix-free implementations must extend
  it with a new method.

  **Inputs**

   * lo: the [`AbstractLO`](@ref) implementation.
   * mesh
   * sbp
   * eqn
   * opts
   * ctx_residual: the ctx for [`physicsRhs`](@ref) or the another right hand
                   side function built on top of it
   * t: current time
   * x: an AbstractVector (although never a PetscVec)

  **Inputs/Outputs**

   * b: vector updated with results (do not overwrite)
"""
function applyLinearOperator(lo::AbstractLO, mesh::AbstractMesh,
                             sbp::AbstractOperator, eqn::AbstractSolutionData,
                             opts::Dict, ctx_residual, t, x::AbstractVector, 
                             b::AbstractVector)


  error("reached AbstractLO applyLinearOperator(), did you forget to define applyLinearOperator() for your AbstractLO implementation?")

end

"""
  Applies the transpose of the linear operator, ie. A.'*x = b
  Similar to [`applyLinearOperator`](@ref), see that function for details.

  Note that not every method supports this.  In particular, Petsc
  matrix-free LinearOperators don't currently expose this 
  (although they could with enough reverse-mode)

"""
function applyLinearOperatorTranspose(lo::AbstractLO, 
                             mesh::AbstractMesh, sbp::AbstractOperator,
                             eqn::AbstractSolutionData, opts::Dict, 
                             ctx_residual, t, x::AbstractVector, 
                             b::AbstractVector)

  error("reached AbstractLO applyLinearOperatorTranspose(), did you forget to define applyLinearOperatorTranspose() for your AbstractLO implementation?")

end


"""
  Similar to [`getBasePC`](@ref) except it gets the underlying linear operator,
  ie. one of [`DenseLO`](@ref), [`SparseDirectLO`](@ref), [`PetscMatLO`](@ref)
  or [`PetscMatFreeLO`](@ref).

  For matrix-explicit methods, this is a good way of getting the underlying
  linear operator object, which contains the matrix in the `A` field (for
  all matrix-explicit linear operators).

  **Inputs**

   * lo: an AbstractLO
"""
function getBaseLO(lo::AbstractLO)

  # this will recurse down to the underlying linear operator
  return getBaseLO(lo.lo_inner)
end

"""
  This function calls either [`getBasePC`](@ref) or [`getBaseLO`](@ref)
  depending on the type of its argument.

  **Inputs**

   * lo: either an `AbstractLO` or `AbstractPC` object

  **Outputs**

   * the base PC or LO object
"""
function getBaseObject(lo::AbstractLO)

  return getBaseLO(lo)
end

function getBaseObject(pc::AbstractPC)
  
  return getBasePC(pc)
end


"""
  Like [`getInnerLO`](@ref), but for linear operators
"""
function getInnerLO(lo::T1, ::Type{T2}) where {T1 <: AbstractLO, T2}
# can't have T2 <: AbstractLO because that would precude Unions like 
# NewtonLinearOperators
  if T1 <: T2
    return lo
  else
    return getInnerLO(lo.lo_inner, T2)
  end

  error("unreachable reached")
end



"""
  Returns `true` if [`calcLinearOperator`](@ref) needs parallel communication
  started before
  it is called, false otherwise.  Defaults to true for non-matrix-free
  linear operators.  
  Users should extend this function with a new method if the default value is
  not correct for a particular linear operator.

  **Inputs**

   * lo:  an `AbstractLO

  **Outputs**

   * Bool
"""
function needParallelData(pc::AbstractLO)
  return true
end

function needParallelData(lo::AbstractPetscMatFreeLO)
  return false
end





"""
  This function frees any memory belonging to external libraries.  Users must
  call this function when they are finished with an AbstractLO
  object.

  Users do not have to define this function for their
  [`AbstractLO`](@ref) types.

  **Inputs**

   * lo: the AbstractLO object
"""
function free(lo::AbstractLO)

  free(getBaseLO(lo))
end

# include implementations
include("pc_none.jl")
include("pc_petscmat.jl")
include("pc_petscmatfree.jl")
include("lo_dense.jl")
include("lo_sparsedirect.jl")
include("lo_petscmat.jl")
include("lo_petscmatfree.jl")
include("ls_standard.jl")
include("utils.jl")

end  # end module
