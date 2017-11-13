# This file provide abstraction around performing linear solves (Ax = b)
# and concrete implementations using Dense matrices, SparseMatrixCSC, and
# Petsc (both matrix-explicit and matrix-free)

module LinearSolvers

using ODLCommonTools
using Utils
using SummationByParts

"""
  Abstract supertype of all linear solvers.  The [`StandardLinearSolver`](@ref)
  implementation should be general enough for everything we do.

  The purpose of this type is to provide a unified interface for managing
  a linear solve, including preconditioning if needed.  Many of the
  operations on this type are delegated to the pc and lo objects.

  The [`AbstractPC`](@ref) and [`AbstractLinearOperator`](@ref) types
  defined in this module are suitible for solving Ax = b when A is the
  Jacobian of the physics.  Other methods (for example, unsteady time marching
  or homotopy methods) should build their own AbstractPC and
  AbstractLinearOperator objects and use them with
  [`StandardLinearSolver`](@ref).

  **Required Fields**

   * pc: an [`AbstractPC`](@ref) object
   * lo: an [`AbstractLinearOperator`](@ref) object
   * shared_mat:  Bool, true if pc and lo share the same matrix object false
                  otherwise (either different matrix objects or matrix-free)

  **Static Parameters**

   * T1: type of pc
   * T2: type of lo
"""
abstract LinearSolver{T1, T2}

"""
  The most commonly used implementation of [`LinearSolver`](@ref).
"""
type StandardLinearSolver{T1, T2} <: LinearSolver{T1, T2}
  pc::T1
  lo::T2
  shared_mat::Bool
end

"""
  Constructor for StandardLinearSolver
"""
function StandardLinearSolver{T1, T2}(pc::T1, lo::T2)

end

"""
  Abstract supertype of all preconditioner types.  Preconditioners can be used
  with iterative methods only.  When using a direct method, [`PCNone`](@ref)
  should be used.

  The purpose of this type is to provide a consistent interface for different
  types of preconditioners.  In particular, it provides control over when the
  preconditioner is recomputed
"""
abstract AbstractPC

# PC interface
"""
  This function calculates the preconditioner.  Every implementation of
  [`AbstractPC`](@ref) should extend this function with a new method

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
function calcPC(pc::AbstractPC, mesh::AbstractMesh, sbp::AbstractSBP,
                eqn::AbstractSolutionData, opts::Dict, ctx_residual, t)

  error("reached AbstractPC calcPC(), did you forget to define calcPC() for your AbstracPC implementation?")
end

"""
  Applies the preconditioner, ie. x = inv(Ap)*b, where Ap is the approximation
  to the matrix A.  Note that Ap itself may not be available for some
  preconditioners, hence there is only an API for applying inv(AP),
  not Ap itself.

  Every implementation of [`AbstractPC`](@ref) should extend this function with
  a new method.

  **Inputs**

   * pc: the [`AbstractPC`](@ref) implementation.
   * b: an AbstractVector (although never a PetscVec)

  **Inputs/Outputs**

   * x: vector overwritten with results
"""
function applyPC(pc::AbstractPC, x::AbstractVector, b::AbstractVector)


  error("reached AbstractPC applyPC(), did you forget to define applyPC() for your AbstracPC implementation?")

end

"""
  Applies the transpose of the preconditioner, ie. x = inv(Ap).'*b.
  Similar to [`applyPC`](@ref), see that function for details.

  Note that not every method supports this.  In particular, Petsc
  matrix-explicit PCs don't currently expose this (although they could)

  **Inputs**

   * pc: the [`AbstractPC`](@ref)
   * b: input vector

  **Inputs/Outputs**

   * x: vector overwritten with result
"""
function applyPCTranspose(pc::AbstractPC, x::AbstractVector, b::AbstractVector)

  
  error("reached AbstractPC applyPCTranspose(), did you forget to define applyPCTranspose() for your AbstracPC implementation?")

end


"""
  Preconditioner type for direct solve.
  Do not use with PetscLinearOperator.
"""
type PCNone <: AbstractPC
end

# LinearOperator interface
"""
  Abstract supertype of all linear operators used for A when solving Ax = b.

  The purpose of this type is to provide a consistent interface for different
  types of linear operators.  This type really combines two notions: what the
  type of the linear operator is, and how it should be solved.

  Any implementation of this type should subtype the appropriate catagory
  of: [`AbstractDenseLO`](@ref), [`AbstractSparseDirectLO`](@ref),
  [`AbstractIterativeMatLO`](@ref), [`AbstractIterativeMatFreeLO`](@ref)

  Note that matrix-explicit implementations can often write a single function
  for all these cases if using an matrix interface functions that are defined
  for all the matrix types.  See [`MatExplicitLO`](@ref)
"""
abstract AbstractLinearOperator

#TODO: doc these

"""
  Linear operator type for Dense matrices.  This is generally used only for
  debugging.
"""
abstract AbstractDenseLO <: AbstractLinearOperator

"""
  Linear operator type for `SparseMatrixCSC` matrices, which use a direct
  solver.
"""
abstract AbstractSparseDirectLO <: AbstractLinearOperator

"""
  Linear operator type for Petsc matrix-explicit.
"""
abstract AbstractIterativeMatLO <: AbstractLinearOperator

"""
  Linear operator type for Petsc matrix-free.
"""
abstract AbstractIterativeMatFreeLO <: AbstractLinearOperator

"""
  Useful union for all the matrix-explicit linear operator types.
  Because matrices have a small set of common interface functions, it is
  sometimes possible to write a single function that works on all the different
  types of matrices.
"""
typealias MatExplicitLO Union{AbstractDenseLO, AbstractSparseDirectLO, AbstractIterativeMatLO}

"""
  This function calculates the linear operator.  Every implementation of
  [`AbstractLinearOperator`](@ref) should extend this function with a new
  method.  For matrix-free operators, this function must exist but need
  not perform any actions.

  **Inputs**

   * lo: the AbstractLinearOperator implementation (updated with new
         matrix).
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
    this function is not called again before the next solve

"""
function calcLinearOperator(lo::AbstractLinearOperator, mesh::AbstractMesh,
                            sbp::AbstractSBP, eqn::AbstractSolutionData,
                            opts::Dict, ctx_residual, t)

  error("reached AbstractLinearOperator calcLinearOperator(), did you forget to define calcLinearOperator() for your AbstracLinearOperator implementation?")
end

"""
  Applies the linear operator, ie. , Ax = b
  
  Matrix-explicit implementations [`AbstractLinearOperator`](@ref) do not have
  to implement this function, though matrix-free implementations must extend
  it with a new method.

  **Inputs**

   * lo: the [`AbstractLinearOperator`](@ref) implementation.
   * x: an AbstractVector (although never a PetscVec)

  **Inputs/Outputs**

   * b: vector overwritten with results
"""
function applyLinearOperator(lo::AbstractLinearOperator, x::AbstractVector, 
                             b::AbstractVector)


  error("reached AbstractLinearOperator applyLinearOperator(), did you forget to define applyLinearOperator() for your AbstracLinearOperator implementation?")

end

"""
  Applies the transpose of the linear operator, ie. A.'*x = b
  Similar to [`applyLinearOperator`](@ref), see that function for details.

  Note that not every method supports this.  In particular, Petsc
  matrix-free LinearOperators don't currently expose this 
  (although they could with enough reverse-mode)

  **Inputs**

   * lo: the [`AbstractLinearOperator`](@ref)
   * x: input vector

  **Inputs/Outputs**

   * b: vector overwritten with result
"""
function applyLinearOperatorTranspose(lo::AbstractLinearOperator, 
                                      x::AbstractVector, b::AbstractVector)

  
  error("reached AbstractLinearOperator applyLinearOperatorTranspose(), did you forget to define applyLinearOperatorTranspose() for your AbstracLinearOperator implementation?")

end



end  # end module
