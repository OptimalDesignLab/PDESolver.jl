# implementation of the StandardLinearSolver
# This is the interface algorithms should use for doing linear solves

"""
  Calculates the preconditioner for the linear solver.  Thsi preconditioner
  will be used for all linear solves until this function is called again.

  For direct solvers, this function calculates the linear operator itself.
  Prefer [`calcPCandLO`](@ref) whenever possible.

  **Inputs**

   * ls: StandardLinearSolver
   * mesh
   * sbp
   * eqn
   * opts
   * ctx_residual: the ctx required by [`physicsRhs`](@ref) like functions
   * t: current time
"""
function calcPC(ls::StandardLinearSolver, mesh::AbstractMesh, sbp::AbstractSBP,
                eqn::AbstractSolutionData, opts::Dict, ctx_residual, t)

  if ls.pc <: PCNONE
    calcLinearOperator(ls.lo, mesh, sbp, eqn, opts, ctx_residual, t)
  else
    calcPC(ls.pc, mesh, sbp, eqn, opts, ctx_residual, t)
  end

  return nothing
end

"""
  Calculates the linear operator.  Use this function only if you want to
  calculate the linear operator and not the preconditioner.
  Prefer [`calcPCandLO`](@ref), which avoids calculating the matrix twice if
  the preconditioner and linear operator share the same matrix

  **Inputs**

   * ls: StandardLinearSolver
   * mesh
   * sbp
   * eqn
   * opts
   * ctx_residual: the ctx required by [`physicsRhs`](@ref) like functions
   * t: current time
"""
function calcLinearOperator(ls::StandardLinearSolver, mesh::AbstractMesh,
                         sbp::AbstractSBP,
                         eqn::AbstractSolutionData, opts::Dict, ctx_residual, t)


  calcLinearOperator(ls.lo, mesh, sbp, eqn, opts, ctx_residual, t)

end


"""
  Calculates both the preconditioner and linear operator.  In the case where
  they share the matrix, the calculation is only performed once.  This function
  should be preferred to callcing [`calcPC`](@ref) and [`calcLO`](@ref) one
  after the other

  **Inputs**

   * ls: StandardLinearSolver
   * mesh
   * sbp
   * eqn
   * opts
   * ctx_residual: the ctx required by [`physicsRhs`](@ref) like functions
   * t: current time
"""
function calcPCandLO(ls::StandardLinearSolver, mesh::AbstractMesh,
                     sbp::AbstractSBP,
                     eqn::AbstractSolutionData, opts::Dict, ctx_residual, t)


  if lo.pc <: PCNONE
    calcLinearOperator(ls.lo, mesh, sbp, eqn, opts, ctx_residual, t)
  elseif ls.is_shared
    calcLinearOperator(ls.lo, mesh, sbp, eqn, opts, ctx_residual, t)
    pc2 = getBasePC(ls.pc)
    pc2.is_setup = false  # we don't call calcPC but this flag still needs to
                          # be set
  else
    calcPC(ls.pc, mesh, sbp, eqn, opts, ctx_residual, t)
    calcLinearOperator(ls.lo, mesh, sbp, eqn, opts, ctx_residual, t)
  end

  return nothing
end

"""
  Apply the preconditioner, ie. x = inv(Ap)*b, where Ap is an approximation
  to the linear operator used for preconditioner.  This function also works
  for matrix-free methods.

  Note that [`calcPC`](@ref) or [`calcPCandLO`](@ref) must be called before
  this function can be used.

  For direct methods, a linear solve is performed because there is no
  preconditioner.
"""
function applyPC(ls::StandardLinearSolver, mesh::AbstractMesh, sbp::AbstractSBP,
                 eqn::AbstractSolutionData, opts::Dict, t, b::AbstractVector, 
                 x::AbstractVector)

  if ls.pc <: PCNone
    # apply the most recent linear operator
    linearSolve(ls, b, x)
  else
    applyPC(ls.pc, mesh, sbp, eqn, opts, t, b, x)
  end

  return nothing
end

"""
  Like [`applyPC`](@ref), but applies the transpose (if possible, otherwise
  throws an error).
"""
function applyPCTranspose(ls::StandardLinearSolver, mesh::AbstractMesh,
                          sbp::AbstractSBP,
                          eqn::AbstractSolutionData, opts::Dict, t,
                          b::AbstractVector, x::AbstractVector)
  if ls.pc <: PCNone
    linearSolvetransposee(ls, b, x)
  else
    applPCTranspose(ls, mesh, sbp, eqn, opts, t, b, x)
  end

  return nothing
end


"""
  Solves the linear system Ax = b for x. This function does not recompute the
  precondtioner or linear operator,  The preconditioner and linear operator
  used are the ones calculated by the most recent call to
  [`calcPC`](@ref), [`calcLinearOperator`](@ref) (or [`calcPCandLO`](@ref).

  For Petsc matrices, this function does the final matrix assembly.

  **Inputs**

   * ls: the [`StandardLinearSolver`](@ref) object.
   * b: the right hand side vector (local process portion only)
  
  **Inputs/Outputs**

   * x: vector overwritten with result (local process portion only)

  Implementation Notes:

  The preconditioner and matrix factorization (for direct solves) might
  be computed during this function, but they will be computed at the state
  correspondong to the last call to the functions that calculate them.
"""
function linearSolve(ls::StandardLinearSolver, b::AbstractVector,
                     x::AbstractVector, verbose=5)

  # call the specific solver

  #TODO: have time_all return a struct
  tmp, t_solve, t_gc, alloc = @time_all _linearSolver(ls, b, x)

  @verbose5 if ls.myrank == 0
    println(BSTDOUT, "matrix solve:")
    print_time_all(BSTDOUT, t_solve, t_gc, alloc)
  end

  return nothing
end

function linearSolveTranspose(ls::StandardLinearSolver, b::AbstractVector,
                     x::AbstractVector, verbose=5)

  # call the specific solver

  #TODO: have time_all return a struct
  tmp, t_solve, t_gc, alloc = @time_all _linearSolver(ls, b, x, trans=true)

  @verbose5 if ls.myrank == 0
    println(BSTDOUT, "matrix transpose solve:")
    print_time_all(BSTDOUT, t_solve, t_gc, alloc)
  end

  return nothing
end

function _linearSolve{Tlo <: AbstractDenseLO, Tpc}(
                      ls::StandardLinearSolver{Tlo, Tpc},
                      b::AbstractVector, x::AbstractVector; trans=false)

  @assert typeof(lo.pc) <: PCNone

  lo2 = getBaseLO(ls.lo)

  # factor is needed
  if !lo2..is_setup
    getrf!(lo2.A. lo2.ipiv)
    lo2.is_setup = true
  end

  tchar = trans ? 'T' : 'N'

  # solve
  copy!(x, b)
  getrs2!(tchar, lo2.A, lo2.ipiv, x)

  return nothing
end


function _linearSolve{Tlo <: AbstractSparseDirectLO, Tpc}(
                      ls::StandardLinearSolver{Tlo, Tpc},
                      b::AbstractVector, x::AbstractVector; trans=false)

  @assert typeof(lo.pc) <: PCNone

  lo2 = getBaseLO(ls.lo)

  # compute factorization if needed
  if !lo2.is_setup
    umfpack_free_numeric(lo2.fac)  # free old factorization
    # note: the matrix stored in the factorization object must alias. lo2.A
    make_zerobased(lo2.A)
    umfpack_numeric!(lo2.fac)
    make_onebased(lo2.A)
    lo2.is_setup = true
  end

  tchar = trans ? UMFPACK_At : UMFPACK_A

  solve_suitesparse(lo2.fac, b, tchar, x)

  return nothing
end

function _linearSolve{Tlo <: PetscLO , Tpc}(
                      ls::StandardLinearSolver{Tlo, Tpc},
                      b::AbstractVector, x::AbstractVector; trans=false)

  myrank = ls.myrank
  pc2 = getBasePC(ls.pc)
  lo2 = getBaseLO(ls.lo)
  @assert !(typeof(pc2) <: PCNone)
  
  if !pc2.is_setup
    setupPC(pc2)
  end

  # prepare the data structures
  assemblePetscData(ls, b)

  # do the solve
  ksp =  ls.ksp
  KSPSetTolerances(ksp, ls.reltol, ls.abstol, ls.dtol, PetscInt(ls.itermax))

  if trans
    KSPSolveTranspose(ksp, lo2.btmp, lo2.xtmp)
  else
    KSPSolve(ksp, lo2.btmp, lo2.xtmp)
  end

  # print convergence info
  @mpi_master begin
    reason = KSPGetConvergedReason(ksp)
    println(BSTDOUT, "KSP converged reason = "< KSPConvergedREasonDict[reason])
    rnorm = KSPGetResidualNorm(ksp)
    @mpi_master println("Linear residual = ", rnorm)
  end
  # copy result back to x
  xtmp, x_ptr = PetscVecGetArrayRead(lo2.xtmp)
  copy!(x, xtmp)
  PetscVecRestoreArrayRead(lo2.xtmp, x_ptr)

  return nothing
end

"""
  Helper function for doing linear solves with Petsc matrices.  It assembles
  A, and Ap, and copies b into the local part of the petsc vector for b.
"""
function assemblePetscData(ls::StandardLinearSolver, b::AbstractVector)

  lo_matfree = isLOMatFree(ls)
  pc_matfree = isPCMatFRee(ls)

  pc2 = getBasePC(ls.pc)
  lo2 = getBaseLO(ls.lo)
  myrank = ls.myrank

  # assemble things
  if !lo_matfree
    PetscMatAssembleBegin(lo2.A, PETSC_MAT_FINAL_ASSEMBLY)
  end

  if !pc_matfree && !lo.shared_mat
    PetscMatAssembleBegin(pc2.Ap, PETSC_MAT_FINAL_ASSEMBLY)
  end

  # copy values into the vector
  btmp, b_ptr = PetscVecGetArray(lo2.btmp)
  copy!(btmp, b)
  PetscVecRestoreArray(lo2.btmp, b_ptr)

  # end assembly
  if !lo_matfree
    PetscMatAssembleEnd(lo2.A, PETSC_MAT_FINAL_ASSEMBLY)
    matinfo = PetscMatGetInfo(lo2.A)
    if matinfo.mallocs > 0.5  # if any mallocs
      println(BSTDERR, "Warning: non-zero number of mallocs for A on process $myrank: $(matinfo.mallocs) mallocs")
    end
  end

  if !pc_matfree && !lo.shared_mat
    PetscMatAssembleEnd(pc2.Ap, PETSC_MAT_FINAL_ASSEMBLY)
    matinfo = PetscMatGetInfo(pc2.Ap)
    if matinfo.mallocs > 0.5  # if any mallocs
      println(BSTDERR, "Warning: non-zero number of mallocs for Ap on process $myrank: $(matinfo.mallocs) mallocs")
    end
  end

  return nothing
end

"""
  Returns true if the linear operator is matrix free, false otherwise
"""
function isLOMatFree(ls::StandardLinearSolver)

  return ls.lo <: AbstractPetscMatFree
end

"""
  Returns true if the preconditioner is matrix free, false otherwise
"""
function isPCMatFree(ls::StandardLinearSolver)

  return ls.pc <: AbstractMatFreePC
end

"""
  Set the tolerances for iterative solves.  Has no effect for direct solves.
  Supplying a negative value results in retaining the original value.

  **Inputs**
  
   * ls: StandardLinearSolver
   * reltol: relative residual tolerance
   * abstol: absolute residual tolerance
   * dtol: divergence tolerance
   * itermax: maximum number of iterations
"""
function setTolerances(ls::StandardLinearSolver, reltol, abstol, dtol, itermax)


  if reltol > 0
    ls.reltol = reltol
  end

  if abstol > 0
    ls.abstol = abstol
  end

  if dtol > 0
    ls.dtol = dtol
  end

  if itermax > 0
    ls.itermax = itermax
  end

  return nothing
end
