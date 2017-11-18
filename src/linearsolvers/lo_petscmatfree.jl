# linear operator implementation for Petsc matrix-free

"""
  Petsc matrix-free linear operator
"""
type PetscMatFreeLO <: AbstractPetscMatFreeLO
  A::PetscMat  # shell matrix
  xtmp::PetscVec  # these are shared with the PC if possible
  btmp::PetscVec
  ctx
  is_finalized::Bool
end

function PetscMatLO(pc::AbstracPetscMatPC, mesh::AbstractMesh,
                    sbp::AbstractSBP, eqn::AbstractSolutionData, opts::Dict)

  A = createPetscMatShell(mesh, sbp, eqn, opts)

  # set operations:
  fptr = cfunction(applyLinearOperator_wrapper, PetscErrorCode, (PetscMat, PetscVec, PetscVec))
  MatShellSetOperation(A, PETSc.MATOP_MULT, fptr)

  fptr = cfunction(applyLinearOperator_wrapperTranspose, PetscErrorCode, (PetscMat, PetscVec, PetscVec))
  MatShellSetOperation(A, PETSc.MATOP_MULT_TRANSPOSE, fptr)


  pc2 = getBasePC(pc)
  if pc2 <: PetscMatPC  # if pc2 has Petsc vectors inside it
    xtmp = pc2.xtmp
    btmp = pc2.btmp
  else
    xtmp = createPetscVec(mesh, sbp, eqn, opts)
    btmp = createPetscVec(mesh, sbp, eqn, opts)
  end

  ctx = C_NULL  # this gets set later
  is_finalized = false

  return new(A, xtmp, btmp, ctx, is_finalized)
end

function free(lo::PetscMatFreeLO)

  if !lo.is_finalized
    if lo.A.pobj != C_NULL
      PetscDestroy(lo.A)
      lo.A.pobj = C_NULL
    end

    if lo.A.xtmp.pobj != C_NULL
      PetscDestroy(lo.xtmp)
      lo.xtmp.pobj = C_NULL
    end

    if lo.btmp.pboj != C_NULL
      PetscDestroy(lo.btmp)
      lo.btmp.pboj = C_NULL
    end
  end

  lo.is_finalized = true

  return nothing
end


"""
  This function sets the ctx for the underlying Petsc matrix object.  The
  end result is that the mesh, sbp, eqn, and opts that are passed into this
  function will be passed to the user-defined [`applyLinearOperator`](@ref)
  during the next linear solve.  See that function for a description of
  the arguments.

  **Inputs**

   * lo: the user defined linear operator
   * mesh
   * sbp
   * eqn
   * opts
   * ctx_residual
   * t
"""
function setLOCtx(lo::AbstractPetscMatFreeLO, mesh::AbstractMesh,
                  sbp::AbstractSBP,
                  eqn::AbstractSolutionData, opts::Dict, ctx_residual, t)

  lo2 = getBaseLO(lo)
  @assert lo2 <: PetscMatFreeLO
  lo2.ctx = (mesh, sbp, eqn, opts, lo, ctx_residual, t)
  MatShellSetContext(lo2.A, pointer(lo2.ctx))

  return nothing
end

function applyLinearOperator_wrapper(A::PetscMat, x::PetscVec, b::PetscVec)

  ctx = MatShellGetContext(A)
  checkLOCtx(ctx)
  mesh = ctx[1]
  sbp = ctx[2]
  eqn = ctx[3]
  opts = ctx[4]
  lo = ctx[5]  # user defined LO
  ctx_residual = ctx[6]
  t = ctx[7]

  btmp, b_ptr = PetscVecGetArray(b)
  xtmp, x_ptr = PetscVecGetArrayRead(x)

  applyLinearOperator(lo, mesh, sbp, eqn, opts, ctx_residual, t, xtmp, btmp)

  PetscVecRestoreArray(b, b_ptr)
  PetscVecRestoreArrayRead(x, x_ptr)

  return nothing
end

function applyLinearOperatorTranspose_wrapper(A::PetscMat, x::PetscVec,
                                              b::PetscVec)

  ctx = MatShellGetContext(A)
  checkLOCtx(ctx)
  mesh = ctx[1]
  sbp = ctx[2]
  eqn = ctx[3]
  opts = ctx[4]
  lo = ctx[5]  # user defined LO
  ctx_residual = ctx[6]
  t = ctx[7]

  btmp, b_ptr = PetscVecGetArray(b)
  xtmp, x_ptr = PetscVecGetArrayRead(x)

  applyLinearOperatorTranspose(lo, mesh, sbp, eqn, opts, ctx_residual, t, 
                               xtmp, btmp)

  PetscVecRestoreArray(b, b_ptr)
  PetscVecRestoreArrayRead(x, x_ptr)

  return nothing
end


"""
  This function verifies the [`PetscMatFreeLO`](@ref) ctx is valid, throwing
  an exception otherwise.

  Currently this only checks whether the members of the tuple have the right
  abstract supertypes

  **Inputs**

   * ctx: the ctx retrieved from the Petsc matrix shell object
"""
function checkLOCtx(ctx)

  mesh = ctx[1]
  sbp = ctx[2]
  eqn = ctx[3]
  opts = ctx[4]
  lo = ctx[5]  # user defined LO
  ctx_residual = ctx[6]
  t = ctx[7]

  lo2 = getBaseLO(lo)
  @assert mesh <: AbstractMesh
  @assert sbp <: AbstractSBP
  @assert eqn <: AbstractSolutionData
  @assert lo <: AbstractPetscMatFreeLO
  @assert lo2 <: PetscMatFreeLO
  @assert t <: Number

  return nothing
end

function getBaseLinearOperator(lo::PetscMatFreeLO)

  # this is the bottom of the recursion tree
  return lo
end
