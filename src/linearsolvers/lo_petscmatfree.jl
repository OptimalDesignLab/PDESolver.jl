# linear operator implementation for Petsc matrix-free

"""
  Petsc matrix-free linear operator.

  **Public Fields**

   * none
"""
type PetscMatFreeLO <: AbstractPetscMatFreeLO
  A::PetscMat  # shell matrix
  xtmp::PetscVec  # these are shared with the PC if possible
  btmp::PetscVec
  ctx
  is_finalized::Bool
 
  nsolves::Int
  ntsolves::Int

  # MPI stuff
  comm::MPI.Comm
  myrank::Int
  commsize::Int
end

"""
  Outer constructor for PetscMatFreeLO

  **Inputs**

   * pc: a Petsc preconditioner of some kind
   * mesh
   * sbp
   * eqn
   * opts
"""
function PetscMatFreeLO(pc::Union{AbstractPetscMatPC, AbstractPetscMatFreePC},
                    mesh::AbstractMesh,
                    sbp::AbstractSBP, eqn::AbstractSolutionData, opts::Dict)

  A = createPetscMatShell(mesh, sbp, eqn, opts)

  # set operations:
  fptr = cfunction(applyLinearOperator_wrapper, PetscErrorCode, (PetscMat, PetscVec, PetscVec))
  MatShellSetOperation(A, PETSc2.MATOP_MULT, fptr)

  fptr = cfunction(applyLinearOperatorTranspose_wrapper, PetscErrorCode, (PetscMat, PetscVec, PetscVec))
  MatShellSetOperation(A, PETSc2.MATOP_MULT_TRANSPOSE, fptr)


  pc2 = getBasePC(pc)
  if typeof(pc2) <: PetscMatPC  # if pc2 has Petsc vectors inside it
    xtmp = pc2.xtmp
    btmp = pc2.btmp
  else
    xtmp = createPetscVec(mesh, sbp, eqn, opts)
    btmp = createPetscVec(mesh, sbp, eqn, opts)
  end

  ctx = C_NULL  # this gets set later
  is_finalized = false
  nsolves = 0
  ntsolves = 0
  comm = eqn.comm
  myrank = eqn.myrank
  commsize = eqn.commsize


  return PetscMatFreeLO(A, xtmp, btmp, ctx, is_finalized, nsolves, ntsolves,
                        comm, myrank, commsize)
end

function free(lo::PetscMatFreeLO)

  if !lo.is_finalized
    if lo.A.pobj != C_NULL
      PetscDestroy(lo.A)
      lo.A.pobj = C_NULL
    end

    if lo.xtmp.pobj != C_NULL
      PetscDestroy(lo.xtmp)
      lo.xtmp.pobj = C_NULL
    end

    if lo.btmp.pobj != C_NULL
      PetscDestroy(lo.btmp)
      lo.btmp.pobj = C_NULL
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
  @assert typeof(lo2) <: PetscMatFreeLO
  lo2.ctx = (mesh, sbp, eqn, opts, lo, ctx_residual, t)
  MatShellSetContext(lo2.A, pointer_from_objref(lo2.ctx))

  return nothing
end

function applyLinearOperator_wrapper(A::PetscMat, x::PetscVec, b::PetscVec)

  ctx_ptr = MatShellGetContext(A)
  @assert ctx_ptr != C_NULL
  ctx = unsafe_pointer_to_objref(ctx_ptr)
  checkLOCtx(ctx)
  mesh = ctx[1]
  sbp = ctx[2]
  eqn = ctx[3]
  opts = ctx[4]
  lo = ctx[5]  # user defined LO
  ctx_residual = ctx[6]
  t = ctx[7]

  btmp = VecGetArray(b)
  xtmp = VecGetArrayRead(x)

  applyLinearOperator(lo, mesh, sbp, eqn, opts, ctx_residual, t, xtmp, btmp)

  VecRestoreArray(b, btmp)
  VecRestoreArrayRead(x, xtmp)

  return PetscErrorCode(0)
end

function applyLinearOperatorTranspose_wrapper(A::PetscMat, x::PetscVec,
                                              b::PetscVec)

  ctx_ptr = MatShellGetContext(A)
  @assert ctx_ptr != C_NULL
  ctx = unsafe_pointer_to_objref(ctx_ptr)
  checkLOCtx(ctx)
  mesh = ctx[1]
  sbp = ctx[2]
  eqn = ctx[3]
  opts = ctx[4]
  lo = ctx[5]  # user defined LO
  ctx_residual = ctx[6]
  t = ctx[7]

  btmp = VecGetArray(b)
  xtmp = VecGetArrayRead(x)

  applyLinearOperatorTranspose(lo, mesh, sbp, eqn, opts, ctx_residual, t, 
                               xtmp, btmp)

  VecRestoreArray(b, btmp)
  VecRestoreArrayRead(x, xtmp)

  return PetscErrorCode(0)
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
  @assert typeof(mesh) <: AbstractMesh
  @assert typeof(sbp) <: AbstractSBP
  @assert typeof(eqn) <: AbstractSolutionData
  @assert typeof(lo) <: AbstractPetscMatFreeLO
  @assert typeof(lo2) <: PetscMatFreeLO
  @assert typeof(t) <: Number

  return nothing
end

function getBaseLO(lo::PetscMatFreeLO)

  # this is the bottom of the recursion tree
  return lo
end
