# Petsc matrix-free implementation

"""
  Type for managing a Petsc matrix-free preconditioner.  Actual preconditioners
  can be built on top of this.  Because matrix-free preconditioners do not have
  a particular form, this type requires the preconditioner to implement the
  entire [`AbstractPC`](@ref) interface.  The benefit of using this type to
  build a preconditioner is that is automatically exposes the preconditioner
  to Petsc.

  The [`calcPC`](@ref) function the user defines must call [`setPCCtx`](@ref).

  **Public Fields**

   * none
"""
mutable struct PetscMatFreePC <: AbstractPetscMatFreePC
  pc::PC
  ctx  # Petsc PC ctx
  is_shared::Bool
  is_finalized::Bool

  nsetups::Int
  napplies::Int
  ntapplies::Int

  # MPI stuff
  comm::MPI.Comm
  myrank::Int
  commsize::Int
end

"""
  Outer constructor

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
"""
function PetscMatFreePC(mesh::AbstractMesh, sbp::AbstractSBP,
                    eqn::AbstractSolutionData, opts::Dict)

  pc = createPetscPC(mesh, sbp, eqn, opts)

  # register wrapper routines with Petsc
  fptr = cfunction(applyPC_wrapper, PetscErrorCode, (PC, PetscVec, PetscVec))
  PCShellSetApply(pc, fptr)

  fptr = cfunction(applyPCTranspose_wrapper, PetscErrorCode, (PC, PetscVec, PetscVec))
  PCShellSetApplyTranspose(pc, fptr)

  fptr = cfunction(calcPC_wrapper, PetscErrorCode, (PC,))
  PCShellSetSetUp(pc, fptr)


  ctx = ()
  is_shared = false
  is_finalized = false

  nsetups = 0
  napplies = 0
  ntapplies = 0

  comm = eqn.comm
  myrank = eqn.myrank
  commsize = eqn.commsize

  # only update the PC when explicitly requested
  PCSetReusePreconditioner(pc, PETSC_TRUE)

  return PetscMatFreePC(pc, ctx, is_shared, is_finalized, nsetups, napplies,
                        ntapplies, comm, myrank, commsize)
end

function free(pc::PetscMatFreePC)

  if !pc.is_finalized
    if pc.pc.pobj != C_NULL
      PetscDestroy(pc.pc)
      pc.pc.pobj = C_NULL
    end
  end

  pc.is_finalized = true

  return nothing
end


"""
  This function should be called by the users [`calcPC`](@ref) function.
  The arguments passed to this function are passed to [`applyPC`](@ref)
  during the next preconditioner application.
  See that function for a description of the arguments

  **Inputs**

   * pc:
   * mesh
   * sbp
   * eqn
   * opts
   * ctx_residual
   * t
"""
function setPCCtx(pc::AbstractPetscMatFreePC, mesh::AbstractMesh,
                  sbp::AbstractSBP,
                  eqn::AbstractSolutionData, opts::Dict, ctx_residual, t)
# users *must* call this function from within their calcPC function

  pc2 = getBasePC(pc)
  @assert typeof(pc2) <: PetscMatFreePC
  pc2.ctx = (mesh, sbp, eqn, opts, pc, ctx_residual, t)
  PCShellSetContext(pc2.pc, pointer_from_objref(pc2.ctx))

  return nothing
end

function calcPC_wrapper(_pc::PC)

  ctx_ptr = PCShellGetContext(_pc)
  @assert ctx_ptr != C_NULL
  ctx = unsafe_pointer_to_objref(ctx_ptr)
  checkPCCtx(ctx)
  mesh = ctx[1]
  sbp = ctx[2]
  eqn = ctx[3]
  opts = ctx[4]
  pc = ctx[5]  # the user PC
  ctx_residual = ctx[6]
  t = ctx[7]

  calcPC(pc, mesh, sbp, eqn, opts, ctx_residual, t)

  return PetscErrorCode(0)
end


# wrapper for Petsc
function applyPC_wrapper(_pc::PC, b::PetscVec, x::PetscVec)

  ctx = unsafe_pointer_to_objref(PCShellGetContext(_pc))
  checkPCCtx(ctx)
  mesh = ctx[1]
  sbp = ctx[2]
  eqn = ctx[3]
  opts = ctx[4]
  pc = ctx[5]  # the user PC
  ctx_residual = ctx[6]
  t = ctx[7]


  # get the local parts of the vectors
  btmp = VecGetArrayRead(b)
  xtmp = VecGetArray(x)


  applyPC(pc, mesh, sbp, eqn, opts, t, btmp, xtmp)

  VecRestoreArrayRead(b, btmp)
  VecRestoreArray(x, xtmp)

  return PetscErrorCode(0)
end

function applyPCTranspose_wrapper(_pc::PC, b::PetscVec, x::PetscVec)

  ctx_ptr = PCShellGetContext(_pc)
  @assert ctx_ptr != C_NULL
  ctx = unsafe_pointer_to_objref(ctx_ptr)
  checkPCCtx(ctx)
  mesh = ctx[1]
  sbp = ctx[2]
  eqn = ctx[3]
  opts = ctx[4]
  pc = ctx[5]  # the user PC
  ctx_residual = ctx[6]
  t = ctx[7]

  # get the local parts of the vectors
  btmp = VecGetArrayRead(b)
  xtmp = VecGetArray(x)


  applyPCTranspose(pc, mesh, sbp, eqn, opts, t, btmp, xtmp)

  VecRestoreArrayRead(b, btmp)
  VecRestoreArray(x, xtmp)

  return PetscErrorCode(0)
end



"""
  This function verifies the [`PetscMatFreePC`](@ref) ctx is valid, throwing
  an exception otherwise.

  Currently this only checks whether the members of the tuple have the right
  abstract supertypes

  **Inputs**

   * ctx: the ctx retrieved from the Petsc matrix-free preconditioner object
"""
function checkPCCtx(ctx)

  mesh = ctx[1]
  sbp = ctx[2]
  eqn = ctx[3]
  opts = ctx[4]
  pc = ctx[5]  # the user PC
#  pc2 = ctx[6]  # PetscMatFreePC
  ctx_residual = ctx[6]
  t = ctx[7]

  pc2 = getBasePC(pc)

  @assert typeof(mesh) <: AbstractMesh
  @assert typeof(sbp) <: AbstractSBP
  @assert typeof(eqn) <: AbstractSolutionData
  @assert typeof(opts) <: Dict
  @assert typeof(pc) <: AbstractPC
  @assert typeof(pc2) <: PetscMatFreePC
  @assert typeof(t) <: Number

  # verify equality with pc2.ctx ?

  return nothing
end

function getBasePC(pc::PetscMatFreePC)

  # this is the bottom of the recursion tree
  return pc
end
