# Petsc matrix-free implementation

"""
  Type for managing a Petsc matrix-free preconditioner.  Actual preconditioners
  can be built on top of this.  Because matrix-free preconditioners do not have
  a particular form, this type requires the preconditioner to implement the
  entire [`AbstractPC`](@ref) interface.  The benefit of using this type to
  build a preconditioner is that is automatically exposes the preconditioner
  to Petsc.

  The [`calcPC`](@ref) function the user defines must call [`setPCCtx`](@ref).

"""
type PetscMatFreePC <: AbstracPetscMatFreePC
  pc::PC
  ctx  # Petsc PC ctx
  is_setup::Bool # needed for consistency with PetscMatPC
end

#TODO: outer constructor

#TODO: rename this to calcPC?
"""
  This function should be called by the users [`calcPC`](@ref) function.
  The arguments passed to this function are passed to [`applyPC`](@ref)
  during the next linear solve.  See that function for a description of
  the arguments

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
  @assert pc2 <: PetscMatFreePC
  pc2.ctx = (mesh, sbp, eqn, opts, pc, ctx_residual, t)
  PCShellSetContex(pc2, pointer(pc2.ctx))

  return nothing
end

function calcPC_wrapper(_pc::PC)

  ctx = unsafe_pointer_to_objref(PCShellGetContex(_pc))
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

  ctx = unsafe_pointer_to_objref(PCShellGetContex(_pc))
  checkPCCtx(ctx)
  mesh = ctx[1]
  sbp = ctx[2]
  eqn = ctx[3]
  opts = ctx[4]
  pc = ctx[5]  # the user PC
  ctx_residual = ctx[6]
  t = ctx[7]

  # get the local parts of the vectors
  btmp, b_ptr = PetscVecGetArrayRead(b)
  xtmp, x_ptr = PetscVecGetArray(x)


  applyPC(pc, mesh, sbp, eqn, opts, t, btmp, xtmp)

  PetscVecRestoreArrayRead(b, b_ptr)
  PetscVecRestoreArray(x, x_ptr)

  return PetscErrorCode(0)
end

function applyPCTranspose_wrapper(_pc::PC, b::PetscVec, x::PetscVec)

  ctx = unsafe_pointer_to_objref(PCShellGetContex(_pc))
  checkPCCtx(ctx)
  mesh = ctx[1]
  sbp = ctx[2]
  eqn = ctx[3]
  opts = ctx[4]
  pc = ctx[5]  # the user PC
  ctx_residual = ctx[6]
  t = ctx[7]

  # get the local parts of the vectors
  btmp, b_ptr = PetscVecGetArrayRead(b)
  xtmp, x_ptr = PetscVecGetArray(x)


  applyPCTranspose(pc, mesh, sbp, eqn, opts, t, btmp, xtmp)

  PetscVecRestoreArrayRead(b, b_ptr)
  PetscVecRestoreArray(x, x_ptr)

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

  @assert mesh <: AbstractMesh
  @assert sbp <: AbstractSBP
  @assert eqn <: AbstractSolutionData
  @assert opts <: Dict
  @assert pc <: AbstracPC
  @assert pc2 <: PetscMatFreePC
  @assert t <: Number

  # verify equality with pc2.ctx ?

  return nothing
end

function getBasePC(pc::PetscMatFreePC)

  # this is the bottom of the recursion tree
  return pc
end
