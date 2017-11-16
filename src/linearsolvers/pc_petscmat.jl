# AbstractPC type of Petsc matrix-explicit

"""
  [`AbstractPC`](@ref) implementation for Petsc matrix-explicit preconditioners.

  Methods in NonlinearSolvers that want to use a Petsc matrix-explicit
  preconditioner are encouraged to use composition, ie. create their own
  [`AbstractPC`](@ref) type and include a PetscMatPC as a field.  This allows
  the nonlinear method to use the functions defined on PetscMatPC as a building.

  **Field**

   * pc: a Petsc PC object
   * Ap: a PetscMat object used to calculate the preconditioner
"""
type PetscMatPC <: AbstractPetscMatPC
  pc::PC  # Petsc PC object
  Ap::PetscMat  # Petsc Mat object
  xtmp::PetscVec  # reusable temporary vector
  btmp::PetscVec
  is_setup::Bool  # is PC already set up
end

#TODO: outer constructor, have a flag to avoid allocating xtmp and btmp

function calcPC(pc::PetscMatPC, mesh::AbstractMesh, sbp::AbstractSBP,
                eqn::AbstractSolutionData, opts::Dict, ctx_residual, t)

  # compute the jacobian here
  physicsJac(mesh, sbp, eqn, opts, pc.Ap, ctx_residual, t)
 
  # don't setup the PC here because PC built on top of this one might
  # modify Ap after calling this function
  pc.is_setup = false

  return nothing
end


function applyPC(pc::AbstractPetscMatPC, mesh::AbstractMesh, sbp::AbstractSBP,
                 eqn::AbstractSolutionData, opts::Dict, t, b::AbstractVector, 
                 x::AbstractVector)
  pc2 = getBasePC(pc) 
  if !pc2.is_setup
    setupPC(pc2)
  end


  # copy into temp vector
  btmp, b_ptr = PetscVecGetArray(pc2.btmp)
  copy!(btmp, b)
  PetscVecRestoreArray(pc2.btmp, b_ptr)

  # call Petsc PCApply
  PetscPCApply(pc2.pc, pc2.btmp, pc2.xtmp)

  # copy back to x
  xtmp, x_ptr = PetscVecGetArrayRead(pc2.xtmp)
  copy!(x, xtmp)
  PetscVecRestoreArrayRead(pc2.xtmp, x_ptr)

  return nothing
end


function applyPCTranspose(pc::AbstracPetscMatPC, mesh::AbstractMesh,
                 sbp::AbstractSBP,
                 eqn::AbstractSolutionData, opts::Dict, t,
                 b::AbstractVector, x::AbstractVector)


  pc2 = getBasePC(pc)
  @assert pc2 <: PetscMatPC

  if !PCApplyTransposeExists(pc2.pc)
    ptype = PCGetType(pc2.pc)
    throw(ErrorException("PCApplyTranspose not defined for PC $ptype"))
  end

  if !pc2.is_setup
    setupPC(pc2)
  end

  # copy into temp vector
  btmp, b_ptr = PetscVecGetArray(pc2.btmp)
  copy!(btmp, b)
  PetscVecRestoreArray(pc2.btmp, b_ptr)

  # call Petsc PCApplyTranspose
  PetscPCApplyTranspose(pc2.pc, pc2.btmp, pc2.xtmp)

  # copy back to x
  xtmp, x_ptr = PetscVecGetArrayRead(pc2.xtmp)
  copy!(x, xtmp)
  PetscVecRestoreArrayRead(pc2.xtmp, x_ptr)

  return nothing
end


"""
  This internal function is used to setup the PC, including setting the flag.
  Users only need to do this if they did not call [`setupPC`](@ref)

  **Inputs**

   * pc: PetscMatPC object
"""
function setupPC(pc::PetscMatPC)

  PCSetUp(pc.pc)
  # this is potentially bad because Petsc will *never* recompute the 
  # preconditioner on its own.  Using higher level functionality like TS
  # or Petsc nonlinear solvers likely won't work in this case
  PCSetReusePreconditioner(pc)
  pc.is_setup = true

  return nothing
end

function getBasePC(pc::PetscMatPC)

  # this is the bottom of the recursion tree
  return pc
end
