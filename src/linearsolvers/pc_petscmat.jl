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
type PetscMatPC <: AbstractPC
  pc::PC  # Petsc PC object
  Ap::PetscMat  # Petsc Mat object
  xtmp::PetscVec  # reusable temporary vector
  btmp::PetscVec
  pc_setup::Bool  # is PC already set up
end

#TODO: outer constructor, have a flag to avoid allocating xtmp and btmp

function calcPC(pc::PetscMatPC, mesh::AbstractMesh, sbp::AbstractSBP,
                eqn::AbstractSolutionData, opts::Dict, ctx_residual, t)

  # compute the jacobian here
  physicsJac(mesh, sbp, eqn, opts, pc.Ap, ctx_residual, t)
 
  # don't setup the PC here because PC built on top of this one might
  # modify Ap after calling this function
  pc.pc_setup = false

  return nothing
end


function applyPC(pc::AbstractPC, mesh::AbstractMesh, sbp::AbstractSBP,
                 eqn::AbstractSolutionData, opts::Dict, t, b::AbstractVector, 
                 x::AbstractVector)
  
  if !pc.pc_setup
    setupPC(pc)
  end


  # copy into temp vector
  btmp, b_ptr = PetscVecGetArray(pc.btmp)
  copy!(btmp, b)
  PetscVecRestoreArray(btmp, b_ptr)

  # call Petsc PCApply
  PetscPCApply(pc.pc, pc.btmp, pc.xtmp)

  # copy back to x
  xtmp, x_ptr = PetscVecGetArrayRead(pc.xtmp)
  copy!(x, xtmp)
  PetscVecRestoreArrayRead(xtmp, x_ptr)

  return nothing
end


function applyPCTranspose(pc::PetscMatPC, mesh::AbstractMesh, sbp::AbstractSBP,
                 eqn::AbstractSolutionData, opts::Dict, t, b::AbstractVector, 
                 x::AbstractVector)


  if !PCApplyTransposeExists(pc)
    throw(ErrorException("PCApplyTranspose not defined for this PC"))
  end

  if !pc.pc_setup
    setupPC(pc)
  end

  # copy into temp vector
  btmp, b_ptr = PetscVecGetArray(pc.btmp)
  copy!(btmp, b)
  PetscVecRestoreArray(btmp, b_ptr)

  # call Petsc PCApplyTranspose
  PetscPCApplyTranspose(pc.pc, pc.btmp, pc.xtmp)

  # copy back to x
  xtmp, x_ptr = PetscVecGetArrayRead(pc.xtmp)
  copy!(x, xtmp)
  PetscVecRestoreArrayRead(xtmp, x_ptr)

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
  pc.pc_setup = true

  return nothing
end
