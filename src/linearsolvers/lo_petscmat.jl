# Linear operator implementation for Petsc matrix-explicit

#TODO: rename the abstract type to AbstractPetscMatLO?
type PetscMatLO <: AbstractIterativeMatLO
  A::PetscMat
  xtmp::PetscVec  # these are shared with the PC if possible
  btmp::PetscVec
end

function calcLinearOperator(lo::AbstractLinearOperator, mesh::AbstractMesh,
                            sbp::AbstractSBP, eqn::AbstractSolutionData,
                            opts::Dict, ctx_residual, t)


  physicsJac(mesh, sbp, eqn, opts, pc.Ap, ctx_residual, t)

  return nothing
end


function applyLinearOperator(lo::PetscMatLO, mesh::AbstractMesh,
                             sbp::AbstractSBP, eqn::AbstractSolutionData,
                             opts::Dict, ctx_residual, t, x::AbstractVector, 
                             b::AbstractVector)

  xtmp, x_ptr = PetscVecGetArray(lo.xtmp)
  copy!(xtmp, x)
  PetscVecRestoreArray(lo.xtmp, x_ptr)

  PetscMatMult(lo.A, lo.xtmp, lo.btmp)

  btmp, b_ptr = PetscVecGetArrayRead(lo.btmp)
  copy!(b, btmp)
  PetscVecRestoreArrayRead(lo.btmp, b_ptr)

  return nothing
end


function applyLinearOperatorTranspose(lo::PetscMatLO, 
                             mesh::AbstractMesh, sbp::AbstractSBP,
                             eqn::AbstractSolutionData, opts::Dict, 
                             ctx_residual, t, x::AbstractVector, 
                             b::AbstractVector)

  xtmp, x_ptr = PetscVecGetArray(lo.xtmp)
  copy!(xtmp, x)
  PetscVecRestoreArray(lo.xtmp, x_ptr)

  PetscMatMultTranspose(lo.A, lo.xtmp, lo.btmp)

  btmp, b_ptr = PetscVecGetArrayRead(lo.btmp)
  copy!(b, btmp)
  PetscVecRestoreArrayRead(lo.btmp, b_ptr)

  return nothing
end


function getBaseLinearOperator(lo::PetscMatLO)

  # this is the bottom of the recursion tree
  return lo
end
