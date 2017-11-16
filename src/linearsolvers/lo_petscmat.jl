# Linear operator implementation for Petsc matrix-explicit

#TODO: rename the abstract type to AbstractPetscMatLO?
type PetscMatLO <: AbstractPetscMatLO
  A::PetscMat
  xtmp::PetscVec  # these are shared with the PC if possible
  btmp::PetscVec
end

function calcLinearOperator(lo::PetscMatLO, mesh::AbstractMesh,
                            sbp::AbstractSBP, eqn::AbstractSolutionData,
                            opts::Dict, ctx_residual, t)


  physicsJac(mesh, sbp, eqn, opts, pc.Ap, ctx_residual, t)

  return nothing
end


function applyLinearOperator(lo::AbstractPetscMatLO, mesh::AbstractMesh,
                             sbp::AbstractSBP, eqn::AbstractSolutionData,
                             opts::Dict, ctx_residual, t, x::AbstractVector, 
                             b::AbstractVector)

  lo2 = getBaseLinearOperator(lo)
  xtmp, x_ptr = PetscVecGetArray(lo2.xtmp)
  copy!(xtmp, x)
  PetscVecRestoreArray(lo2.xtmp, x_ptr)

  PetscMatMult(lo2.A, lo2.xtmp, lo2.btmp)

  btmp, b_ptr = PetscVecGetArrayRead(lo2.btmp)
  copy!(b, btmp)
  PetscVecRestoreArrayRead(lo2.btmp, b_ptr)

  return nothing
end


function applyLinearOperatorTranspose(lo::AbstractPetscMatLO, 
                             mesh::AbstractMesh, sbp::AbstractSBP,
                             eqn::AbstractSolutionData, opts::Dict, 
                             ctx_residual, t, x::AbstractVector, 
                             b::AbstractVector)
  lo2 = getBaseLinearOperator(lo)
  xtmp, x_ptr = PetscVecGetArray(lo2.xtmp)
  copy!(xtmp, x)
  PetscVecRestoreArray(lo2.xtmp, x_ptr)

  PetscMatMultTranspose(lo2.A, lo2.xtmp, lo2.btmp)

  btmp, b_ptr = PetscVecGetArrayRead(lo2.btmp)
  copy!(b, btmp)
  PetscVecRestoreArrayRead(lo2.btmp, b_ptr)

  return nothing
end


function getBaseLinearOperator(lo::PetscMatLO)

  # this is the bottom of the recursion tree
  return lo
end
