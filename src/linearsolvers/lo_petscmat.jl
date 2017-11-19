# Linear operator implementation for Petsc matrix-explicit

"""
  Petsc matrix-explicit linear operator.
"""
type PetscMatLO <: AbstractPetscMatLO
  A::PetscMat
  xtmp::PetscVec  # these are shared with the PC if possible
  btmp::PetscVec
  is_finalized::Bool
end

function PetscMatLO(pc::AbstractPetscMatPC, mesh::AbstractMesh,
                    sbp::AbstractSBP, eqn::AbstractSolutionData, opts::Dict)

  pc2 = getBasePC(pc)
  if !opts["use_jac_precon"] && !(typeof(pc) <: PetscMatFreeLO) # share matrix if possible
    A = pc2.Ap
  else
    A = createPetscMat(mesh, sbp, eqn, opts)
  end

  if typeof(pc2) <: PetscMatPC  # if pc2 has Petsc vectors inside it
    xtmp = pc2.xtmp
    btmp = pc2.btmp
  else
    xtmp = createPetscVec(mesh, sbp, eqn, opts)
    btmp = createPetscVec(mesh, sbp, eqn, opts)
  end

  return new(A, xtmp, btmp, false)
end

function free(lo::PetscMatLO)

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


function calcLinearOperator(lo::PetscMatLO, mesh::AbstractMesh,
                            sbp::AbstractSBP, eqn::AbstractSolutionData,
                            opts::Dict, ctx_residual, t)


#  physicsJac(mesh, sbp, eqn, opts, pc.Ap, ctx_residual, t)

  return nothing
end


function applyLinearOperator(lo::AbstractPetscMatLO, mesh::AbstractMesh,
                             sbp::AbstractSBP, eqn::AbstractSolutionData,
                             opts::Dict, ctx_residual, t, x::AbstractVector, 
                             b::AbstractVector)

  lo2 = getBaseLO(lo)
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
  lo2 = getBaseLO(lo)
  xtmp, x_ptr = PetscVecGetArray(lo2.xtmp)
  copy!(xtmp, x)
  PetscVecRestoreArray(lo2.xtmp, x_ptr)

  PetscMatMultTranspose(lo2.A, lo2.xtmp, lo2.btmp)

  btmp, b_ptr = PetscVecGetArrayRead(lo2.btmp)
  copy!(b, btmp)
  PetscVecRestoreArrayRead(lo2.btmp, b_ptr)

  return nothing
end


function getBaseLO(lo::PetscMatLO)

  # this is the bottom of the recursion tree
  return lo
end
