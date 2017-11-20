# Linear operator implementation for Petsc matrix-explicit

"""
  Petsc matrix-explicit linear operator.
"""
type PetscMatLO <: AbstractPetscMatLO
  A::PetscMat
  xtmp::PetscVec  # these are shared with the PC if possible
  btmp::PetscVec
  is_setup::Array{Bool, 1}  # is matrix assembled (shared with PC.is_assembled)
                            # if A is shared
  is_finalized::Bool

  # MPI stuff
  comm::MPI.Comm
  myrank::Int
  commsize::Int
end

function PetscMatLO(pc::AbstractPetscMatPC, mesh::AbstractMesh,
                    sbp::AbstractSBP, eqn::AbstractSolutionData, opts::Dict)

  pc2 = getBasePC(pc)
  if !opts["use_jac_precond"] && !(typeof(pc) <: PetscMatFreeLO) # share matrix if possible
    A = pc2.Ap
    is_setup = pc2.is_assembled  # share the indicator array
  else
    A = createPetscMat(mesh, sbp, eqn, opts)
    is_setup = Bool[false]
  end

  if typeof(pc2) <: PetscMatPC  # if pc2 has Petsc vectors inside it
    xtmp = pc2.xtmp
    btmp = pc2.btmp
  else
    xtmp = createPetscVec(mesh, sbp, eqn, opts)
    btmp = createPetscVec(mesh, sbp, eqn, opts)
  end

  is_finalized = false
  comm = eqn.comm
  myrank = eqn.myrank
  commsize = eqn.commsize


  return PetscMatLO(A, xtmp, btmp, is_setup, is_finalized, comm, myrank,
                    commsize)
end


function free(lo::PetscMatLO)

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


function calcLinearOperator(lo::PetscMatLO, mesh::AbstractMesh,
                            sbp::AbstractSBP, eqn::AbstractSolutionData,
                            opts::Dict, ctx_residual, t)


#  physicsJac(mesh, sbp, eqn, opts, pc.Ap, ctx_residual, t)
  setIsSetup(lo, false)

  return nothing
end


function applyLinearOperator(lo::AbstractPetscMatLO, mesh::AbstractMesh,
                             sbp::AbstractSBP, eqn::AbstractSolutionData,
                             opts::Dict, ctx_residual, t, x::AbstractVector, 
                             b::AbstractVector)

  lo2 = getBaseLO(lo)
  if !getIsSetup(lo2)
    assemblePetscData(lo, x, lo.xtmp)
  end

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

  if !getIsSetup(lo2)
    assemblePetscData(lo, x, lo.xtmp)
  end

  PetscMatMultTranspose(lo2.A, lo2.xtmp, lo2.btmp)

  btmp, b_ptr = PetscVecGetArrayRead(lo2.btmp)
  copy!(b, btmp)
  PetscVecRestoreArrayRead(lo2.btmp, b_ptr)

  return nothing
end


function assemblePetscData(lo::PetscMatLO, b::AbstractVector,
                           b_petsc::PetscVec)

  myrank = lo.myrank

  # assemble things
  if !getIsSetup(lo)
    PetscMatAssemblyBegin(lo.A, PETSC_MAT_FINAL_ASSEMBLY)
  end

  # copy values into the vector
  btmp, b_ptr = PetscVecGetArray(b_petsc)
  copy!(btmp, b)
  PetscVecRestoreArray(b_petsc, b_ptr)

  if !getIsSetup(lo)
    PetscMatAssemblyEnd(lo.A, PETSC_MAT_FINAL_ASSEMBLY)
    setIsSetup(lo, true)
    matinfo = PetscMatGetInfo(lo.A, PETSc.MAT_LOCAL)
    if matinfo.mallocs > 0.5  # if any mallocs
      println(BSTDERR, "Warning: non-zero number of mallocs for A on process $myrank: $(matinfo.mallocs) mallocs")
    end
  end

  return nothing
end

# this function works for all linear operators
function getIsSetup(lo::AbstractLinearOperator)

  lo2 = getBaseLO(lo)
  return lo2.is_setup[1]
end

function setIsSetup(lo::AbstractLinearOperator, val::Bool)
  
  lo2 = getBaseLO(lo)
  lo2.is_setup[1] = val

  return nothing
end



function getBaseLO(lo::PetscMatLO)

  # this is the bottom of the recursion tree
  return lo
end
