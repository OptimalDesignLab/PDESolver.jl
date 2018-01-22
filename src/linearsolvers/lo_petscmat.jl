# Linear operator implementation for Petsc matrix-explicit

"""
  Petsc matrix-explicit linear operator.

  Note that if an interface common to both Petsc and regular matrices is used
  when accessing the `Ap` field, it is possible to write functions that
  operate on both regular and Petsc linear operators.

  **Public Fields**

   * A: a PetscMat

"""
type PetscMatLO <: AbstractPetscMatLO
  A::PetscMat
  xtmp::PetscVec  # these are shared with the PC if possible
  btmp::PetscVec
  is_setup::Array{Bool, 1}  # is matrix assembled (shared with PC.is_assembled)
                            # if A is shared
  is_shared::Bool
  is_finalized::Bool
  nassemblies::Array{Int, 1}
  nsolves::Int
  ntsolves::Int

  # MPI stuff
  comm::MPI.Comm
  myrank::Int
  commsize::Int
end

"""
  Outer constructor for PetscMatLO

  **Inputs**

   * pc: a Petsc preconditioner of some kind
   * mesh
   * sbp
   * eqn
   * opts
"""
function PetscMatLO(pc::AbstractPetscPC, mesh::AbstractMesh,
                    sbp::AbstractSBP, eqn::AbstractSolutionData, opts::Dict)

  pc2 = getBasePC(pc)

  # share matrix if possible
  if !opts["use_jac_precond"] && !(typeof(pc) <: PetscMatFreeLO) 
    A = pc2.A
    is_setup = pc2.is_assembled  # share the indicator array
    is_shared = true
    nassemblies = pc2.nassemblies
  else
    A = createPetscMat(mesh, sbp, eqn, opts)
    is_setup = Bool[false]
    is_shared = true
    nassemblies = Int[0]
  end

  if typeof(pc2) <: PetscMatPC  # if pc2 has Petsc vectors inside it
    xtmp = pc2.xtmp
    btmp = pc2.btmp
  else
    xtmp = createPetscVec(mesh, sbp, eqn, opts)
    btmp = createPetscVec(mesh, sbp, eqn, opts)
  end

  is_finalized = false
  nsolves = 0
  ntsolves = 0
  comm = eqn.comm
  myrank = eqn.myrank
  commsize = eqn.commsize


  return PetscMatLO(A, xtmp, btmp, is_setup, is_shared, is_finalized,
                    nassemblies, nsolves, ntsolves, comm, myrank, commsize)
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


#  physicsJac(mesh, sbp, eqn, opts, pc.A, ctx_residual, t)
  setIsSetup(lo, false)

  return nothing
end


function applyLinearOperator(lo::AbstractPetscMatLO, mesh::AbstractMesh,
                             sbp::AbstractSBP, eqn::AbstractSolutionData,
                             opts::Dict, ctx_residual, t, x::AbstractVector, 
                             b::AbstractVector)

  lo2 = getBaseLO(lo)

  # assemble Ap (if needed)
  assemblePetscData(lo2, x, lo2.xtmp)

  MatMult(lo2.A, lo2.xtmp, lo2.btmp)

  btmp = VecGetArrayRead(lo2.btmp)
  copy!(b, btmp)
  VecRestoreArrayRead(lo2.btmp, btmp)

  return nothing
end


function applyLinearOperatorTranspose(lo::AbstractPetscMatLO, 
                             mesh::AbstractMesh, sbp::AbstractSBP,
                             eqn::AbstractSolutionData, opts::Dict, 
                             ctx_residual, t, x::AbstractVector, 
                             b::AbstractVector)
  lo2 = getBaseLO(lo)

  # assemble Ap (if needed)
  assemblePetscData(lo2, x, lo.xtmp)


  MatMultTranspose(lo2.A, lo2.xtmp, lo2.btmp)

  btmp = VecGetArrayRead(lo2.btmp)
  copy!(b, btmp)
  VecRestoreArrayRead(lo2.btmp, btmp)

  return nothing
end


function assemblePetscData(lo::PetscMatLO, b::AbstractVector,
                           b_petsc::PetscVec)

  myrank = lo.myrank

  # assemble things
  if !getIsSetup(lo)
    MatAssemblyBegin(lo.A, MAT_FINAL_ASSEMBLY)
  end

  # copy values into the vector
  btmp = VecGetArray(b_petsc)
  copy!(btmp, b)
  VecRestoreArray(b_petsc, btmp)

  if !getIsSetup(lo)
    MatAssemblyEnd(lo.A, MAT_FINAL_ASSEMBLY)
    setIsSetup(lo, true)
    lo.nassemblies[1] += 1
    matinfo = MatGetInfo(lo.A, PETSc2.MAT_LOCAL)
    if matinfo.mallocs > 0.5  # if any mallocs
      println(BSTDERR, "Warning: non-zero number of mallocs for A on process $myrank: $(matinfo.mallocs) mallocs")
    end
  end

  return nothing
end

# this function works for all linear operators
function getIsSetup(lo::AbstractLO)

  lo2 = getBaseLO(lo)
  return lo2.is_setup[1]
end

function setIsSetup(lo::AbstractLO, val::Bool)
  
  lo2 = getBaseLO(lo)
  lo2.is_setup[1] = val

  return nothing
end



function getBaseLO(lo::PetscMatLO)

  # this is the bottom of the recursion tree
  return lo
end
