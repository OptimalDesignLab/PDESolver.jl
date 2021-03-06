# AbstractPC type of Petsc matrix-explicit

"""
  [`AbstractPC`](@ref) implementation for Petsc matrix-explicit preconditioners.
 
  **Public Fields**

   * A: a PetscMat object used to calculate the preconditioner
"""
mutable struct PetscMatPC <: AbstractPetscMatPC
  pc::PC  # Petsc PC object
  A::PetscMat  # Petsc Mat object
  xtmp::PetscVec  # reusable temporary vector
  btmp::PetscVec
  is_assembled::Array{Bool, 1}  # is A assembled
  is_setup::Bool  # is PC already set up
  is_shared::Bool
  is_finalized::Bool
  nassemblies::Array{Int, 1}  # number of matrix assemblies, shared with lo
                              # if A is shared
  nsetups::Int  # number of PC setups
  napplies::Int  # number of PC applications
  ntapplies::Int  # number of transpose applies

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
function PetscMatPC(mesh::AbstractMesh, sbp::AbstractOperator,
                    eqn::AbstractSolutionData, opts::Dict)

  pc = createPetscPC(mesh, sbp, eqn, opts)
  A = createPetscMat(mesh, sbp, eqn, opts)
  xtmp = createPetscVec(mesh, sbp, eqn, opts)
  btmp = createPetscVec(mesh, sbp, eqn, opts)
  is_assembled = Bool[false]
  is_setup = false

  # we can't reliably detect whether the LO (which is constructed *after* the
  # PC) will share the jacobian or not.  We make a best guess here and correct
  # the is_shared flag in the LO constructor if needed
  if opts["use_jac_precond"] || opts["jac_type"] == 4  # this must match the condition in LOPetscMat
    is_shared = false
  else
    is_shared = true
  end

  is_finalized = false
  nassemblies = Int[0]
  nsetups = 0
  napplies = 0
  ntapplies = 0

  comm = eqn.comm
  myrank = eqn.myrank
  commsize = eqn.commsize

  return PetscMatPC(pc, A, xtmp, btmp, is_assembled, is_setup, is_shared,
                    is_finalized, nassemblies, nsetups, napplies, ntapplies,
                    comm, myrank, commsize)
end

function free(pc::PetscMatPC)

  if !pc.is_finalized
    if pc.pc.pobj != C_NULL
      PetscDestroy(pc.pc)
      pc.pc.pobj = C_NULL
    end

    if pc.A.pobj != C_NULL
      PetscDestroy(pc.A)
      pc.A.pobj = C_NULL
    end

    if pc.xtmp.pobj != C_NULL
      PetscDestroy(pc.xtmp)
      pc.xtmp.pobj = C_NULL
    end

    if pc.btmp.pobj != C_NULL
      PetscDestroy(pc.btmp)
      pc.btmp.pobj = C_NULL
    end
  end

  pc.is_finalized = true

  return nothing
end




function calcPC(pc::PetscMatPC, mesh::AbstractMesh, sbp::AbstractOperator,
                eqn::AbstractSolutionData, opts::Dict, ctx_residual, t)

  # compute the jacobian here
#  physicsJac(mesh, sbp, eqn, opts, pc.A, ctx_residual, t)
 
  # don't setup the PC here because PC built on top of this one might
  # modify A after calling this function
  setIsAssembled(pc, false)
  pc.is_setup = false

  return nothing
end


function applyPC(pc::AbstractPetscMatPC, mesh::AbstractMesh, sbp::AbstractOperator,
                 eqn::AbstractSolutionData, opts::Dict, t, b::AbstractVector, 
                 x::AbstractVector)

  pc2 = getBasePC(pc)

  # assemble matrix (if needed)
  assemblePetscData(pc2, b, pc2.btmp)

  if !pc2.is_setup
    setupPC(pc2)
  end

  # call Petsc PCApply
  PCApply(pc2.pc, pc2.btmp, pc2.xtmp)

  # copy back to x
  xtmp = VecGetArrayRead(pc2.xtmp)
  copy!(x, xtmp)
  VecRestoreArrayRead(pc2.xtmp, xtmp)

  pc2.napplies += 1

  return nothing
end


function applyPCTranspose(pc::AbstractPetscMatPC, mesh::AbstractMesh,
                 sbp::AbstractOperator,
                 eqn::AbstractSolutionData, opts::Dict, t,
                 b::AbstractVector, x::AbstractVector)


  pc2 = getBasePC(pc)
  @assert typeof(pc2) <: PetscMatPC

  if !PCApplyTransposeExists(pc2.pc)
    ptype = PCGetType(pc2.pc)
    throw(ErrorException("PCApplyTranspose not defined for PC $ptype"))
  end

  # assemble matrix (if needed)
  assemblePetscData(pc2, b, pc2.btmp)

  if !pc2.is_setup
    setupPC(pc2)
  end

  # call Petsc PCApplyTranspose
  PCApplyTranspose(pc2.pc, pc2.btmp, pc2.xtmp)

  # copy back to x
  xtmp = VecGetArrayRead(pc2.xtmp)
  copy!(x, xtmp)
  VecRestoreArrayRead(pc2.xtmp, xtmp)

  pc2.ntapplies += 1
  return nothing
end


"""
  This internal function is used to setup the PC, including setting the flag.
  This function doesn't actually compute the preconditioner, but it tells Petsc
  to do it the next time it is needed

  **Inputs**

   * pc: PetscMatPC object
"""
function setupPC(pc::PetscMatPC)

  PCSetReusePreconditioner(pc.pc, PETSC_FALSE)
#  tsetup = @elapsed PCSetUp(pc.pc)
#  tsetup = @elapsed PCSetUp(pc.pc)
#  println(BSTDOUT, "PCSetup time = ", tsetup)
  # this is potentially bad because Petsc will *never* recompute the 
  # preconditioner on its own.  Using higher level functionality like TS
  # or Petsc nonlinear solvers likely won't work in this case
#  PCSetReusePreconditioner(pc.pc, PETSC_TRUE)
  pc.nsetups += 1
  pc.is_setup = true

  return nothing
end

function assemblePetscData(pc::PetscMatPC, b::AbstractVector, b_petsc::PetscVec)

  myrank = pc.myrank

  if !getIsAssembled(pc)
    MatAssemblyBegin(pc.A, MAT_FINAL_ASSEMBLY)
  end

  # copy values into the vector
  btmp = VecGetArray(b_petsc)
  copy!(btmp, b)
  VecRestoreArray(b_petsc, btmp)

  if !getIsAssembled(pc)
    MatAssemblyEnd(pc.A, MAT_FINAL_ASSEMBLY)
    setIsAssembled(pc, true)
    pc.nassemblies[1] += 1
    matinfo = MatGetInfo(pc.A, PETSc2.MAT_LOCAL)
    if matinfo.mallocs > 0.5  # if any mallocs
      println(BSTDERR, "Warning: non-zero number of mallocs for A on process $myrank: $(matinfo.mallocs) mallocs")
    end
  end

  return nothing
end


"""
  Internal function for recording whether pc.A is assembled or not
"""
function setIsAssembled(pc::PetscMatPC, val::Bool)

  pc.is_assembled[1] = val

  return nothing
end

"""
  Internal function for retrieving whether pc.A is assembled or not.
"""
function getIsAssembled(pc::PetscMatPC)

  return pc.is_assembled[1]
end


function getBasePC(pc::PetscMatPC)

  # this is the bottom of the recursion tree
  return pc
end
