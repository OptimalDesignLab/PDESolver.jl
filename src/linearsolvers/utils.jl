# helper functions

"""
  Creates and preallocates an explicit Petsc matrix to be used as the Jacobian

  Petsc must already be initialized and have its options set.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts:
                 
"""
function createPetscMat(mesh::AbstractMesh, sbp::AbstractSBP,
                        eqn::AbstractSolutionData, opts)

  const mattype = PETSc.MATMPIAIJ # should this be BAIJ?
  numDofPerNode = mesh.numDofPerNode

  comm = eqn.comm
  myrank = mesh.myrank
  obj_size = PetscInt(mesh.numDof)  # length of vectors, side length of matrices

  # create the matrix
  A = PetscMat(comm)
  PetscMatSetFromOptions(A)
  PetscMatSetType(A, mattype)
  PetscMatSetSizes(A, obj_size, obj_size, PETSC_DECIDE, PETSC_DECIDE)
  if mesh.isDG
    MatSetOption(A, PETSc.MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_TRUE)
  end
  MatSetOption(A, PETSc.MAT_ROW_ORIENTED, PETSC_FALSE)

  dnnz = zeros(PetscInt, mesh.numDof)  # diagonal non zeros per row
  onnz = zeros(PetscInt, mesh.numDof)
  bs = PetscInt(mesh.numDofPerNode)  # block size

  # calculate number of non zeros per row for A
  for i=1:mesh.numDof
    dnnz[i] = mesh.sparsity_counts[1, i]
    onnz[i] = mesh.sparsity_counts[2, i]
  end

  # preallocate A
  @mpi_master println(BSTDOUT, "preallocating A")

  #TODO: set block size?
  PetscMatMPIAIJSetPreallocation(A, PetscInt(0),  dnnz, PetscInt(0), onnz)
  #PetscMatXAIJSetPreallocation(Ap, bs, dnnz, onnz, PetscIntNullArray, PetscIntNullArray)
  PetscMatZeroEntries(A)
  matinfo = PetscMatGetInfo(A, Int32(1))
  @mpi_master println(BSTDOUT, "A block size = ", matinfo.block_size)

  return A
end

"""
  Create a Petsc shell matrix

  **Inputs**
"""
function createPetscMatShell(mesh::AbstractMesh, sbp::AbstractSBP,
                             eqn::AbstractSolutionData, opts)

  obj_size = PetscInt(mesh.numDof)  # length of vectors, side length of matrices
  comm = eqn.comm

  ctx_ptr = C_NULL # this gets set later
  A = MatCreateShell(comm, obj_size, obj_size, PETSC_DECIDE, PETSC_DECIDE, ctx_ptr)
  PetscMatSetFromOptions(A)  # necessary?
  PetscSetUp(A)

  return A
end

"""
  Create a parallel Petsc vector of size mesh.numdof

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
"""
function createPetscVec(mesh::AbstractMesh, sbp::AbstractSBP,
                        eqn::AbstractSolutionData, opts)

  if PetscInitialized() == 0 # PETSc Not initialized before
    PetscInitialize()
  end
  comm = eqn.comm
  const vectype = VECMPI
  obj_size = PetscInt(mesh.numDof)  # length of vectors, side length of matrices
  b = PetscVec(comm)
  PetscVecSetType(b, vectype)
  PetscVecSetSizes(b, obj_size, PETSC_DECIDE)

  return b
end

"""
  Create the Petsc PC object.

  Petsc must already be initialized and have its options set.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts

"""
function createPetscPC(mesh::AbstractMesh, sbp::AbstractSBP,
                        eqn::AbstractSolutionData, opts)

  if PetscInitialized() == 0 # PETSc Not initialized before
    PetscInitialize()
  end
  comm = eqn.comm

  pc = PCCreate(comm)
  PCSetFromOptions(pc)

  return pc
end



"""
  Create the KSP object

  Petsc must already be initialized and have its options set.

  **Inputs**

   * pc: a fully initailized Petsc preconditioner, either [`PetscMatPC`](@ref)
         or [`PetscMatFreePC`](@ref)
   * lo: a fully initialized Petsc linear operator, either [`PetscMatLO`](@ref)
         or [`PetscMatFreeLO`](@ref)
   * comm: the MPI communicator used to defined the pc and lo

"""
function createKSP(pc::Union{PetscMatPC, PetscMatFreePC},
                   lo::Union{PetscMatLO, PetscMatFreeLO}, comm::MPI.Comm)

  ksp = KSP(comm)
  KSPSetFromOptions(ksp)
  KSPSetPC(ksp, pc.pc)

  if pc <: PetscMatFreePC
    Ap = lo.A  # Ap is never used, so set it to whatever
  else  # matrix-explicit PC
    Ap = pc.Ap
  end

  KSPSetOperators(ksp, A, Ap)

  return ksp
end
