# helper functions

"""
  Creates and preallocates an explicit Petsc matrix to be used as the Jacobian.
  Matrices created by this routine require the arguments to MatSetValues to
  be column-major.

  Petsc must already be initialized and have its options set.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts:
                 
"""
function createPetscMat(mesh::AbstractMesh, sbp::AbstractSBP,
                        eqn::AbstractSolutionData, opts)

#  const mattype = PETSc2.MATMPIAIJ # should this be BAIJ?
  const mattype = PETSc2.MATMPIBAIJ # should this be BAIJ?
  numDofPerNode = mesh.numDofPerNode

  comm = eqn.comm
  myrank = mesh.myrank
  obj_size = PetscInt(mesh.numDof)  # length of vectors, side length of matrices

  # create the matrix
  A = PetscMat(comm)
  SetFromOptions(A)
  MatSetType(A, mattype)
  MatSetSizes(A, obj_size, obj_size, PETSC_DECIDE, PETSC_DECIDE)
  if mesh.isDG
    MatSetOption(A, PETSc2.MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_TRUE)
  end

  if mattype == PETSc2.MATMPIAIJ
    bs = 1
  else
    bs = PetscInt(mesh.numDofPerNode)  # block size
  end

  @assert mesh.numDof % bs == 0
  nblocks = div(mesh.numDof, bs)
  dnnz = zeros(PetscInt, nblocks)  # diagonal non zeros per row
  onnz = zeros(PetscInt, nblocks)

  # calculate number of non zeros per row for A
  if bs == 1
    for i=1:mesh.numDof
      dnnz[i] = mesh.sparsity_counts[1, i]
      onnz[i] = mesh.sparsity_counts[2, i]
    end
  else
    
#    if opts["calc_jac_explicit"]
      disctype = INVISCID  # TODO: update this when merging with viscous
      face_type = getFaceType(mesh.sbpface)
      _dnnz, _onnz = getBlockSparsityCounts(mesh, mesh.sbpface, disctype, face_type)
      dnnz2 = convert(Vector{PetscInt}, _dnnz)
      onnz2 = convert(Vector{PetscInt}, _onnz)
    
#    else  # use coloring
      
      for i=1:mesh.numDof
        # this writes the information bs times to each entry
        block_i = div(i - 1, bs) + 1
        @assert mesh.sparsity_counts[1, i] % bs == 0
        @assert mesh.sparsity_counts[2, i] % bs == 0
        dnnz[block_i] = div(mesh.sparsity_counts[1, i], bs)
        onnz[block_i] = div(mesh.sparsity_counts[2, i], bs)
      end
#    end  # end if calc_jac_explicit
   
    
    println("dnnz comparison = ", hcat(dnnz, dnnz2))
    
  end  # end if bs == 1


  # preallocate A
  @mpi_master println(BSTDOUT, "preallocating A")

  if mattype == PETSc2.MATMPIAIJ
    MatMPIAIJSetPreallocation(A, PetscInt(0),  dnnz, PetscInt(0), onnz)
  else
    if opts["calc_jac_explicit"]
      MatXAIJSetPreallocation(A, bs, dnnz2, onnz2, PetscIntNullArray, PetscIntNullArray)
    else
      MatXAIJSetPreallocation(A, bs, dnnz, onnz, PetscIntNullArray, PetscIntNullArray)
    end
  end
  MatZeroEntries(A)
  matinfo = MatGetInfo(A, PETSc2.MAT_LOCAL)
  @mpi_master println(BSTDOUT, "A block size = ", matinfo.block_size)
  @mpi_master println(BSTDOUT, "A info = ", matinfo)

  # Petsc objects if this comes before preallocation
  MatSetOption(A, PETSc2.MAT_ROW_ORIENTED, PETSC_FALSE)


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
  SetFromOptions(A)  # necessary?
  SetUp(A)

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

  @assert PetscInitialized()

  comm = eqn.comm
  const vectype = VECMPI
  obj_size = PetscInt(mesh.numDof)  # length of vectors, side length of matrices
  b = PetscVec(comm)
  VecSetType(b, vectype)
  VecSetSizes(b, obj_size, PETSC_DECIDE)

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

  @assert PetscInitialized()
  comm = eqn.comm

  pc = PC(comm)
  SetFromOptions(pc)

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
  SetFromOptions(ksp)
  KSPSetPC(ksp, pc.pc)

  if typeof(pc) <: PetscMatFreePC
    Ap = lo.A  # Ap is never used, so set it to whatever
  else  # matrix-explicit PC
    Ap = pc.A
  end

  SetOperators(ksp, lo.A, Ap)

  return ksp
end

"""
  Get the integer representing the sbpface type
"""
function getFaceType(sbpface::AbstractFace)
  if typeof(sbpface) <: DenseFace
    face_type = 1
  elseif typeof(sbpface) <: SparseFace
    face_type = 2
  else
    error("unrecogized sbpface type $(typeof(sbpface))")
  end

  return face_type
end

