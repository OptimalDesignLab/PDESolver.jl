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
function createPetscMat(mesh::AbstractMesh, sbp::AbstractOperator,
                        eqn::AbstractSolutionData, opts)

#  const mattype = PETSc2.MATMPIAIJ # should this be BAIJ?
#  const mattype = PETSc2.MATMPIBAIJ # should this be BAIJ?
  const mattype = PETSc2.MATBAIJ
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
    
    disctype = getSparsityPattern(mesh, sbp, eqn, opts)
    face_type = getFaceType(mesh.sbpface)
    _dnnz, _onnz = getBlockSparsityCounts(mesh, mesh.sbpface, disctype, face_type)
    dnnz = convert(Vector{PetscInt}, _dnnz)
    onnz = convert(Vector{PetscInt}, _onnz)
  end  # end if bs == 1


  # preallocate A
  @mpi_master println(BSTDOUT, "preallocating A")

#  if mattype == PETSc2.MATMPIAIJ
#    MatMPIAIJSetPreallocation(A, PetscInt(0),  dnnz, PetscInt(0), onnz)
#  else
    MatXAIJSetPreallocation(A, bs, dnnz, onnz, PetscIntNullArray, PetscIntNullArray)
#  end
  MatZeroEntries(A)

  # Petsc objects if this comes before preallocation
  MatSetOption(A, PETSc2.MAT_ROW_ORIENTED, PETSC_FALSE)

  if opts["addShockCapturing"]
    writeSparsityPattern(A, mesh, sbp, eqn, opts, disctype)
  end

  return A
end

"""
  Create a Petsc shell matrix

  **Inputs**
"""
function createPetscMatShell(mesh::AbstractMesh, sbp::AbstractOperator,
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
function createPetscVec(mesh::AbstractMesh, sbp::AbstractOperator,
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
function createPetscPC(mesh::AbstractMesh, sbp::AbstractOperator,
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

"""
  This function writes all zeros to the matrix.  This should only be used
  immediately after the matrix is created, to write the row and column
  indices to the matrix.  Otherwise use `MatZeroEntries`

  For most cases this is unnecessary, but for shock capturing it is required,
  because the physics module will only write the diffusion term entries to
  the elements of the matrix that have a shock in them.  On the first assemply,
  Petsc will remove any non-written to entries from the matrix structure.
  Thus if the shock moves position after the first Jacobian assembly, there
  won't be enough space in the matrix.  This function solves the problem by
  writing to all the entries (although this result in the Jacobian using more
  memory than is strictly necessary.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * A: the Petsc matrix
   * disc_type: the sparsity pattern enum
"""
function writeSparsityPattern(mesh, sbp, eqn, opts, A::PetscMat, disc_type)

  assem = _AssembleElement(A, mesh, sbp, eqn, opts)

  jac = zeros(mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerElement,
              mesh.numNodesPerElement)

  # volume terms
  for i=1:mesh.numEl
    assembleElement(assem, mesh, i, jac)
  end

  # interface terms
  for i=1:mesh.numInterfaces
    iface_i = mesh.interfaces[i]
    if disc_type == INVISCID
      assembleInterface(assem, mesh.sbpface, mesh, iface_i, jac, jac, jac, jac)
    elseif disc_type == VISCOUSTIGHT
      assembleIinterfaceVisc(assem, mesh.sbpface, mesh, iface_i, jac, jac, jac, jac)
    elseif disc_type == COLORING
      assembleInterfaceFull(assem, mesh, iface_i, jac, jac, jac, jac)
    end  # else disc_type == VISCOUS: don't currently have assembly functions
         # for this, and I'm not sure the sparsity pattern is correct either.
  end

  # shared face terms
  for peer=1:mesh.npeers
    ifaces = mesh.shared_interfaces[peer]
    for i=1:length(ifaces)
      iface_i = ifaces[i]

      if disc_type == INVISCID
        assembleSharedFAce(assem, mesh.sbpface, mesh, iface_i, jac, jac)
      elseif disc_type == VISCOUSTIGHT
        assembleSharedFaceVisc(assem, mesh.sbpface, mesh, iface_i, jac, jac)
      elseif disc_type == COLORING
        assembleSharedInterfaceFull(assem, mesh, iface_i, jac, jac)
      end
    end
  end

  return nothing
end
