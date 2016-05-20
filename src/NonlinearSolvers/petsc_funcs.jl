# petc_funcs.jl: function used to create, solve, and destroy petsc object

@doc """
### NonlinearSolver.createPetscData

  This function creates all the objects needed to use the Petsc solver.

  Inputs:
    mesh: AbstractMesh object
    pmesh: AbstractMesh object used for preconditioning
    sbp: SBP operator
    eqn: AbstractSolutionData object
    opts: options dictonary
    newton_data:  newton iteration object
    func: residual evaluation function

  Outputs:
    A: PetscMat for the Jacobian
    Ap: PetscMat for the preconditioning Jacobian
    x: PetscVec for the solution of Ax=b
    b: Petscvec for the b in Ax=b
    ksp: KSP context
    ctx: the context for the KSP context.  Is a tuple of all the objects other
         than A, x, and b needed to do jacobian vector products

    Aliasing restrictions: none

"""->
function createPetscData(mesh::AbstractMesh, pmesh::AbstractMesh, sbp, eqn::AbstractSolutionData, opts, newton_data::NewtonData, func)
# initialize petsc and create Jacobian matrix A, and vectors x, b, and the
# ksp context
# serial only
# func residual evaluation function

jac_type = opts["jac_type"]::Int
myrank = newton_data.myrank

const vectype = VECMPI
const mattype = PETSc.MATMPIAIJ # should this be BAIJ?
# initialize Petsc
#PetscInitialize(["-malloc", "-malloc_debug", "-malloc_dump", "-ksp_monitor", "-pc_type", "ilu", "-pc_factor_levels", "4" ])

#PetscInitialize(["-malloc", "-malloc_debug", "-malloc_dump", "-sub_pc_factor_levels", "4", "ksp_gmres_modifiedgramschmidt", "-ksp_pc_side", "right", "-ksp_gmres_restart", "1000" ])
numDofPerNode = mesh.numDofPerNode
PetscInitialize(["-malloc", "-malloc_debug", "-ksp_monitor", "kspout",  "-pc_type", "bjacobi", "-sub_pc_type", "ilu", "-sub_pc_factor_levels", "4", "ksp_gmres_modifiedgramschmidt", "-ksp_pc_side", "right", "-ksp_gmres_restart", "30", "-log_summary" ])
comm = MPI.COMM_WORLD

obj_size = PetscInt(mesh.numDof)  # length of vectors, side length of matrices
@mpi_master println("creating b")
b = PetscVec(comm)
PetscVecSetType(b, vectype)
PetscVecSetSizes(b, obj_size, PETSC_DECIDE)

@mpi_master println("creating x")
x = PetscVec(comm)
PetscVecSetType(x, vectype)
PetscVecSetSizes(x, obj_size, PETSC_DECIDE)

# tuple of all the objects needed to evaluate the residual
# only used by Mat Shell
ctx = (mesh, sbp, eqn, opts, newton_data, func)  # tuple of all the objects needed to evalute the
                                    # residual
if jac_type == 3  # explicit sparse jacobian
  @mpi_master println("creating A")
  A = PetscMat(comm)
  PetscMatSetFromOptions(A)
  PetscMatSetType(A, mattype)
  PetscMatSetSizes(A, obj_size, obj_size, PETSC_DECIDE, PETSC_DECIDE)

elseif jac_type == 4  # jacobian-vector product
  # create matrix shell
  ctx_ptr = pointer_from_objref(ctx)  # make a pointer from the tuple
  A = MatCreateShell(comm, obj_size, obj_size, obj_size, PETSC_DECIDE, PETSC_DECIDE)
  PetscMatSetFromOptions(A)  # necessary?
  PetscSetUp(A)

  # give PETSc the function pointer of Jacobian vector product function
  fptr = cfunction(calcJacVecProd_wrapper, PetscErrorCode, (PetscMat, PetscVec, PetscVec))
  MatShellSetOperation(A, PETSc.MATOP_MULT, fptr)

else
  println(STDERR, "Unsupported jacobian type requested")
  println(STDERR, "jac_type = ", jac_type)
end


if jac_type == 4 || opts["use_jac_precond"]
  @mpi_master println("creating Ap")  # used for preconditioner
  Ap = PetscMat(comm)
  PetscMatSetFromOptions(Ap)
  PetscMatSetType(Ap, mattype)
  PetscMatSetSizes(Ap, PetscInt(mesh.numDof), PetscInt(mesh.numDof), PetscInt(mesh.numDof), PetscInt(mesh.numDof))

  @mpi_master println("type of Ap = ", MatGetType(Ap))
else
  Ap = A
end

@mpi_master println("preallocating Petsc Matrices")
# prellocate matrix
dnnz = zeros(PetscInt, mesh.numDof)  # diagonal non zeros per row
onnz = zeros(PetscInt, mesh.numDof)
#dnnzu = zeros(PetscInt, 1)  # only needed for symmetric matrices
#onnzu = zeros(PetscInt, 1)  # only needed for symmetric matrices
bs = PetscInt(mesh.numDofPerNode)  # block size

# calculate number of non zeros per row for A
for i=1:mesh.numDof
#  max_dof = mesh.sparsity_nodebnds[2, i]
#  min_dof = mesh.sparsity_nodebnds[1, i]
#  nnz_i = max_dof - min_dof + 1
  dnnz[i] = mesh.sparsity_counts[1, i]
  onnz[i] = mesh.sparsity_counts[2, i]
#  println("row ", i," has ", nnz_i, " non zero entries")
end
myrank = mesh.myrank


# preallocate A
if jac_type == 3

  PetscMatMPIAIJSetPreallocation(A, PetscInt(0),  dnnz, PetscInt(0), onnz)

  MatSetOption(A, PETSc.MAT_ROW_ORIENTED, PETSC_FALSE)
  PetscMatZeroEntries(A)
  matinfo = PetscMatGetInfo(A, Int32(1))
  @mpi_master println("A block size = ", matinfo.block_size)

#  PetscMatAssemblyBegin(A, PETSC_MAT_FLUSH_ASSEMBLY)
#  PetscMatAssemblyEnd(A, PETSC_MAT_FLUSH_ASSEMBLY)
end

# preallocate Ap
if jac_type == 4 || opts["use_jac_precond"]
  # calculate number of nonzeros per row for A[
  for i=1:mesh.numNodes
#    max_dof = pmesh.sparsity_nodebnds[2, i]
#    min_dof = pmesh.sparsity_nodebnds[1, i]
#    nnz_i = max_dof - min_dof + 1
#    dnnz[i] = nnz_i
    dnnz[i] = mesh.sparsity_counts_node[1, i]
    onnz[i] = mesh.sparsity_counts_node[2, i]

  #  println("row ", i," has ", nnz_i, " non zero entries")
  end

  PetscMatXAIJSetPreallocation(Ap, bs, dnnz, onnz, dnnzu, onnzu)

  MatSetOption(Ap, PETSc.MAT_ROW_ORIENTED, PETSC_FALSE)
  matinfo = PetscMatGetInfo(Ap, Int32(1))
  @mpi_master println("Ap block size = ", matinfo.block_size)

  # zero initialize the matrix just in case
  PetscMatZeroEntries(Ap)

  PetscMatAssemblyBegin(Ap, PETSC_MAT_FLUSH_ASSEMBLY)
  PetscMatAssemblyEnd(Ap, PETSC_MAT_FLUSH_ASSEMBLY)
end


# set some options
# MatSetValuesBlocked will interpret the array of values as being column
# major

# create KSP contex
ksp = KSP(comm)

KSPSetFromOptions(ksp)
KSPSetOperators(ksp, A, Ap)  # this was A, Ap

# set: rtol, abstol, dtol, maxits
#KSPSetTolerances(ksp, 1e-2, 1e-12, 1e5, PetscInt(1000))
#KSPSetUp(ksp)


pc = KSPGetPC(ksp)
pc_type = PCGetType(pc)
@mpi_master println("pc_type = ", pc_type)
#=
if pc_type == "bjacobi"
  n_local, first_local, ksp_arr = PCBJacobiGetSubKSP(pc)
  println("n_local = ", n_local, ", first_local = ", first_local)
  println("length(ksp_arr) = ", length(ksp_arr))

  sub_ksp = ksp_arr[first_local + 1]
  sub_pc = KSPGetPC(sub_ksp)
  pc_subtype = PCGetType(sub_pc)
  println("pc_subtype = ", pc_subtype)

  fill_level = PCFactorGetLevels(sub_pc)
  println("preconditioner using fill level = ", fill_level)
end

if pc_type == "ilu"
  fill_level = PCFactorGetLevels(sub_pc)
  println("preconditioner using fill level = ", fill_level)
end
=#


return A, Ap, x, b, ksp, ctx

end

@doc """
### NonlinearSolvers.destroyPetsc

  Destroy all the objects when finished with Petsc

  Inputs: see createPetscData

"""->
function destroyPetsc(A::PetscMat, Ap::PetscMat, x::PetscVec,  b::PetscVec, ksp::KSP)
# destory Petsc data structures, finalize Petsc

PetscDestroy(A)
if A.pobj != Ap.pobj
  PetscDestroy(Ap)
end
PetscDestroy(x)
PetscDestroy(b)
PetscDestroy(ksp)
PetscFinalize()

return nothing

end

@doc """
### NonlinearSolver.petscSolve

  This function performs the petsc solve x = inv(A)*b. More specifically, it 
  copies he right hand side into the PetscVec b, assembles them, and performs
  the solve.

  Inputs:
    newton_data: Newton's method object
    A:  PetscMat that is the Jacobian
    Ap: PetscMat that is the preconditioner
    x: PetscVec to store the solution in
    b: PetscVec to put the right hand side in
    ksp:  KSP context
    opts: options dictonary
    res_0: vector to use for the right hand side

  Inputs/Outputs:
    delta_res_vec: julia vector ot return the solution x in

  Aliasing restrictions: none

"""->
function petscSolve(newton_data::NewtonData, A::PetscMat, Ap::PetscMat, x::PetscVec, b::PetscVec, ksp::KSP, opts, res_0::AbstractVector, delta_res_vec::AbstractVector, dof_offset::Integer=0 )

  # solve the system for the newton step, write it to delta_res_vec
  # writing it to delta_res_vec is an unecessary copy, becasue we could
  # write it directly to eqn.q, but for consistency we do it anyways

#  return nothing
  jac_type = opts["jac_type"]::Int
  myrank = newton_data.myrank
  # copy res_0 into b
  # create the index array
  numDof = length(delta_res_vec)
#  idx = zeros(PetscInt, numDof)
#  for i=1:numDof
#    idx[i] = i - 1 + dof_offset
#  end

  b_petsc, b_ptr = PetscVecGetArray(b)
  for i=1:numDof
    b_petsc[i] = res_0[i]
  end
  PetscVecRestoreArray(b, b_ptr)

  # copy into Petsc and assemble
  # should do this by direct array access
#  PetscVecSetValues(b, idx, res_0, PETSC_INSERT_VALUES)
#  PetscVecAssemblyBegin(b)
#  PetscVecAssemblyEnd(b)

  #=
  PetscVecSetValues(x, idx, res_0, PETSC_INSERT_VALUES)
  PetscVecAssemblyBegin(x)
  PetscVecAssemblyEnd(x)
  =#


  # assemble matrix
  PetscMatAssemblyBegin(Ap, PETSC_MAT_FINAL_ASSEMBLY)
  PetscMatAssemblyEnd(Ap, PETSC_MAT_FINAL_ASSEMBLY)


  # should this be overlapped with the vector assembly?
  if jac_type != 4  && !(A === Ap)
    PetscMatAssemblyBegin(A, PETSC_MAT_FINAL_ASSEMBLY)
    PetscMatAssemblyEnd(A, PETSC_MAT_FINAL_ASSEMBLY)
  end

  if jac_type != 4
    matinfo = PetscMatGetInfo(A, MAT_LOCAL)
    if matinfo.mallocs > 0.5  # if any mallocs
      println("Caution: non zero number of mallocs for A")
      println("  number of mallocs = ", matinfo.mallocs)
    end
  end

    matinfo = PetscMatGetInfo(Ap, MAT_LOCAL)

    if matinfo.mallocs > 0.5  # if any mallocs
      println("Caution: non zero number of mallocs for Ap")
      println("  number of mallocs = ", matinfo.mallocs)
    end
 
  return nothing
  # only call this first time?
  # what happens when A and Ap change?
  # this is not truely necessary for the common case, because KSPSolve
  # calls it if needed
  # it is necessary to call KSPSetUp before getting the preconditioner
  # context in some cases

  KSPSetTolerances(ksp, newton_data.reltol, newton_data.abstol, 
                   newton_data.dtol, PetscInt(newton_data.itermax))

  KSPSetUp(ksp)

  nx = PetscVecGetSize(x)
  nb = PetscVecGetSize(b)
  KSPSolve(ksp, b, x)
  reason = KSPGetConvergedReason(ksp)
  @mpi_master println("KSP converged reason = ", KSPConvergedReasonDict[reason])
  rnorm = KSPGetResidualNorm(ksp)
  @mpi_master println("Linear residual = ", rnorm)

  # copy solution from x into delta_res_vec
  arr, ptr_arr = PetscVecGetArray(x)
  for i=1:numDof
    delta_res_vec[i] = arr[i]
  end

  PetscVecRestoreArray(x, ptr_arr)

  # zero out the Jacobian for next use
  if opts["jac_type"] != 4
    PetscMatZeroEntries(A)
  end
  PetscMatZeroEntries(Ap)
  return nothing
end


