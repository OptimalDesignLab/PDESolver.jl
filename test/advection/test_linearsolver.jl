# tests for the LinearSolver module

function test_linearsolver()


  fname = "input_vals_linear.jl"
  mesh, sbp, eqn, opts = run_advection(fname)

  facts("----- Testing Linear Solver -----") do
    test_dense(mesh, sbp, eqn, opts)
    test_sparsedirect(mesh, sbp, eqn, opts)
    test_petscmat(mesh, sbp, eqn, opts)
    test_petscmatfree(mesh, sbp, eqn, opts)
    test_petscmat_matfree(mesh, sbp, eqn, opts)
  end


end

function test_dense(mesh, sbp, eqn, opts)

  pc = PCNone(mesh, sbp, eqn, opts)
  lo = DenseLO(pc, mesh, sbp, eqn, opts)
  ls = StandardLinearSolver(pc, lo, eqn.comm)

  ctx_residual = (evalResidual,)
  t = 0.0
  calcLinearOperator(ls, mesh, sbp, eqn, opts, ctx_residual, t)
  rand!(lo.A)
  for i=1:mesh.numDof
    lo.A[i, i] = 5  # make non-singular
  end

  A_orig = copy(lo.A)
  b = rand(mesh.numDof)

  # test product
  c = lo.A*b
  c2 = zeros(c)
  applyLinearOperator(ls.lo, mesh, sbp, eqn, opts, ctx_residual, t, b, c2)
  @fact norm(c - c2) --> roughly(0.0, atol=1e-12)

  # test transpose product
  c = lo.A.'*b
  applyLinearOperatorTranspose(ls.lo, mesh, sbp, eqn, opts, ctx_residual, t, b, c2)
  @fact norm(c - c2) --> roughly(0.0, atol=1e-12)

  # test solve
  x = lo.A\b
  x2 = zeros(x)
  linearSolve(ls, b, x2)
  @fact norm(x2 - x) --> roughly(0.0, atol=1e-12)

  # test PC
  x3 = zeros(x2)
  applyPC(ls, mesh, sbp, eqn, opts, t, b, x3)
  @fact norm(x3 - x2) --> roughly(0.0, atol=1e-12)

  # test transpose solve
  x = A_orig.'\b
  x2 = zeros(x)
  linearSolveTranspose(ls, b, x2)
  @fact norm(x2 - x) --> roughly(0.0, atol=1e-12)

  # test PC
  x3 = zeros(x2)
  applyPCTranspose(ls, mesh, sbp, eqn, opts, t, b, x3)
  @fact norm(x3 - x2) --> roughly(0.0, atol=1e-12)


  
  # test product (again)
  c = A_orig*b
  c2 = zeros(c)
  @fact_throws applyLinearOperator(ls.lo, mesh, sbp, eqn, opts, ctx_residual, t, b, c2)
  #=
  @fact norm(c - c2) --> roughly(0.0, atol=1e-12)
  println("c = \n", c)
  println("c2 = \n", c2)
  println("diff = \n", c2 - c)
  println("ipiv = \n", ls.lo.ipiv)
  =#
  

  # transpose product throws
  @fact_throws applyLinearOperatorTranspose(ls.lo, mesh, sbp, eqn, opts, ctx_residual, t, b, c2)

  # test that the factorization is not performed again

  x = A_orig\b
  x3 = zeros(x2)
  linearSolve(ls, b, x3)
  @fact norm(x3 - x) --> roughly(0.0, atol=1e-12)
  @fact ls.lo.nfactorizations --> 1
  @fact ls.lo.nsolves --> 3

  free(ls)

  return nothing
end


function test_sparsedirect(mesh, sbp, eqn, opts)

  # test sparse direct solve
  println("testing sparse direct linear solver")
  pc = PCNone(mesh, sbp, eqn, opts)
  lo = SparseDirectLO(pc, mesh, sbp, eqn, opts)
  ls = StandardLinearSolver(pc, lo, eqn.comm)

  ctx_residual = (evalResidual,)
  t = 0.0
  calcLinearOperator(ls, mesh, sbp, eqn, opts, ctx_residual, t)
  rand!(lo.A.nzval)
  for i=1:mesh.numDof
    lo.A[i, i] = 5  # make non-singular
  end

  A_orig = copy(lo.A)
  b = rand(mesh.numDof)

  # test product
  c = lo.A*b
  c2 = zeros(c)
  applyLinearOperator(ls.lo, mesh, sbp, eqn, opts, ctx_residual, t, b, c2)
  @fact norm(c - c2) --> roughly(0.0, atol=1e-12)

  # test transpose product
  c = A_orig.'*b
  applyLinearOperatorTranspose(ls.lo, mesh, sbp, eqn, opts, ctx_residual, t, b, c2)
  @fact norm(c - c2) --> roughly(0.0, atol=1e-12)

  # test solve
  x = factorize(lo.A)\b
  x2 = zeros(x)
  linearSolve(ls, b, x2)
  @fact norm(x2 - x) --> roughly(0.0, atol=1e-12)

  # test PC
  x3 = zeros(x2)
  applyPC(ls, mesh, sbp, eqn, opts, t, b, x3)
  @fact norm(x3 - x2) --> roughly(0.0, atol=1e-12)

  # test transpose solve
  x = A_orig.'\b
  x2 = zeros(x)
  linearSolveTranspose(ls, b, x2)
  @fact norm(x2 - x) --> roughly(0.0, atol=1e-12)

  # test PC
  x3 = zeros(x2)
  applyPCTranspose(ls, mesh, sbp, eqn, opts, t, b, x3)
  @fact norm(x3 - x2) --> roughly(0.0, atol=1e-12)




  
  # test product (again)
  c = lo.A*b
  c2 = zeros(c)
  applyLinearOperator(ls.lo, mesh, sbp, eqn, opts, ctx_residual, t, b, c2)
  @fact norm(c - c2) --> roughly(0.0, atol=1e-12)

  # test transpose product (again)
  c = A_orig.'*b
  applyLinearOperatorTranspose(ls.lo, mesh, sbp, eqn, opts, ctx_residual, t, b, c2)
  @fact norm(c - c2) --> roughly(0.0, atol=1e-12)

  # test solve (again)
  x = factorize(lo.A)\b
  x2 = zeros(x)
  linearSolve(ls, b, x2)
  @fact norm(x2 - x) --> roughly(0.0, atol=1e-12)
  @fact ls.lo.nfactorizations --> 1
  @fact ls.lo.nsolves --> 3

  # test transpose solve
  x = lo.A.'\b
  x2 = zeros(x)
  linearSolveTranspose(ls, b, x2)
  @fact norm(x2 - x) --> roughly(0.0, atol=1e-12)
  @fact ls.lo.nfactorizations --> 1
  @fact ls.lo.ntsolves --> 3


  free(ls)

  return nothing
end


function test_petscmat(mesh, sbp, eqn, opts)

  println("Testing Petsc matrix-explicit linear operator")
  if PetscInitialized() == 0
    PetscInitialize()
  end
  PetscSetOptions(opts["petsc_options"])

  vals = rand(mesh.numDof, mesh.numDof)
  for i=1:mesh.numDof
    vals[i, i] = 5  # make nearly diagonally dominant
  end

  pc = PetscMatPC(mesh, sbp, eqn, opts)
  lo = PetscMatLO(pc, mesh, sbp, eqn, opts)
  ls = StandardLinearSolver(pc, lo, eqn.comm)
  ctx_residual = (evalResidual,)
  t = 0.0

  # for testing
  MatSetOption(lo.A, PETSc.MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE)

  calcLinearOperator(ls, mesh, sbp, eqn, opts, ctx_residual, t)
  idx = collect(PetscInt, 1:mesh.numDof)
  idy = copy(idx)
  set_values1!(lo.A, idx, idy, vals, PETSC_ADD_VALUES)

  b = rand(mesh.numDof)

  # test product
  c = vals*b
  c2 = zeros(c)
  applyLinearOperator(ls.lo, mesh, sbp, eqn, opts, ctx_residual, t, b, c2)
  @fact norm(c - c2) --> roughly(0.0, atol=1e-12)

  # test transpose product
  c = vals.'*b
  applyLinearOperatorTranspose(ls.lo, mesh, sbp, eqn, opts, ctx_residual, t, b, c2)
  @fact norm(c - c2) --> roughly(0.0, atol=1e-12)

  setTolerances(ls, 1e-16, 1e-16, -1, -1)
  # test solve
  x = vals\b
  x2 = zeros(x)
  linearSolve(ls, b, x2)
  @fact norm(x2 - x) --> roughly(0.0, atol=1e-12)

  # test transpose solve
  x = vals.'\b
  x2 = zeros(x)
  linearSolveTranspose(ls, b, x2)
  @fact norm(x2 - x) --> roughly(0.0, atol=1e-12)

  # test product (again)
  c = vals*b
  c2 = zeros(c)
  applyLinearOperator(ls.lo, mesh, sbp, eqn, opts, ctx_residual, t, b, c2)
  @fact norm(c - c2) --> roughly(0.0, atol=1e-12)
  @fact lo.nassemblies[1] --> 1

  # test transpose product (again)
  c = vals.'*b
  applyLinearOperatorTranspose(ls.lo, mesh, sbp, eqn, opts, ctx_residual, t, b, c2)
  @fact norm(c - c2) --> roughly(0.0, atol=1e-12)
  @fact lo.nassemblies[1] --> 1

  # test solve (again)
  x = vals\b
  x2 = zeros(x)
  linearSolve(ls, b, x2)
  @fact norm(x2 - x) --> roughly(0.0, atol=1e-12)
  @fact lo.nassemblies[1] --> 1
  @fact lo.nsolves --> 2

  # test transpose solve (again)
  x = vals.'\b
  x2 = zeros(x)
  linearSolveTranspose(ls, b, x2)
  @fact norm(x2 - x) --> roughly(0.0, atol=1e-12)
  @fact lo.nassemblies[1] --> 1
  @fact lo.ntsolves --> 2
  @fact pc.nsetups --> 1



  free(ls)

  return nothing
end


function test_petscmatfree(mesh, sbp, eqn, opts)

  println("testing matrix-free operator")

  vals = rand(mesh.numDof, mesh.numDof)
  for i=1:mesh.numDof
    vals[i, i] = 20  # make really diagonally dominant
  end

  # set petsc options
  petsc_opts = opts["petsc_options"]
  petsc_opts["-pc_type"] = "shell"
  PetscSetOptions(petsc_opts)

  pc = TestMatFreePC(mesh, sbp, eqn, opts)
  lo = TestMatFreeLO(pc, mesh, sbp, eqn, opts)
  ls = StandardLinearSolver(pc, lo, eqn.comm)
  ctx_residual = (evalResidual, vals)
  t = 0.0

  setTolerances(ls, 1e-16, 1e-16, -1, -1)

  calcPC(ls, mesh ,sbp, eqn, opts, ctx_residual, t)
  calcLinearOperator(ls, mesh, sbp, eqn, opts, ctx_residual, t)

  b = rand(mesh.numDof)

  # test product
  c = vals*b
  c2 = zeros(c)
  applyLinearOperator(ls.lo, mesh, sbp, eqn, opts, ctx_residual, t, b, c2)
  @fact norm(c - c2) --> roughly(0.0, atol=1e-12)

  # test transpose product
  c = vals.'*b
  applyLinearOperatorTranspose(ls.lo, mesh, sbp, eqn, opts, ctx_residual, t, b, c2)
  @fact norm(c - c2) --> roughly(0.0, atol=1e-12)

  # test PC product
  x = inv(diagm(diag(vals)))*b
  x2 = zeros(x)
  applyPC(ls, mesh, sbp, eqn, opts, t, b, x2)
  @fact norm(x2 - x) --> roughly(0.0, atol=1e-12)

  # test PC transpose product
  x = inv(diagm(diag(vals)))*b
  x2 = zeros(x)
  applyPCTranspose(ls, mesh, sbp, eqn, opts, t, b, x2)
  @fact norm(x2 - x) --> roughly(0.0, atol=1e-12)

  # test solve
  x = vals\b
  x2 = zeros(x)
  linearSolve(ls, b, x2)
  @fact norm(x2 - x) --> roughly(0.0, atol=1e-12)

  # test transpose solve
  x = vals.'\b
  x2 = zeros(x)
  linearSolveTranspose(ls, b, x2)
  @fact norm(x2 - x) --> roughly(0.0, atol=1e-12)

  ### Test everything again

  # test product
  c = vals*b
  c2 = zeros(c)
  applyLinearOperator(ls.lo, mesh, sbp, eqn, opts, ctx_residual, t, b, c2)
  @fact norm(c - c2) --> roughly(0.0, atol=1e-12)

  # test transpose product
  c = vals.'*b
  applyLinearOperatorTranspose(ls.lo, mesh, sbp, eqn, opts, ctx_residual, t, b, c2)
  @fact norm(c - c2) --> roughly(0.0, atol=1e-12)

  # test PC product
  x = inv(diagm(diag(vals)))*b
  x2 = zeros(x)
  applyPC(ls, mesh, sbp, eqn, opts, t, b, x2)
  @fact norm(x2 - x) --> roughly(0.0, atol=1e-12)

  # test PC transpose product
  x = inv(diagm(diag(vals)))*b
  x2 = zeros(x)
  applyPCTranspose(ls, mesh, sbp, eqn, opts, t, b, x2)
  @fact norm(x2 - x) --> roughly(0.0, atol=1e-12)

  # test solve
  x = vals\b
  x2 = zeros(x)
  linearSolve(ls, b, x2)
  @fact norm(x2 - x) --> roughly(0.0, atol=1e-12)
  @fact lo.lo_inner.nsolves --> 2

  # test transpose solve
  x = vals.'\b
  x2 = zeros(x)
  linearSolveTranspose(ls, b, x2)
  @fact norm(x2 - x) --> roughly(0.0, atol=1e-12)
  @fact lo.lo_inner.ntsolves --> 2


  free(ls)

  return nothing
end


function test_petscmat_matfree(mesh, sbp, eqn, opts)
# test matrix-free linear operator with a matrix-free PC
  println("testing matrix-free operator with matrix-explicit PC")

  vals = rand(mesh.numDof, mesh.numDof)
  for i=1:mesh.numDof
    vals[i, i] = 20  # make really diagonally dominant
  end

  # set petsc options
  petsc_opts = opts["petsc_options"]
  petsc_opts["-pc_type"] = "jacobi"
  PetscSetOptions(petsc_opts)

  pc = TestMatPC(mesh, sbp, eqn, opts)
  lo = TestMatFreeLO(pc, mesh, sbp, eqn, opts)
  ls = StandardLinearSolver(pc, lo, eqn.comm)
  ctx_residual = (evalResidual, vals)
  t = 0.0

  setTolerances(ls, 1e-16, 1e-16, -1, -1)

  # for testing
  MatSetOption(getBasePC(pc).Ap, PETSc.MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE)


  calcPC(ls, mesh ,sbp, eqn, opts, ctx_residual, t)
  calcLinearOperator(ls, mesh, sbp, eqn, opts, ctx_residual, t)

  b = rand(mesh.numDof)

  # test product
  c = vals*b
  c2 = zeros(c)
  applyLinearOperator(ls.lo, mesh, sbp, eqn, opts, ctx_residual, t, b, c2)
  @fact norm(c - c2) --> roughly(0.0, atol=1e-12)

  # test transpose product
  c = vals.'*b
  applyLinearOperatorTranspose(ls.lo, mesh, sbp, eqn, opts, ctx_residual, t, b, c2)
  @fact norm(c - c2) --> roughly(0.0, atol=1e-12)

  # test PC product
  x = inv(diagm(diag(vals)))*b
  x2 = zeros(x)
  applyPC(ls, mesh, sbp, eqn, opts, t, b, x2)
  @fact norm(x2 - x) --> roughly(0.0, atol=1e-12)

  # test PC transpose product
  x = inv(diagm(diag(vals)))*b
  x2 = zeros(x)
  applyPCTranspose(ls, mesh, sbp, eqn, opts, t, b, x2)
  @fact norm(x2 - x) --> roughly(0.0, atol=1e-12)

  # test solve
  x = vals\b
  x2 = zeros(x)
  linearSolve(ls, b, x2)
  @fact norm(x2 - x) --> roughly(0.0, atol=1e-12)

  # test transpose solve
  x = vals.'\b
  x2 = zeros(x)
  linearSolveTranspose(ls, b, x2)
  @fact norm(x2 - x) --> roughly(0.0, atol=1e-12)


  free(ls)


  return nothing
end

#------------------------------------------------------------------------------
# define PC and LO
import LinearSolvers: calcPC, applyPC, applyPCTranspose, calcLinearOperator,
                      applyLinearOperator, applyLinearOperatorTranspose

# define matrix-free precondtioner for testing
# diagonal preconditioning
type TestMatFreePC <: AbstractPetscMatFreePC
  pc_inner::PetscMatFreePC
  diag::Array{Float64, 1}  # inverse of diagonal of matrix
end

function TestMatFreePC(mesh::AbstractMesh, sbp::AbstractSBP,
                       eqn::AbstractSolutionData, opts)

  pc_inner = PetscMatFreePC(mesh, sbp, eqn, opts)
  diag = Array(Float64, mesh.numDof)

  return TestMatFreePC(pc_inner, diag)
end

function calcPC(pc::TestMatFreePC, mesh::AbstractMesh, sbp::AbstractSBP,
                eqn::AbstractSolutionData, opts::Dict, ctx_residual, t)

  f = ctx_residual[1]
  vals = ctx_residual[2]  # for testing only

  for i=1:size(vals, 1)
    pc.diag[i] = 1/vals[i, i]
  end

  setPCCtx(pc, mesh, sbp, eqn, opts, ctx_residual, t)
  pc2 = getBasePC(pc)
  pc2.nsetups += 1

  return nothing
end


function applyPC(pc::TestMatFreePC, mesh::AbstractMesh, sbp::AbstractSBP,
                eqn::AbstractSolutionData, opts::Dict, t,
                b::AbstractVector, x::AbstractVector)

  for i=1:length(pc.diag)
    x[i] = pc.diag[i]*b[i]
  end

  pc2 = getBasePC(pc)
  pc2.napplies += 1

  return nothing
end

function applyPCTranspose(pc::TestMatFreePC, mesh::AbstractMesh,
                sbp::AbstractSBP,
                eqn::AbstractSolutionData, opts::Dict, t,
                b::AbstractVector, x::AbstractVector)

  for i=1:length(pc.diag)
    x[i] = pc.diag[i]*b[i]
  end

  pc2 = getBasePC(pc)
  pc2.ntapplies += 1

  return nothing
end

#------------------------------------------------------------------------------
# define matrix-explicit PC
type TestMatPC <: AbstractPetscMatPC
  pc_inner::PetscMatPC
end

function TestMatPC(mesh::AbstractMesh, sbp::AbstractSBP,
                       eqn::AbstractSolutionData, opts)

  pc_inner = PetscMatPC(mesh, sbp, eqn, opts)

  return TestMatPC(pc_inner)
end

function calcPC(pc::TestMatPC, mesh::AbstractMesh, sbp::AbstractSBP,
                eqn::AbstractSolutionData, opts::Dict, ctx_residual, t)

  f = ctx_residual[1]
  vals = ctx_residual[2]  # for testing only

  pc2 = getBasePC(pc)
  calcPC(pc2, mesh, sbp, eqn, opts, ctx_residual, t)

  idx = collect(PetscInt, 1:mesh.numDof)
  idy = copy(idx)
  set_values1!(pc2.Ap, idx, idy, vals, PETSC_ADD_VALUES)

  return nothing
end

# other functions do not need to be defined for matrix-explicit PC


#------------------------------------------------------------------------------
# define LO

# define fake matrix-free linear operator for testing
type TestMatFreeLO <: AbstractPetscMatFreeLO
  lo_inner::PetscMatFreeLO
  vals::Array{Float64, 2}  # the matrix
end

function TestMatFreeLO(pc::AbstractPC, mesh::AbstractMesh, sbp::AbstractSBP,
                       eqn::AbstractSolutionData, opts)


  lo_inner = PetscMatFreeLO(pc, mesh, sbp, eqn, opts)
  vals = zeros(mesh.numDof, mesh.numDof)

  return TestMatFreeLO(lo_inner, vals)
end


function calcLinearOperator(lo::TestMatFreeLO, mesh::AbstractMesh,
                sbp::AbstractSBP,
                eqn::AbstractSolutionData, opts::Dict, ctx_residual, t)

  f = ctx_residual[1]
  vals = ctx_residual[2]
  copy!(lo.vals, vals)

  setLOCtx(lo, mesh, sbp, eqn, opts, ctx_residual, t)

  return nothing
end

function applyLinearOperator(lo::TestMatFreeLO, mesh::AbstractMesh,
                sbp::AbstractSBP,
                eqn::AbstractSolutionData, opts::Dict, ctx_residual, t,
                x::AbstractVector, b::AbstractVector)

  smallmatvec!(lo.vals, x, b)

  return nothing
end

function applyLinearOperatorTranspose(lo::TestMatFreeLO, mesh::AbstractMesh,
                sbp::AbstractSBP,
                eqn::AbstractSolutionData, opts::Dict, ctx_residual, t,
                x::AbstractVector, b::AbstractVector)

  smallmatTvec!(lo.vals, x, b)

  return nothing
end





add_func1!(AdvectionTests, test_linearsolver, [TAG_TMP, TAG_SHORTTEST])


