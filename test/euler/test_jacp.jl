"""
  Does basic testing of explicit jacobian calculation in parallel
"""
function test_jac_parallel()

  # SBPOmega
  fname = "input_vals_jac3dp.jl"
  fname2 = "input_vals_jac_tmp.jl"
  mesh, sbp, eqn, opts = run_solver(fname)

  MPI.Barrier(mesh.comm)
  test_jac_parallel_inner(mesh, sbp, eqn, opts)
  opts["preallocate_jacobian_coloring"] = true
  test_jac_parallel_inner(mesh, sbp, eqn, opts, is_prealloc_exact=true, set_prealloc=false)

  return nothing
end


add_func1!(EulerTests, test_jac_parallel, [TAG_SHORTTEST, TAG_JAC]) 


"""
  Does more thorough testing of jacobian calculation in parallel
"""
function test_jac_parallel_long()

  @testset "----- Testing jacobian assembly long -----" begin
    fname = "input_vals_jac3dp.jl"
    fname2 = "input_vals_jac_tmp.jl"

    myrank = MPI.Comm_rank(MPI.COMM_WORLD)

    # SBPGamma
    if myrank == 0
      fname4 = "input_vals_jac_tmp.jl"
      opts_tmp = read_input_file(fname)
      opts_tmp["operator_type"] = "SBPGamma"
      make_input(opts_tmp, fname2)
    end
    MPI.Barrier(MPI.COMM_WORLD)
    mesh4, sbp4, eqn4, opts4 = run_solver(fname2)

    # SBPDiagonalE
    if myrank == 0
      opts_tmp = read_input_file(fname)
      opts_tmp["operator_type"] = "SBPDiagonalE"
      make_input(opts_tmp, fname2)
    end
    MPI.Barrier(MPI.COMM_WORLD)
    mesh5, sbp5, eqn5, opts5 = run_solver(fname2)

    if myrank == 0
      opts_tmp = read_input_file(fname)
      opts_tmp["operator_type"] = "SBPDiagonalE"
      opts_tmp["use_Minv"] = true
      make_input(opts_tmp, fname2)
    end
    MPI.Barrier(MPI.COMM_WORLD)
    mesh6, sbp6, eqn6, opts6 = run_solver(fname2)

    opts4_tmp = copy(opts4)
    test_jac_parallel_inner(mesh4, sbp4, eqn4, opts4)
    test_jac_homotopy(mesh4, sbp4, eqn4, opts4_tmp)


    test_jac_parallel_inner(mesh5, sbp5, eqn5, opts5)

    # run test twice to make sure arrays are zeroed out correctly
    test_jac_parallel_inner(mesh5, sbp5, eqn5, opts5)

    test_jac_parallel_inner(mesh6, sbp6, eqn6, opts6)

  end

  return nothing
end

add_func1!(EulerTests, test_jac_parallel_long, [TAG_LONGTEST, TAG_JAC]) 



function test_jac_parallel_inner(mesh, sbp, eqn, opts; is_prealloc_exact=true, set_prealloc=true)

  # use a spatially varying solution
  icfunc = EulerEquationMod.ICDict["ICExp"]
  icfunc(mesh, sbp, eqn, opts, eqn.q_vec)
  array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)

  # get the correct differentiated flux function (this is needed because the
  # input file set calc_jac_explicit = false
  eqn.flux_func_diff = EulerEquationMod.FluxDict_diff[opts["Flux_name"]]


  startSolutionExchange(mesh, sbp, eqn, opts)

  opts["calc_jac_explicit"] = false
  pc1, lo1 = NonlinearSolvers.getNewtonPCandLO(mesh, sbp, eqn, opts)
  opts["calc_jac_explicit"] = true
  val_orig = opts["preallocate_jacobian_coloring"]
  if set_prealloc
    opts["preallocate_jacobian_coloring"] = false
  end
  pc2, lo2 = NonlinearSolvers.getNewtonPCandLO(mesh, sbp, eqn, opts)
  opts["preallocate_jacobian_coloring"] = val_orig

  jac1 = getBaseLO(lo1).A
  jac2 = getBaseLO(lo2).A

  assembler = NonlinearSolvers._AssembleElementData(getBaseLO(lo2).A, mesh, sbp, eqn, opts)

  # compute jacobian via coloring
  opts["calc_jac_explicit"] = false
  ctx_residual = (evalResidual,)
  NonlinearSolvers.physicsJac(mesh, sbp, eqn, opts, jac1, ctx_residual)

  # compute jacobian explicitly
  opts["calc_jac_explicit"] = true
  evalJacobian(mesh, sbp, eqn, opts, assembler)

  assembly_begin(jac1, MAT_FINAL_ASSEMBLY)
  assembly_begin(jac2, MAT_FINAL_ASSEMBLY)

  # multiply against a random vector to make sure the jacobian is
  # the same
  for i=1:10
    x = rand(PetscScalar, mesh.numDof)
    b1 = zeros(PetscScalar, mesh.numDof)
    b2 = zeros(PetscScalar, mesh.numDof)

    t = 0.0
    applyLinearOperator(lo1, mesh, sbp, eqn, opts, ctx_residual, t, x, b1)
    applyLinearOperator(lo2, mesh, sbp, eqn, opts, ctx_residual, t, x, b2)

    @test isapprox( norm(b1 - b2), 0.0) atol=1e-12
  end

  A = getBaseLO(lo2).A
  if typeof(A) <: PetscMat
    matinfo = MatGetInfo(A, PETSc2.MAT_LOCAL)
    if is_prealloc_exact
      @test ( matinfo.nz_unneeded )== 0
    else
      @test  matinfo.nz_unneeded  > 0
    end
  end


  free(lo1)
  free(lo2)
  free(pc1)
  free(pc2)

  return nothing
end


function test_jac_homotopy(mesh, sbp, eqn, opts)

  println("\nTesting homotopy jacobian")

  # use a spatially varying solution
  icfunc = EulerEquationMod.ICDict["ICExp"]
  icfunc(mesh, sbp, eqn, opts, eqn.q_vec)
  array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)

  # get the correct differentiated flux function (this is needed because the
  # input file set calc_jac_explicit = false
  eqn.flux_func_diff = EulerEquationMod.FluxDict_diff["RoeFlux"]

  res1 = zeros(eqn.res)
  res2 = zeros(eqn.res)
  println("\ncomputing regular homotopy dissipation")
#  EulerEquationMod.calcHomotopyDiss(mesh, sbp, eqn, opts, res1)
  println("\ncomputing new homotopy dissipation")
  h = 1e-20
  pert = Complex128(0, h)
  eqn.q[1] += pert
  EulerEquationMod.calcHomotopyDiss(mesh, sbp, eqn, opts, res2)
  eqn.q[1] -= pert
#=
  println("diffnorm = ", vecnorm(res1 - res2))
  println("res1 = \n", res1)
  println("res2 = \n", res2)
  println("diff = \n", res1 - res2)
  @assert vecnorm(res1 - res2) < 1e-13
=#
  startSolutionExchange(mesh, sbp, eqn, opts, wait=true)

  println("constructing first operator")
  opts["calc_jac_explicit"] = false
  pc1, lo1 = NonlinearSolvers.getHomotopyPCandLO(mesh, sbp, eqn, opts)

  println("constructing second operator")
  opts["calc_jac_explicit"] = true
  pc2, lo2 = NonlinearSolvers.getHomotopyPCandLO(mesh, sbp, eqn, opts)

  jac1 = getBaseLO(lo1).A
  jac2 = getBaseLO(lo2).A

  assembler = NonlinearSolvers._AssembleElementData(getBaseLO(lo2).A, mesh, sbp, eqn, opts)

  function _evalHomotopy(mesh, sbp, eqn, opts, t)
    evalHomotopy(mesh, sbp, eqn, opts, eqn.res, t)
  end

  ctx_residual = (_evalHomotopy,)
  println("\nevaluating jacobians")

  opts["calc_jac_explicit"] = false
  NonlinearSolvers.physicsJac(mesh, sbp, eqn, opts, jac1, ctx_residual)

  # compute jacobian explicitly
  opts["calc_jac_explicit"] = true

  
  evalHomotopyJacobian(mesh, sbp, eqn, opts, assembler, lo2.lambda)

  assembly_begin(jac1, MAT_FINAL_ASSEMBLY)
  assembly_begin(jac2, MAT_FINAL_ASSEMBLY)

  # multiply against a random vector to make sure the jacobian is
  # the same
  for i=1:10
    x = rand(PetscScalar, mesh.numDof)
    b1 = zeros(PetscScalar, mesh.numDof)
    b2 = zeros(PetscScalar, mesh.numDof)

    t = 0.0
    applyLinearOperator(lo1, mesh, sbp, eqn, opts, ctx_residual, t, x, b1)
    applyLinearOperator(lo2, mesh, sbp, eqn, opts, ctx_residual, t, x, b2)

    @test isapprox( norm(b1 - b2), 0.0) atol=1e-12
  end

  free(lo1)
  free(lo2)
  free(pc1)
  free(pc2)

  return nothing
end

