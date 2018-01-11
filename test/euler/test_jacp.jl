function test_jac_parallel()

  # SBPOmega
  fname = "input_vals_jac3dp.jl"
  fname2 = "input_vals_jac_tmp.jl"
  mesh, sbp, eqn, opts = run_solver(fname)

  # SBPGamma
  if mesh.myrank == 0
    fname4 = "input_vals_jac_tmp.jl"
    opts_tmp = read_input_file(fname)
    opts_tmp["operator_type"] = "SBPGamma"
    make_input(opts_tmp, fname2)
  end
  MPI.Barrier(mesh.comm)
  mesh4, sbp4, eqn4, opts4 = run_solver(fname2)

  # SBPDiagonalE
  if mesh.myrank == 0
    opts_tmp = read_input_file(fname)
    opts_tmp["operator_type"] = "SBPDiagonalE"
    make_input(opts_tmp, fname2)
  end
  MPI.Barrier(mesh.comm)
  mesh5, sbp5, eqn5, opts5 = run_solver(fname2)

  if mesh.myrank == 0
    opts_tmp = read_input_file(fname)
    opts_tmp["operator_type"] = "SBPDiagonalE"
    opts_tmp["use_Minv"] = true
    make_input(opts_tmp, fname2)
  end
  MPI.Barrier(mesh.comm)
  mesh6, sbp6, eqn6, opts6 = run_solver(fname2)

  MPI.Barrier(mesh.comm)
  test_jac_parallel_inner(mesh, sbp, eqn, opts)

  test_jac_parallel_inner(mesh4, sbp4, eqn4, opts4)

  test_jac_parallel_inner(mesh5, sbp5, eqn5, opts5)

  # run test twice to make sure arrays are zeroed out correctly
  test_jac_parallel_inner(mesh5, sbp5, eqn5, opts5)

  test_jac_parallel_inner(mesh6, sbp6, eqn6, opts6)

  return nothing
end


add_func1!(EulerTests, test_jac_parallel, [TAG_SHORTTEST, TAG_JAC]) 


function test_jac_parallel_inner(mesh, sbp, eqn, opts)

  # use a spatially varying solution
  icfunc = EulerEquationMod.ICDict["ICExp"]
  icfunc(mesh, sbp, eqn, opts, eqn.q_vec)
  disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)

  # get the correct differentiated flux function (this is needed because the
  # input file set calc_jac_explicit = false
  eqn.flux_func_diff = EulerEquationMod.FluxDict_diff[opts["Flux_name"]]


  startSolutionExchange(mesh, sbp, eqn, opts)

  pc1, lo1 = NonlinearSolvers.getNewtonPCandLO(mesh, sbp, eqn, opts)
  pc2, lo2 = NonlinearSolvers.getNewtonPCandLO(mesh, sbp, eqn, opts)

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

    @fact norm(b1 - b2) --> roughly(0.0, atol=1e-12)
  end

  free(lo1)
  free(lo2)
  free(pc1)
  free(pc2)

  return nothing
end




