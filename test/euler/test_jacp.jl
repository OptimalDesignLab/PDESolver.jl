function test_jac_parallel()

  mesh, sbp, eqn, opts = run_solver("input_vals_jac3dp.jl")

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

#  opts["addBoundaryIntegrals"] = false #TODO

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

    @fact norm(b1 - b2) --> roughly(0.0, atol=1e-14)
  end

  free(lo1)
  free(lo2)
  free(pc1)
  free(pc2)

  return nothing
end


add_func1!(EulerTests, test_jac_parallel, [TAG_SHORTTEST, TAG_JAC]) 



