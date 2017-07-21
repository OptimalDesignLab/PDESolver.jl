# test staggered grid things

function test_staggered()

  fname = "input_vals_staggered.jl"
  mesh, sbp, eqn, opts = run_solver(fname)

  # evalute the residual on the original, second mesh, make sure they are
  # the same
  mesh2 = mesh.mesh2
  sbp2 = mesh.sbp2

  println("testing primary mesh")
  fill!(eqn.res, 0.0)
  evalResidual(mesh, sbp, eqn, opts)
  res_orig = copy(eqn.res)

  println("testing secondary mesh")
  fill!(eqn.res, 0.0)
  evalResidual(mesh2, sbp2, eqn, opts)

  @fact vecnorm(res_orig - eqn.res) --> roughly(0.0, atol=1e-13)

  # test a non-uniform condition
  EulerEquationMod.ICExp(mesh, sbp, eqn, opts, eqn.q_vec)
  disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)

  fill!(eqn.res, 0.0)
  evalResidual(mesh, sbp, eqn, opts)
  res_orig = copy(eqn.res)

  fill!(eqn.res, 0.0)
  evalResidual(mesh2, sbp2, eqn, opts)

  @fact vecnorm(res_orig - eqn.res) --> roughly(0.0, atol=1e-13)

  return nothing
end

add_func1!(EulerTests, test_staggered, [TAG_SHORTTEST, TAG_STAGGERED])


