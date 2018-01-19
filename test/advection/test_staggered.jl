# tests for staggered grid

function test_staggered()

  fname = "input_vals_staggered.jl"
  mesh, sbp, eqn, opts = run_solver(fname)

  # test volume integrals when both grids are same
  opts["use_staggered_grid"] = false
  fill!(eqn.res, 0.0)
  AdvectionEquationMod.evalVolumeIntegrals(mesh, sbp, eqn, opts)

  res_orig = copy(eqn.res)
  fill!(eqn.res, 0.0)
  AdvectionEquationMod.calcVolumeIntegralsStaggered(mesh, mesh.mesh2, sbp, mesh.sbp2, eqn, opts)

  @fact vecnorm(eqn.res - res_orig) --> roughly(0.0, atol=1e-13)

  # test face integrals
  fill!(eqn.res, 0.0)
  AdvectionEquationMod.evalFaceIntegrals(mesh, sbp, eqn, opts)

  res_orig = copy(eqn.res)
  fill!(eqn.res, 0.0)
  AdvectionEquationMod.calcFaceIntegralsStaggered_nopre(mesh, mesh.mesh2, sbp, mesh.sbp2, eqn, opts, eqn.flux_func)

  @fact vecnorm(eqn.res - res_orig) --> roughly(0.0, atol=1e-13)



  # use different operators on the grids
  delete!(opts, "use_staggered_grid")
  opts["operator_type2"] = "SBPDiagonalE"
  opts["order2"] = 3
  opts["IC_name"] = "ICp1"
  opts["order"] = 2

  fname = make_input(opts, "input_vals_staggered2")

  mesh, sbp, eqn, opts = run_solver(fname)

  # check if the D operators are the same
  println("checking D")
  for d=1:2
    sbp2 = mesh.sbp2
    Hinv = inv(diagm(sbp2.w))
    D_tilde = mesh.I_F2S*Hinv*sbp2.Q[:, :, d]*mesh.I_S2F

    Hinv = inv(diagm(sbp.w))
    D_exact = Hinv*sbp.Q[:, :, d]
  end

  println("checking different operators")
  # because the operators are exact in this case, the residuals should be the
  # same
  opts["use_staggered_grid"] = false
  fill!(eqn.res, 0.0)
  AdvectionEquationMod.evalVolumeIntegrals(mesh, sbp, eqn, opts)

  res_orig = copy(eqn.res)
  fill!(eqn.res, 0.0)
  AdvectionEquationMod.calcVolumeIntegralsStaggered(mesh, mesh.mesh2, sbp, mesh.sbp2, eqn, opts)

  @fact vecnorm(eqn.res - res_orig) --> roughly(0.0, atol=1e-13)

  # test face integrals
  fill!(eqn.res, 0.0)
  AdvectionEquationMod.evalFaceIntegrals(mesh, sbp, eqn, opts)

  res_orig = copy(eqn.res)
  fill!(eqn.res, 0.0)
  AdvectionEquationMod.calcFaceIntegralsStaggered_nopre(mesh, mesh.mesh2, sbp, mesh.sbp2, eqn, opts, eqn.flux_func)

  @fact vecnorm(eqn.res - res_orig) --> roughly(0.0, atol=1e-13)




  return nothing
end

add_func1!(AdvectionTests, test_staggered, [TAG_SHORTTEST, TAG_STAGGERED])
