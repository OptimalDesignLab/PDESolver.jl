# test staggered grid things

function test_staggered()

  fname = "input_vals_staggered.jl"
  mesh, sbp, eqn, opts = run_solver(fname)

  test_staggered_same(mesh, sbp, eqn, opts)


  
  delete!(opts, "use_staggered_grid")
  opts["operator_type"] = "SBPOmega"
  opts["operator_type2"] = "SBPGamma"
  opts["order"] = 2
  opts["order2"] = 3

  fname = make_input(opts, "input_vals_staggered2")

  mesh, sbp, eqn, opts = run_solver(fname)
 
  test_staggered_different(mesh, sbp, eqn, opts)


  return nothing
end

"""
  Test staggered grid functions when operators on both meshes are the samea

  This function is called from test_staggered()
"""
function test_staggered_same(mesh, sbp, eqn, opts)

  # force it to not use the staggered grid
  opts["use_staggered_grid"] = false
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


  # test staggered grid volume integrals
  println("testing volume integrals")
  fill!(eqn.res, 0.0)
  EulerEquationMod.calcVolumeIntegralsSplitFormCurvilinear(mesh, sbp, eqn, opts,
                                                           eqn.flux_func)

  println("testing staggered volume integrals")
  res_orig = copy(eqn.res)
  fill!(eqn.res, 0.0)
  EulerEquationMod.calcVolumeIntegralsSplitFormCurvilinear(mesh, mesh2, sbp, sbp2, eqn, opts, eqn.flux_func)

  @fact vecnorm(eqn.res - res_orig) --> roughly(0.0, atol=1e-13)


  fill!(eqn.res, 0.0)
  EulerEquationMod.getFaceElementIntegral(mesh, sbp, eqn, eqn.face_element_integral_func, eqn.flux_func, mesh.sbpface, mesh.interfaces)

  res_orig = copy(eqn.res)
  fill!(eqn.res, 0.0)
  EulerEquationMod.getFaceElementIntegral(mesh, mesh2, sbp, sbp2, eqn, eqn.face_element_integral_func, eqn.flux_func, mesh2.sbpface, mesh.interfaces)

  @fact vecnorm(eqn.res - res_orig) --> roughly(0.0, atol=1e-13)

  fill!(eqn.res, 0.0)

  return nothing
end

"""
  Test staggered grid functions when the meshes are different

  This function is called from test_staggered()
"""
function test_staggered_different(mesh, sbp, eqn, opts)

  mesh2 = mesh.mesh2
  sbp2 = mesh.sbp2

  # test interpolation
  # the operators are 
  q_ex = zeros(5, mesh.numNodesPerElement)
  q = zeros(q_ex)
  q2 = zeros(5, mesh.mesh2.numNodesPerElement)
  q2_ex = zeros(q2)

  # sbp_s is 2nd order, so make sure a 2nd order polynomial is interpolated
  # exactly
  for i=1:mesh.numNodesPerElement
    x = mesh.coords[1, i, 1]
    y = mesh.coords[2, i, 1]

    q_ex[1, i] = 1
    q_ex[2, i] = x + 1
    q_ex[3, i] = y + 1
    q_ex[4, i] = x + y + 1
    q_ex[5, i] = x*x + y*y + x + y + 1
  end
 
  for i=1:mesh2.numNodesPerElement
    x = mesh2.coords[1, i, 1]
    y = mesh2.coords[2, i, 1]

    q2_ex[1, i] = 1
    q2_ex[2, i] = x + 1
    q2_ex[3, i] = y + 1
    q2_ex[4, i] = x + y + 1
    q2_ex[5, i] = x*x + y*y + x + y + 1
  end

  smallmatmat!(q_ex, mesh.I_S2FT, q2)
  @fact norm(q2 - q2_ex) --> roughly(0.0, atol=1e-13)

  smallmatmat!(q2_ex, mesh.I_F2ST, q)
  @fact norm(q - q_ex) --> roughly(0.0, atol=1e-13)
 
  #=  
  # test staggered grid volume integrals
  # S is skew-symmetric, so c.' S c = 0
  # this is basically the test of conservation
  EulerEquationMod.ICRhoe1E2U3(mesh, sbp, eqn, opts, eqn.q_vec)
  disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)

  q_node = eqn.q[:, 1, 1]

  res_orig = copy(eqn.res)
  =#
  fill!(eqn.res, 0.0)
  EulerEquationMod.calcVolumeIntegralsSplitFormCurvilinear(mesh, mesh2, sbp, sbp2, eqn, opts, eqn.flux_func)
  @fact sum(eqn.res) --> roughly(0.0, atol=1e-13)


  fill!(eqn.res, 0.0)
  applyPoly(mesh, sbp, eqn, opts, sbp.degree)
  penalty_functor = EulerEquationMod.FaceElementDict["ELFPenaltyFaceIntegral"]

  EulerEquationMod.getFaceElementIntegral(mesh, mesh2, sbp, sbp2, eqn, penalty_functor, eqn.flux_func, mesh2.sbpface, mesh.interfaces)
  
  @fact vecnorm(eqn.res) --> roughly(0.0, atol=1e-13)

  return nothing
end

add_func1!(EulerTests, test_staggered, [TAG_SHORTTEST, TAG_STAGGERED])


