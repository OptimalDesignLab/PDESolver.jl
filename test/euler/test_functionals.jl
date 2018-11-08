# test functional calculation

"""
  Tests for functional computations.  The input file loads a 2 element mesh,
  square, -1 to 1,
  with a uniform flow (ICRho1E2U3).  The functional calculation should be exact
  in this case, so we can test against an analytical value.  The entire
  boundary is on BC 1
"""
function test_functionals()

  @testset "Testing functionals" begin


  # test all functional derivatives
  opts = Dict{String, Any}(
    "physics" => "Euler",
    "operator_type" => "SBPOmega",
    "dimensions" => 2,
    "run_type" => 5,
    "jac_method" => 2,
    "jac_type" => 2,
    "order" => 1,
    "IC_name" => "ICIsentropicVortex",
    "use_DG" => true,
    "volume_integral_type" => 2,
    "Volume_flux_name" => "IRFlux",
    "face_integral_type" => 2,
    "FaceElementIntegral_name" => "ESLFFaceIntegral",
    "Flux_name" => "IRFlux",
    "numBC" => 3,
    "BC1" => [0],
    "BC1_name" => "isentropicVortexBC",  # outlet
    "BC2" => [2],
    "BC2_name" => "isentropicVortexBC", # inlet
    "BC3" => [1, 3],
    "BC3_name" => "noPenetrationBC",  # was noPenetrationBC
    "aoa" => 0.0,
    "smb_name" => "SRCMESHES/vortex_3x3_.smb",
    "dmg_name" => ".null",
    "itermax" => 20,
    "res_abstol" => 1e-9,
    "res_reltol" => 1e-9,
    "do_postproc" => true,
    "exact_soln_func" => "ICIsentropicVortex",
    "force_solution_complex" => true,
    "force_mesh_complex" => true,
    "need_adjoint" => true,
    )

  mesh, sbp, eqn, opts = solvePDE(opts)

  testEntropyDissFunctional(mesh, sbp, eqn, opts)

  # test derivative of all functionals
  for funcname in keys(EulerEquationMod.FunctionalDict)
    println("testing functional", funcname)
    obj = createFunctional(mesh, sbp, eqn, opts, funcname, [1, 3])
    test_functional_deriv(mesh, sbp, eqn, opts, obj)
  end


  end  # end testset

  return nothing
end


"""
  Test the entropy dissipation functional.  This function is defined over all
  faces, rather than a boundary, so it is a bit special
"""
function testEntropyDissFunctional(mesh, sbp, eqn, opts)

  opts2 = read_input_file("input_vals_channel.jl")
  opts2["Flux_name"] = "RoeFlux"
  opts2["use_DG"] = true
  opts2["solve"] = false

  mesh2, sbp2, eqn2, opts2 = solvePDE(opts2)

  obj = createFunctional(mesh2, sbp2, eqn2, opts2, "entropyflux", [1])
  val = evalFunctional(mesh2, sbp2, eqn2, opts2, obj)
  # the functional is u_i * U dot n_i, so for a constant field around a closed
  # curve it is zero
  @test isapprox(val, 0.0) atol=1e-13


  # compute the functional the regular way
  func = createFunctional(mesh, sbp, eqn, opts, "entropydissipation", Int[])
  J1 = evalFunctional(mesh, sbp, eqn, opts, func)

  # compute the functional using evalResidual
  opts["addVolumeIntegrals"] = false
  opts["addBoundaryIntegrals"] = false
  eqn.face_element_integral_func = EulerEquationMod.ELFPenaltyFaceIntegral(mesh, eqn)
  array3DTo1D(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
  w_vec = copy(eqn.q_vec)
  EulerEquationMod.convertToIR(mesh, sbp, eqn, opts, w_vec)

  J2 = dot(w_vec, eqn.res_vec)

  @test isapprox(J1, J2) atol=1e-13

  return nothing
end


"""
  Test the derivative of a boundary functional wrt to q

  The eqn object must have been created with Tsol = Complex128 for this to work
"""
function test_functional_deriv(mesh, sbp, eqn, opts, func)

  h = 1e-20
  pert = Complex128(0, h)

  # use a spatially varying solution
  icfunc = EulerEquationMod.ICDict["ICRho1E2U3"]
  eqn.q_vec .+= 0.1*rand(length(eqn.q_vec))
  icfunc(mesh, sbp, eqn, opts, eqn.q_vec)
  array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)

  q_dot = rand(size(eqn.q))
#  q_dot = zeros(eqn.q)
#  q_dot[1, 1, 1] = 100
  q_bar = zeros(eqn.q)

  evalFunctionalDeriv_q(mesh, sbp, eqn, opts, func, q_bar)
  val = sum(q_bar .* q_dot)


  eqn.q_vec .+= pert*vec(q_dot)
  f = evalFunctional(mesh, sbp, eqn, opts, func)
  val2 = imag(f)/h
  eqn.q_vec .-= pert*vec(q_dot)

  @test abs(val - val2) < 1e-12

  return nothing
end

  


add_func1!(EulerTests, test_functionals, [TAG_FUNCTIONAL, TAG_SHORTTEST, TAG_TMP])
