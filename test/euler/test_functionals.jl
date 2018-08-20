# test functional calculation

global const test_functionals_inputfile = "input_vals_channel.jl"
global const test_functionals_moddict = Dict{String, Any}(
  "Flux_name" => "RoeFlux",
  "use_DG" => true,
  "new_fname" => "input_vals_test_functionals.jl"
)

"""
  Tests for functional computations.  The input file loads a 2 element mesh,
  square, -1 to 1,
  with a uniform flow (ICRho1E2U3).  The functional calculation should be exact
  in this case, so we can test against an analytical value.  The entire
  boundary is on BC 1
"""
function test_functionals(mesh, sbp, eqn, opts)

  @testset "Testing functionals" begin
    obj = createFunctional(mesh, sbp, eqn, opts, "entropyflux", [1])
    val = evalFunctional(mesh, sbp, eqn, opts, obj)
    # the functional is u_i * U dot n_i, so for a constant field around a closed
    # curve it is zero
    @test isapprox(val, 0.0) atol=1e-13

    println("\n\n\nAbout to test entropy dissipation functional")
    testEntropyDissFunctional()

  end  # end testset

  return nothing
end


function testEntropyDissFunctional()

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
    )

  mesh, sbp, eqn, opts = solvePDE(opts)

  # compute the functional the regular way
  func = createFunctional(mesh, sbp, eqn, opts, "entropydissipation", Int[])
  J1 = evalFunctional(mesh, sbp, eqn, opts, func)

  # compute the functional using evalResidual
  opts["addVolumeIntegrals"] = false
  opts["addBoundaryIntegrals"] = false
  eqn.face_element_integral_func = EulerEquationMod.ELFPenaltyFaceIntegral()
  array3DTo1D(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
  w_vec = copy(eqn.q_vec)
  EulerEquationMod.convertToIR(mesh, sbp, eqn, opts, w_vec)

  J2 = dot(w_vec, eqn.res_vec)

  @test isapprox(J1, J2) atol=1e-13

  return nothing
end







add_func3!(EulerTests, test_functionals, test_functionals_inputfile, test_functionals_moddict,  [TAG_FUNCTIONAL, TAG_SHORTTEST, TAG_TMP])
