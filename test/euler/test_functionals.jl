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

  end

  return nothing
end




add_func3!(EulerTests, test_functionals, test_functionals_inputfile, test_functionals_moddict,  [TAG_FUNCTIONAL, TAG_SHORTTEST])
