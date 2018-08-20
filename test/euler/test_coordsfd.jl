# test finite difference derivative with respect to the mesh vertex coordinates

"""
  Test complex stepping the derivative of the residual with repsect to the
  mesh coordinate nodes
"""
function test_coordsfd()

  println("testing coordinates finite difference")

  @testset "Testing Coordinate Field Derivative" begin
    fname = "input_vals_coordsfd.jl"
    mesh, sbp, eqn, opts = solvePDE(fname)
    func = createFunctional(mesh, sbp, eqn, opts, "lift",  [1, 2, 3])

    testMeshDerivatives(mesh, sbp, opts)
    testResidualDerivatives(mesh, sbp, eqn, opts)
    testDJDx(mesh, sbp, eqn, opts, func)

    func = createFunctional(mesh, sbp, eqn, opts, "entropydissipation", Int[])
    testDJDx(mesh, sbp, eqn, opts, func)
  end

end


add_func1!(EulerTests, test_coordsfd, [TAG_SHORTTEST, TAG_COMPLEX, TAG_ADJOINT, TAG_TMP])
