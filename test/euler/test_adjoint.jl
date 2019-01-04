#Test functional Integrate and adjoint for euler equation.

@doc """
Euler Equation -- test_adjoint

The function tests for the correctness of objective and non-objective functional
computation and then tests if the adjoint vector is being computed correctly.
Adjoint vector computation is based on an objective function. This is a serial
test and uses the input file called

`input_vals_vortex_adjoint_DG.jl`

"""->

function test_adjoint()

  @testset "Checking complete derivative of a functional using adjoint vector" begin
    # println("testing adjoint functions\n")
    #fname = "input_vals_vortex_adjoint_DG.jl"
    fname = "input_vals_airfoil.jl"
    mesh, sbp, eqn, opts = solvePDE(fname)
    @assert opts["aoa"] == 2.0

    lift = createFunctional(mesh, sbp, eqn, opts, 1)

    h = 1e-20
    pert = Complex128(0, h)

    pc, lo = getNewtonPCandLO(mesh, sbp, eqn, opts)
    ls = StandardLinearSolver(pc, lo, eqn.comm, opts)
    adjoint_vec = zeros(Complex128, mesh.numDof)
    calcAdjoint(mesh, sbp, eqn, opts, ls, lift, adjoint_vec, recalc_jac=true, recalc_pc=true)

    # compute dR/daoa and dJ/daoa
    eqn.params.aoa += pert
    evalResidual(mesh, sbp, eqn, opts)
    array3DTo1D(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
    dRdaoa = imag(eqn.res_vec)/h

    val = evalFunctional(mesh, sbp, eqn, opts, lift)
    dJdaoa = imag(val)/h
    eqn.params.aoa -= pert

    DJDaoa = dJdaoa + dot(adjoint_vec, dRdaoa)

    # save value for parallel tests
    fname = "BoundaryForceData_values.dat"
    file = open(fname, "w")
    println(file, val)
    close(file)

    # Check complete derivatives w.r.t alpha using finite difference
    pert = 1e-6
    eqn.params.aoa += pert
    solvePDE(mesh, sbp, eqn, opts, mesh)
    val2 = evalFunctional(mesh, sbp, eqn, opts, lift)
    eqn.params.aoa -= pert

    DJDaoa_fd = (val2 - real(val))/pert
    err_val = norm(DJDaoa - DJDaoa_fd, 2)

    @test isapprox( err_val, 0.0) atol= 1e-6

  end # End  testset("Checking complete derivative of a functional using adjoint vector")

end # End function test_adjoint

add_func1!(EulerTests, test_adjoint, [TAG_ADJOINT, TAG_LONGTEST])
