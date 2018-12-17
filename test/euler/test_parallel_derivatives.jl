# test_parallel_derivatives.jl

function test_parallel_derivatives()

  @testset "--- Checking Complete Derivative w.r.t angle of attack for a 2D airfoil ---" begin

    fname = "./input_vals_airfoil_parallel.jl"

    # Get the adjoint vector
    mesh, sbp, eqn, opts = solvePDE(fname)
    lift = createFunctional(mesh, sbp, eqn, opts, 1)
    val = EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, lift)

    @testset "Checking functional J in parallel" begin

      serial_vals = readdlm("BoundaryForceData_values.dat")
      @test isapprox( val, serial_vals[1]) atol= 1e-14
    end # End  testset("Checking ∂J/∂aoa")

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

    term2 = dot(adjoint_vec, dRdaoa)
    term2 = MPI.Allreduce(term2, MPI.SUM, eqn.comm)
    DJDaoa = dJdaoa + term2

    # Check complete derivatives w.r.t alpha using finite difference
    pert = 1e-6
    eqn.params.aoa += pert
    solvePDE(mesh, sbp, eqn, opts, mesh)
    val2 = evalFunctional(mesh, sbp, eqn, opts, lift)
    eqn.params.aoa -= pert

    DJDaoa_fd = (val2 - real(val))/pert
    err_val = norm(DJDaoa - DJDaoa_fd, 2)

    @test isapprox( err_val, 0.0) atol= 1e-6






#=
    lift = lift.lift_val
    pc, lo = getNewtonPCandLO(mesh, sbp, eqn, opts)
    ls = StandardLinearSolver(pc, lo, eqn.comm, opts)
    adjoint_vec = zeros(Complex128, mesh.numDof)


    calcAdjoint(mesh, sbp, eqn, opts, ls, lift, adjoint_vec, recalc_jac=true, recalc_pc=true)
    dJdaoa = EulerEquationMod.eval_dJdaoa(mesh, sbp, eqn, opts, lift, "lift", adjoint_vec)

    # Check complete derivatives w.r.t alpha using finite difference
    pert = 1e-6
    eqn.params.aoa += pert
    solvePDE(mesh, sbp, eqn, opts, mesh)
    EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, lift)
    lift_pert = lift.lift_val
    drag_pert = lift.drag_val
    eqn.params.aoa -= pert

    dJdaoa_fd = (lift_pert - lift)/pert

    @testset "Check complete derivative dJ/daoa" begin
      deriv_err = norm(dJdaoa_fd - dJdaoa, 2)
      @test isapprox( deriv_err, 0.0) atol=1e-6
    end
=#
  end # End Facts Checking Complete Derivative w.r.t angle of attack for a 2D airfoil

  return nothing
end

add_func1!(EulerTests, test_parallel_derivatives, [TAG_PARALLEL_DERIVATIVES, TAG_LONGTEST, TAG_TMP])
