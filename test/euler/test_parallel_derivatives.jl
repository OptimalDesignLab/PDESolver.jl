# test_parallel_derivatives.jl

function test_parallel_derivatives()

  facts("--- Checking Complete Derivative w.r.t angle of attack for a 2D airfoil ---") do

    ARGS[1] = "./input_vals_airfoil_parallel.jl"

    # Get the adjoint vector
    include("../../src/solver/euler/startup.jl")
    objective = EulerEquationMod.createObjectiveFunctionalData(mesh, sbp, eqn, opts)
    EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, objective)

    context("Checking functional J and ∂J/∂aoa in parallel") do

      serial_vals = readdlm("BoundaryForceData_values.dat")
      @fact objective.lift_val --> roughly(serial_vals[1], atol = 1e-14)
      @fact objective.drag_val --> roughly(serial_vals[2], atol = 1e-14)
      @fact objective.dLiftdaoa --> roughly(serial_vals[3], atol = 1e-14)
      @fact objective.dDragdaoa --> roughly(serial_vals[4], atol = 1e-14)

      # @fact objective.lift_val --> roughly(0.008165584931809718, atol = 1e-14)
      # @fact objective.drag_val --> roughly(-0.0005471055309820266, atol = 1e-14)
      # @fact objective.dLiftdaoa --> roughly(0.0005471055309820266, atol = 1e-14)
      # @fact objective.dDragdaoa --> roughly(0.008165584931809718, atol = 1e-14)
    end # End context("Checking ∂J/∂aoa")

    lift = objective.lift_val
    pc, lo = getNewtonPCandLO(mesh, sbp, eqn, opts)
    ls = StandardLinearSolver(pc, lo, eqn.comm, opts)
    adjoint_vec = zeros(Complex128, mesh.numDof)


    calcAdjoint(mesh, sbp, eqn, opts, ls, objective, adjoint_vec, recalc_jac=true, recalc_pc=true)
    dJdaoa = EulerEquationMod.eval_dJdaoa(mesh, sbp, eqn, opts, objective, "lift", adjoint_vec)

    # Check complete derivatives w.r.t alpha using finite difference
    pert = 1e-6
    eqn.params.aoa += pert
    EulerEquationMod.solve_euler(mesh, sbp, eqn, opts, mesh)
    EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, objective)
    lift_pert = objective.lift_val
    drag_pert = objective.drag_val
    eqn.params.aoa -= pert

    dJdaoa_fd = (lift_pert - lift)/pert

    context("Check complete derivative dJ/daoa") do
      deriv_err = norm(dJdaoa_fd - dJdaoa, 2)
      @fact deriv_err --> roughly(0.0, atol=1e-6)
    end

  end # End Facts Checking Complete Derivative w.r.t angle of attack for a 2D airfoil

  return nothing
end

add_func1!(EulerTests, test_parallel_derivatives, [TAG_PARALLEL_DERIVATIVES, TAG_LONGTEST])
