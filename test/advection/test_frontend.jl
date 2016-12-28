function test_frontend()

  facts("----- Testing PDESolver Module frontend -----") do
    physics_name = "Advection"
    @fact haskey(PDESolver.PhysicsModDict, physics_name) --> true
    mod, func = retrieve_physics(physics_name)
    # FactCheck doesn't compare functions correctly, so use == instead
    @fact mod == AdvectionEquationMod --> true
    @fact func == run_advection --> true


    # test run_solver()
    ARGS[1] = "input_vals_channel.jl"
    mesh, sbp, eqn, opts = run_advection(ARGS[1])
    mesh, sbp, eqn2, opts = run_solver(ARGS[1])

    @fact eqn.q_vec --> roughly(eqn2.q_vec, atol=1e-13)

  end  # end facts block

  return nothing
    
end

add_func1!(AdvectionTests, test_frontend, [TAG_FRONTEND])
