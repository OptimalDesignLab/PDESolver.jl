function test_frontend()

  facts("----- Testing PDESolver Module frontend -----") do
    physics_name = "Advection"
    @fact haskey(PDESolver.PhysicsModDict, physics_name) --> true
    mod, func = retrieve_physics(physics_name)
    # FactCheck doesn't compare functions correctly, so use == instead
    @fact mod == AdvectionEquationMod --> true
    @fact func == run_advection --> true
  end  # end facts block

  return nothing
    
end

add_func1!(AdvectionTests, test_frontend)
