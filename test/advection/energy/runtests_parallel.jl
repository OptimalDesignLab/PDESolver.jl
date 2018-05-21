function test_energy(mesh, sbp, eqn, opts)
  energy_final = calcNorm(eqn, eqn.q_vec)
  q_initial = zeros(eqn.q_vec)
  ICfunc = AdvectionEquationMod.ICDict[opts["IC_name"]]
  ICfunc(mesh, sbp, eqn, opts, q_initial)
  energy_initial = calcNorm(eqn, q_initial)

  @fact abs(energy_initial - energy_final) --> less_than(1e-12)
end


function test_energy_parallel()
  facts("----- Testing Energy Stability -----") do

    start_dir = pwd()

    ARGS[1] = "input_vals_periodic.jl"
    cd("./2dp4")
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])
    test_energy(mesh, sbp, eqn, opts)

    cd("../3dp4")
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])
    test_energy(mesh, sbp, eqn, opts)


  end

  return nothing
end

test_energy_parallel()
