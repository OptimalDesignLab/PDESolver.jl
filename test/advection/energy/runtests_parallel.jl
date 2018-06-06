function test_energy(mesh, sbp, eqn, opts)
  energy_final = calcNorm(eqn, eqn.q_vec)
  q_initial = zeros(eqn.q_vec)
  ICfunc = AdvectionEquationMod.ICDict[opts["IC_name"]]
  ICfunc(mesh, sbp, eqn, opts, q_initial)
  energy_initial = calcNorm(eqn, q_initial)

  @test  abs(energy_initial - energy_final)  < 1e-12
end


function test_energy_parallel()
  @testset "----- Testing Energy Stability -----" begin

    start_dir = pwd()

    fname = "input_vals_periodic.jl"
    cd("./2dp4")
    mesh, sbp, eqn, opts = solvePDE(fname)
    test_energy(mesh, sbp, eqn, opts)

    cd("../3dp4")
    mesh, sbp, eqn, opts = solvePDE(fname)
    test_energy(mesh, sbp, eqn, opts)


  end

  return nothing
end

test_energy_parallel()
