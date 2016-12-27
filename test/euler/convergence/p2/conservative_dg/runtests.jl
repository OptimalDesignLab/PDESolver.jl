#global const STARTUP_PATH = joinpath(Pkg.dir("PDESolver"), "src/solver/euler/startup.jl")
#using FactCheck

#run(`cptest.sh`)
#run(`cperr.sh`)

function test_convergence_p2_dg()
  facts("---- P2 Conservative DG Convergence Tests -----") do
    start_dir = pwd()

    resize!(ARGS, 1)

    cd("./m1")
    ARGS[1] = "input_vals_vortex3.jl"
    mesh, sbp, eqn, opts = run_euler(ARGS[1])
    ARGS[1] = "input_vals_vortex4.jl"
    mesh, sbp, eqn, opts = run_euler(ARGS[1])

    cd("../m2")
    ARGS[1] = "input_vals_vortex3.jl"
    mesh, sbp, eqn, opts = run_euler(ARGS[1])
    ARGS[1] = "input_vals_vortex4.jl"
    mesh, sbp, eqn, opts = run_euler(ARGS[1])

    cd("..")
    include("calc_line.jl")

    slope = calc_line()
    println("slope = ", slope)

    data = readdlm("err_data.dat")
    err_vals = data[:, 2]
    #println("err_vals = ", err_vals)

    slope_val = 3.00
    slope_margin = 0.1

    @fact slope => greater_than(slope_val - slope_margin)
    @fact slope => less_than(slope_val + slope_margin)

    err_val = 0.000877
    slope_fac = 1.25
    println("err_vals[1] = ", err_vals[1])
    @fact err_vals[1] => greater_than(err_val/slope_fac)
    @fact err_vals[1] => less_than(err_val*slope_fac)

  end

  return nothing
end

test_convergence_p2_dg()
