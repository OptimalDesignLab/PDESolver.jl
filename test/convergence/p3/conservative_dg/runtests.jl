global const STARTUP_PATH = joinpath(Pkg.dir("PDESolver"), "src/solver/euler/startup.jl")
using FactCheck

#run(`cptest.sh`)
#run(`cperr.sh`)
facts("---- P1 Conservative Convergence Tests -----") do
start_dir = pwd()

resize!(ARGS, 1)

cd("./m1")
ARGS[1] = "input_vals_vortex3.jl"
include(STARTUP_PATH)
ARGS[1] = "input_vals_vortex4.jl"
include(STARTUP_PATH)

cd("../m2")
ARGS[1] = "input_vals_vortex3.jl"
include(STARTUP_PATH)
ARGS[1] = "input_vals_vortex4.jl"
include(STARTUP_PATH)

cd("..")
include("calc_line.jl")

slope = calc_line()
println("slope = ", slope)

data = readdlm("err_data.dat")
err_vals = data[:, 2]
#println("err_vals = ", err_vals)

slope_val = 3.84
slope_margin = 0.1

@fact slope => greater_than(slope_val - slope_margin)
@fact slope => less_than(slope_val + slope_margin)

err_val = 6.87e-5
slope_fac = 1.25
println("err_vals[1] = ", err_vals[1])
@fact err_vals[1] => greater_than(err_val/slope_fac)
@fact err_vals[1] => less_than(err_val*slope_fac)

end
