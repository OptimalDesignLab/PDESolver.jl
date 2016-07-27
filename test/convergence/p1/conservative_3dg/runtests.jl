global const STARTUP_PATH = joinpath(Pkg.dir("PDESolver"), "src/solver/euler/startup.jl")
using FactCheck

#run(`cptest.sh`)
#run(`cperr.sh`)
facts("---- P1 Conservative DG Convergence Tests -----") do
start_dir = pwd()

resize!(ARGS, 1)

cd("./m1")
ARGS[1] = "input_vals_vortex3.jl"
include(STARTUP_PATH)

cd("../m2")
ARGS[1] = "input_vals_vortex3.jl"
include(STARTUP_PATH)

cd("..")
include("calc_line.jl")

slope = calc_line()
println("slope = ", slope)

data = readdlm("err_data.dat")
err_vals = data[:, 2]
#println("err_vals = ", err_vals)

slope_val = 1.75
slope_margin = 0.1

@fact slope --> greater_than(slope_val - slope_margin)
@fact slope --> less_than(slope_val + slope_margin)

end
