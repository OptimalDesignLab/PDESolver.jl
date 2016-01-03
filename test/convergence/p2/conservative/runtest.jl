#global const STARTUP_PATH = joinpath(Pkg.dir("PDESolver"), "src/solver/euler/startup.jl")

#run(`cptest.sh`)
#run(`cperr.sh`)

start_dir = pwd()

resize!(ARGS, 1)

cd("./m4")
ARGS[1] = "input_vals_vortex3.jl"
include(STARTUP_PATH)
ARGS[1] = "input_vals_vortex4.jl"
include(STARTUP_PATH)

cd("../m5")
ARGS[1] = "input_vals_vortex3.jl"
include(STARTUP_PATH)
ARGS[1] = "input_vals_vortex4.jl"
include(STARTUP_PATH)

cd("..")
include("calc_line.jl")

slope = calc_line()
#println("slope = ", slope)

data = readdlm("err_data.dat")
err_vals = data[:, 2]
#println("err_vals = ", err_vals)

slope_val = 3.12
slope_margin = 0.1

@fact slope => greater_than(slope_val - slope_margin)
@fact slope => less_than(slope_val + slope_margin)

err_val = 0.000793
slope_fac = 1.25

@fact err_vals[1] => greater_than(err_val/slope_fac)
@fact err_vals[1] => less_than(err_vals*slope_fac)
