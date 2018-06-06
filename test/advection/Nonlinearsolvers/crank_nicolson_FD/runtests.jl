@testset "---- Crank-Nicolson Convergence Tests, Finite Difference Jacobian -----" begin
start_dir = pwd()

cd("./m1")
fname = "input_vals1.jl"
mesh, sbp, eqn, opts = solvePDE(fname)

cd("../m2")
fname = "input_vals1.jl"
mesh, sbp, eqn, opts = solvePDE(fname)

cd("../m3")
fname = "input_vals1.jl"
mesh, sbp, eqn, opts = solvePDE(fname)

cd("..")
include("calc_line.jl")

slope = calc_line()
# println("slope = ", slope)

data = readdlm("err_data.dat")
err_vals = data[:, 2]
#println("err_vals = ", err_vals)

slope_val = 2.00
slope_margin = 0.1

@test  slope  > slope_val - slope_margin
@test  slope  < slope_val + slope_margin

err_val = 0.09095728504176116 
slope_fac = 1.25
# println("err_vals[1] = ", err_vals[1])
@test  err_vals[1]  > err_val/slope_fac
@test  err_vals[1]  < err_val*slope_fac

end
