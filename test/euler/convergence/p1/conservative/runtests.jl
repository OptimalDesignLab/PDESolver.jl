function run_convergence_p1_conservative()
  @testset "---- P1 Conservative Convergence Tests -----" begin
    start_dir = pwd()
    println("pwd = ", pwd())

    resize!(ARGS, 1)

    cd("./m1")
    ARGS[1] = "input_vals_vortex3.jl"
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])
    ARGS[1] = "input_vals_vortex4.jl"
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])

    cd("../m2")
    ARGS[1] = "input_vals_vortex3.jl"
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])
    ARGS[1] = "input_vals_vortex4.jl"
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])

    cd("..")
    include("calc_line.jl")

    slope = calc_line()
    println("slope = ", slope)

    data = readdlm("err_data.dat")
    err_vals = data[:, 2]
    #println("err_vals = ", err_vals)

    slope_val = 1.80
    slope_margin = 0.1

    @test  slope  > slope_val - slope_margin
    @test  slope  < slope_val + slope_margin

    err_val = 0.0522
    slope_fac = 1.25
    println("err_vals[1] = ", err_vals[1])
    @test  err_vals[1]  > err_val/slope_fac
    @test  err_vals[1]  < err_val*slope_fac

  end  # end facts block

  return nothing
end
