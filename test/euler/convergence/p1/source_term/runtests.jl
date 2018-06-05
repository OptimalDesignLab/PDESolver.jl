function test_convergence_p1_source()
  @testset "---- P1 Conservative DG source -----" begin
    start_dir = pwd()

    resize!(ARGS, 1)

    cd("./m1")
    ARGS[1] = "input_vals_vortex3.jl"
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])

    cd("../m2")
    ARGS[1] = "input_vals_vortex3.jl"
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])

    cd("..")

    slope = calc_line()
    println("slope = ", slope)

    data = readdlm("err_data.dat")
    err_vals = data[:, 2]
    #println("err_vals = ", err_vals)

    slope_val = 2.30
    slope_margin = 0.1

    @test  slope  > slope_val - slope_margin
    @test  slope  < slope_val + slope_margin

  end  # end facts block

  return nothing
end
