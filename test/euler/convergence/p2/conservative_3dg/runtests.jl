function test_convergence_p2_3dg()
  @testset "---- P2 Conservative 3DG Convergence Tests -----" begin
    start_dir = pwd()

    cd("./m1")
    fname = "input_vals_vortex3.jl"
    mesh, sbp, eqn, opts = solvePDE(fname)

    cd("../m2")
    fname = "input_vals_vortex3.jl"
    mesh, sbp, eqn, opts = solvePDE(fname)

    cd("..")

    slope = calc_line2()
    println("slope = ", slope)

    data = readdlm("err_data.dat")
    err_vals = data[:, 2]
    #println("err_vals = ", err_vals)

    slope_val = 3.08
    slope_margin = 0.1

    @test  slope  > slope_val - slope_margin
    @test  slope  < slope_val + slope_margin

  end

  return nothing
end
