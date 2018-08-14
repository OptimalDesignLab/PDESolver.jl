function test_convergence_p3_entropy()
  @testset "---- P3 Entropy Convergence Tests -----" begin
  start_dir = pwd()

  cd("./m1")
  fname = "input_vals_vortex3.jl"
  mesh, sbp, eqn, opts = solvePDE(fname)
  fname = "input_vals_vortex4.jl"
  mesh, sbp, eqn, opts = solvePDE(fname)

  cd("../m2")
  fname = "input_vals_vortex3.jl"
  mesh, sbp, eqn, opts = solvePDE(fname)
  fname = "input_vals_vortex4.jl"
  mesh, sbp, eqn, opts = solvePDE(fname)

  cd("..")

  slope = calc_line()
  println("slope = ", slope)

  data = readdlm("err_data.dat")
  err_vals = data[:, 2]
  #println("err_vals = ", err_vals)

  slope_val = 3.64
  slope_margin = 0.1

  @test  slope  > slope_val - slope_margin
  @test  slope  < slope_val + slope_margin

  err_val = 5.76e-5
  slope_fac = 1.25
  println("err_vals[1] = ", err_vals[1])
  @test  err_vals[1]  > err_val/slope_fac
  @test  err_vals[1]  < err_val*slope_fac

  end

  return nothing
end
