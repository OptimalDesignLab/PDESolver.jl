function test_CN_parallel()
    
  start_dir = pwd()
  facts("---- Crank-Nicolson Convergence Tests, PETSc + CS Jacobian -----") do

    cd(dirname(@__FILE__))
    cd("./m1")
    ARGS[1] = "input_vals1.jl"
    mesh, sbp, eqn, opts = run_advection(ARGS[1])

    cd("../m2")
    ARGS[1] = "input_vals1.jl"
    mesh, sbp, eqn, opts = run_advection(ARGS[1])

#    cd("../m3")
#    ARGS[1] = "input_vals1.jl"
#    mesh, sbp, eqn, opts = run_advection(ARGS[1])

    cd("..")
    include(joinpath(pwd(), "calc_line.jl"))

    MPI.Barrier(mesh.comm)
    slope = calc_line()
    # println("slope = ", slope)

    MPI.Barrier(mesh.comm)
    
    data = readdlm("err_data.dat")
    err_vals = data[:, 2]
    #println("err_vals = ", err_vals)

    slope_val = 2.00
    slope_margin = 0.1

    @fact slope --> greater_than(slope_val - slope_margin)
    @fact slope --> less_than(slope_val + slope_margin)

    err_val = 0.09095728504176116 
    slope_fac = 1.25
    # println("err_vals[1] = ", err_vals[1])
    @fact err_vals[1] --> greater_than(err_val/slope_fac)
    @fact err_vals[1] --> less_than(err_val*slope_fac)

  end

  facts("---- Testing restart -----") do
    cd("./m1")
    mesh, sbp, eqn, opts = run_advection("input_vals_restart")
    MPI.Barrier(mesh.comm)
    data = readdlm("error_calc.dat")
    datas = readdlm("../../crank_nicolson_PETSc_serial/m1/error_calc.dat")

    @fact data[1] --> roughly(datas[1], atol=1e-13)
    @fact data[2] --> roughly(datas[2], atol=1e-13)

    cd("..")


  end

  cd(start_dir)

  return nothing
end

add_func1!(AdvectionTests, test_CN_parallel, [TAG_SHORTTEST])
