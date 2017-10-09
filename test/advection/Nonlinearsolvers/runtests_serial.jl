include(joinpath(dirname(@__FILE__), "crank_nicolson_CS", "calc_line.jl"))

function test_CN()
  facts("---- Crank-Nicolson Convergence Tests, Complex Step Jacobian -----") do
    start_dir = pwd()

    cd(dirname(@__FILE__))
    resize!(ARGS, 1)

    # =================== CN, CS tests ===================
    cd("./crank_nicolson_CS/")
    println("======",pwd())

    cd("./m1")
    println("======", pwd())
    ARGS[1] = "input_vals1.jl"
    mesh, sbp, eqn, opts = run_advection(ARGS[1])

    cd("../m2")
    println("======", pwd())
    ARGS[1] = "input_vals1.jl"
    mesh, sbp, eqn, opts = run_advection(ARGS[1])

#    cd("../m3")
#    println("======", pwd())
#    ARGS[1] = "input_vals1.jl"
#    mesh, sbp, eqn, opts = run_advection(ARGS[1])

    cd("..")
    println("======", pwd())
    ARGS[1] = "calc_line.jl"  #???
#    mesh, sbp, eqn, opts = run_advection(ARGS[1])

    slope = calc_line()
    # println("slope = ", slope)

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

    cd("../")
    # =================== CN, FD tests ===================
    cd("./crank_nicolson_FD/")

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
    ARGS[1] = "calc_line.jl"
#    mesh, sbp, eqn, opts = run_advection(ARGS[1])  #???

    slope = calc_line()
    # println("slope = ", slope)

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

    cd("../")
    # =================== CN, PETSc CS tests =================== 
    cd("./crank_nicolson_PETSc_serial/")

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
    ARGS[1] = "calc_line.jl"
#    mesh, sbp, eqn, opts = run_advection(ARGS[1]) #???

    slope = calc_line()
    # println("slope = ", slope)

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

    cd(start_dir)
  end  # end facts block

end

add_func1!(AdvectionTests, test_CN, [TAG_SHORTTEST, TAG_CN])
