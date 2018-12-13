include("calc_line.jl")

function test_CN_parallel()
    
  start_dir = pwd()
  @testset "---- Crank-Nicolson Convergence Tests, PETSc + CS Jacobian -----" begin

    cd(dirname(@__FILE__))
    cd("./m1")
    fname = "input_vals1.jl"
    mesh, sbp, eqn, opts = solvePDE(fname)

    cd("../m2")
    fname = "input_vals1.jl"
    mesh, sbp, eqn, opts = solvePDE(fname)

#    cd("../m3")
#    fname = "input_vals1.jl"
#    mesh, sbp, eqn, opts = solvePDE(fname)

    cd("..")
#    include(joinpath(pwd(), "calc_line.jl"))

    MPI.Barrier(mesh.comm)
    slope = calc_line()
    # println("slope = ", slope)

    MPI.Barrier(mesh.comm)
    
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

  @testset "---- Testing restart -----" begin
    cd("./m1")
    mesh, sbp, eqn, opts = solvePDE("input_vals_restart")
    MPI.Barrier(mesh.comm)
    data = readdlm("error_calc.dat")
    datas = readdlm("../../crank_nicolson_PETSc_serial/m1/error_calc.dat")

    @test isapprox( data[1], datas[1]) atol=1e-13
    @test isapprox( data[2], datas[2]) atol=1e-13

    cd("..")


  end

  @testset "--- Testing matrix-free Petsc -----" begin
    cd("./m1")
    opts = read_input("input_vals1.jl")
    
    # clear out the old options
#    PetscClearOptions(opts["petsc_options"])
    opts["petsc_options"] = Dict{AbstractString, AbstractString}(
                                "-pc_type" => "none",
                                "-malloc" => "",
                                "-malloc_debug" => "",
                                "-ksp_monitor" => "",
                                "-ksp_gmres_modifiedgramschmidt" => "",
                                "-ksp_gmres_restart" => "30",
                                )

    mesh, sbp, eqn, opts = solvePDE("input_vals_restart")
    MPI.Barrier(mesh.comm)
    data = readdlm("error_calc.dat")
    datas = readdlm("../../crank_nicolson_PETSc_serial/m1/error_calc.dat")

    @test isapprox( data[1], datas[1]) atol=1e-13
    @test isapprox( data[2], datas[2]) atol=1e-13

    cd("..")


  end

  cd(start_dir)

  return nothing
end

add_func1!(AdvectionTests, test_CN_parallel, [TAG_SHORTTEST])
