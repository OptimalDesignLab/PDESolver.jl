# run tests in parallel with 2 processes

push!(LOAD_PATH, abspath(joinpath(pwd(), "..")))

using PDESolver
#using Base.Test
using Base.Test
using ODLCommonTools
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
using Utils
using EulerEquationMod
using LinearSolvers
using NonlinearSolvers   # non-linear solvers
using OptimizationInterface
using ArrayViews
import ArrayViews.view
import MPI
using Input

#------------------------------------------------------------------------------
# define test list
using TestSystem
include("../tags.jl")

global const EulerTests = TestList()
# define global const tags here

# NOTE: For adding parallel tests to TICON, Write them in a separate file and
# include them below `global const EulerTests = TestList()`.

include("test_parallel_derivatives.jl")

"""
  Run the parallel tests and compare against serial results run as part of
  the serial tests
"""
function test_parallel2()
  @testset "----- Testing Parallel -----" begin

    start_dir = pwd()

    # test rk4
    cd("./rk4/parallel")
    fname = "input_vals_parallel.jl"
    mesh, sbp, eqn, opts = solvePDE(fname)

    datas = readdlm("../serial/error_calc.dat")
    datap = readdlm("error_calc.dat")

    @test isapprox( datas[1], datap[1]) atol=1e-13
    @test isapprox( datas[2], datap[2]) atol=1e-13

    # test staggered_parallel
    cd("../staggered_parallel")
    fname = "input_vals_parallel.jl"
    mesh, sbp, eqn, opts = solvePDE(fname)

    datas = readdlm("../staggered_serial/error_calc.dat")
    datap = readdlm("error_calc.dat")

    @test isapprox( datas[1], datap[1]) atol=1e-13
    @test isapprox( datas[2], datap[2]) atol=1e-13


    cd("../../")

    cd("./lserk/parallel")
    fname = "input_vals_parallel.jl"
    mesh, sbp, eqn, opts = solvePDE(fname)

    datas = readdlm("../serial/error_calc.dat")
    datap = readdlm("error_calc.dat")

    @test isapprox( datas[1], datap[1]) atol=1e-13
    @test isapprox( datas[2], datap[2]) atol=1e-13

    cd("../../")

    # test newton
    cd("./newton/parallel")
    fname = "input_vals_parallel.jl"
    mesh, sbp, eqn, opts = solvePDE(fname)

    datas = readdlm("../serial/error_calc.dat")
    datap = readdlm("./error_calc.dat")
    @test isapprox( datas[1], datap[1]) atol=1e-13

    fill!(eqn.q, 1.0)
    fill!(eqn.res, 0.0)

    for i=1:mesh.npeers
      fill!(eqn.shared_data[i].q_send, 1.0)
      fill!(eqn.shared_data[i].q_recv, 1.0)
    end
    evalHomotopy(mesh, sbp, eqn, opts, eqn.res)

    @test isapprox( vecnorm(eqn.res), 0.0) atol=1e-13

    cd(start_dir)

  end

  return nothing
end

function test_parallel_nopre()

  @testset "----- Testing parallel nopre -----" begin
    start_dir = pwd()

    # test rk4
    cd("./rk4/parallel")
    fname = "input_vals_parallel.jl"

    mesh, sbp, eqn, opts = solvePDE(fname)

    EulerEquationMod.evalResidual(mesh, sbp, eqn, opts)
    res_orig = copy(eqn.res)

    opts["precompute_face_flux"] = false
    fill!(eqn.res, 0.0)

    EulerEquationMod.evalResidual(mesh, sbp, eqn, opts)
    
    @test isapprox( norm(vec(eqn.res - res_orig)), 0.0) atol=1e-13

    cd(start_dir)

    # test_newton
    cd("./newton/parallel")
    fname = "input_vals_parallel.jl"
    mesh, sbp, eqn, opts = solvePDE(fname)

    fill!(eqn.res, 0.0)
    EulerEquationMod.evalResidual(mesh, sbp, eqn, opts)
    res_orig = copy(eqn.res)

    opts["precompute_face_flux"] = false
    fill!(eqn.res, 0.0)

    EulerEquationMod.evalResidual(mesh, sbp, eqn, opts)
    
    @test isapprox( norm(vec(eqn.res - res_orig)), 0.0) atol=1e-13

    cd(start_dir)

  end  # end facts block

  return nothing
end

function test_restart()
  @testset "----- Testing restart -----" begin
    start_dir = pwd()

    # test rk4
    cd("./rk4/parallel")
    fname = "input_vals_restart"
    mesh, sbp, eqn, opts = solvePDE(fname)

    datas = readdlm("../serial/error_calc.dat")
    datap = readdlm("error_calc.dat")

    @test isapprox( datas[1], datap[1]) atol=1e-13
    @test isapprox( datas[2], datap[2]) atol=1e-13

    cd("../../lserk/parallel")
    fname = "input_vals_restart"
    mesh, sbp, eqn, opts = solvePDE(fname)

    datas = readdlm("../serial/error_calc.dat")
    datap = readdlm("error_calc.dat")

    @test isapprox( datas[1], datap[1]) atol=1e-13
    @test isapprox( datas[2], datap[2]) atol=1e-13


    cd(start_dir)
  end

  return nothing
end


#test_parallel2()
add_func1!(EulerTests, test_parallel2, [TAG_SHORTTEST])
add_func1!(EulerTests, test_parallel_nopre, [TAG_SHORTTEST])
add_func1!(EulerTests, test_restart, [TAG_SHORTTEST])

#------------------------------------------------------------------------------
# run tests
@testset "----- Running Euler 2 process tests -----" begin
  runTestSystem(EulerTests, solvePDE, ARGS)
end
