push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/advection"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))
include(joinpath(Pkg.dir("PDESolver"), "src/input/make_input.jl"))

using PDESolver
#using Base.Test
using FactCheck
using ODLCommonTools
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
using AdvectionEquationMod
using ForwardDiff
using NonlinearSolvers   # non-linear solvers
using ArrayViews

global const STARTUP_PATH = joinpath(Pkg.dir("PDESolver"), "src/solver/advection/startup_advection.jl")
# insert a command line argument
resize!(ARGS, 1)

facts("----- Testing Parallel -----") do

  start_dir = pwd()
  cd ("./rk4/parallel")
  ARGS[1] = "input_vals_parallel_runp.jl"
  include(STARTUP_PATH)

  datas = readdlm("../serial/error_calc.dat")
  datap = readdlm("error_calc.dat")

  @fact datas[1] --> roughly(datap[1], atol=1e-13)
  cd("../../")

  cd("./newton/parallel")
  ARGS[1] = "input_vals_parallel.jl"
  include(STARTUP_PATH)

  datas = readdlm("../serial/error_calc.dat")
  datap = readdlm("./error_calc.dat")
  @fact datap[1] --> roughly(datap[1], atol=1e-13)

  cd(start_dir)

end


