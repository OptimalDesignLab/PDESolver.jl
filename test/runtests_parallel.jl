# run tests in parallel

push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/euler"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))


using PDESolver
#using Base.Test
using FactCheck
using ODLCommonTools
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
using EulerEquationMod
using ForwardDiff
using NonlinearSolvers   # non-linear solvers
using ArrayViews
include( joinpath(Pkg.dir("PDESolver"), "src/solver/euler/complexify.jl"))
include( joinpath(Pkg.dir("PDESolver"), "src/input/make_input.jl"))
global const STARTUP_PATH = joinpath(Pkg.dir("PDESolver"), "src/solver/euler/startup.jl")


start_dir = pwd()
resize!(ARGS, 1)
ARGS[1] = "input_vals_parallel.jl"

facts("----- Testing Parallel -----") do

  start_dir = pwd()
  cd ("./rk4/parallel")
  ARGS[1] = "input_vals_parallel.jl"
  include(STARTUP_PATH)

  datas = readdlm("../serial/error_calc.dat")
  datap = readdlm("error_calc.dat")

  @fact datas[1] --> roughly(datap[1], atol=1e-13)
  @fact datas[2] --> roughly(datap[2], atol=1e-13)
  cd("../../")
#=
  cd("./newton/parallel")
  ARGS[1] = "input_vals_parallel.jl"
  include(STARTUP_PATH)

  datas = readdlm("../serial/error_calc.dat")
  datap = readdlm("./error_calc.dat")
  @fact datap[1] --> roughly(datap[1], atol=1e-13)

  cd(start_dir)
=#
end


