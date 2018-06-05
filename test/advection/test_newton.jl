# Run advection tests

include(joinpath(Pkg.dir("PDESolver"), "src/input/make_input.jl"))

using PDESolver
#using Base.Test
using FactCheck
using ODLCommonTools
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
using AdvectionEquationMod
using NonlinearSolvers   # non-linear solvers
using ArrayViews

global const STARTUP_PATH = joinpath(Pkg.dir("PDESolver"), "src/solver/advection/startup.jl")
# insert a command line argument
resize!(ARGS, 1)
ARGS[1] = "input_vals_newton.jl"

@testset "----- Testing Newtons Method: Finite Difference -----" begin
  mesh, sbp, eqn, opts = solvePDE(ARGS[1])
  @test  calcNorm(eqn, eqn.res_vec, strongres=true)  < opts["res_abstol"]
end

@testset "----- Testing Newtons Method: Complex Step -----" begin
  arg_dict["jac_method"] = 2
  make_input(arg_dict, "input_vals_newton2")
  ARGS[1] = "input_vals_newton2.jl"
  mesh, sbp, eqn, opts = solvePDE(ARGS[1])
  @test  calcNorm(eqn, eqn.res_vec, strongres=true)  < opts["res_abstol"]
end

@testset "----- Testing Newtons Method: Petsc -----" begin
  arg_dict["jac_type"] = 3
  make_input(arg_dict, "input_vals_newton3")
  ARGS[1] = "input_vals_newton3.jl"
  mesh, sbp, eqn, opts = solvePDE(ARGS[1])
  @test  calcNorm(eqn, eqn.res_vec, strongres=true)  < opts["res_abstol"]
end
