# Run advection tests

include(joinpath(Pkg.dir("PDESolver"), "src/input/make_input.jl"))

using PDESolver
#using Base.Test
using FactCheck
using ODLCommonTools
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
using AdvectionEquationMod
#using ForwardDiff
using NonlinearSolvers   # non-linear solvers
using ArrayViews

global const STARTUP_PATH = joinpath(Pkg.dir("PDESolver"), "src/solver/advection/startup.jl")
# insert a command line argument
resize!(ARGS, 1)
ARGS[1] = "input_vals_newton.jl"

facts("----- Testing Newtons Method: Finite Difference -----") do
  mesh, sbp, eqn, opts = solvePDE(ARGS[1])
  @fact calcNorm(eqn, eqn.res_vec, strongres=true) --> less_than(opts["res_abstol"])
end

facts("----- Testing Newtons Method: Complex Step -----") do
  arg_dict["jac_method"] = 2
  make_input(arg_dict, "input_vals_newton2")
  ARGS[1] = "input_vals_newton2.jl"
  mesh, sbp, eqn, opts = solvePDE(ARGS[1])
  @fact calcNorm(eqn, eqn.res_vec, strongres=true) --> less_than(opts["res_abstol"])
end

facts("----- Testing Newtons Method: Petsc -----") do
  arg_dict["jac_type"] = 3
  make_input(arg_dict, "input_vals_newton3")
  ARGS[1] = "input_vals_newton3.jl"
  mesh, sbp, eqn, opts = solvePDE(ARGS[1])
  @fact calcNorm(eqn, eqn.res_vec, strongres=true) --> less_than(opts["res_abstol"])
end
