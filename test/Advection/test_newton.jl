# Run advection tests

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
ARGS[1] = "input_vals_newton.jl"

facts("----- Testing Newtons Method: Finite Difference -----") do
  include(STARTUP_PATH)
  @fact calcNorm(eqn, eqn.res_vec, strongres=true) --> less_than(opts["res_abstol"])
end

facts("----- Testing Newtons Method: Complex Step -----") do
  arg_dict["run_type"] = 5
  make_input(arg_dict, "input_vals_newton2")
  ARGS[1] = "input_vals_newton2.jl"
  include(STARTUP_PATH)
  @fact calcNorm(eqn, eqn.res_vec, strongres=true) --> less_than(opts["res_abstol"])
end

facts("----- Testing Newtons Method: Petsc -----") do
  arg_dict["jac_type"] = 3
  make_input(arg_dict, "input_vals_newton3")
  ARGS[1] = "input_vals_newton3.jl"
  include(STARTUP_PATH)
  @fact calcNorm(eqn, eqn.res_vec, strongres=true) --> less_than(opts["res_abstol"])
end

# run a serial case to compare parallel against later
cd ("./newton/serial")
arg_dict["smb_name"] = "src/mesh_files/serial2.smb"
arg_dict["dmg_name"] = "src/mesh_files/serial2.dmg"
make_input(arg_dict, "input_vals_serial")
ARGS[1] = "input_vals_serial.jl"
include(STARTUP_PATH)

# make parallel file
cd("../parallel")
arg_dict["smb_name"] = "src/mesh_files/psquare2.smb"
arg_dict["dmg_name"] = "src/mesh_files/psquare2.dmg"
arg_dict["krylov_abstol"] = 1e-12
make_input(arg_dict, "input_vals_parallel")
cd("../../")


