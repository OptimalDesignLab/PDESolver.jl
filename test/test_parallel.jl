# run tests related to parallel functionality that are not actually parallel
#=
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
=#

resize!(ARGS, 1)
ARGS[1] = "input_vals_vortex3.jl"

cd("./rk4/serial")
include(STARTUP_PATH)
cd("../parallel")

arg_dict["smb_name"] = "src/mesh_files/psquare2.smb"
make_input(arg_dict, "input_vals_parallel")

