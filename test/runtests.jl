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


# insert a command line argument
resize!(ARGS, 1)
ARGS[1] = "input_vals_channel.jl"
include("test_empty.jl")
include("test_lowlevel.jl")
#include("test_simplemesh.jl")
include("test_modes.jl")
FactCheck.exitstatus()


# write your own tests here
# @test 1 == 1


# using SummationByParts
# #using Base.Test
# using FactCheck
# 
# include("test_orthopoly.jl")
# include("test_symcubatures.jl")
# include("test_cubature.jl")
# include("test_SummationByParts.jl")
# 
# FactCheck.exitstatus()
