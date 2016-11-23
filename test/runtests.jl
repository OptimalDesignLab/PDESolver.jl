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
# insert a command line argument
resize!(ARGS, 1)
ARGS[1] = "input_vals_channel.jl"
include("test_empty.jl")
include("test_input.jl")
include("test_complexify.jl")
include("test_lowlevel.jl")
include("test_dg.jl")
#include("test_simplemesh.jl")
include("test_GLS3.jl")
# TODO: uncomment when SBP is fixed
#include("test_modes.jl")
include("test_3d.jl")
include("test_adjoint.jl")
include("test_flux.jl")
include("test_ESS.jl")
include("test_rk4.jl")


cd("./convergence")
include(joinpath(pwd(), "runtests.jl"))
cd("..")

include("Utils.jl")
include("test_parallel.jl")

if MPI.Initialized()
  MPI.Finalize()
end
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
