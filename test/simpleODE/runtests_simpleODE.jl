# Run advection tests

push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/simpleODE"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))
include(joinpath(Pkg.dir("PDESolver"), "src/input/make_input.jl"))

using PDESolver
#using Base.Test
using FactCheck
using ODLCommonTools
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
using SimpleODEMod
using ForwardDiff
using NonlinearSolvers   # non-linear solvers
using ArrayViews

function clean_dict(collection)
  for i in keys(collection)
    delete!(collection, i)
  end
end

global const STARTUP_PATH = joinpath(Pkg.dir("PDESolver"), "src/solver/simpleODE/startup_simpleODE.jl")
# insert a command line argument
resize!(ARGS, 1)

facts("---- testing SimpleODE ----") do
  cd("./eqn4/")
  ARGS[1] = "input_vals_simpleODE.jl"
  include(STARTUP_PATH)

  for i = 1:length(eqn.q_vec)
    @fact eqn.q_vec[i] --> roughly(4.0, atol=1e-10)
  end
end

#=
include("test_empty.jl")
#include("test_input.jl")
include("test_lowlevel.jl")
include("test_3d.jl")
include("test_gamma.jl")
include("test_mms.jl")
include("test_jac.jl")
include("test_GLS2.jl")
include("test_dg.jl")
include("test_functional_integrate.jl")
include("test_parallel.jl")

start_dir = pwd()
cd("./energy")
include( joinpath(pwd(), "runtests.jl"))
cd(start_dir)

cd("./Nonlinearsolvers/")
include(joinpath(pwd(), "runtests_serial.jl"))
cd("../")
=#

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