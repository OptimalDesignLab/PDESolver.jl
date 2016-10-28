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

facts("----- Testing Parallel4 -----") do

  include("test_ESS_parallel.jl")

end

if MPI.Initialized()
  MPI.Finalize()
end

FactCheck.exitstatus()


