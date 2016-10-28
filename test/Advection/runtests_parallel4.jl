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

facts("----- Testing Parallel 4 -----") do

  start_dir = pwd()
  cd("./energy")
  include(joinpath(pwd(), "runtests_parallel.jl"))
  cd(start_dir)
 

end

if MPI.Initialized()
  MPI.Finalize()
end
FactCheck.exitstatus()

