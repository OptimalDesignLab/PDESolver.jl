# run 4 processor tests

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

#------------------------------------------------------------------------------
# define tests and tags

include("../TestSystem.jl")
# define tags that will be used

# test list
global const AdvectionTests = TestList()



function runtests_parallel4()
  facts("----- Testing Parallel 4 -----") do

    start_dir = pwd()
    cd("./energy")
    include(joinpath(pwd(), "runtests_parallel.jl"))
    cd(start_dir)
  end

  return nothing
end

#runtests_parallel4()
add_func1!(AdvectionTests, runtests_parallel4)

#------------------------------------------------------------------------------
# run tests
facts("----- Running Advection 2 processor tests -----") do
  nargs = length(ARGS)
  if nargs == 0
    tags = ASCIIString[TAG_DEFAULT]
  else
    tags = Array(ASCIIString, nargs)
    copy!(tags, ARGS)
  end

  resize!(ARGS, 1)
  ARGS[1] = ""
  run_testlist(AdvectionTests, tags)
end

#------------------------------------------------------------------------------
# cleanup


if MPI.Initialized()
  MPI.Finalize()
end
FactCheck.exitstatus()

