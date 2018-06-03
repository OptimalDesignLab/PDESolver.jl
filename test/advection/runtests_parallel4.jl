# run 4 processor tests

push!(LOAD_PATH, abspath(joinpath(pwd(), "..")))

using PDESolver
#using Base.Test
using FactCheck
using ODLCommonTools
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
using AdvectionEquationMod
using LinearSolvers
using NonlinearSolvers   # non-linear solvers
using OptimizationInterface
using ArrayViews
import ArrayViews.view
using Utils
using Input
using PETSc2


#------------------------------------------------------------------------------
# define tests and tags

#include("../TestSystem.jl")
using TestSystem
# define tags that will be used
include("../tags.jl")

# test list
global const AdvectionTests = TestList()

include("Nonlinearsolvers/crank_nicolson_PETSc_parallel/runtests.jl")
include("test_parallel2.jl")
"""
  Test energy stability in parallel
"""
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
add_func1!(AdvectionTests, runtests_parallel4, [TAG_SHORTTEST])

#------------------------------------------------------------------------------
# run tests
facts("----- Running Advection 4 processor tests -----") do
  nargs = length(ARGS)
  if nargs == 0
    tags = String[TAG_DEFAULT]
  else
    tags = Array{String}(nargs)
    copy!(tags, ARGS)
  end

  resize!(ARGS, 1)
  ARGS[1] = ""
  run_testlist(AdvectionTests, solvePDE, tags)
end

#------------------------------------------------------------------------------
# cleanup

FactCheck.exitstatus()

