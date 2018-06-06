# run 4 processor tests

push!(LOAD_PATH, abspath(joinpath(pwd(), "..")))

using PDESolver
#using Base.Test
using Base.Test
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
include("test_parallel4.jl")
"""
  Test energy stability in parallel
"""
function runtests_parallel4()
  @testset "----- Testing Parallel 4 -----" begin

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
@testset "----- Running Advection 4 processor tests -----" begin
  runTestSystem(AdvectionTests, solvePDE, ARGS)
end

#------------------------------------------------------------------------------
# cleanup

