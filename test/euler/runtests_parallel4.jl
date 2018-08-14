# run tests in parallel

push!(LOAD_PATH, abspath(joinpath(pwd(), "..")))

using PDESolver
#using Base.Test
using Base.Test
using ODLCommonTools
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
using Utils
using EulerEquationMod
using LinearSolvers
using NonlinearSolvers   # non-linear solvers
using OptimizationInterface
using ArrayViews
import ArrayViews.view
import MPI
using Input
using PETSc2

#------------------------------------------------------------------------------
# define test list
#include("../TestSystem.jl")
using TestSystem
include("../tags.jl")

global const EulerTests = TestList()
# define global const tags here


include("test_ESS_parallel.jl")
include("test_jacp.jl")

#------------------------------------------------------------------------------
# run tests
@testset "----- Running Euler 4 process tests -----" begin
  runTestSystem(EulerTests, solvePDE, ARGS)
end

#------------------------------------------------------------------------------
# cleanup
