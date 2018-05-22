# run tests in parallel

push!(LOAD_PATH, abspath(joinpath(pwd(), "..")))

using PDESolver
#using Base.Test
using FactCheck
using ODLCommonTools
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
using Utils
using EulerEquationMod
using LinearSolvers
using NonlinearSolvers   # non-linear solvers
using OptimizationInterface
using ArrayViews
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
facts("----- Running Euler 4 process tests -----") do

  nargs = length(ARGS)
  if nargs == 0
    tags = ASCIIString[TAG_DEFAULT]
  else
    tags = Array(ASCIIString, nargs)
    copy!(tags, ARGS)
  end

  resize!(ARGS, 1)
  ARGS[1] = ""
  run_testlist(EulerTests, solvePDE, tags)
end

#------------------------------------------------------------------------------
# cleanup

FactCheck.exitstatus()


