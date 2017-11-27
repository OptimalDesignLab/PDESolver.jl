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
using ForwardDiff
using NonlinearSolvers   # non-linear solvers
using ArrayViews
import MPI
using Input

#------------------------------------------------------------------------------
# define test list
#include("../TestSystem.jl")
using TestSystem
include("../tags.jl")

global const EulerTests = TestList()
# define global const tags here


include("test_ESS_parallel.jl")

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

  tags = ASCIIString[TAG_DEFAULT]
  resize!(ARGS, 1)
  ARGS[1] = ""
  run_testlist(EulerTests, run_euler, tags)
end

#------------------------------------------------------------------------------
# cleanup

FactCheck.exitstatus()


