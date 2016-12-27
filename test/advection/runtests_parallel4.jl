# run 4 processor tests

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
using Utils

global const STARTUP_PATH = joinpath(Pkg.dir("PDESolver"), "src/solver/advection/startup.jl")

#------------------------------------------------------------------------------
# define tests and tags

include("../TestSystem.jl")
# define tags that will be used

# test list
global const AdvectionTests = TestList()

include("Nonlinearsolvers/crank_nicolson_PETSc_parallel/runtests.jl")

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
add_func1!(AdvectionTests, runtests_parallel4)

#------------------------------------------------------------------------------
# run tests
facts("----- Running Advection 4 processor tests -----") do
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

