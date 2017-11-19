# Run advection tests

push!(LOAD_PATH, abspath(joinpath(pwd(), "..")))

using PDESolver
#using Base.Test
using FactCheck
using ODLCommonTools
using Utils
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
using AdvectionEquationMod
using ForwardDiff
using NonlinearSolvers   # non-linear solvers
using ArrayViews
import ODLCommonTools.sview
using Input
using LinearSolvers

function clean_dict(collection)
  for i in keys(collection)
    delete!(collection, i)
  end
end

#------------------------------------------------------------------------------
# define tests and tags
using TestSystem
# define tags that will be used
include("../tags.jl")

# test list
global const AdvectionTests = TestList()

include("test_frontend.jl")
include("test_linearsolver.jl")
include("test_lowlevel.jl")
include("test_eqn_deepcopy.jl")
include("test_3d.jl")
include("test_utils.jl")
include("test_gamma.jl")
include("test_mms.jl")
include("test_jac.jl")
include("test_GLS2.jl")
# test_newton.jl?
include("test_dg.jl")
include("test_staggered.jl")
include("test_adjoint.jl")
include("test_parallel.jl")
include( "./energy/runtests.jl")
#cd("./Nonlinearsolvers/")
include(joinpath(pwd(), "Nonlinearsolvers", "runtests_serial.jl"))
#cd("../")


#------------------------------------------------------------------------------
# run tests
facts("----- Running Advection tests -----") do
  nargs = length(ARGS)
  if nargs == 0
    tags = ASCIIString[TAG_DEFAULT]
  else
    tags = Array(ASCIIString, nargs)
    copy!(tags, ARGS)
  end

  resize!(ARGS, 1)
  ARGS[1] = ""
  run_testlist(AdvectionTests, run_advection, tags)
end


#------------------------------------------------------------------------------
# cleanup

# define global variable if needed
# this trick allows running the test files for multiple physics in the same
# session without finalizing MPI too soon
if !isdefined(:TestFinalizeMPI)
  TestFinalizeMPI = true
end

if MPI.Initialized() && TestFinalizeMPI
  MPI.Finalize()
end

FactCheck.exitstatus()
