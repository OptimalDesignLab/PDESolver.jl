using FactCheck
using ODLCommonTools
import ODLCommonTools.sview
using SummationByParts  # SBP operators

push!(LOAD_PATH, abspath(joinpath(pwd(), "..")))

using PDESolver
#using Base.Test
using PdePumiInterface  # common mesh interface - pumi
using EulerEquationMod
using ForwardDiff
using NonlinearSolvers   # non-linear solvers
using ArrayViews
using EulerEquationMod
using Utils
using MPI
using Input

if !MPI.Initialized()
  MPI.Init()
end

#------------------------------------------------------------------------------
# define tests and tags

using TestSystem
# define tags that will be used
include("../tags.jl")

# test list
global const EulerTests = TestList()


include("test_empty.jl")
include("test_input.jl")
include("test_complexify.jl")
include("test_lowlevel.jl")

include("test_dg.jl")
#include("test_simplemesh.jl")
include("test_GLS3.jl")
# TODO: uncomment when SBP is fixed
#include("test_modes.jl")
include("test_3d.jl")

include("test_adjoint.jl")
include("test_flux.jl")
include("test_ESS.jl")
include("test_rk4.jl")
include(joinpath("./convergence/runtests.jl"))
include("Utils.jl")
include("test_parallel.jl")

#------------------------------------------------------------------------------
# run tests
facts("----- Running Euler tests -----") do
  nargs = length(ARGS)
  if nargs == 0
    tags = ASCIIString[TAG_DEFAULT]
  else
    tags = Array(ASCIIString, nargs)
    copy!(tags, ARGS)
  end

  resize!(ARGS, 1)
  ARGS[1] = ""
  run_testlist(EulerTests, run_euler, tags)
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

