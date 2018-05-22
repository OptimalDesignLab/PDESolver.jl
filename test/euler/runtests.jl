using FactCheck
using ODLCommonTools
import ODLCommonTools.sview
using SummationByParts  # SBP operators

push!(LOAD_PATH, abspath(joinpath(pwd(), "..")))

using PDESolver
#using Base.Test
using PdePumiInterface  # common mesh interface - pumi
using EulerEquationMod
#using ForwardDiff
using LinearSolvers
using NonlinearSolvers   # non-linear solvers
using OptimizationInterface
using ArrayViews
using EulerEquationMod
using Utils
using MPI
using Input
using PETSc2

#------------------------------------------------------------------------------
# define tests and tags

using TestSystem
# define tags that will be used
include("../tags.jl")

# test list
global const EulerTests = TestList()


include("test_eqn_deepcopy.jl")     # note: eqn gets written random values to it, so anything that 
include("test_empty.jl")
include("test_input.jl")
include("test_complexify.jl")
include("test_lowlevel.jl")

include("test_dg.jl")
#include("test_simplemesh.jl")
include("test_GLS3.jl")
include("test_modes.jl")
include("test_3d.jl")
include("test_bc.jl")
include("test_interp.jl")
include("test_jac.jl")
include("test_adjoint.jl")
include("test_reversemode.jl")
include("test_flux.jl")
include("test_ESS.jl")
include("test_curvilinear.jl")
include("test_rk4.jl")
include(joinpath("./convergence/runtests.jl"))
include("Utils.jl")
include("test_parallel.jl")
include("test_homotopy.jl")
include("test_staggered.jl")
include("test_checkpoint.jl")

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
  run_testlist(EulerTests, solvePDE, tags)
end


#------------------------------------------------------------------------------
# cleanup

FactCheck.exitstatus()
