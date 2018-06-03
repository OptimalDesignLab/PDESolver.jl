using Base.Test
using ODLCommonTools
import ODLCommonTools.sview
using SummationByParts  # SBP operators

push!(LOAD_PATH, abspath(joinpath(pwd(), "..")))

using PDESolver
#using Base.Test
using PdePumiInterface  # common mesh interface - pumi
using EllipticEquationMod
using LinearSolvers
using NonlinearSolvers   # non-linear solvers
using OptimizationInterface
using ArrayViews
import ArrayViews.view
using EllipticEquationMod
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
global const EllipticTests = TestList()

include("test_conv_rate.jl")

#------------------------------------------------------------------------------
# run tests
@testset "----- Running Elliptic tests -----" begin
  nargs = length(ARGS)
  if nargs == 0
    tags = String[TAG_DEFAULT]
  else
    tags = Array{String}(nargs)
    copy!(tags, ARGS)
  end

  resize!(ARGS, 1)
  ARGS[1] = ""
  run_testlist(EllipticTests, solvePDE, tags)
end


#------------------------------------------------------------------------------
# cleanup
