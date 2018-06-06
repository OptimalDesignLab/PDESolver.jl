# Run advection tests

push!(LOAD_PATH, abspath(joinpath(pwd(), "..")))

using MPI
using PDESolver
#using Base.Test
using Base.Test
using ODLCommonTools
using Utils
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
using AdvectionEquationMod
using NonlinearSolvers   # non-linear solversa
using OptimizationInterface
using ArrayViews
import ODLCommonTools.sview
import ArrayViews.view
using Input
using LinearSolvers
using PETSc2

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
include("test_parallel_serialpart.jl")
include( "./energy/runtests.jl")
#cd("./Nonlinearsolvers/")
include(joinpath(pwd(), "Nonlinearsolvers", "runtests_serial.jl"))
#cd("../")


#------------------------------------------------------------------------------
# run tests
@testset "----- Running Advection tests -----" begin
  runTestSystem(AdvectionTests, solvePDE, ARGS)
end

println("finished running tests")

#------------------------------------------------------------------------------
# cleanup
