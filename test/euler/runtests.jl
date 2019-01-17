using Base.Test
using ODLCommonTools
import ODLCommonTools.sview
using SummationByParts  # SBP operators

push!(LOAD_PATH, abspath(joinpath(pwd(), "..")))  # get TestSystem
push!(LOAD_PATH, abspath(joinpath(pwd(), "../common")))  # get TestCommon

using PDESolver
#using Base.Test
using TestCommon
using PdePumiInterface  # common mesh interface - pumi
using EulerEquationMod
using LinearSolvers
using NonlinearSolvers   # non-linear solvers
using OptimizationInterface
using ArrayViews
import ArrayViews.view
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
include("test_functionals.jl")
include("test_adjoint.jl")
include("test_flux.jl")
include("test_ESS.jl")
include("test_curvilinear.jl")
include("test_rk4.jl")
include(joinpath("./convergence/runtests.jl"))
include("Utils.jl")
include("test_parallel_serialpart.jl")
include("test_homotopy.jl")
include("test_staggered.jl")
include("test_checkpoint.jl")
include("test_coordsfd.jl")
include("test_shock_capturing.jl")

#------------------------------------------------------------------------------
# run tests
@testset "----- Running Euler tests -----" begin

  runTestSystem(EulerTests, solvePDE, ARGS)
end


#------------------------------------------------------------------------------
# cleanup
