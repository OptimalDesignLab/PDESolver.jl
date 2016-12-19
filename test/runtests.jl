using FactCheck
using ArrayViews
using ODLCommonTools
import ODLCommonTools.sview
using SummationByParts  # SBP operators
include( joinpath(Pkg.dir("PDESolver"), "src/solver/euler/complexify.jl"))
include( joinpath(Pkg.dir("PDESolver"), "src/input/make_input.jl"))

push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/euler"))

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

if !MPI.Initialized()
  MPI.Init()
end

global const STARTUP_PATH = joinpath(Pkg.dir("PDESolver"), "src/solver/euler/startup.jl")


#------------------------------------------------------------------------------
# define tests and tags

include("./TestSystem.jl")
# define tags that will be used
global const TAG_COMPLEX = "tag_complex"
global const TAG_BC = "tag_bc"
global const TAG_FLUX = "tag_flux"
global const TAG_ENTROPYVARS = "tag_entropyvars"
global const TAG_VOLUMEINTEGRALS = "tag_volumeintegral"
global const TAG_CONVERGENCE = "tag_convergence"
global const TAG_NLSOVLERS = "tag_nlsolvers"
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
  run_testlist(EulerTests, tags)
end


#------------------------------------------------------------------------------
# cleanup
if MPI.Initialized()
  MPI.Finalize()
end
FactCheck.exitstatus()

