# Run advection tests

push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/advection"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))
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

function clean_dict(collection)
  for i in keys(collection)
    delete!(collection, i)
  end
end

global const STARTUP_PATH = joinpath(Pkg.dir("PDESolver"), "src/solver/advection/startup_advection.jl")

#------------------------------------------------------------------------------
# define tests and tags

include("../TestSystem.jl")
# define tags that will be used
global const TAG_COMPLEX = "tag_complex"
global const TAG_BC = "tag_bc"
global const TAG_FLUX = "tag_flux"
global const TAG_VOLUMEINTEGRALS = "tag_volumeintegral"
global const TAG_CONVERGENCE = "tag_convergence"

# test list
global const AdvectionTests = TestList()

include("test_lowlevel.jl")
include("test_3d.jl")
include("test_gamma.jl")
include("test_mms.jl")
include("test_jac.jl")
include("test_GLS2.jl")
include("test_dg.jl")
include("test_functional_integrate.jl")
include("test_parallel.jl")
include( "./energy/runtests.jl")

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
  run_testlist(AdvectionTests, tags)
end



if MPI.Initialized()
  MPI.Finalize()
end

FactCheck.exitstatus()


# write your own tests here
# @test 1 == 1


# using SummationByParts
# #using Base.Test
# using FactCheck
# 
# include("test_orthopoly.jl")
# include("test_symcubatures.jl")
# include("test_cubature.jl")
# include("test_SummationByParts.jl")
# 
# FactCheck.exitstatus()
