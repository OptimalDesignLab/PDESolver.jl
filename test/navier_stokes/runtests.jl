using PDESolver
using NavierStokesMod
using Base.Test
using ODLCommonTools
import ODLCommonTools.sview

push!(LOAD_PATH, abspath(joinpath(pwd(), "..")))  # get TestSystem

#------------------------------------------------------------------------------
# define tests and tags

using TestSystem
# define tags that will be used
include("../tags.jl")

# test list
global const NSTests = TestList()

# include test files
include("test_viscous.jl")


#------------------------------------------------------------------------------
# run tests
@testset "----- Running Euler tests -----" begin

  runTestSystem(NSTests, solvePDE, ARGS)
end


