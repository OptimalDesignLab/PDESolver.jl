# Run advection tests

push!(LOAD_PATH, abspath(joinpath(pwd(), "..")))

using PDESolver
#using Base.Test
using Base.Test
using ODLCommonTools
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
using Utils
using SimpleODEMod
using NonlinearSolvers   # non-linear solvers
using ArrayViews
import ArrayViews.view
using Input

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
global const SimpleODETests = TestList()


function test_eq4()
  @testset "---- testing SimpleODE ----" begin
    start_dir = pwd()
    cd("./eqn4/")
    fname = "input_vals_simpleODE.jl"
    mesh, sbp, eqn, opts = solvePDE(fname)

    for i = 1:length(eqn.q_vec)
      @test isapprox( eqn.q_vec[i], 4.0) atol=1e-10
    end

    cd(start_dir)
  end  # end facts block

  return nothing
end

add_func1!(SimpleODETests, test_eq4, [TAG_SHORTTEST])

#------------------------------------------------------------------------------
# run tests
@testset "----- Running SimpleODE tests -----" begin
  nargs = length(ARGS)
  if nargs == 0
    tags = String[TAG_DEFAULT]
  else
    tags = Array{String}(nargs)
    copy!(tags, ARGS)
  end

  resize!(ARGS, 1)
  ARGS[1] = ""
  run_testlist(SimpleODETests, solvePDE, tags)
end


#------------------------------------------------------------------------------
# cleanup
