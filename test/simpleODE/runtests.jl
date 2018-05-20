# Run advection tests

push!(LOAD_PATH, abspath(joinpath(pwd(), "..")))

using PDESolver
#using Base.Test
using FactCheck
using ODLCommonTools
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
using Utils
using SimpleODEMod
using ForwardDiff
using NonlinearSolvers   # non-linear solvers
using ArrayViews
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
  facts("---- testing SimpleODE ----") do
    start_dir = pwd()
    cd("./eqn4/")
    ARGS[1] = "input_vals_simpleODE.jl"
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])

    for i = 1:length(eqn.q_vec)
      @fact eqn.q_vec[i] --> roughly(4.0, atol=1e-10)
    end

    cd(start_dir)
  end  # end facts block

  return nothing
end

add_func1!(SimpleODETests, test_eq4, [TAG_SHORTTEST])

#------------------------------------------------------------------------------
# run tests
facts("----- Running SimpleODE tests -----") do
  nargs = length(ARGS)
  if nargs == 0
    tags = ASCIIString[TAG_DEFAULT]
  else
    tags = Array(ASCIIString, nargs)
    copy!(tags, ARGS)
  end

  resize!(ARGS, 1)
  ARGS[1] = ""
  run_testlist(SimpleODETests, solvePDE, tags)
end


#------------------------------------------------------------------------------
# cleanup

FactCheck.exitstatus()

