# Run advection tests

push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/simpleODE"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))
include(joinpath(Pkg.dir("PDESolver"), "src/input/make_input.jl"))

using PDESolver
#using Base.Test
using FactCheck
using ODLCommonTools
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
using SimpleODEMod
using ForwardDiff
using NonlinearSolvers   # non-linear solvers
using ArrayViews

function clean_dict(collection)
  for i in keys(collection)
    delete!(collection, i)
  end
end

global const STARTUP_PATH = joinpath(Pkg.dir("PDESolver"), "src/solver/simpleODE/startup_simpleODE.jl")

#------------------------------------------------------------------------------
# define tests and tags


include("../TestSystem.jl")
# define tags that will be used

# test list
global const SimpleODETests = TestList()


function test_eq4()
  facts("---- testing SimpleODE ----") do
    start_dir = pwd()
    cd("./eqn4/")
    ARGS[1] = "input_vals_simpleODE.jl"
    include(STARTUP_PATH)

    for i = 1:length(eqn.q_vec)
      @fact eqn.q_vec[i] --> roughly(4.0, atol=1e-10)
    end

    cd(start_dir)
  end  # end facts block

  return nothing
end

add_func1!(SimpleODETests, test_eq4)

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
  run_testlist(SimpleODETests, tags)
end


#------------------------------------------------------------------------------
# cleanup
if MPI.Initialized()
  MPI.Finalize()
end

FactCheck.exitstatus()

