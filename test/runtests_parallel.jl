# run tests in parallel with 2 processes

push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/euler"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))


using PDESolver
#using Base.Test
using FactCheck
using ODLCommonTools
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
using EulerEquationMod
using ForwardDiff
using NonlinearSolvers   # non-linear solvers
using ArrayViews
include( joinpath(Pkg.dir("PDESolver"), "src/solver/euler/complexify.jl"))
include( joinpath(Pkg.dir("PDESolver"), "src/input/make_input.jl"))
global const STARTUP_PATH = joinpath(Pkg.dir("PDESolver"), "src/solver/euler/startup.jl")

#------------------------------------------------------------------------------
# define test list
include("./TestSystem.jl")
global const EulerTests = TestList()
# define global const tags here


function test_parallel2()
  facts("----- Testing Parallel -----") do

    start_dir = pwd()
    cd ("./rk4/parallel")
    ARGS[1] = "input_vals_parallel.jl"
    include(STARTUP_PATH)

    datas = readdlm("../serial/error_calc.dat")
    datap = readdlm("error_calc.dat")

    @fact datas[1] --> roughly(datap[1], atol=1e-13)
    @fact datas[2] --> roughly(datap[2], atol=1e-13)
    cd("../../")

    cd("./newton/parallel")
    ARGS[1] = "input_vals_parallel.jl"
    include(STARTUP_PATH)

    datas = readdlm("../serial/error_calc.dat")
    datap = readdlm("./error_calc.dat")
    @fact datas[1] --> roughly(datap[1], atol=1e-13)

    cd(start_dir)

  end

  return nothing
end

#test_parallel2()
add_func1!(EulerTests, test_parallel2)

#------------------------------------------------------------------------------
# run tests
facts("----- Running Euler 2 process tests -----") do
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

if MPI.Initialized()
  MPI.Finalize()
end



FactCheck.exitstatus()


