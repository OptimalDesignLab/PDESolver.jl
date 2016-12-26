# run 2 processor tests
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

global const STARTUP_PATH = joinpath(Pkg.dir("PDESolver"), "src/solver/advection/startup.jl")

#------------------------------------------------------------------------------
# define tests and tags

include("../TestSystem.jl")
# define tags that will be used

# test list
global const AdvectionTests = TestList()


"""
  Run parallel tests and compare to serial results calculated as part of
  serial tests.
"""
function runtests_parallel()
  facts("----- Testing Parallel -----") do

    start_dir = pwd()
    cd ("./rk4/parallel")
    ARGS[1] = "input_vals_parallel_runp.jl"
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

    cd("./rk4_3d/parallel")
    ARGS[1] = "input_vals_parallel.jl"
    include(STARTUP_PATH)

    datas = readdlm("../serial/error_calc.dat")
    datap = readdlm("error_calc.dat")
    @fact datas[1] --> roughly(datap[1], atol=1e-13)

    cd(start_dir)

    cd("./newton_3d/parallel")
    ARGS[1] = "input_vals_parallel.jl"
    include(STARTUP_PATH)
    datas = readdlm("../serial/error_calc.dat")
    datap = readdlm("error_calc.dat")
    @fact datas[1] --> roughly(datap[1], atol=1e-13)
  end  # end facts block

  return nothing
end

#runtests_parallel()
add_func1!(AdvectionTests, runtests_parallel)

#------------------------------------------------------------------------------
# run tests
facts("----- Running Advection 2 processor tests -----") do
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

#------------------------------------------------------------------------------
# cleanup

if MPI.Initialized()
  MPI.Finalize()
end
FactCheck.exitstatus()

