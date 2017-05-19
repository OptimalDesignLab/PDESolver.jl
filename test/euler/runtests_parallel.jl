# run tests in parallel with 2 processes

push!(LOAD_PATH, abspath(joinpath(pwd(), "..")))

using PDESolver
#using Base.Test
using FactCheck
using ODLCommonTools
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
using Utils
using EulerEquationMod
using ForwardDiff
using NonlinearSolvers   # non-linear solvers
using ArrayViews
import MPI
using Input

#------------------------------------------------------------------------------
# define test list
using TestSystem
include("../tags.jl")

global const EulerTests = TestList()
# define global const tags here

# NOTE: For adding parallel tests to TICON, Write them in a separate file and
# include them below `global const EulerTests = TestList()`.

"""
  Run the parallel tests and compare against serial results run as part of
  the serial tests
"""
function test_parallel2()
  facts("----- Testing Parallel -----") do

    start_dir = pwd()
    cd ("./rk4/parallel")
    ARGS[1] = "input_vals_parallel.jl"
    mesh, sbp, eqn, opts = run_euler(ARGS[1])

    datas = readdlm("../serial/error_calc.dat")
    datap = readdlm("error_calc.dat")

    @fact datas[1] --> roughly(datap[1], atol=1e-13)
    @fact datas[2] --> roughly(datap[2], atol=1e-13)
    cd("../../")

    cd("./newton/parallel")
    ARGS[1] = "input_vals_parallel.jl"
    mesh, sbp, eqn, opts = run_euler(ARGS[1])

    datas = readdlm("../serial/error_calc.dat")
    datap = readdlm("./error_calc.dat")
    @fact datas[1] --> roughly(datap[1], atol=1e-13)

    cd(start_dir)

  end

  return nothing
end

#test_parallel2()
add_func1!(EulerTests, test_parallel2, [TAG_SHORTTEST])

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
  run_testlist(EulerTests, run_euler, tags)
end

#=
facts("--- Testing Functional Computation On a Boundary in Parallel ---") do

  ARGS[1] = "input_vals_vortex_adjoint_DG_parallel.jl"
  include("../src/solver/euler/startup.jl")

  @fact mesh.isDG --> true
  @fact opts["calc_functional"] --> true
  @fact opts["functional_error"] --> true
  @fact opts["functional_name1"] --> "drag"
  @fact opts["analytical_functional_val"] --> roughly(-1/1.4, atol = 1e-13)
  @fact opts["geom_edges_functional1"] --> [3]

  fname = "./functional_error1.dat"
  relative_error = readdlm(fname)

  @fact relative_error[1] --> roughly(0.000177342284, atol = 1e-6)

  # rm("./functional_error1.dat") # Delete the file


end  # End do
=#

# define global variable if needed
# this trick allows running the test files for multiple physics in the same
# session without finalizing MPI too soon
if !isdefined(:TestFinalizeMPI)
  TestFinalizeMPI = true
end


if MPI.Initialized() && TestFinalizeMPI
  MPI.Finalize()
end

FactCheck.exitstatus()


