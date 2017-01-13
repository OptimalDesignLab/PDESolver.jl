# run 2 processor tests

push!(LOAD_PATH, abspath(joinpath(pwd(), "..")))

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
using Utils
using Input

#------------------------------------------------------------------------------
# define tests and tags

#include("../TestSystem.jl")
using TestSystem
# define tags that will be used
include("../tags.jl")

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
    mesh, sbp, eqn, opts = run_advection(ARGS[1])

    datas = readdlm("../serial/error_calc.dat")
    datap = readdlm("error_calc.dat")

    @fact datas[1] --> roughly(datap[1], atol=1e-13)
    @fact datas[2] --> roughly(datap[2], atol=1e-13)
    cd("../../")

    cd("./newton/parallel")
    ARGS[1] = "input_vals_parallel.jl"
    mesh, sbp, eqn, opts = run_advection(ARGS[1])

    datas = readdlm("../serial/error_calc.dat")
    datap = readdlm("./error_calc.dat")
    @fact datas[1] --> roughly(datap[1], atol=1e-13)

    cd(start_dir)

    cd("./rk4_3d/parallel")
    ARGS[1] = "input_vals_parallel.jl"
    mesh, sbp, eqn, opts = run_advection(ARGS[1])

    datas = readdlm("../serial/error_calc.dat")
    datap = readdlm("error_calc.dat")
    @fact datas[1] --> roughly(datap[1], atol=1e-13)

    cd(start_dir)

    cd("./newton_3d/parallel")
    ARGS[1] = "input_vals_parallel.jl"
    mesh, sbp, eqn, opts = run_advection(ARGS[1])
    datas = readdlm("../serial/error_calc.dat")
    datap = readdlm("error_calc.dat")
    @fact datas[1] --> roughly(datap[1], atol=1e-13)

    cd(start_dir)
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
  run_testlist(AdvectionTests, run_advection, tags)
end

#------------------------------------------------------------------------------
# cleanup

#=
facts("----- Testing Functional Computation On Boundary In Parallel -----") do

  resize!(ARGS, 1)
  ARGS[1] = "input_vals_functional_DG_parallel.jl"
  include(STARTUP_PATH)

  @fact mesh.isDG --> true
  @fact opts["functional_name1"] --> "qflux"
  @fact opts["functional_error"] --> true
  @fact opts["smb_name"] --> "src/mesh_files/gsquare2np2.smb"
  @fact opts["analytical_functional_val"] --> roughly(2*(exp(1) - 1), atol=1e-12)
  @fact opts["geom_edges_functional1"] --> [1,2]

  fname = "./functional_error1.dat"
  error = readdlm(fname)

  @fact error[1] --> roughly(0.00681567877682826, atol=1e-6)

end
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

