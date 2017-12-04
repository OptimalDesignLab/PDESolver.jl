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
using OptimizationInterface
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
    cd("./rk4/parallel")
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
add_func1!(AdvectionTests, runtests_parallel, [TAG_SHORTTEST])

function test_precompute()
  facts("----- testing non-precompute functions -----") do
    start_dir = pwd()

    # test rk4
    cd("./rk4/parallel")
    ARGS[1] = "input_vals_parallel_runp.jl"
    #TODO: set opts["solve"] = false before doing this
    mesh, sbp, eqn, opts = run_advection(ARGS[1])

    fill!(eqn.res, 0.0)
    evalResidual(mesh, sbp, eqn, opts)

    res_orig = copy(eqn.res)

    opts["precompute_face_flux"] = false
    evalResidual(mesh, sbp, eqn, opts)

    @fact norm(vec(eqn.res - res_orig)) --> roughly(0.0, atol=1e-13)

    # test newton
    cd(start_dir)
    cd("./newton/parallel")
    ARGS[1] = "input_vals_parallel.jl"
    mesh, sbp, eqn, opts = run_advection(ARGS[1])

    fill!(eqn.res, 0.0)
    evalResidual(mesh, sbp, eqn, opts)
    res_orig = copy(eqn.res)

    opts["precompute_face_flux"] = false
    evalResidual(mesh, sbp, eqn, opts)

    @fact norm(vec(eqn.res - res_orig)) --> roughly(0.0, atol=1e-13)

    cd(start_dir)
  end


  return nothing
end

add_func1!(AdvectionTests, test_precompute, [TAG_SHORTTEST])

function test_adjoint_parallel()

  facts("--- Testing Adjoint Computation on a Geometric Boundary ---") do

    resize!(ARGS, 1)
    ARGS[1] = "input_vals_functional_DG_parallel.jl"
    # include(STARTUP_PATH)
    mesh, sbp, eqn, opts = run_advection(ARGS[1])

    objective = AdvectionEquationMod.createObjectiveFunctionalData(mesh, sbp, eqn, opts)
    AdvectionEquationMod.evalFunctional(mesh, sbp, eqn, opts, objective)
    adjoint_vec = zeros(Complex{Float64}, mesh.numDof)
    AdvectionEquationMod.calcAdjoint(mesh, sbp, eqn, opts, objective, adjoint_vec)

    for i = 1:length(adjoint_vec)
      @fact real(adjoint_vec[i]) --> roughly(1.0 , atol=1e-10)
    end
   
  end # End facts("--- Testing Functional Computation on a Geometric Boundary ---")

  return nothing
end

add_func1!(AdvectionTests, test_adjoint_parallel, [TAG_ADJOINT, TAG_LONGTEST])

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

FactCheck.exitstatus()
