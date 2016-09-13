# run tests in parallel

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


start_dir = pwd()
resize!(ARGS, 1)

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

if MPI.Initialized()
  MPI.Finalize()
end

FactCheck.exitstatus()


