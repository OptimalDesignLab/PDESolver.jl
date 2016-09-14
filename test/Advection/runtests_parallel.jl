push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/advection"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))
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

global const STARTUP_PATH = joinpath(Pkg.dir("PDESolver"), "src/solver/advection/startup_advection.jl")
# insert a command line argument
resize!(ARGS, 1)

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

end

cd(start_dir)

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

if MPI.Initialized()
  MPI.Finalize()
end
FactCheck.exitstatus()

