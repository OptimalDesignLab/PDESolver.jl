
#=
push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/advection"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))

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
=#

function make_input_mms(degree; dg=false)
  arg_dict["order"] = degree
  arg_dict["IC_name"] = "ICp$degree"
  arg_dict["BC1_name"] = string("p", degree, "BC")
  arg_dict["SRCname"] = "SRCp$degree"

  if dg
    arg_dict["use_DG"] = true
    arg_dict["Flux_name"] = "LFFlux"
  else
    arg_dict["use_DG"] = false
  end

  fname = "input_vals_mms$degree.jl"
  rmfile(fname)
  f = open(fname, "w")
  print(f, "arg_dict = ")
  println(f, arg_dict)

  return fname
end

facts("----- Testing using manufactured polynomials -----") do

  println("  -----testing degree 1 polynomial -----")
  resize!(ARGS, 1)
  ARGS[1] = "input_vals_mms.jl"
  include(STARTUP_PATH)
  fill!(eqn.res, 0.0)
  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
  println("eqn.q_vec = ", eqn.q_vec)
  AdvectionEquationMod.evalAdvection(mesh, sbp, eqn, opts)
  eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
  @fact eqn.res_vec --> roughly(zeros(mesh.numDof), atol=1e-12)

  fname = make_input_mms(1, dg=true)
  include(STARTUP_PATH)
  fill!(eqn.res, 0.0)
  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
  println("eqn.q_vec = ", eqn.q_vec)
  AdvectionEquationMod.evalAdvection(mesh, sbp, eqn, opts)
  eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
  @fact eqn.res_vec --> roughly(zeros(mesh.numDof), atol=1e-12)




  println("  -----testing degree 2 polynomial -----")
  fname = make_input_mms(2)
  ARGS[1] = fname
  include(STARTUP_PATH)
  fill!(eqn.res, 0.0)
  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
  println("eqn.q_vec = ", eqn.q_vec)
  AdvectionEquationMod.evalAdvection(mesh, sbp, eqn, opts)
  eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
  @fact eqn.res_vec --> roughly(zeros(mesh.numDof), atol=1e-12)

  fname = make_input_mms(2, dg=true)
  include(STARTUP_PATH)
  fill!(eqn.res, 0.0)
  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
  println("eqn.q_vec = ", eqn.q_vec)
  AdvectionEquationMod.evalAdvection(mesh, sbp, eqn, opts)
  eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
  @fact eqn.res_vec --> roughly(zeros(mesh.numDof), atol=1e-12)



  println("  -----testing degree 3 polynomial -----")
  fname = make_input_mms(3)
  ARGS[1] = fname
  include(STARTUP_PATH)
  fill!(eqn.res, 0.0)
  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
  println("eqn.q_vec = ", eqn.q_vec)
  AdvectionEquationMod.evalAdvection(mesh, sbp, eqn, opts)
  eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
  @fact eqn.res_vec --> roughly(zeros(mesh.numDof), atol=1e-12)

  fname = make_input_mms(3, dg=true)
  include(STARTUP_PATH)
  fill!(eqn.res, 0.0)
  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
  println("eqn.q_vec = ", eqn.q_vec)
  AdvectionEquationMod.evalAdvection(mesh, sbp, eqn, opts)
  eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
  @fact eqn.res_vec --> roughly(zeros(mesh.numDof), atol=1e-12)



  println("  -----testing degree 4 polynomial -----")
  fname = make_input_mms(4)
  ARGS[1] = fname
  include(STARTUP_PATH)
  fill!(eqn.res, 0.0)
  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
#  println("eqn.q_vec = ", eqn.q_vec)
  AdvectionEquationMod.evalAdvection(mesh, sbp, eqn, opts)
  eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
  @fact eqn.res_vec --> roughly(zeros(mesh.numDof), atol=1e-12)

  fname = make_input_mms(4, dg=true)
  include(STARTUP_PATH)
  fill!(eqn.res, 0.0)
  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
  println("eqn.q_vec = ", eqn.q_vec)
  AdvectionEquationMod.evalAdvection(mesh, sbp, eqn, opts)
  eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
  @fact eqn.res_vec --> roughly(zeros(mesh.numDof), atol=1e-12)



end



