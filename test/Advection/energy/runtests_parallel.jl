#=
push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/advection"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))
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

resize!(ARGS, 1)

global const STARTUP_PATH = joinpath(Pkg.dir("PDESolver"), "src/solver/advection/startup_advection.jl")
=#

facts("----- Testing Energy Stability -----") do

  function test_energy()
    energy_final = calcNorm(eqn, eqn.q_vec)
    q_initial = zeros(eqn.q_vec)
    ICfunc(mesh, sbp, eqn, opts, q_initial)
    energy_initial = calcNorm(eqn, q_initial)

    @fact abs(energy_initial - energy_final) --> less_than(1e-12)
  end


  start_dir = pwd()

  ARGS[1] = "input_vals_periodic.jl"
  cd("./2dp4")
  include(STARTUP_PATH)
  test_energy()

  cd("../3dp4")
  include(STARTUP_PATH)
  test_energy()


end
