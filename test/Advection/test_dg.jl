# test the basic DG functions

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

ARGS[1] = "input_vals_channel.jl"
include(STARTUP_PATH)

arg_dict["use_DG"] = true
arg_dict["Flux_name"] = "LFFlux"

make_input(arg_dict, "input_vals_channelDG")
ARGS[1] = "input_vals_channelDG.jl"
include(STARTUP_PATH)

facts("----- Testing DG Flux ------") do

  dxidx1 = mesh.dxidx_face[:, :, 1, 1]
  println("dxidx = ", dxidx1)
  nrm = view(sbp.facenormal, :, mesh.interfaces[1].faceL)
  println("nrm = ", nrm)
  alpha = [eqn.alpha_x, eqn.alpha_y]
  alpha_n = sum((dxidx1*alpha).*nrm)
  println("alpha_n = ", alpha_n)
  qL = 1.0
  qR = 2.0
  flux_test = alpha_n*(qL + qR)/2

  flux_func = AdvectionEquationMod.FluxDict["LFFlux"]
  flux_code = flux_func(qL, qR, eqn.alpha_x, eqn.alpha_y, dxidx1, nrm, eqn.params)

  @fact flux_code --> roughly(flux_test, atol=1e-13)

  eqn.q_bndry[1, 1, :, 1] = 1.0
  eqn.q_bndry[1, 2, :, 1] = 2.0

  AdvectionEquationMod.calcFaceFlux(mesh, sbp, eqn, eqn.flux_func, mesh.interfaces, eqn.flux_face)

  for i=1:mesh.sbpface.numnodes
    @fact eqn.flux_face[1, i, 1] --> roughly(-flux_test, atol=1e-13)
  end


end

facts("----- Testing DG Boundary Condition -----") do

  for i=1:mesh.sbpface.numnodes
    eqn.q_face[1, i, 1] = 2.0
  end

  AdvectionEquationMod.calcBoundaryFlux(mesh, sbp, eqn, mesh.bndry_funcs[1], mesh.bndryfaces, eqn.bndryflux)

  val_code = dot(mesh.sbpface.wface, eqn.bndryflux[1, :, 1])
  val_test = 4*eqn.q_face[1,1,1]*mesh.alpha_x
  @fact val_code --> roughly(val_test, atol=1e-13)


end
