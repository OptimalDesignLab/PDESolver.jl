# test the basic DG functions
#=
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
=#
# insert a command line argument
resize!(ARGS, 1)

ARGS[1] = "input_vals_channel.jl"
include(STARTUP_PATH)

arg_dict["use_DG"] = true
arg_dict["Flux_name"] = "LFFlux"
arg_dict["BC1_name"] = "p1BC"

make_input(arg_dict, "input_vals_channelDG")
ARGS[1] = "input_vals_channelDG.jl"
include(STARTUP_PATH)

facts("----- Testing DG Flux ------") do
  eqn.params.LFalpha = 1.0
  dxidx1 = mesh.dxidx_face[:, :, 1, 1]
  nrm = view(sbp.facenormal, :, mesh.interfaces[1].faceL)
  alpha = [eqn.alpha_x, eqn.alpha_y]
  alpha_n = sum((dxidx1*alpha).*nrm)
  qL = 1.0
  qR = 2.0
  flux_test = alpha_n*(qL + qR)/2

  flux_func = AdvectionEquationMod.FluxDict["LFFlux"]
  flux_code = flux_func(qL, qR, eqn.alpha_x, eqn.alpha_y, dxidx1, nrm, eqn.params)

  @fact flux_code --> roughly(flux_test, atol=1e-13)

  eqn.q_face[1, 1, :, 1] = 1.0
  eqn.q_face[1, 2, :, 1] = 2.0

  AdvectionEquationMod.calcFaceFlux(mesh, sbp, eqn, eqn.flux_func, mesh.interfaces, eqn.flux_face)

  for i=1:mesh.sbpface.numnodes
    @fact eqn.flux_face[1, i, 1] --> roughly(-flux_test, atol=1e-13)
  end


end

facts("\n----- Testing DG Boundary Condition -----") do

  eqn.params.LFalpha = 1.0

  for i=1:mesh.sbpface.numnodes
    eqn.q_bndry[1, i, :] = 2.0
  end

  # test use of eqn.q_bndry for BC
  eqn.alpha_x = -1.0
  eqn.alpha_y = -1.0
  range_idx = 1:mesh.numBoundaryEdges
  AdvectionEquationMod.calcBoundaryFlux(mesh, sbp, eqn, mesh.bndry_funcs[1], range_idx, mesh.bndryfaces, eqn.bndryflux)

  val_code = 0.0
  for i=1:mesh.sbpface.numnodes
    val_code += mesh.sbpface.wface[i]*eqn.bndryflux[1, i, 1]
  end
  val_test = 4*eqn.q_bndry[1,1,1]*eqn.alpha_x
  @fact val_code --> roughly(val_test, atol=1e-13)


  # test use of the boundary condition value
  eqn.alpha_x = 1.0
  eqn.alpha_y = 1.0
  bndry_coords = mesh.coords_bndry[:, :, 1]

  AdvectionEquationMod.calcBoundaryFlux(mesh, sbp, eqn, mesh.bndry_funcs[1], range_idx, mesh.bndryfaces, eqn.bndryflux)
  val_code = 0.0
  for i=1:mesh.sbpface.numnodes
    val_code += mesh.sbpface.wface[i]*eqn.bndryflux[1, i, 1]
  end
  val_test = 12.0

  @fact val_code --> roughly(val_test, atol=1e-13)


  # check that the interpolation and coordinates match
  fill!(eqn.q_bndry, 0.0)
  AdvectionEquationMod.ICp1(mesh, sbp, eqn, opts, eqn.q_vec)
  mesh.bndry_funcs[1] = AdvectionEquationMod.BCDict["p1BC"]
  AdvectionEquationMod.evalBndry(mesh, sbp, eqn)

  for i=1:mesh.numBoundaryEdges
    for j=1:mesh.sbpface.numnodes
      coords = mesh.coords_bndry[:, j, i]
      q_test = AdvectionEquationMod.calc_p1(coords, eqn.alpha_x, eqn.alpha_y, 0.0)
      q_code = eqn.q_bndry[1, j, i]
      @fact q_code --> roughly(q_test, atol=1e-13)
    end
  end



end
