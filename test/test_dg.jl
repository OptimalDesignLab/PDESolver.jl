# tests for DG functionality

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
# insert a command line argument
resize!(ARGS, 1)
ARGS[1] = "input_vals_channel.jl"
include(STARTUP_PATH)

arg_dict["Flux_name"] = "RoeFlux"
arg_dict["use_DG"] = true

make_input(arg_dict, "input_vals_channel_dg")
ARGS[1] = "input_vals_channel_dg.jl"
include(STARTUP_PATH)


facts("----- Testing DG flux -----") do

  # test the Roe Flux
  uL = [1.0, 2.0, 3.0, 7.0]
  uR = copy(uL)
  flux_roe = zeros(4)
  flux_euler = zeros(4)
  func = EulerEquationMod.FluxDict["RoeFlux"]
  for i=1:mesh.numInterfaces
    iface = mesh.interfaces[i]
    for j=1:mesh.sbpface.numnodes
      dxidx = mesh.dxidx_face[:, :, j, i]
      eqn.aux_vars_bndry[1, j, i] = EulerEquationMod.calcPressure(eqn.params, uL)
      aux_vars = eqn.aux_vars_face[:, j, i]
      nrm = sbp.facenormal[:, iface.faceL]


      nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
      ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]
      nrm_scaled = [nx, ny]
      EulerEquationMod.calcEulerFlux(eqn.params, uL, aux_vars, nrm_scaled, flux_euler)

      func(uL, uR, aux_vars, dxidx, nrm, flux_roe, eqn.params)

      @fact flux_roe --> roughly(-flux_euler, atol=1e-13)
    end
  end

  # now test calcFaceFlux
  fill!(eqn.flux_face, 0.0)
  for i=1:mesh.numInterfaces
    for j=1:mesh.sbpface.numnodes
      eqn.q_face[:, 1, j, i] = uL
      eqn.q_face[:, 2, j, i] = uL
      eqn.aux_vars_face[1, j, i] = EulerEquationMod.calcPressure(eqn.params, uL)
    end
  end

  EulerEquationMod.calcFaceFlux(mesh, sbp, eqn, func, mesh.interfaces, eqn.flux_face)
  for i=1:mesh.numInterfaces
    iface = mesh.interfaces[i]
    for j=1:mesh.sbpface.numnodes
      dxidx = mesh.dxidx_face[:, :, j, i]
      aux_vars = eqn.aux_vars_face[:, j, i]
      nrm = sbp.facenormal[:, iface.faceL]


      nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
      ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]
      nrm_scaled = [nx, ny]
      EulerEquationMod.calcEulerFlux(eqn.params, uL, aux_vars, nrm_scaled, flux_euler)

      @fact eqn.flux_face[:, j, i] --> roughly(flux_euler, atol=1e-13)
    end
  end



end


facts("----- Testing DG Boundary -----") do

  EulerEquationMod.ICRho1E2U3(mesh, sbp, eqn, opts, eqn.q_vec)
  EulerEquationMod.interpolateBoundary(mesh, sbp, eqn, opts, eqn.q, eqn.q_bndry)
  mesh.bndry_funcs[1:end] = EulerEquationMod.BCDict["Rho1E2U3BC"]

  # check that the interpolation worked
  for i=1:mesh.numBoundaryFaces
    for j=1:mesh.sbpface.numnodes
      @fact eqn.q_bndry[:, j, i] --> roughly( [1.0, 0.35355, 0.35355, 2.0], atol=1e-13)
    end
  end

  uL = eqn.q_bndry[:, 1, 1]
  flux_euler = zeros(4)
  EulerEquationMod.getBCFluxes(mesh, sbp, eqn, opts)

  for i=1:mesh.numBoundaryFaces
    bndry_i = mesh.bndryfaces[i]
    for j=1:mesh.sbpface.numnodes
      dxidx = mesh.dxidx_bndry[:, :, j, i]
      eqn.aux_vars_bndry[1, j, i] = EulerEquationMod.calcPressure(eqn.params, eqn.q_bndry[:, j, i])
      aux_vars = eqn.aux_vars_bndry[:, j, i]
      nrm = sbp.facenormal[:, bndry_i.face]

      nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
      ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]
      nrm_scaled = [nx, ny]
      EulerEquationMod.calcEulerFlux(eqn.params, uL, aux_vars, nrm_scaled, flux_euler)

      @fact eqn.bndryflux[:, j, i] --> roughly(-flux_euler, atol=1e-13)
    end
  end


end

# reset eqn

include(STARTUP_PATH)

facts("----- Testing Uniform Channel -----") do

  calcResidual(mesh, sbp, eqn, opts, evalEuler)

  for i=1:mesh.numDof
    @fact eqn.res_vec[i] --> roughly(0.0, atol=1e-13)
  end

end
