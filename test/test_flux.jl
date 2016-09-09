
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
facts("----- Testing Numerical Fluxes -----") do

  ARGS[1] = "input_vals_channel_dg.jl"
  include(STARTUP_PATH)
  println("testing 2d")

  qL = [1.0, 2.0, 3.0, 7.0]
  qR = qL + 1

  F_euler = zeros(qL)
  F_num = zeros(F_euler)
  F_num2 = zeros(F_euler)


  aux_vars = zeros(1)
  aux_vars[1] = EulerEquationMod.calcPressure(eqn.params, qL)

  nrm = [1.0, 1]

  # get the euler flux
  EulerEquationMod.calcEulerFlux(eqn.params, qL, aux_vars, nrm, F_euler)

  function test_symmetric_flux(functor, F_num, F_num2)

    # test symmetry
    functor(eqn.params, qL, qR, aux_vars, nrm, F_num)
    functor(eqn.params, qR, qL, aux_vars, nrm, F_num2)
    for i=1:length(F_num)
      @fact F_num[i] --> roughly(F_num2[i], atol=1e-12)
    end

    # test consistency
    functor(eqn.params, qL, qL, aux_vars, nrm, F_num)
    for i=1:length(F_num)
      @fact F_num[i] --> roughly(F_euler[i])
    end

  end

  functor = EulerEquationMod.FluxDict["StandardFlux"]
  println("testing StandardFlux")
  test_symmetric_flux(functor, F_num, F_num2)

  println("testing DucrosFlux")
  functor = EulerEquationMod.FluxDict["DucrosFlux"]
  test_symmetric_flux(functor, F_num, F_num2)

  println("testing IRFlux")
  functor = EulerEquationMod.FluxDict["IRFlux"]
  test_symmetric_flux(functor, F_num, F_num2)

  # test calculating -Q*f + Eij*f_star_ij = Q^T_ij f_star_ij
  fill!(eqn.res, 0.0)
  EulerEquationMod.calcVolumeIntegralsSplitForm(mesh, sbp, eqn, opts, eqn.volume_flux_func)
  res_split = copy(eqn.res)
  fill!(eqn.res, 0.0)

  opts["Q_transpose"] = false
  EulerEquationMod.getEulerFlux(mesh, sbp, eqn, opts)
  EulerEquationMod.evalVolumeIntegrals(mesh, sbp, eqn, opts)
  opts["Q_transpose"] = true

  E = zeros(sbp.Q)
  for dim=1:2
    E[:, :, dim] = sbp.Q[:, :, dim] + sbp.Q[:, :, dim].'
  end
  F_tmp = zeros(4)
  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      q_j = sview(eqn.q, :, j, i)
      aux_vars = sview(eqn.aux_vars, :, j, i)
      for k=1:mesh.numNodesPerElement
        q_k = sview(eqn.q, :, k, i)
        for dim=1:2
          nrm = mesh.dxidx[dim, :, j, i]
          EulerEquationMod.calcEulerFlux_standard(eqn.params, q_j, q_k, aux_vars, nrm, F_tmp)
          eqn.res[:, j, i] += 2*E[j, k, dim]*F_tmp
        end
      end
    end
  end

  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      for k=1:size(res_split, 1)
        @fact res_split[k, j, i] --> roughly(eqn.res[k, j, i], atol=1e-12)
      end
    end
  end



  # test 3D
  ARGS[1] = "input_vals_3d.jl"
  include(STARTUP_PATH)

  println("testing 3d")
  qL =  [1., 2, 3, 4, 15]
  qR =  qL + 1
  F_euler = zeros(qL)
  F_num = zeros(F_euler)
  F_num2 = zeros(F_euler)
  nrm = [1., 1, 1]

  # get the euler flux
  EulerEquationMod.calcEulerFlux(eqn.params, qL, aux_vars, nrm, F_euler)
  aux_vars[1] = EulerEquationMod.calcPressure(eqn.params, qL)
  println("testing StandardFlux")
  functor = EulerEquationMod.FluxDict["StandardFlux"]
  test_symmetric_flux(functor, F_num, F_num2)

  println("testing DucrosFlux")
  functor = EulerEquationMod.FluxDict["DucrosFlux"]
  test_symmetric_flux(functor, F_num, F_num2)

  println("testing IRFlux")
  functor = EulerEquationMod.FluxDict["IRFlux"]
  test_symmetric_flux(functor, F_num, F_num2)


  # test calculating -Q*f + Eij*f_star_ij = Q^T_ij f_star_ij
  fill!(eqn.res, 0.0)
  EulerEquationMod.calcVolumeIntegralsSplitForm(mesh, sbp, eqn, opts, eqn.volume_flux_func)
  res_split = copy(eqn.res)
  fill!(eqn.res, 0.0)

  opts["Q_transpose"] = false
  EulerEquationMod.getEulerFlux(mesh, sbp, eqn, opts)
  EulerEquationMod.evalVolumeIntegrals(mesh, sbp, eqn, opts)
  opts["Q_transpose"] = true

  E = zeros(sbp.Q)
  for dim=1:3
    E[:, :, dim] = sbp.Q[:, :, dim] + sbp.Q[:, :, dim].'
  end
  F_tmp = zeros(5)
  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      q_j = sview(eqn.q, :, j, i)
      aux_vars = sview(eqn.aux_vars, :, j, i)
      for k=1:mesh.numNodesPerElement
        q_k = sview(eqn.q, :, k, i)
        for dim=1:3
          nrm = mesh.dxidx[dim, :, j, i]
          EulerEquationMod.calcEulerFlux_standard(eqn.params, q_j, q_k, aux_vars, nrm, F_tmp)
          eqn.res[:, j, i] += 2*E[j, k, dim]*F_tmp
        end
      end
    end
  end

  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      for k=1:size(res_split, 1)
        @fact res_split[k, j, i] --> roughly(eqn.res[k, j, i], atol=1e-12)
      end
    end
  end





end
