# test_flux.jl: test flux calculation

# logarithmic mean used by IR flux
function logmean(aL, aR)
  xi = aL/aR
  f = (xi - 1)/(xi + 1)
  u = f*f
  eps = 1e-2
  if ( u < eps)
    F = 1.0 + u/3.0 + u*u/5.0 + u*u*u/7.0
  else
    F = log(xi)/2.0/f
  end

  return (aL + aR)/(2*F)
end

"""
  Test that a flux is symmetric.  This is not a test function itself, but it
  is called by test functions
"""
function test_symmetric_flux(functor, params, qL, qR, aux_vars, nrm, F_num, F_num2, F_euler)

  # test symmetry
  functor(params, qL, qR, aux_vars, nrm, F_num)
  functor(params, qR, qL, aux_vars, nrm, F_num2)
  for i=1:length(F_num)
    @fact F_num[i] --> roughly(F_num2[i], atol=1e-12)
  end

  # test consistency
  functor(params, qL, qL, aux_vars, nrm, F_num)
  for i=1:length(F_num)
    @fact F_num[i] --> roughly(F_euler[i])
  end

end

"""
  An inefficient version of the IR flux, used for comparison to the efficient
  version.
"""
function ir_flux(params, qL, qR, nrm)

  pL = EulerEquationMod.calcPressure(params, qL)
  pR = EulerEquationMod.calcPressure(params, qR)
  z5_ln = logmean(sqrt(qL[1]*pL), sqrt(qR[1]*pR))
  rho_hat = 0.5*(sqrt(qL[1]/pL) + sqrt(qR[1]/pR))*z5_ln

  z1L = sqrt(qL[1]/pL); z1R = sqrt(qR[1]/pR)
  z1_avg = 0.5*( z1L + z1R )

  u_hat = 0.5*(z1L*qL[2]/qL[1] + z1R*qR[2]/qR[1])/z1_avg
  v_hat = 0.5*(z1L*qL[3]/qL[1] + z1R*qR[3]/qR[1])/z1_avg

  p1_hat = 0.5*( sqrt(qL[1]*pL) + sqrt(qR[1]*pR) )/z1_avg
  z1_ln = logmean(sqrt(qL[1]/pL), sqrt(qR[1]/pR))

  p2_hat = (params.gamma + 1)*z5_ln/(2*params.gamma*z1_ln) + params.gamma_1*0.5*(sqrt(qL[1]*pL) + sqrt(qR[1]*pR))/(2*params.gamma*z1_avg)
  h_hat = params.gamma*p2_hat/(rho_hat*params.gamma_1) + 0.5*(u_hat*u_hat + v_hat*v_hat)

  fx = [rho_hat*u_hat, rho_hat*u_hat*u_hat + p1_hat, rho_hat*u_hat*v_hat, rho_hat*u_hat*h_hat]
  fy = [rho_hat*v_hat, rho_hat*u_hat*v_hat, rho_hat*v_hat*v_hat + p1_hat, rho_hat*v_hat*h_hat]

  return fx*nrm[1] + fy*nrm[2]

end


"""
  Test calculation of numerical flux functions in 2D
"""
function test_flux_2d()
  ARGS[1] = "input_vals_channel_dg.jl"
  mesh, sbp, eqn, opts = run_euler(ARGS[1])

  facts("----- Testing 2D Numerical Fluxes -----") do
 
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

    functor = EulerEquationMod.FluxDict["StandardFlux"]
    test_symmetric_flux(functor, eqn.params, qL, qR, aux_vars, nrm, F_num, F_num2, F_euler)

    functor = EulerEquationMod.FluxDict["DucrosFlux"]
    test_symmetric_flux(functor, eqn.params, qL, qR, aux_vars, nrm, F_num, F_num2, F_euler)

    functor = EulerEquationMod.FluxDict["IRFlux"]
    test_symmetric_flux(functor, eqn.params, qL, qR, aux_vars, nrm, F_num, F_num2, F_euler)

    # test against calculated solution
    nrm = [1., 2.0]
    flux_test = ir_flux(eqn.params, qL, qR, nrm)

    F_code = zeros(4)
    functor(eqn.params, qL, qR, aux_vars, nrm, F_code)
    @fact F_code --> roughly(flux_test, atol=1e-12)

    # test calculating -Q*f = -(2*S_ij f_star_ij + Eij*f_star_ij)
    # set eqn.q to something interesting
    ic_func = ICDict["ICExp"]
    ic_func(mesh, sbp, eqn, opts, eqn.q_vec)

    disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    fill!(eqn.res, 0.0)
    EulerEquationMod.calcVolumeIntegralsSplitForm(mesh, sbp, eqn, opts, eqn.volume_flux_func)
    res_split = copy(eqn.res)
    fill!(eqn.res, 0.0)

    # check that 1^T * volume integrals using S * F_Star * 1 == 0, because 
    # S is skew symmetric and F_star is symmetric
    for i=1:mesh.numEl
      val = sum(res_split[:, :, i])
      @fact val --> roughly(0.0, atol=1e-13)
    end

    opts["Q_transpose"] = false
    EulerEquationMod.getEulerFlux(mesh, sbp, eqn, opts)
    EulerEquationMod.evalVolumeIntegrals(mesh, sbp, eqn, opts)
    opts["Q_transpose"] = true

    @fact eqn.volume_flux_func --> EulerEquationMod.FluxDict["StandardFlux"]
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
            res_split[:, j, i] -= E[j, k, dim]*F_tmp
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

    # test that constant field -> 0 residual
    ic_func = EulerEquationMod.ICDict[opts["IC_name"]]
    ic_func(mesh, sbp, eqn, opts, eqn.q_vec)
    disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)

    fill!(eqn.res, 0.0)
    opts["volume_integral_type"] = 2
    opts["Volume_flux_name"] = "IRFlux"
    EulerEquationMod.init(mesh, sbp, eqn, opts)
    EulerEquationMod.evalResidual(mesh, sbp, eqn, opts)
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        for k=1:mesh.numDofPerNode
          @fact eqn.res[k, j, i] --> roughly(0.0, atol=1e-12)
        end
      end
    end

  end  # end facts block

  return nothing
end  # end function

#test_flux_2d()
add_func1!(EulerTests, test_flux_2d, [TAG_VOLUMEINTEGRALS, TAG_FLUX])

"""
  Test calculation of numerical flux functions in 3D
"""
function test_flux_3d()
  # test 3D
  ARGS[1] = "input_vals_3d.jl"
  mesh, sbp, eqn, opts = run_euler(ARGS[1])
  facts("----- testing 3D  Numerical Fluxes -----") do

    qL =  [1., 2, 3, 4, 15]
    qR =  qL + 1
    F_euler = zeros(qL)
    F_num = zeros(F_euler)
    F_num2 = zeros(F_euler)
    nrm = [1., 1, 1]

    aux_vars = zeros(1)
    aux_vars[1] = EulerEquationMod.calcPressure(eqn.params, qL)


    # get the euler flux
    EulerEquationMod.calcEulerFlux(eqn.params, qL, aux_vars, nrm, F_euler)
    aux_vars[1] = EulerEquationMod.calcPressure(eqn.params, qL)
    functor = EulerEquationMod.FluxDict["StandardFlux"]
    test_symmetric_flux(functor, eqn.params, qL, qR, aux_vars, nrm, F_num, F_num2, F_euler)

    functor = EulerEquationMod.FluxDict["DucrosFlux"]
    test_symmetric_flux(functor, eqn.params, qL, qR, aux_vars, nrm, F_num, F_num2, F_euler)

    functor = EulerEquationMod.FluxDict["IRFlux"]
    test_symmetric_flux(functor, eqn.params, qL, qR, aux_vars, nrm, F_num, F_num2, F_euler)


    # test calculating -Q*f = -(2*S_ij f_star_ij + Eij*f_star_ij)
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
            res_split[:, j, i] -= E[j, k, dim]*F_tmp
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

  end  # end facts block

  return nothing
end


#test_flux_3d()
add_func1!(EulerTests, test_flux_3d, [TAG_VOLUMEINTEGRALS, TAG_FLUX])
