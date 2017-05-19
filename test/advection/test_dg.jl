# test advection discontinuous-Galerkin functions

global const test_dg_inputfile = "input_vals_channelDG.jl"

"""
  Test face flux for DG.
"""
function test_dg_flux(mesh, sbp, eqn, opts)
  facts("----- Testing DG Flux ------") do
    eqn.params.LFalpha = 1.0
    dxidx1 = mesh.dxidx_face[:, :, 1, 1]
    nrm = sview(sbp.facenormal, :, mesh.interfaces[1].faceL)
    alpha = [eqn.params.alpha_x, eqn.params.alpha_y]
    alpha_n = sum((dxidx1*alpha).*nrm)
    qL = 1.0
    qR = 2.0
    flux_test = alpha_n*(qL + qR)/2

    flux_func = AdvectionEquationMod.FluxDict["LFFlux"]
    flux_code = flux_func(qL, qR, dxidx1, nrm, eqn.params)

    @fact flux_code --> roughly(flux_test, atol=1e-13)

    eqn.q_face[1, 1, :, 1] = 1.0
    eqn.q_face[1, 2, :, 1] = 2.0

    AdvectionEquationMod.calcFaceFlux(mesh, sbp, eqn, eqn.flux_func, mesh.interfaces, eqn.flux_face)

    for i=1:mesh.sbpface.numnodes
      @fact eqn.flux_face[1, i, 1] --> roughly(-flux_test, atol=1e-13)
    end

  end  # end facts block

  return nothing
end

#test_dg_flux(mesh, sbp, eqn, opts)
add_func2!(AdvectionTests, test_dg_flux, test_dg_inputfile, [TAG_FLUX, TAG_SHORTTEST])

"""
  Test boundary conditions for DG.
"""
function test_dg_bc(mesh, sbp, eqn, opts)
  facts("\n----- Testing DG Boundary Condition -----") do

    eqn.params.LFalpha = 1.0

    for i=1:mesh.sbpface.numnodes
      eqn.q_bndry[1, i, :] = 2.0
    end

    # test use of eqn.q_bndry for BC
    eqn.params.alpha_x = -1.0
    eqn.params.alpha_y = -1.0
    range_idx = 1:mesh.numBoundaryFaces
    AdvectionEquationMod.calcBoundaryFlux(mesh, sbp, eqn, mesh.bndry_funcs[1], range_idx, mesh.bndryfaces, eqn.bndryflux)

    val_code = 0.0
    for i=1:mesh.sbpface.numnodes
      val_code += mesh.sbpface.wface[i]*eqn.bndryflux[1, i, 1]
    end
    val_test = 4*eqn.q_bndry[1,1,1]*eqn.params.alpha_x
    @fact val_code --> roughly(val_test, atol=1e-13)


    # test use of the boundary condition value
    eqn.params.alpha_x = 1.0
    eqn.params.alpha_y = 1.0
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
    AdvectionEquationMod.evalBoundaryIntegrals(mesh, sbp, eqn)

    for i=1:mesh.numBoundaryFaces
      for j=1:mesh.sbpface.numnodes
        coords = mesh.coords_bndry[:, j, i]
        q_test = AdvectionEquationMod.calc_p1(coords, eqn.params, 0.0)
        q_code = eqn.q_bndry[1, j, i]
        @fact q_code --> roughly(q_test, atol=1e-13)
      end
    end

  end  # end facts block

  return nothing
end  # end function

#test_dg_bc(mesh, sbp, eqn, opts)
add_func2!(AdvectionTests, test_dg_flux, test_dg_inputfile, [TAG_BC, TAG_SHORTTEST])
