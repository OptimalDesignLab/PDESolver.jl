# tests for DG functionality

# input file to modify, same for all test functions
const test_dg_inputfile = "input_vals_channel_dg.jl"
#const test_dg_moddict = Dict{String, Any}("Flux_name" => "RoeFlux", "use_DG" => true, "new_fname" => "input_vals_channel_dg")

"""
  This functino tests tests calculating fluxes as used for the DG face
  integrals.
"""
function test_dg_flux(mesh, sbp, eqn, opts)
  @testset "----- Testing DG flux -----" begin

    # test the Roe Flux
    uL = [1.0, 2.0, 3.0, 7.0]
    uR = copy(uL)
    flux_roe = zeros(4)
    flux_euler = zeros(4)
    func = EulerEquationMod.FluxDict["RoeFlux"]
    for i=1:mesh.numInterfaces
      iface = mesh.interfaces[i]
      for j=1:mesh.sbpface.numnodes
        eqn.aux_vars_bndry[1, j, i] = EulerEquationMod.calcPressure(eqn.params, uL)
        aux_vars = eqn.aux_vars_face[:, j, i]

        nrm_scaled = sview(mesh.nrm_face, :, j, i)
        EulerEquationMod.calcEulerFlux(eqn.params, uL, aux_vars, nrm_scaled, flux_euler)

        func(eqn.params, uL, uR, aux_vars, nrm_scaled, flux_roe)

        @test isapprox( flux_roe, flux_euler) atol=1e-13
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
        aux_vars = eqn.aux_vars_face[:, j, i]
        nrm = mesh.sbpface.normal[:, iface.faceL]

        nrm_scaled = sview(mesh.nrm_face, :, j, i)
        EulerEquationMod.calcEulerFlux(eqn.params, uL, aux_vars, nrm_scaled, flux_euler)

        @test isapprox( eqn.flux_face[:, j, i], flux_euler) atol=1e-13
      end
    end

  end  # end facts block

  return nothing
end  # end function

#test_dg_flux(mesh, sbp, eqn, opts)
add_func2!(EulerTests, test_dg_flux, test_dg_inputfile, [TAG_FLUX, TAG_SHORTTEST])

"""
  This function tests DG boundary integrals, including interpolation to
  interfaces.
"""
function test_dg_boundary(mesh, sbp, eqn, opts)

  @testset "----- Testing DG Boundary -----" begin

    EulerEquationMod.ICRho1E2U3(mesh, sbp, eqn, opts, eqn.q_vec)
    EulerEquationMod.interpolateBoundary(mesh, sbp, eqn, opts, eqn.q, eqn.q_bndry, eqn.aux_vars_bndry)
    mesh.bndry_funcs[1:end] = EulerEquationMod.BCDict["Rho1E2U3BC"]

    # check that the interpolation worked
    for i=1:mesh.numBoundaryFaces
      for j=1:mesh.sbpface.numnodes
        @test isapprox( eqn.q_bndry[:, j, i], [1.0, 0.35355, 0.35355, 2.0]) atol=1e-13
      end
    end

    uL = eqn.q_bndry[:, 1, 1]
    flux_euler = zeros(4)
    EulerEquationMod.getBCFluxes(mesh, sbp, eqn, opts)

    for i=1:mesh.numBoundaryFaces
      bndry_i = mesh.bndryfaces[i]
      for j=1:mesh.sbpface.numnodes
        eqn.aux_vars_bndry[1, j, i] = EulerEquationMod.calcPressure(eqn.params, eqn.q_bndry[:, j, i])
        aux_vars = eqn.aux_vars_bndry[:, j, i]
        nrm_scaled = sview(mesh.nrm_bndry, :, j, i)

        EulerEquationMod.calcEulerFlux(eqn.params, uL, aux_vars, nrm_scaled, flux_euler)

        @test isapprox( eqn.bndryflux[:, j, i], flux_euler) atol=1e-13
      end
    end


  end  # end facts block

  return nothing
end  # end function

#test_dg_boundary(mesh, sbp, eqn, opts)
add_func2!(EulerTests, test_dg_boundary, test_dg_inputfile, [TAG_BC, TAG_SHORTTEST])

"""
  This functions tests that a uniform flow gives zero residual
"""
function test_dg_uniform(mesh, sbp, eqn, opts)

  # reset eqn
  mesh, sbp, eqn, opts = solvePDE(ARGS[1])

  @testset "----- Testing Uniform Channel -----" begin

    physicsRhs(mesh, sbp, eqn, opts, eqn.res_vec, (evalResidual,))

    for i=1:mesh.numDof
      @test isapprox( eqn.res_vec[i], 0.0) atol=1e-13
    end

  end  # end facts block

  return nothing
end  # end function

#test_dg_uniform(mesh, sbp, eqn, opts)
add_func2!(EulerTests, test_dg_uniform, test_dg_inputfile, [TAG_SHORTTEST])
