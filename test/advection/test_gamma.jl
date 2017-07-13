# test SBPGamma operators

global const test_gamma_inputfile = "input_vals_3d_gamma.jl"

"""
  Test weakdifferentiate and that uniform flow goes to zero residual
"""
function test_gamma_sbp(mesh, sbp, eqn, opts)
  facts("----- Testing SummationByParts -----") do
    q = ones(1, mesh.numNodesPerElement, 2)
    res = zeros(q)
    weakdifferentiate!(sbp, 1, q, res, trans=false)

    for i=1:2
      for j=1:mesh.numNodesPerElement
        @fact res[1, j, i] --> roughly(0.0, atol=1e-13)
      end
    end

    # test constant IC -> zero residual
    opts["use_src_term"] = false
    opts["BC1_name"] = "constantBC"
    fill!(eqn.res, 0.0)
    fill!(eqn.q_vec, 2.0)
    fill!(eqn.q, 2.0)
    AdvectionEquationMod.ICConstant(mesh, sbp, eqn, opts, eqn.q_vec)
    disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    AdvectionEquationMod.getBCFunctors(mesh, sbp, eqn, opts)
    evalResidual(mesh, sbp, eqn, opts)

    for i=1:length(eqn.res)
      @fact eqn.res[i] --> roughly(0.0, atol=1e-13)
    end

  end  # end facts block

  return nothing
end

#test_gamma_sbp(mesh, sbp, eqn, opts)
add_func2!(AdvectionTests, test_gamma_sbp, test_gamma_inputfile, [TAG_SHORTTEST])

"""
  Test boundary conditions.  This is not a test function, but is used
  by test functions.
"""
function test_bc(flux_exp, mesh, sbp, eqn, opts)
# test summing boundary condition
 fill!(eqn.bndryflux, 0.0)
 AdvectionEquationMod.evalBoundaryIntegrals(mesh, sbp, eqn, opts)

 flux_in = zero(eltype(eqn.bndryflux))
 flux_out = zero(flux_in)
 for i=1:mesh.numBoundaryFaces
   net_flux_face = zero(flux_in)
   for j=1:mesh.numNodesPerFace
     net_flux_face += mesh.sbpface.wface[j]*eqn.bndryflux[1, j, i]
   end
   if net_flux_face > 0
     flux_in += net_flux_face
   else
     flux_out += net_flux_face
   end
 end

 # flux = area of face * solution * alpha
 @fact abs(flux_in) --> roughly(flux_exp, atol=1e-13)
 @fact abs(flux_out) --> roughly(abs(flux_in), atol=1e-13)
end

"""
  Test the Roe solver used for boundary flux calculations.
"""
function test_gamma_bcsolver(mesh, sbp, eqn, opts)
  facts("----- Testing BCSolver -----") do
    # check that the solver produces the regular flux when qL = qR
    q2 = rand(1, mesh.numNodesPerElement, mesh.numEl)
    q1 = eqn.q
    eqn.q = q2
    fill!(eqn.res, 0.0)
    AdvectionEquationMod.evalVolumeIntegrals(mesh, sbp, eqn, opts)
    nrm = [1., 1, 1]  # arbirary normal vector
    nrm2 = zeros(nrm)
    alphas_xy = [eqn.params.alpha_x, eqn.params.alpha_y, eqn.params.alpha_z]
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        u = eqn.q[1, j, i]
        dxidx = sview(mesh.dxidx, :, :, j, i)
        calcBCNormal(eqn.params, dxidx, nrm, nrm2)
        # calculate the flux directly
        alphas_parametric = dxidx*alphas_xy
        flux_parametric = alphas_parametric*u
        net_flux = sum(flux_parametric.*nrm)

        # calculate boundary flux
        bndryflux_calc = AdvectionEquationMod.RoeSolver(eqn.params, u, u, nrm2)
        # calculate flux from evalVolumeIntegrals
        bndryflux_weak = zero(eltype(eqn.flux_parametric))
        for d=1:3
           bndryflux_weak += nrm[d]*eqn.flux_parametric[1, j, i, d]
        end

        for d=1:3
          @fact eqn.flux_parametric[1, j, i, d] --> roughly(flux_parametric[d], atol=1e-13)
        end

        @fact bndryflux_calc --> roughly(net_flux, atol=1e-13)
      end
    end
    eqn.q = q1
  end  # end facts block

  return nothing
end

#test_gamma_bcsolver(mesh, sbp, eqn, opts)
add_func2!(AdvectionTests, test_gamma_bcsolver, test_gamma_inputfile, [TAG_BC, TAG_FLUX, TAG_SHORTTEST])

"""
  Test boundary flux calculation in all 3 directions
"""
function test_gamma_bcflux(mesh, sbp, eqn, opts)
   facts("----- Testing boundary flux calculation -----") do
     # set alpha_x = 1. all others zero, q = constant, check flux
     fill!(eqn.q, 2.0)
     fill!(eqn.res, 0.0)
     eqn.params.alpha_x = 1.0
     eqn.params.alpha_y = 0.0
     eqn.params.alpha_z = 0.0

     opts["BC1_name"] = "constantBC"
     opts["numBC"] = 1
     AdvectionEquationMod.getBCFunctors(mesh, sbp, eqn, opts)

     test_bc(8.0, mesh, sbp, eqn, opts)
     eqn.params.alpha_x = 0.0
     eqn.params.alpha_y = 1.0
     eqn.params.alpha_z = 0.0

     test_bc(8.0, mesh, sbp, eqn, opts)
     eqn.params.alpha_x = 0.0
     eqn.params.alpha_y = 0.0
     eqn.params.alpha_z = 1.0

     test_bc(8.0, mesh, sbp, eqn, opts)
     eqn.params.alpha_x = 1.0
     eqn.params.alpha_y = 1.0
     eqn.params.alpha_z = 1.0
     test_bc(3*8.0, mesh, sbp, eqn, opts)
   end  # end facts block

  return nothing
end

#test_gamma_bcflux(mesh, sbp, eqn, opts)
add_func2!(AdvectionTests, test_gamma_bcflux, test_gamma_inputfile, [TAG_BC, TAG_SHORTTEST])

"""
  Test face flux calculation
"""
function test_gamma_faceflux(mesh, sbp, eqn, opts)
  facts("----- Testing face flux -----") do
    # the interpolation should be exact for this case
    AdvectionEquationMod.ICp1(mesh, sbp, eqn, opts, eqn.q_vec)
    fill!(eqn.q_face, 0.0)
    fill!(eqn.flux_face, 0.0)
    fill!(eqn.res, 0.0)
    AdvectionEquationMod.evalFaceIntegrals(mesh, sbp, eqn, opts)
    sbpface = mesh.sbpface
    for i=1:mesh.numInterfaces
      iface_i = mesh.interfaces[i]
      el_i = iface_i.elementL
      face_i = iface_i.faceL
      coords_el = mesh.vert_coords[:, :, el_i]
      vertmap = mesh.topo.face_verts[:, face_i]
      coords_faceverts = coords_el[:, vertmap].'
      coords_face = SummationByParts.SymCubatures.calcnodes(sbpface.cub, coords_faceverts)
      for j=1:mesh.numNodesPerFace
        coords_j = coords_face[:, j]
        q_exp = AdvectionEquationMod.calc_p1(eqn.params, coords_j, 0.0)
        q_calc = eqn.q_face[1, 1, j, i]
        q_calc2 = eqn.q_face[1, 2, j, i]
        @fact q_calc --> roughly(q_exp, atol=1e-13)
        @fact q_calc2 --> roughly(q_exp, atol=1e-13)
      end
    end
  end  # end facts block

  return nothing
end

#test_gamma_faceflux(mesh, sbp, eqn, opts)
add_func2!(AdvectionTests, test_gamma_faceflux, test_gamma_inputfile, [TAG_BC, TAG_SHORTTEST])
