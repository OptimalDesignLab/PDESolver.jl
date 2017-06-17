# Tests for all reverse mode derivatives
@doc """
Euler Equation -- test_reversemode

Tests all the functions pertaining to reverse mode  for the Euler Equation
module

"""->

function test_reversemode()

  resize!(ARGS, 1)
  ARGS[1] = "input_vals_vortex_reversemode.jl"
  include("../../src/solver/euler/startup.jl")

  facts("--- Testing Pressure derivative in reverse mode ---") do

    press_bar = complex(rand(Float64),0)
    q_bar = zeros(Complex128, mesh.numDofPerNode)
    q_bar_complex = zeros(Complex128, mesh.numDofPerNode)
    pert = complex(0,1e-20)
    # println("type of eqn.q = ", typeof(eqn.q))
    for i = 1:mesh.numEl
      for j = 1:mesh.numNodesPerElement
        q_vals = sview(eqn.q, :,j,i)
        fill!(q_bar, 0.0)
        EulerEquationMod.calcPressure_revq(eqn.params, q_vals, press_bar, q_bar)

        # Check agains complex step
        for k = 1:mesh.numDofPerNode
          q_vals[k] += pert
          press_complex = EulerEquationMod.calcPressure(eqn.params, q_vals)
          press_complex = imag(press_complex)/imag(pert)
          q_bar_complex[k] = press_complex*press_bar
          q_vals[k] -= pert
          error = norm(q_bar_complex[k] - q_bar[k], 2)
          @fact error --> roughly(0.0, atol=1e-10)
        end
      end # End for j = 1:mesh.numNodesPerElement
    end   # End for i = 1:mesh.numEl

  end # facts("--- Testing Pressure derivative in reverse mode ---")

  facts("--- Testing Euler Flux derivative in Reverse mode ---") do

    Tdim = mesh.dim

    context("Checking reversemode derivative w.r.t mesh metrics") do
      # Create a random vector
      F_bar = rand(mesh.numDofPerNode) # For 2D
      flux = zeros(Complex128, mesh.numDofPerNode) # For complex step
      dir_bar = zeros(Float64, Tdim) # In 2D
      dir_bar_complex = zeros(Complex128, Tdim)
      nrm = zeros(Complex128, Tdim)
      pert = complex(0, 1e-20)

      ctr = 0
      for i = 1:mesh.numEl
        for j = 1:mesh.numNodesPerElement
          q_vals = sview(eqn.q, :,j,i)
          aux_vars = sview(eqn.aux_vars, :, j, i)
          for k=1:Tdim  # loop over dimensions
            for p=1:Tdim
              nrm[p] = mesh.dxidx[k, p, j, i]
            end
            fill!(dir_bar,0.0)
            EulerEquationMod.calcEulerFlux_revm(eqn.params, q_vals, aux_vars,
                               nrm, F_bar, dir_bar)

            # Do the complex step in normal
            fill!(dir_bar_complex, 0.0)
            for p = 1:Tdim
              nrm[p] += pert
              fill!(flux, 0.0)
              EulerEquationMod.calcEulerFlux(eqn.params, q_vals, aux_vars, nrm, flux)
              flux[:] = imag(flux[:])/imag(pert)
              nrm[p] -= pert
              dir_bar_complex[p] = dot(F_bar, flux)
              error = norm(dir_bar_complex[p] - dir_bar[p], 2)
              @fact error --> roughly(0.0, atol=1e-10)
            end # End for p = 1:Tdim
          end # End for k=1:Tdim
        end # End for j = 1:mesh.numNodesPerElement
      end   # End for i = 1:mesh.numEl

    end # End context("Checking reversemode derivative w.r.t mesh metrics")

    context("Checking reverse mode derivative w.r.t solution q") do

      # Create a random vector
      F_bar = rand(mesh.numDofPerNode) # For 2D
      flux = zeros(Complex128, mesh.numDofPerNode) # For complex step
      nrm = zeros(Complex128, Tdim)
      q_bar = zeros(mesh.numDofPerNode)
      q_bar_complex = zeros(Complex128, mesh.numDofPerNode)
      pert = complex(0, 1e-20)

      ctr = 0
      for i = 1:mesh.numEl
        for j = 1:mesh.numNodesPerElement
          fill!(q_bar, 0.0)
          q_vals = sview(eqn.q, :,j,i)
          aux_vars = sview(eqn.aux_vars, :, j, i)
          for k=1:Tdim  # loop over dimensions
            for p=1:Tdim
              nrm[p] = mesh.dxidx[k, p, j, i]
            end
            # Reverse mode w.r.t q
            fill!(q_bar, 0.0)
            EulerEquationMod.calcEulerFlux_revq(eqn.params, q_vals, aux_vars,
                                                nrm, F_bar, q_bar)

            # Do complex step in q
            for p = 1:mesh.numDofPerNode
              q_vals[p] += pert
              fill!(flux, 0.0)
              EulerEquationMod.calcEulerFlux(eqn.params, q_vals, aux_vars, nrm, flux)
              flux[:] = imag(flux[:])/imag(pert)
              q_bar_complex[p] = dot(F_bar, flux)
              error = norm(q_bar_complex[p] - q_bar[p],2)
              @fact error --> roughly(0.0, atol=1e-10)
              q_vals[p] -= pert
            end # End for p = 1:mesh.numDofPerNode
          end   # End for k = 1:Tdim
        end     # End for j = 1:mesh.numNodesPerElement
      end       # End for i = 1:mesh.numEl

    end # End context("Cehcking reverse mode derivative w.r.t solution q") do

  end # End facts("--- Testing Euler Flux computation in Reverse mode ---")

  facts("--- Testing Boundary Functional In Reverse Mode ---") do

    # Create functional object
    drag = EulerEquationMod.createObjectiveFunctionalData(mesh, sbp, eqn, opts)
    EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, drag)

    context("Checking Boundary Functional Integrand w.r.t nrm") do

      # Uses conservative variables
      Tdim = mesh.dim
      val_bar = rand(Tdim) # Random seed
      nxny_bar = zeros(Float64, 2)
      nxny_bar_complex = zeros(Complex128, 2)
      pert = complex(0, 1e-20) # Complex step perturbation

      # Test on geometric edge 3 (0 based indexing) with no penetration BC
      start_index = mesh.bndry_offsets[4]
      end_index = mesh.bndry_offsets[5]
      idx_range = start_index:(end_index-1)
      bndry_facenums = sview(mesh.bndryfaces, idx_range) # faces on geometric edge i

      nfaces = length(bndry_facenums)
      phys_nrm = zeros(Complex128, Tdim)
      boundary_integrand = zeros(Complex128, drag.ndof)

      ctr = 0
      for i = 1:nfaces
        bndry_i = bndry_facenums[i]
        global_facenum = idx_range[i]
        for j = 1:mesh.sbpface.numnodes
          q = sview(eqn.q_bndry, :, j, global_facenum)
          aux_vars = sview(eqn.aux_vars_bndry, :, j, global_facenum)
          x = sview(mesh.coords_bndry, :, j, global_facenum)
          copy!(phys_nrm, mesh.nrm_face[:, j, i])
          phys_nrm = sview(mesh.nrm_bndry, :, j, global_facenum)
#          dxidx = sview(mesh.dxidx_bndry, :, :, j, global_facenum)
#          nrm = sview(sbp.facenormal, :, bndry_i.face)
#          for k = 1:Tdim
#            phys_nrm[k] = dxidx[1,k]*nrm[1] + dxidx[2,k]*nrm[2]
#          end # End for k = 1:Tdim
          node_info = Int[1,j,i]
          fill!(boundary_integrand, 0.0)

          # Reverse mode
          fill!(nxny_bar, 0.0)
          EulerEquationMod.calcBoundaryFunctionalIntegrand_revm(eqn.params, q,
                       aux_vars, phys_nrm, node_info, drag, nxny_bar, val_bar)

          # Do Complex step
          for k = 1:Tdim
            phys_nrm[k] += pert
            EulerEquationMod.calcBoundaryFunctionalIntegrand(eqn.params, q, aux_vars, phys_nrm,
                                        node_info, drag, boundary_integrand)
            boundary_integrand[:] = imag(boundary_integrand[:])/imag(pert)
            nxny_bar_complex[k] = dot(boundary_integrand, val_bar)
            phys_nrm[k] -= pert
            error = norm(nxny_bar_complex[k] - nxny_bar[k], 2)
            @fact error --> roughly(0.0, atol=1e-10)
          end # End for k = 1:Tdim
        end # End for j = 1:mesh.sbpface.numnodes
      end   # End for i = 1:nfaces
    end # End context("Checking Boundary Functional Integrand")

    context("Checking Complete boundary functional lift in reverse") do

      fill!(mesh.nrm_bndry_bar, 0.0)
      EulerEquationMod.evalFunctional_revm(mesh, sbp, eqn, opts, drag, "lift")

      # Check lift values against complex step
      pert = complex(0, 1e-20) # Complex step perturbation
      local_functional_val = zeros(Complex128, drag.ndof) # Local processor share
      bndry_force = drag.bndry_force
      fill!(bndry_force, 0.0)
      functional_edges = drag.geom_faces_functional
      phys_nrm = zeros(Complex128, 2)

      # Get bndry_offsets for the functional edge concerned
      for itr = 1:length(functional_edges)
        g_edge_number = functional_edges[itr] # Extract geometric edge number
        itr2 = 0
        for itr2 = 1:mesh.numBC
          if findfirst(mesh.bndry_geo_nums[itr2],g_edge_number) > 0
            break
          end
        end
        start_index = mesh.bndry_offsets[itr2]
        end_index = mesh.bndry_offsets[itr2+1]
        idx_range = start_index:(end_index-1)
        bndry_facenums = sview(mesh.bndryfaces, idx_range) # faces on geometric edge i
        nfaces = length(bndry_facenums)
        for i = 1:nfaces
          bndry_i = bndry_facenums[i]
          global_facenum = idx_range[i]
          for j = 1:mesh.sbpface.numnodes
#            dxidx_bar = sview(mesh.dxidx_bndry_bar, :, :, j, global_facenum)
#            dxidx = sview(mesh.dxidx_bndry, :, :, j, global_facenum)
             nrm = sview(mesh.nrm_bndry, :, j, global_facenum)
             nrm_bar = sview(mesh.nrm_bndry_bar, :, j, global_facenum)
            # Perturb dxidx_bndry
            for k = 1:length(nrm)
              nrm[k] += pert
              EulerEquationMod.calcBndryFunctional(mesh, sbp, eqn, opts, drag)
              dlift_dm = imag(drag.lift_val)/imag(pert)
              error = norm(nrm_bar[k] - dlift_dm, 2)
              @fact error --> roughly(0.0, atol=1e-12)
              nrm[k] -= pert
            end
          end  # End for j = 1:mesh.sbpface.numnodes
        end    # End for i = 1:nfaces
      end # End for itr = 1:length(functional_edges)


    end # End context("Checking Complete boundary functional lift in reverse")

    context("Checking Complete boundary functional drag in reverse") do
      fill!(mesh.nrm_bndry_bar, 0.0)
#      fill!(mesh.dxidx_bndry_bar, 0.0) # Clear out for this test
      EulerEquationMod.evalFunctional_revm(mesh, sbp, eqn, opts, drag, "drag")

      # Check lift values against complex step
      pert = complex(0, 1e-20) # Complex step perturbation
      local_functional_val = zeros(Complex128, drag.ndof) # Local processor share
      bndry_force = drag.bndry_force
      fill!(bndry_force, 0.0)
      functional_edges = drag.geom_faces_functional
      phys_nrm = zeros(Complex128, 2)

      # Get bndry_offsets for the functional edge concerned
      for itr = 1:length(functional_edges)
        g_edge_number = functional_edges[itr] # Extract geometric edge number
        itr2 = 0
        for itr2 = 1:mesh.numBC
          if findfirst(mesh.bndry_geo_nums[itr2],g_edge_number) > 0
            break
          end
        end
        start_index = mesh.bndry_offsets[itr2]
        end_index = mesh.bndry_offsets[itr2+1]
        idx_range = start_index:(end_index-1)
        bndry_facenums = sview(mesh.bndryfaces, idx_range) # faces on geometric edge i
        nfaces = length(bndry_facenums)
        for i = 1:nfaces
          bndry_i = bndry_facenums[i]
          global_facenum = idx_range[i]
          for j = 1:mesh.sbpface.numnodes
#            dxidx_bar = sview(mesh.dxidx_bndry_bar, :, :, j, global_facenum)
#            dxidx = sview(mesh.dxidx_bndry, :, :, j, global_facenum)
             nrm = sview(mesh.nrm_bndry, :, j, global_facenum)
             nrm_bar = sview(mesh.nrm_bndry_bar, :, j, global_facenum)
            # Perturb dxidx_bndry
            for k = 1:length(nrm)
              nrm[k] += pert
              EulerEquationMod.calcBndryFunctional(mesh, sbp, eqn, opts, drag)
              ddrag_dm = imag(drag.drag_val)/imag(pert)
              error = norm(nrm_bar[k] - ddrag_dm, 2)
              @fact error --> roughly(0.0, atol=1e-12)
              nrm[k] -= pert
            end
          end  # End for j = 1:mesh.sbpface.numnodes
        end    # End for i = 1:nfaces
      end # End for itr = 1:length(functional_edges)


    end # End context("Checking Complete boundary functional drag in reverse")
  end # End facts("--- Testing Boundary Functional In Reverse Mode ---")


  facts("--- Testing SAT terms in Reverse Mode ---") do

    q = Complex128[2.0043681897362733,0.040161434857338515,-1.3465473815098652,2.241635694978014]
    nrm2 = Complex128[-0.07115741913664266,-0.005089279059529922]
    vel = Complex128[0.020036954818477136,-0.671806401840292]
    H = 1.4753802296828 + 0im
    u = q[2]/q[1]
    v = q[3]/q[1]
    vel_bar = zeros(Complex128, mesh.dim)
    psi = rand(Float64, 4) + zeros(Complex128, 4)
    dq_bar = zeros(Complex128, 4)
    nrm2_bar = zeros(Complex128, mesh.dim)
    sat = zeros(Complex128, 4)
    pert = complex(0, 1e-20) # Complex step perturbation
    EulerEquationMod.calcSAT_revm(eqn.params, nrm2, q, [u,v], H, psi,
            nrm2_bar)
    for k = 1:length(nrm2)
      nrm2[k] += pert
      EulerEquationMod.calcSAT(eqn.params, nrm2, q, sat, u, v, H)
      dSat = imag(sat[:])/imag(pert)
      complex_valbar_SAT = dot(psi, dSat)
      nrm2[k] -= pert
      error = norm(complex_valbar_SAT - nrm2_bar[k], 2)
      @fact error --> roughly(0.0, atol=1e-10)
    end # End for k = 1:length(nrm2)

  end # End facts("--- Testing SAT terms in Reverse Mode ---")

  facts("--- Testing Roe Solver in Reverse Mode ---") do

    fill!(mesh.nrm_bndry_bar, 0.0)
    # EulerEquationMod.dataPrep(mesh, sbp, eqn, opts)
    params = eqn.params
    Tdim = mesh.dim
    val_bar = rand(Tdim) # Random seed
    pert = complex(0, 1e-20) # Complex step perturbation
    complex_flux = zeros(Complex128, 4)
    psi = rand(Float64, 4) + zeros(Complex128, 4)
    qg = ones(Complex128, 4)

    # Test on geometric edge 3 (0 based indexing) with no penetration BC
    start_index = mesh.bndry_offsets[4]
    end_index = mesh.bndry_offsets[5]
    idx_range = start_index:(end_index-1)
    bndry_facenums = sview(mesh.bndryfaces, idx_range) # faces on geometric edge i

    nfaces = length(bndry_facenums)
    phys_nrm = zeros(Complex128, Tdim)
#    nrm = zeros(, Tdim)

    for i = 1:nfaces
      bndry_i = bndry_facenums[i]
      global_facenum = idx_range[i]
      for j = 1:mesh.sbpface.numnodes
        q = sview(eqn.q_bndry, :, j, global_facenum)
        aux_vars = sview(eqn.aux_vars_bndry, :, j, global_facenum)
        x = sview(mesh.coords_bndry, :, j, global_facenum)
#        dxidx = sview(mesh.dxidx_bndry, :, :, j, global_facenum)
#        nrm[:] = sbp.facenormal[:, bndry_i.face]
#        nrm = sview(sbp.facenormal, :, bndry_i.face)
        # println("nrm = $(real(nrm))")
#        dxidx_bar = zeros(Complex128, 2,2)
        nrm = sview(mesh.nrm_bndry, :, j, global_facenum)
        nrm_bar = sview(mesh.nrm_bndry_bar, :, j, global_facenum)
        EulerEquationMod.RoeSolver_revm(params, q, qg, aux_vars, nrm, psi, nrm_bar)
        for k = 1:length(nrm)
          nrm[k] += pert
          EulerEquationMod.RoeSolver(params, q, qg, aux_vars, nrm, complex_flux)
          dRoeFlux = imag(complex_flux)/imag(pert)
          complex_psi_dRoeFlux = dot(psi, dRoeFlux)
          error = norm(nrm_bar[k] - complex_psi_dRoeFlux, 2)
          @fact error --> roughly(0.0, atol = 1e-12)
          nrm[k] -= pert
        end # End for k = 1:Tdim
      end
    end

  end # End facts ("--- Testing Roe Solver in Reverse Mode ---")

  facts("--- Testing reverse mode for BC functors ---") do

    # Populate eqn.bndryflux_bar
    for i = 1:length(eqn.bndryflux_bar)
      eqn.bndryflux_bar[i] = 1.0 + 0.0im # randn() + zero(Complex128) # 1.0 + 0.0im
    end
    EulerEquationMod.init_revm(mesh, sbp, eqn, opts)
    fill!(mesh.nrm_bndry_bar, 0.0)

    context("Checking reverse mode for noPenetrationBC") do

      # EulerEquationMod.dataPrep(mesh, sbp, eqn, opts)
      functor_rev = mesh.bndry_funcs_revm[4]
      functor = mesh.bndry_funcs[4]
#      fill!(mesh.dxidx_bndry_bar, 0.0)
      fill!(mesh.nrm_bndry_bar, 0.0)

      start_index = mesh.bndry_offsets[4]
      end_index = mesh.bndry_offsets[5]
      idx_range = start_index:(end_index-1)
      bndry_facenums = sview(mesh.bndryfaces, idx_range) # faces on geometric edge i
      bndryflux_bar = sview(eqn.bndryflux_bar, :, :, idx_range)

      EulerEquationMod.calcBoundaryFlux_revm(mesh, sbp, eqn, functor_rev,
                                      idx_range, bndry_facenums, bndryflux_bar)

      pert = complex(0, 1e-20) # Complex step perturbation
      nfaces = length(bndry_facenums)
#      nrm = zeros(Complex128, 2)
      dBndryfluxdm = zeros(Complex128, 4)
      for i = 1:nfaces
        bndry_i = bndry_facenums[i]
        global_facenum = idx_range[i]
        for j = 1:mesh.sbpface.numnodes
          q = sview(eqn.q_bndry, :, j, global_facenum)
          aux_vars = sview(eqn.aux_vars_bndry, :, j, global_facenum)
          x = sview(mesh.coords_bndry, :, j, global_facenum)
#          dxidx = sview(mesh.dxidx_bndry, :, :, j, global_facenum)
#          nrm[:] = sbp.facenormal[:, bndry_i.face]
#          nrm = sview(sbp.facenormal, :, bndry_i.face)
          nrm = sview(mesh.nrm_bndry, :, j, global_facenum)
          nrm_bar = sview(mesh.nrm_bndry_bar, :, j, global_facenum)
          tmpflux = zeros(Complex128, 4)
          bndryflux_bar_i = sview(bndryflux_bar, :, j, i)
          # println("bndryflux_bar_i = $bndryflux_bar_i")
#          dxidx_bndry_bar = sview(mesh.dxidx_bndry_bar, :, :, j, global_facenum)
          for k = 1:length(nrm)
            nrm[k] += pert
            functor(q, aux_vars, x, nrm, tmpflux, eqn.params)
            tmpflux[:] = imag(tmpflux[:])/imag(pert)
            dot_product = dot(real(bndryflux_bar_i), real(tmpflux))
            error = norm(nrm_bar[k] - dot_product, 2)
            # println("error = $error")
            @fact error --> roughly(0.0, atol = 1e-12)
            # println("dot_product = $dot_product, dxidx_bndry_bar[$k] = $(real(dxidx_bndry_bar[k]))")
            nrm[k] -= pert
          end # End for k = 1:Tdim
        end
      end

    end # End context("Checking noPenetrationBC_revm")

    context("Checking reverse mode for isentropicVortexBC") do

      # EulerEquationMod.dataPrep(mesh, sbp, eqn, opts)
      functor_rev = mesh.bndry_funcs_revm[2]
      functor = mesh.bndry_funcs[2]

      start_index = mesh.bndry_offsets[2]
      end_index = mesh.bndry_offsets[3]
      idx_range = start_index:(end_index-1)
      bndry_facenums = sview(mesh.bndryfaces, idx_range) # faces on geometric edge i
      bndryflux_bar = sview(eqn.bndryflux_bar, :, :, idx_range)

      EulerEquationMod.calcBoundaryFlux_revm(mesh, sbp, eqn, functor_rev,
                                      idx_range, bndry_facenums, bndryflux_bar)
      # println("bndryflux_bar[:,1,1] = $(bndryflux_bar[:,1,1])")
      pert = complex(0, 1e-20) # Complex step perturbation
      nfaces = length(bndry_facenums)
#      nrm = zeros(Complex128, 2)
      dBndryfluxdm = zeros(Complex128, 4)
      for i = 1:nfaces
        bndry_i = bndry_facenums[i]
        global_facenum = idx_range[i]
        for j = 1:mesh.sbpface.numnodes
          q = sview(eqn.q_bndry, :, j, global_facenum)
          aux_vars = sview(eqn.aux_vars_bndry, :, j, global_facenum)
          x = sview(mesh.coords_bndry, :, j, global_facenum)
#          dxidx = sview(mesh.dxidx_bndry, :, :, j, global_facenum)
#          nrm[:] = sbp.facenormal[:, bndry_i.face]
#          nrm = sview(sbp.facenormal, :, bndry_i.face)
          nrm = sview(mesh.nrm_bndry, :, j, global_facenum)
          nrm_bar = sview(mesh.nrm_bndry_bar, :, j, global_facenum)
          tmpflux = zeros(Complex128, 4)
          bndryflux_bar_i = sview(bndryflux_bar, :, j, i)
          # println("bndryflux_bar_i = $bndryflux_bar_i")
#          dxidx_bndry_bar = sview(mesh.dxidx_bndry_bar, :, :, j, global_facenum)
          for k = 1:length(nrm)
            nrm[k] += pert
            functor(q, aux_vars, x, nrm, tmpflux, eqn.params)
            tmpflux[:] = imag(tmpflux[:])/imag(pert)
            dot_product = dot(real(bndryflux_bar_i), real(tmpflux))
            error = norm(nrm_bar[k] - dot_product, 2)
            # println("error = $error")
            @fact error --> roughly(0.0, atol = 1e-12)
            # println("dot_product = $dot_product, dxidx_bndry_bar[$k] = $(real(dxidx_bndry_bar[k]))")
            nrm[k] -= pert
          end # End for k = 1:Tdim
        end
      end

    end # End context("Checking reverse mode for isentropicVortexBC")

  end # Endfacts("--- Testing reverse mode for BC functors ---")

  facts("--- Testing reverse mode for face fluxes w.r.t mesh metrics ---") do
    EulerEquationMod.init_revm(mesh, sbp, eqn, opts)
    fill!(mesh.nrm_face_bar, 0.0)
    for i = 1:length(eqn.flux_face_bar)
      eqn.flux_face_bar[:] = 1.0 + 0.0im # randn() + 0.0im
    end

    @assert opts["face_integral_type"] == 1
    EulerEquationMod.calcFaceFlux_revm(mesh, sbp, eqn, eqn.flux_func_bar,
                                       mesh.interfaces, eqn.flux_face_bar)

    # Check against complex step
    pert = complex(0, 1e-20) # Complex step perturbation
    tmpflux = zeros(Complex128, 4)

    functor = eqn.flux_func
    nfaces = length(mesh.interfaces)
#    nrm = zeros(Complex128, size(sbp.facenormal,1))
    for i=1:nfaces  # loop over faces
      interface_i = mesh.interfaces[i]
      for j = 1:mesh.numNodesPerFace
        eL = interface_i.elementL
        fL = interface_i.faceL

        # get components
        qL = sview(eqn.q_face, :, 1, j, i)
        qR = sview(eqn.q_face, :, 2, j, i)
#        dxidx = sview(mesh.dxidx_face, :, :, j, i)
#        dxidx_bar = sview(mesh.dxidx_face_bar, :, :, j, i)
        nrm = sview(mesh.nrm_face, :, j, i)
        nrm_bar = sview(mesh.nrm_face_bar, :, j, i)
        aux_vars = sview(eqn.aux_vars_face, :, j, i)
#        nrm = sview(sbp.facenormal, :, fL)
#        nrm[:] = sbp.facenormal[:,fL]

        flux_j_bar = sview(eqn.flux_face_bar, :, j, i)
        # println("dxidx_bar = $(real(dxidx_bar))")
        for k = 1:length(nrm)
          nrm[k] += pert
          functor(eqn.params, qL, qR, aux_vars, nrm, tmpflux)
          tmpflux[:] = imag(tmpflux[:])/imag(pert)
          dot_product = dot(real(flux_j_bar),real(tmpflux))
          error = norm(nrm_bar[k] - dot_product, 2)
          @fact error --> roughly(0.0, atol = 1e-12)
          nrm[k] -= pert
        end # End for k = 1:length(dxidx)
      end
    end

  end # End facts("--- Testing reverse mode for face fluxes w.r.t mesh metrics ---")


  facts("--- Testing evalrevm_transposeproduct ---") do

    adjoint_vec = zeros(Complex128, mesh.numDof)
    for i = 1:mesh.numDof
      adjoint_vec[i] = randn() + 0.0im
    end

    fill!(mesh.nrm_face_bar, 0.0)
    fill!(mesh.nrm_bndry_bar, 0.0)
    EulerEquationMod.init_revm(mesh, sbp, eqn, opts) # initialize reversemode functors
    EulerEquationMod.evalrevm_transposeproduct(mesh, sbp, eqn, opts, adjoint_vec)

    # Copy flux_bar values for testing against complex step
    vol_flux_bar = deepcopy(eqn.flux_parametric_bar)
    f = open("vol_flux_bar.dat", "w")
    println(f, vol_flux_bar)
    close(f)

    bndryflux_bar_copy = deepcopy(eqn.bndryflux_bar)
    f = open("bndryflux_bar.dat", "w")
    println(f, bndryflux_bar_copy)
    close(f)

    flux_face_bar_copy = deepcopy(eqn.flux_face_bar)
    f = open("flux_face_bar.dat", "w")
    println(f, flux_face_bar_copy)
    close(f)

    pert = complex(0, 1e-20)
    dir_bar_complex = zeros(Complex128, 2)
    flux = zeros(Complex128, 4)
    Tdim = 2
    nrm = zeros(Complex128, 2)

    # Check reverse mode for volume flux bar components
    for i = 1:mesh.numEl
      for j = 1:mesh.numNodesPerElement
        q_vals = sview(eqn.q, :,j,i)
        aux_vars = sview(eqn.aux_vars, :, j, i)
        for k=1:Tdim  # loop over dimensions
          for p=1:Tdim
            nrm[p] = mesh.dxidx[k, p, j, i]
          end
          dxidx_bar = sview(mesh.dxidx_bar, k, :, j, i)
          F_bar = sview(vol_flux_bar, :, j, i, k)

          # Do the complex step in normal
          fill!(dir_bar_complex, 0.0)
          for p = 1:Tdim
            nrm[p] += pert
            fill!(flux, 0.0)
            EulerEquationMod.calcEulerFlux(eqn.params, q_vals, aux_vars, nrm, flux)
            flux[:] = imag(flux[:])/imag(pert)
            nrm[p] -= pert
            dir_bar_complex[p] = dot(F_bar, flux)
            error = norm(dir_bar_complex[p] - dxidx_bar[p], 2)
            @fact error --> roughly(0.0, atol = 1e-12)
          end # End for p = 1:Tdim
        end # End for k=1:Tdim
      end # End for j = 1:mesh.numNodesPerElement
    end   # End for i = 1:mesh.numEl

    # Check reverse mode for face fluxes
    @assert opts["face_integral_type"] == 1
    tmpflux = zeros(Complex128, 4)

    functor = eqn.flux_func
    nfaces = length(mesh.interfaces)
#    nrm = zeros(Complex128, size(sbp.facenormal,1))
    for i=1:nfaces  # loop over faces
      interface_i = mesh.interfaces[i]
      for j = 1:mesh.numNodesPerFace
        eL = interface_i.elementL
        fL = interface_i.faceL

        # get components
        qL = sview(eqn.q_face, :, 1, j, i)
        qR = sview(eqn.q_face, :, 2, j, i)
#        dxidx = sview(mesh.dxidx_face, :, :, j, i)
#        dxidx_bar = sview(mesh.dxidx_face_bar, :, :, j, i)
        nrm = sview(mesh.nrm_face, :, j, i)
        nrm_bar = sview(mesh.nrm_face_bar, :, j, i)
        aux_vars = sview(eqn.aux_vars_face, :, j, i)
 #       nrm[:] = sbp.facenormal[:,fL]
#        nrm = sview(sbp.facenormal, :, fL)

        flux_j_bar = sview(flux_face_bar_copy, :, j, i)
        for k = 1:length(nrm)
          nrm[k] += pert
          functor(eqn.params, qL, qR, aux_vars, nrm, tmpflux)
          tmpflux[:] = imag(tmpflux[:])/imag(pert)
          dot_product = dot(real(flux_j_bar),real(tmpflux))
          error = norm(nrm_bar[k] - dot_product, 2)
          @fact error --> roughly(0.0, atol = 1e-12)
          nrm[k] -= pert
        end # End for k = 1:length(dxidx)
      end
    end

    # Check reverse mode for boundary fluxes
#    fill!(mesh.nrm_bndry_bar, 0.0)
    for ctr=1:mesh.numBC
      functor_i = mesh.bndry_funcs[ctr]
      start_index = mesh.bndry_offsets[ctr]
      end_index = mesh.bndry_offsets[ctr+1]
      idx_range = start_index:end_index
      bndry_facenums_i = sview(mesh.bndryfaces, start_index:(end_index - 1))
      bndryflux_bar = sview(bndryflux_bar_copy, :, :, start_index:(end_index - 1))

      nfaces = length(bndry_facenums_i)
#      nrm = zeros(Complex128, 2)
      dBndryfluxdm = zeros(Complex128, 4)
      for i = 1:nfaces
        bndry_i = bndry_facenums_i[i]
        global_facenum = idx_range[i]
        for j = 1:mesh.sbpface.numnodes
          q = sview(eqn.q_bndry, :, j, global_facenum)
          aux_vars = sview(eqn.aux_vars_bndry, :, j, global_facenum)
          x = sview(mesh.coords_bndry, :, j, global_facenum)
#          dxidx = sview(mesh.dxidx_bndry, :, :, j, global_facenum)
 #         nrm[:] = sbp.facenormal[:, bndry_i.face]
#          nrm = sview(sbp.facenormal, :, bndry_i.face)
          nrm = sview(mesh.nrm_bndry, :, j, global_facenum)
          nrm_bar = sview(mesh.nrm_bndry_bar, :, j, global_facenum)
          tmpflux = zeros(Complex128, 4)
          bndryflux_bar_i = sview(bndryflux_bar, :, j, i)
#          dxidx_bndry_bar = sview(mesh.dxidx_bndry_bar, :, :, j, global_facenum)
          for k = 1:length(nrm)
            nrm[k] += pert
            functor_i(q, aux_vars, x, nrm, tmpflux, eqn.params)
            tmpflux[:] = imag(tmpflux[:])/imag(pert)
            dot_product = dot(real(bndryflux_bar_i), real(tmpflux))
            error = norm(nrm_bar[k] - dot_product, 2)
            @fact error --> roughly(0.0, atol = 1e-12)
            nrm[k] -= pert
          end # End for k = 1:Tdim
        end
      end
    end

  end # End facts("--- Testing evalrevm_transposeproduct ---")


  return nothing
end # End function test_reversemode

add_func1!(EulerTests, test_reversemode, [TAG_REVERSEMODE, TAG_SHORTTEST])
