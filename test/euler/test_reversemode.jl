# Tests for all reverse mode derivatives
@doc """
Euler Equation -- test_reversemode

Tests all the functions pertaining to reverse mode  for the Euler Equation
module

"""->

function test_reversemode()

  resize!(ARGS, 1)
  ARGS[1] = "input_vals_vortex_adjoint_DG.jl"
  include("../../src/solver/euler/startup.jl")

  facts("--- Testing Pressure derivative in reverse mode ---") do

    press_bar = complex(rand(Float64),0)
    q_bar = zeros(Complex128, mesh.numDofPerNode)
    q_bar_complex = zeros(Complex128, mesh.numDofPerNode)
    pert = complex(0,1e-20)

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

    context("Checking reversemode derivative w.r.t mesh metrics") do
      # Create a random vector
      Tdim = mesh.dim
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
      Tdim = mesh.dim
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
          dxidx = sview(mesh.dxidx_bndry, :, :, j, global_facenum)
          nrm = sview(sbp.facenormal, :, bndry_i.face)
          for k = 1:Tdim
            phys_nrm[k] = dxidx[1,k]*nrm[1] + dxidx[2,k]*nrm[2]
          end # End for k = 1:Tdim
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

  end # End facts("--- Testing Boundary Functional In Reverse Mode ---")

  return nothing
end # End function test_reversemode

add_func1!(EulerTests, test_reversemode, [TAG_REVERSEMODE])
