# Tests for all reverse mode derivatives
@doc """
Euler Equation -- test_reversemode

Tests all the functions pertaining to reverse mode  for the Euler Equation
module

"""->

function test_reversemode()
#=
  facts("--- Testing Euler Flux derivative in Reverse mode ---") do
    resize!(ARGS, 1)
    ARGS[1] = "input_vals_vortex_adjoint_DG.jl"
    include("../../src/solver/euler/startup.jl")

    # Create a random vector
    Tdim = mesh.dim
    F_bar = rand(Tdim+2) # For 2D
    flux = zeros(Complex128, Tdim+2) # For complex step
    dir_bar = zeros(Float64, Tdim) # In 2D
    dir_bar_complex = zeros(Complex128, Tdim)
    nrm = zeros(Complex128, Tdim)
    q_bar = zeros(Tdim + 2)
    pert = complex(0, 1e-20)


    for i = 1:mesh.numEl
      for j = 1:mesh.numNodesPerElement
        fill!(q_bar, 0.0)
        q_vals = sview(eqn.q, :,j,i)
        aux_vars = sview(eqn.aux_vars, :, j, i)
        for k=1:Tdim  # loop over dimensions
          for p=1:Tdim
            nrm[p] = mesh.dxidx[k, p, j, i]
          end
          # println("size of q_vals = $(size(q_vals))")
          fill!(dir_bar,0.0)
          EulerEquationMod.calcEulerFlux_revm(eqn.params, q_vals, aux_vars,
                             nrm, F_bar, q_bar, dir_bar)

          # Now do the complex step
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

  end # End facts("--- Testing Euler Flux computation in Reverse mode ---")
=#
  facts("--- Testing Boundary Functional In Reverse Mode ---") do

    resize!(ARGS, 1)
    ARGS[1] = "input_vals_vortex_adjoint_DG.jl"
    include("../../src/solver/euler/startup.jl")
    drag = EulerEquationMod.createObjectiveFunctionalData(mesh, sbp, eqn, opts)
    EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, drag)

    context("Checking Boundary Functional Integrand w.r.t nrm") do

      # Uses conservative variables
      Tdim = mesh.dim
      val_bar = rand(Tdim) # Random seed
      nxny_bar = zeros(Float64, 2)
      pert = complex(0, 1e-20) # Complex step perturbation

      # Test on geometric edge 3 (0 based indexing) with no penetration BC
      start_index = mesh.bndry_offsets[4]
      end_index = mesh.bndry_offsets[5]
      idx_range = start_index:(end_index-1)
      bndry_facenums = sview(mesh.bndryfaces, idx_range) # faces on geometric edge i

      nfaces = length(bndry_facenums)
      boundary_integrand = zeros(Complex128, drag.ndof, mesh.sbpface.numnodes, nfaces)
      phys_nrm = zeros(Complex128, Tdim)

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
            # nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
            # ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]
            phys_nrm[k] = dxidx[1,k]*nrm[1] + dxidx[2,k]*nrm[2]
          end # End for k = 1:Tdim
          node_info = Int[1,j,i]
          b_integrand_ji = sview(boundary_integrand,:,j,i)

          # Reverse mode
          fill!(nxny_bar, 0.0)
          EulerEquationMod.calcBoundaryFunctionalIntegrand_revm(eqn.params, q, aux_vars, phys_nrm,
                                             node_info, drag, nxny_bar, val_bar)

          # Do Complex step
          for k = 1:Tdim
            phys_nrm(k)
          end # End for k = 1:Tdim

        end # End for j = 1:mesh.sbpface.numnodes
      end   # End for i = 1:nfaces

    end # End context("Checking Boundary Functional Integrand")


  end # End facts("--- Testing Boundary Functional In Reverse Mode ---")

  return nothing
end

# test_reverseMode()
add_func1!(EulerTests, test_reversemode, [TAG_REVERSEMODE])
