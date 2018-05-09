# Tests for all forward mode derivatives
@doc """
Euler Equation -- test_forwardmode

Tests all the functions pertaining to forward mode  for the Euler Equation
module

"""->

function test_forwardmode()

  resize!(ARGS, 1)
  ARGS[1] = "input_vals_vortex_forwardmode.jl"
  include("../../src/solver/euler/startup.jl")
  Tmsh, Tsol, Tres = EulerEquationMod.getTypeParameters(mesh, eqn)

  facts("--- Testing Boundary Functional In Forward Mode ---") do

    # Create functional object
    drag = EulerEquationMod.createObjectiveFunctionalData(mesh, sbp, eqn, opts)
    EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, drag)

    context("Checking Boundary Functional Integrand w.r.t q") do

      # Uses conservative variables
      Tdim = mesh.dim
      val_bar = rand(Tdim) +  zeros(Complex128, Tdim)# Random seed
      pert = complex(0, 1e-20) # Complex step perturbation

      # Test on geometric edge 3 (0 based indexing) with no penetration BC
      start_index = mesh.bndry_offsets[4]
      end_index = mesh.bndry_offsets[5]
      idx_range = start_index:(end_index-1)


      bndry_facenums = sview(mesh.bndryfaces, idx_range) # faces on geometric edge i

      nfaces = length(bndry_facenums)
      phys_nrm = zeros(Complex128, Tdim)
      boundary_integrand = zeros(Complex128, drag.ndof)
      # the fwd mode quantity we're chking
      boundary_integrand_diff = zeros(Complex128, length(boundary_integrand), mesh.numDofPerNode)   # size (2,4)


      for i = 1:nfaces
        bndry_i = bndry_facenums[i]
        global_facenum = idx_range[i]

        q = zeros(Tsol, mesh.numDofPerNode)   # new
        for j = 1:mesh.sbpface.numnodes
          # Yes, we need to be looping over bndry nodes, not volume nodes
          q = sview(eqn.q_bndry, :, j, global_facenum)                # q at this node. all dof. should be size (4,)
          aux_vars = sview(eqn.aux_vars_bndry, :, j, global_facenum)
          x = sview(mesh.coords_bndry, :, j, global_facenum)
          phys_nrm = ro_sview(mesh.nrm_bndry, :, j, global_facenum)
          node_info = Int[1,j,i]        # only one BC, so first element 'itr' is 1
          fill!(boundary_integrand_diff, 0.0)

          # Forward mode
          EulerEquationMod.calcBoundaryFunctionalIntegrand_diff(eqn.params, q, aux_vars,
                                                                phys_nrm, node_info,
                                                                drag, boundary_integrand_diff)

          # Do Complex step on q
          for k = 1:mesh.numDofPerNode
            fill!(boundary_integrand, 0.0)
            q[k] += pert
            EulerEquationMod.calcBoundaryFunctionalIntegrand(eqn.params, q, aux_vars,
                                                             phys_nrm, node_info,
                                                             drag, boundary_integrand)
            boundary_integrand[:] = imag(boundary_integrand[:])/imag(pert)
            q[k] -= pert
            # checking the col of b_i_d that corresponds to this CS'd boundary integrand
            error = norm(boundary_integrand_diff[:,k] - boundary_integrand[:], 2)
            @fact error --> roughly(0.0, atol=1e-10)
          end # End for k = 1:length(q)
        end # End for j = 1:mesh.sbpface.numnodes
      end   # End for i = 1:nfaces
    end # End context("Checking Boundary Functional Integrand w.r.t q")
  end # End facts("--- Testing Boundary Functional In Reverse Mode ---")

  return nothing
end # End function test_forwardmode

add_func1!(EulerTests, test_forwardmode, [TAG_REVERSEMODE, TAG_SHORTTEST])
