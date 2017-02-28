# Tests for all reverse mode derivatives
@doc """
Euler Equation -- test_reversemode

Tests all the functions pertaining to reverse mode  for the Euler Equation
module

"""->

function test_reversemode()

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

    err_ctr = 0

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
          # println("dir_bar = $(dir_bar)")

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

    println("error count = $err_ctr")
  end # End facts("--- Testing Euler Flux computation in Reverse mode ---")

  return nothing
end

# test_reverseMode()
add_func1!(EulerTests, test_reversemode, [TAG_REVERSEMODE])
