# Test functional Integrate and adjoint for euler equation.

function test_adjoint(mesh, sbp, eqn, opts)
  facts("--- Testing Functional Computation On a Boundary ----") do

    @fact mesh.isDG --> true
    @fact opts["functional_name1"] --> "drag"
    @fact opts["analytical_functional_val"] --> roughly(-1/1.4, atol = 1e-13)
    @fact opts["geom_edges_functional1"] --> [4]

    fname = "./functional_error1.dat"
    relative_error = readdlm(fname)

    @fact relative_error[1] --> roughly(0.00018849571, atol = 1e-6)

    rm("./functional_error1.dat") # Delete the file

  end  # End do

  return nothing
end  # end function

#test_adjoint(mesh, sbp, eqn, opts)
add_func2!(EulerTests, test_adjoint, "input_vals_vortex_adjoint_DG.jl")
