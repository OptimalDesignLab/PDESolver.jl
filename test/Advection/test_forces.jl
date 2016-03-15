# test boundary forces

ARGS[1] = "input_vals_force.jl"
include("../../src/solver/advection/startup_advection.jl")  # initialization and construction\

facts("--- Testing Boundary Force Computation ---") do
  geometric_edge_number = 1
  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
  force = AdvectionEquationMod.calcBndryforces(mesh, sbp, eqn, opts, 
          geometric_edge_number)
  @fact force --> roughly(-47.33907493660998, atol=1e-12)
end