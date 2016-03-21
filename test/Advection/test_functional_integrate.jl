# test boundary forces

ARGS[1] = "input_vals_functional_CG.jl"
include("../../src/solver/advection/startup_advection.jl")  # initialization and construction

facts("--- Testing Boundary Functional Computation on CG mesh ---") do
  geometric_edge_number = 1
  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
  functional_val = AdvectionEquationMod.calcBndryfunctional(mesh, sbp, eqn,
  	               opts, geometric_edge_number)
  @fact functional_val --> roughly(-47.33907493660998, atol=1e-12)
end

ARGS[1] = "input_vals_functional_DG.jl"
include("../../src/solver/advection/startup_advection.jl")  # initialization and construction

facts("--- Testing Boundary Functional Computation on DG mesh ---") do
  @fact mesh.isDG --> true
  geometric_edge_number = 1
  boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
  functional_val = AdvectionEquationMod.calcBndryfunctional(mesh, sbp, eqn, opts, 
          geometric_edge_number)
  @fact functional_val --> roughly(-47.18786881412863, atol=1e-12)
end
