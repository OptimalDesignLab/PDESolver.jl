# Test functional Integrate and adjoint for euler equation.

facts("--- Testing Functional Computation On a Boundary ---") do
  ARGS[1] = "input_vals_vortex_adjoint_DG.jl"
  include("../src/solver/euler/startup.jl")  # initialization and construction

  @fact mesh.isDG --> true
  @fact opts["calc_functional"] --> true
  @fact opts["functional_error"] --> true
  @fact opts["functional_name1"] --> "drag"
  @fact opts["analytical_functional_val"] --> roughly(-1/1.4, atol = 1e-13)
  @fact opts["geom_edges_functional1"] --> [3]

  fname = "./functional_error1.dat"
  relative_error = readdlm(fname)

  @fact relative_error[1] --> roughly(0.00018849571, atol = 1e-6)

  rm("./functional_error1.dat") # Delete the file

end  # End do

facts("--- Testing Objective Function Computation On a Boundary ---") do

  include("./input_vals_vortex_adjoint_DG.jl")
  arg_dict["calc_functional"] = false
  arg_dict["objective_function"] = "drag"
  arg_dict["geom_faces_objective"] = [3]

  f = open("input_vals_vortex_objective_computation_DG.jl", "w")
  println(f, "arg_dict = ")
  println(f, arg_dict)
  close(f)

  ARGS[1] = "input_vals_vortex_objective_computation_DG.jl"
  include("../src/solver/euler/startup.jl")

  # Assert basic facts
  @assert mesh.isDG == true
  @assert opts["calc_functional"] == false

  # Check facts
  @fact opts["objective_function"] --> "drag"
  @fact opts["geom_faces_objective"] --> [3]

  drag = EulerEquationMod.OptimizationData{Tsol}(mesh, sbp, opts)
  EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, drag)

  @fact drag.is_objective_fn --> true
  analytical_val = -1/1.4
  drag_err = norm(drag.val - analytical_val)
  @fact drag_err --> roughly(0.0001346397, atol = 1e-6)

end # End facts("--- Testing Objective Function Computation On a Boundary ---")


#=

facts("--- Testing Functional Computation On a Boundary ---") do
  include("./input_vals_vortex_adjoint_DG.jl")
  arg_dict["smb_name"] = "src/mesh_files/gvortex1np2.smb"
  arg_dict["run_type"] = 1
  arg_dict["jac_type"] = 3
  arg_dict["newton_globalize_euler"] = true
  f = open("input_vals_vortex_adjoint_DG_parallel.jl", "w")
  println(f, "arg_dict = ")
  println(f, arg_dict)
  close(f)

  ARGS[1] = "input_vals_vortex_adjoint_DG_parallel.jl"
  include("../src/solver/euler/startup.jl")

  @fact mesh.isDG --> true
  @fact opts["calc_functional"] --> true
  @fact opts["functional_error"] --> true
  @fact opts["functional_name1"] --> "drag"
  @fact opts["analytical_functional_val"] --> roughly(-1/1.4, atol = 1e-13)
  @fact opts["geom_edges_functional1"] --> [3]

  fname = "./functional_error1.dat"
  relative_error = readdlm(fname)

  @fact relative_error[1] --> roughly(0.000177342284, atol = 1e-6)

  rm("./functional_error1.dat") # Delete the file
  rm("./input_vals_vortex_adjoint_DG_parallel.jl")


end  # End do
=#
