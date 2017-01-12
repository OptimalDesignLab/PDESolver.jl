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

facts("--- Tesing adjoint computation on the boundary for DG Meshes---") do

  resize!(ARGS, 1)
  ARGS[1] = "input_vals_airfoil.jl"
  include("../src/solver/euler/startup.jl")

  context("--- Checking partial dR/dq Calculation") do
    
    # Copy all the original values
    disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    orig_q_vec = deepcopy(eqn.q_vec)
    original_res_vec = copy(eqn.res_vec)
    
    rand_vec = rand(length(eqn.q_vec))
    fill!(eqn.res, 0.0)
    fill!(eqn.res_vec, 0.0)
    res_norm = NonlinearSolvers.calcResidual(mesh, sbp, eqn, opts, evalEuler)
    boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
    res_jac, jacData = EulerEquationMod.calcResidualJacobian(mesh, sbp, eqn, opts)
    contract_vec = res_jac*rand_vec

    # Check against FD
    copy!(eqn.q_vec, orig_q_vec)
    for i = 1:length(q_vec)
      eqn.q_vec[i] += 1e-6*rand_vec[i]
    end
    disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    fill!(eqn.res, 0.0)
    fill!(eqn.res_vec, 0.0)
    res_norm = NonlinearSolvers.calcResidual(mesh, sbp, eqn, opts, evalEuler)
    assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
    partialRpartialu = (eqn.res_vec - original_res_vec)/1e-6

    for i = 1:length(partialRpartialu)
      @fact abs(real(contract_vec[i] - partialRpartialu[i])) --> roughly(0.0, atol = 1e-6)
      # println(f,real(contract_vec[i] - partialRpartialu[i]))
    end

    for i = 1:length(q_vec)
      eqn.q_vec[i] = orig_q_vec[i]
    end

  end # End Checking dR/dq
  
  context("--- Checking partial dJ/dq Calculation") do
    
    func_deriv_arr = zeros(eqn.q)
    func_deriv = zeros(eqn.q_vec)
    functional_edges = opts["geom_faces_objective"]
    functional_name = EulerEquationMod.FunctionalDict[opts["objective_function"]]

    boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
    EulerEquationMod.calcFunctionalDeriv(mesh, sbp, eqn, opts, functional_name, functional_edges,
                          objective, func_deriv_arr)  # populate df_dq_bndry
    assembleSolution(mesh, sbp, eqn, opts, func_deriv_arr, func_deriv)

    rand_vec = rand(length(eqn.q_vec))
    contract_val = dot(rand_vec,func_deriv)

    # Check with finite difference
    eqn.q_vec += 1e-6*rand_vec
    disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, objective)
    dJdu_fd = (objective.val-orig_Ju)/1e-6

    @fact norm(dJdu_fd - contract_val, 2) --> roughly(0.0, atol = 1e-8)

  end # End dJ/dq

end # End do
=#
