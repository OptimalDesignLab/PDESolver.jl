#Test functional Integrate and adjoint for euler equation.
function test_adjoint()

  facts("--- Testing Functional Computation On a Boundary ---") do

    ARGS[1] = "input_vals_vortex_adjoint_DG.jl"
    mesh, sbp, eqn, opts, pmesh = EulerEquationMod.createObjects(ARGS[1])
    @assert mesh.isDG == true
    @assert opts["jac_method"] == 2
    @assert opts["run_type"] == 5

    context("Checking Functional Object Creation") do

      lift = EulerEquationMod.createFunctionalData(mesh, sbp, eqn, opts, opts["num_functionals"])
      @fact lift.is_objective_fn --> false
      @fact lift.geom_faces_functional --> [3]
      @fact lift.ndof --> 2
      @fact lift.bndry_force --> Complex{Float64}[0.0, 0.0]
      @fact lift.lift_val --> zero(Complex{Float64})
      @fact lift.drag_val --> zero(Complex{Float64})
      @fact lift.dLiftdAlpha --> zero(Complex{Float64})
      @fact lift.dDragdAlpha --> zero(Complex{Float64})

    end # End context("Checking Functional Object Creation")

    drag = EulerEquationMod.createObjectiveFunctionalData(mesh, sbp, eqn, opts)

    context("Checking Objective Functional Object Creation") do

      @fact drag.is_objective_fn --> true
      @fact drag.geom_faces_functional --> [3]
      @fact drag.ndof --> 2
      @fact drag.bndry_force --> Complex{Float64}[0.0, 0.0]
      @fact drag.lift_val --> zero(Complex{Float64})
      @fact drag.drag_val --> zero(Complex{Float64})
      @fact drag.dLiftdAlpha --> zero(Complex{Float64})
      @fact drag.dDragdAlpha --> zero(Complex{Float64})

    end # context("Checking Objective Functional Object Creation")

    EulerEquationMod.solve_euler(mesh, sbp, eqn, opts, pmesh)
    EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, drag)

    context("Checking Functional Computation") do

      analytical_drag = -1/1.4
      drag_error = norm(drag.drag_val - analytical_drag, 2)
      @fact drag_error --> roughly(0.0, atol = 1e-2)

    end # End context("Checking Functional Computation")

  end # End facts("--- Testing Functional Computation On a Boundary ---")

  facts("--- Tesing adjoint computation on the boundary for DG Meshes---") do
    # println("testing adjoint functions\n")
    resize!(ARGS, 1)
    ARGS[1] = "input_vals_airfoil.jl"
    include("../../src/solver/euler/startup.jl")
    @assert opts["aoa"] == 2.0*pi/180

    lift = EulerEquationMod.createObjectiveFunctionalData(mesh, sbp, eqn, opts)
    EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, lift)

    disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    orig_Ju = deepcopy(lift.lift_val) # Copy the original objective value
    orig_q_vec = deepcopy(eqn.q_vec)
    original_res_vec = deepcopy(eqn.res_vec)

    context("Checking partial dR/dq Calculation") do

      # Copy all the original values
      disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
      rand_vec = rand(length(eqn.q_vec))
      fill!(eqn.res, 0.0)
      fill!(eqn.res_vec, 0.0)
      res_norm = NonlinearSolvers.calcResidual(mesh, sbp, eqn, opts, evalResidual)
      boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
      res_jac, jacData = EulerEquationMod.calcResidualJacobian(mesh, sbp, eqn, opts)
      contract_vec = res_jac*rand_vec

      # Check against FD
      copy!(eqn.q_vec, orig_q_vec)
      for i = 1:length(eqn.q_vec)
        eqn.q_vec[i] += 1e-6*rand_vec[i]
      end
      disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
      fill!(eqn.res, 0.0)
      fill!(eqn.res_vec, 0.0)
      res_norm = NonlinearSolvers.calcResidual(mesh, sbp, eqn, opts, evalResidual)
      assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
      partialRpartialu = (eqn.res_vec - original_res_vec)/1e-6

      for i = 1:length(partialRpartialu)
        @fact abs(real(contract_vec[i] - partialRpartialu[i])) --> roughly(0.0,
                                                                    atol = 1e-5)
        # println(f,real(contract_vec[i] - partialRpartialu[i]))
      end

      for i = 1:length(eqn.q_vec)
        eqn.q_vec[i] = orig_q_vec[i]
      end

    end # End context("--- Checking partial dR/dq Calculation")

    context("Checking partial dJ/dq Calculation") do

      func_deriv_arr = zeros(eqn.q)
      func_deriv = zeros(eqn.q_vec)
      functional_edges = opts["geom_faces_objective"]
      functional_name = EulerEquationMod.FunctionalDict[opts["objective_function"]]

      boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
      EulerEquationMod.calcFunctionalDeriv(mesh, sbp, eqn, opts, lift,
                                        func_deriv_arr)  # populate df_dq_bndry
      assembleSolution(mesh, sbp, eqn, opts, func_deriv_arr, func_deriv)

      rand_vec = rand(length(eqn.q_vec))
      contract_val = dot(rand_vec,func_deriv)

      # Check with finite difference
      eqn.q_vec += 1e-6*rand_vec
      disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
      EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, lift)
      dJdu_fd = (lift.lift_val-orig_Ju)/1e-6

      @fact norm(dJdu_fd - contract_val, 2) --> roughly(0.0, atol = 1e-8)

      # Restore q_vec
      for i = 1:length(eqn.q_vec)
        eqn.q_vec[i] = orig_q_vec[i]
      end

    end # End context("--- Checking partial dJ/dq Calculation")


    context("checking derivative computation using adjoint vector") do

      @assert opts["aoa"] == 2.0*pi/180

      EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, lift)
      dJdAlpha = lift.dLiftdAlpha

      adjoint_vec = zeros(Complex128, mesh.numDof)
      EulerEquationMod.calcAdjoint(mesh, sbp, eqn, opts, lift, adjoint_vec)

      # Check dJdALpha against the complex step method
      @assert opts["epsilon"] == 1e-20
      eqn.params.aoa += opts["epsilon"]*im
      EulerEquationMod.calcBndryFunctional(mesh, sbp, eqn, opts, lift)
      dJdAlpha_comp = imag(lift.lift_val)/opts["epsilon"]

      # Get the partial derivative of the residual vector w.r.t aoa
      eqn.params.aoa = opts["aoa"]
      eqn.params.aoa += opts["epsilon"]*im # Imaginary perturbation
      fill!(eqn.res_vec, 0.0)
      fill!(eqn.res, 0.0)
      res_norm = NonlinearSolvers.calcResidual(mesh, sbp, eqn, opts, evalResidual)
      dRdAlpha = imag(eqn.res_vec)/opts["epsilon"]

      dLdx_adjoint = dJdAlpha + dot(adjoint_vec, dRdAlpha)
      eqn.params.aoa = opts["aoa"]

      #----- Finite Differencing -----#
      pert = 1e-6 # FD perturbation
      eqn.params.aoa = opts["aoa"] + pert
      fill!(eqn.q, 0.0)
      fill!(eqn.q_vec, 0.0)
      fill!(eqn.res, 0.0)
      fill!(eqn.res_vec, 0.0)

      # Rerun with the perturbed value
      ICfunc_name = opts["IC_name"]
      ICfunc = EulerEquationMod.ICDict[ICfunc_name]
      ICfunc(mesh, sbp, eqn, opts, eqn.q_vec)
      pmesh = mesh
      EulerEquationMod.init(mesh, sbp, eqn, opts, pmesh)
      call_nlsolver(mesh, sbp, eqn, opts, pmesh)

      # lift.lift_val = 0.0
      EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, lift)
      dLdx = (real(lift.lift_val) - orig_Ju)/pert

      errfd_norm = norm(dLdx - dLdx_adjoint,2)
      println("errfd_norm = $errfd_norm")
      @fact errfd_norm --> roughly(0.0, atol = 1e-6)

    end # End context("checking derivative computation using adjoint vector")

  end # End facts("--- Tesing adjoint computation on the boundary for DG Meshes---")

  #=
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
    include("../../src/solver/euler/startup.jl")

    # Assert basic facts
    @assert mesh.isDG == true
    @assert opts["jac_method"] == 2
    @assert opts["calc_functional"] == false

    # Check facts
    @fact opts["objective_function"] --> "drag"
    @fact opts["geom_faces_objective"] --> [3]

    drag = EulerEquationMod.OptimizationData{Complex128}(mesh, sbp, opts)
    EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, drag)

    @fact drag.is_objective_fn --> true
    analytical_val = -1/1.4
    println("drag.val = ", drag.val)
    drag_err = norm(drag.val - analytical_val)
    @fact drag_err --> roughly(0.0001, atol = 1e-4)

  end # End facts("--- Testing Objective Function Computation On a Boundary ---")



  =#
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


  end  # End facts("--- Testing Functional Computation On a Boundary ---") do
  =#
end # End function test_adjoint

add_func1!(EulerTests, test_adjoint, [TAG_ADJOINT])
