#Test functional Integrate and adjoint for euler equation.

@doc """
Euler Equation -- test_adjoint

The function tests for the correctness of objective and non-objective functional
computation and then tests if the adjoint vector is being computed correctly.
Adjoint vector computation is based on an objective function. This is a serial
test and uses the input file called

`input_vals_vortex_adjoint_DG.jl`

"""->

function test_adjoint()

  facts("--- Testing Functional Computation On a Boundary ---") do

    ARGS[1] = "input_vals_vortex_adjoint_DG.jl"
    mesh, sbp, eqn, opts, pmesh = EulerEquationMod.createObjects(ARGS[1])
    @assert mesh.isDG == true
    @assert opts["jac_method"] == 2
    @assert opts["run_type"] == 5

    context("Checking Functional Object Creation") do

      lift = EulerEquationMod.createFunctionalData(mesh, sbp, eqn, opts,
                                                   opts["num_functionals"])
      @fact lift.bcnums --> [4]
      @fact lift.ndof --> 2
      @fact lift.bndry_force --> Complex{Float64}[0.0, 0.0]
      @fact lift.lift_val --> zero(Complex{Float64})
      @fact lift.drag_val --> zero(Complex{Float64})
      @fact lift.dLiftdaoa --> zero(Complex{Float64})
      @fact lift.dDragdaoa --> zero(Complex{Float64})

    end # End context("Checking Functional Object Creation")

    drag = EulerEquationMod.createObjectiveFunctionalData(mesh, sbp, eqn, opts)

    context("Checking Objective Functional Object Creation") do

      @fact drag.bcnums --> [4]
      @fact drag.ndof --> 2
      @fact drag.bndry_force --> Complex{Float64}[0.0, 0.0]
      @fact drag.lift_val --> zero(Complex{Float64})
      @fact drag.drag_val --> zero(Complex{Float64})
      @fact drag.dLiftdaoa --> zero(Complex{Float64})
      @fact drag.dDragdaoa --> zero(Complex{Float64})


    end # context("Checking Objective Functional Object Creation")

    context("Checking Functional Computation Before Solve") do
      massflow = EulerEquationMod.MassFlowDataConstructor(Complex128, mesh, sbp, eqn, opts, [3])
      EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, massflow)
      @fact massflow.val --> roughly(0.0, atol=1e-13)

      # check derivative
      func_deriv = zeros(eqn.res)
      func_deriv2 = zeros(eqn.res) 
      h = 1e-20
      pert = Complex128(0, h)
      for i=1:length(eqn.q)
        eqn.q[i] += pert
        boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
        EulerEquationMod.calcBndryFunctional(mesh, sbp, eqn, opts, massflow)
        func_deriv[i] = imag(massflow.val)/h
        eqn.q[i] -= pert
      end

      EulerEquationMod.calcFunctionalDeriv(mesh, sbp, eqn, opts, massflow, func_deriv2)

      @fact vecnorm(func_deriv - func_deriv2) --> roughly(0.0, atol=1e-13)
    end  # context Checking Functional Computation Before Solve


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
    include("../../src/solver/euler/startup.jl")  #TODO: use run_euler
    @assert opts["aoa"] == 2.0

    lift = EulerEquationMod.createObjectiveFunctionalData(mesh, sbp, eqn, opts)
    EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, lift)

    context("Checking functional derivative w.r.t angle of attack") do
      dLiftdaoa = lift.dLiftdaoa
      dDragdaoa = lift.dDragdaoa
      pert = complex(0, 1e-20)
      eqn.params.aoa += pert
      EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, lift)
      dLiftdaoa_c = imag(lift.lift_val)/imag(pert)
      dDragdaoa_c = imag(lift.drag_val)/imag(pert)
      error_lift = norm(dLiftdaoa - dLiftdaoa_c)
      error_drag = norm(dDragdaoa - dDragdaoa_c)
      eqn.params.aoa -= pert
    end # End context("Checking functional derivative w.r.t angle of attack")

    EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, lift)
    # Write the values of BoundaryForceData to file
    fname = "BoundaryForceData_values.dat"
    file = open(fname, "w")
    println(file, lift.lift_val)
    println(file, lift.drag_val)
    println(file, lift.dLiftdaoa)
    println(file, lift.dDragdaoa)
    close(file)

    disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    orig_Ju = deepcopy(lift.lift_val) # Copy the original objective value
    orig_q_vec = deepcopy(eqn.q_vec)
    original_res_vec = deepcopy(eqn.res_vec)

    context("Checking partial dJ/dq Calculation") do

      func_deriv_arr = zeros(eqn.q)
      func_deriv = zeros(eqn.q_vec)

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

    context("Checking complete derivative of a functional using adjoint vector") do

      EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, lift)
      lift_val = lift.lift_val

      pc, lo = getNewtonPCandLO(mesh, sbp, eqn, opts)
      ls = StandardLinearSolver(pc, lo, eqn.comm, opts)
      adjoint_vec = zeros(Complex128, mesh.numDof)
      calcAdjoint(mesh, sbp, eqn, opts, ls, lift, adjoint_vec, recalc_jac=true, recalc_pc=true)

      # Get the complete derivative of the function
      dJdaoa = EulerEquationMod.eval_dJdaoa(mesh, sbp, eqn, opts, lift, "lift", adjoint_vec)

      # Check complete derivatives w.r.t alpha using finite difference
      pert = 1e-6
      eqn.params.aoa += pert
      EulerEquationMod.solve_euler(mesh, sbp, eqn, opts, mesh)
      EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, lift)
      lift_pert = lift.lift_val
      eqn.params.aoa -= pert

      dJdaoa_fd = (lift_pert - lift_val)/pert
      err_val = norm(dJdaoa - dJdaoa_fd, 2)

      @fact err_val --> roughly(0.0, atol = 1e-6)

    end # End context("Checking complete derivative of a functional using adjoint vector")

  end # End facts("--- Tesing adjoint computation on the boundary for DG Meshes---")

end # End function test_adjoint

add_func1!(EulerTests, test_adjoint, [TAG_ADJOINT, TAG_LONGTEST])
