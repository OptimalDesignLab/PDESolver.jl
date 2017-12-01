# test boundary forces
@doc """
Advection Equation -- test_adjoint

The function tests for the correctness of objective and non-objective functional
computation and then tests if the adjoint vector is being computed correctly.
Adjoint vector computation is based on an objective function. This is a serial
test and uses the input file called

`input_vals_functional_DG.jl`

"""->
function test_adjoint()
  #= TODO: Uncomment when CG meshes support mesh.bndry_geo_nums

  facts("--- Testing Boundary Functional Computation on CG Mesh ---") do
    clean_dict(arg_dict)
    ARGS[1] = "input_vals_functional_CG.jl"
    include("../../src/solver/advection/startup_advection.jl")  # initialization and construction
    println("use_DG = ", arg_dict["use_DG"])


    @fact mesh.isDG --> false
    @fact opts["functional_name1"] --> "qflux"
    @fact opts["analytical_functional_val"] --> roughly(2*(exp(1) - 1), atol=1e-12)
    @fact opts["geom_edges_functional1"] --> [2,3]

    fname = "./functional_error1.dat"
    error = readdlm(fname)

    @fact error[1] --> roughly(0.0060826244541961885, atol=1e-6)

    rm("./functional_error1.dat") # Delete the file
  end
  =#

  facts("--- Testing Boundary Functional & Adjoint Computation On DG Mesh ---") do

    ARGS[1] = "input_vals_functional_DG.jl"
    mesh, sbp, eqn, opts, pmesh = AdvectionEquationMod.createObjects(ARGS[1])

    @assert mesh.isDG == true
    @assert opts["jac_method"] == 2
    @assert opts["run_type"] == 5

    context("Checking Functional Object Creation") do

      functional = AdvectionEquationMod.createFunctionalData(mesh, sbp, eqn,
                   opts, opts["num_functionals"])

      @fact functional.bcnums --> [2,3]
      @fact functional.val --> zero(Complex{Float64})
      @fact functional.target_qflux --> zero(Complex{Float64})

    end # End context("Checking Functional Object Creation")

    objective = AdvectionEquationMod.createObjectiveFunctionalData(mesh, sbp, eqn, opts)

    context("Checking Objective Functional Object Creation") do

      @fact objective.bcnums --> [2,3]
      @fact objective.val --> zero(Complex{Float64})
      @fact objective.target_qflux --> zero(Complex{Float64})

    end # End context("Checking Objective Functional Object Creation")

    AdvectionEquationMod.solve_advection(mesh, sbp, eqn, opts, pmesh)
    AdvectionEquationMod.evalFunctional(mesh, sbp, eqn, opts, objective)

    context("Checking Functional Computation") do

      analytical_val = 3.0
      functional_error = norm(real(objective.val) - analytical_val,2)
      @fact functional_error --> roughly(0.0, atol=1e-12)

      # test another functional
      func = AdvectionEquationMod.IntegralQDataConstructor(Complex128, mesh, sbp, eqn, opts, opts["functional_bcs1"], Int[])
      AdvectionEquationMod.calcBndryFunctional(mesh, sbp, eqn, opts, func)
      @fact func.val --> roughly(analytical_val, atol=1e-13)


    end # End context("Checking Functional Computation")

    context("Checking Adjoint Computation on DG mesh") do

      adjoint_vec = zeros(Complex{Float64}, mesh.numDof)
      AdvectionEquationMod.calcAdjoint(mesh, sbp, eqn, opts, objective, adjoint_vec)

      for i = 1:length(adjoint_vec)
        @fact real(adjoint_vec[i]) --> roughly(1.0 , atol=1e-10)
      end

    end # End context("Checking Adjoint Computation on DG mesh")

  end # End facts("--- Testing Boundary Functional Computation On DG Mesh ---")

end # End function test_adjoint

add_func1!(AdvectionTests, test_adjoint, [TAG_ADJOINT, TAG_SHORTTEST])
