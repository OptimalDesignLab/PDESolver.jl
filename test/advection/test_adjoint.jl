# test boundary forces
@doc """
Advection Equation -- test_adjoint

The function tests for the correctness functional
computation and then tests if the adjoint vector is being computed correctly.
This is a serial test and uses the input file called

`input_vals_functional_DG.jl`

"""->
function test_adjoint()
  #= TODO: Uncomment when CG meshes support mesh.bndry_geo_nums

  @testset "--- Testing Boundary Functional Computation on CG Mesh ---" begin
    clean_dict(arg_dict)
    fname = "input_vals_functional_CG.jl"
    mesh, sbp, eqn, opts = solvePDE(fname)
    println("use_DG = ", arg_dict["use_DG"])


    @test ( mesh.isDG )== false
    @test ( opts["functional_name1"] )== "qflux"
    @test isapprox( opts["analytical_functional_val"], 2*(exp(1) - 1)) atol=1e-12
    @test ( opts["geom_edges_functional1"] )== [2,3]

    fname = "./functional_error1.dat"
    error = readdlm(fname)

    @test isapprox( error[1], 0.0060826244541961885) atol=1e-6

    rm("./functional_error1.dat") # Delete the file
  end
  =#

  @testset "--- Testing Boundary Functional & Adjoint Computation On DG Mesh ---" begin

    fname = "input_vals_functional_DG.jl"
    mesh, sbp, eqn, opts, pmesh = createObjects(fname)

    @assert mesh.isDG == true
    @assert opts["jac_method"] == 2
    @assert opts["run_type"] == 5

    functional = createFunctional(mesh, sbp, eqn,
                                  opts, opts["num_functionals"])


    @testset "Checking Functional Object Creation" begin

      @test ( functional.bcnums )== [2,3]
      @test ( functional.val )== zero(Complex{Float64})
      @test ( functional.target_qflux )== zero(Complex{Float64})

    end # End context("Checking Functional Object Creation")

    solvePDE(mesh, sbp, eqn, opts, pmesh)
    evalFunctional(mesh, sbp, eqn, opts, functional)

    @testset "Checking Functional Computation" begin

      analytical_val = 3.0
      functional_error = norm(real(functional.val) - analytical_val,2)
      @test isapprox( functional_error, 0.0) atol=1e-12

      # test another functional
      func = AdvectionEquationMod.IntegralQDataConstructor(Complex128, mesh, sbp, eqn, opts, opts["functional_bcs1"])
      AdvectionEquationMod.calcBndryFunctional(mesh, sbp, eqn, opts, func)
      @test isapprox( func.val, analytical_val) atol=1e-13


    end # End context("Checking Functional Computation")

    @testset "Checking Adjoint Computation on DG mesh" begin

      adjoint_vec = zeros(Complex{Float64}, mesh.numDof)
      pc, lo = getNewtonPCandLO(mesh, sbp, eqn, opts)
      ls = StandardLinearSolver(pc, lo, eqn.comm, opts)
      calcAdjoint(mesh, sbp, eqn, opts, ls, functional, adjoint_vec, recalc_jac=true, recalc_pc=true)

      for i = 1:length(adjoint_vec)
        @test isapprox( real(adjoint_vec[i]), 1.0) atol=1e-10
      end

    end # End context("Checking Adjoint Computation on DG mesh")

  end # End testset("--- Testing Boundary Functional Computation On DG Mesh ---")

end # End function test_adjoint

add_func1!(AdvectionTests, test_adjoint, [TAG_ADJOINT, TAG_SHORTTEST])
