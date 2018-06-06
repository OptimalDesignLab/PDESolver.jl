function test_modes()
  @testset "--- Testing Sparse/Dense Jacobian ---" begin

    fname = "input_vals_vortex1.jl"
    println("\n\ntesting ", fname)
  #  fname = "input_vals_vortex5.jl"
    mesh, sbp, eqn, opts = run_solver(fname)

    @test  calcNorm(eqn, eqn.res_vec)  < 1e-9


    println("\n\n----- Testing jacobian vector product -----")
    # test jacobian vector product
    # load a new initial condition
    ICFunc = EulerEquationMod.ICDict["ICIsentropicVortex"]
    ICFunc(mesh, sbp, eqn, opts, eqn.q_vec)
    NonlinearSolvers.array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec)
#=  # unneeded after LinearOperator introduction
    # compute jacobian
    Tsol = eltype(eqn.q)
    Tres = eltype(eqn.res)
    jac = SparseMatrixCSC(mesh, typeof(real(eqn.res[1])))
    epsilon = 1e-20
    pert = complex(0, epsilon)
    newton_data = NonlinearSolvers.NewtonData{Tsol, Tres}(mesh, sbp, eqn, opts)
    NonlinearSolvers.calcJacobianSparse(mesh, sbp, eqn, opts, EulerEquationMod.evalResidual, Array{Float64}( 0,0,0), pert, jac)

    # to mat-vec product
    newton_data = NonlinearSolvers.NewtonData{Tsol, Tres}(mesh, sbp, eqn, opts)
    v = ones(mesh.numDof)  # vector to multiply jacobian against
    result1 = jac*v
    result2 = zeros(mesh.numDof)

    NonlinearSolvers.calcJacVecProd(mesh, sbp, eqn, opts, pert, physicsRhs, (EulerEquationMod.evalResidual,), v, result2)

    # check the two products are equal
    for i=1:mesh.numDof
      @test isapprox( result1[i], result2[i]) atol=1e-13
    end
=#
  #  results_all = hcat(result1, result2)
  #  println("results_all = \n", results_all)

  #  results_diff = result1 - result2
  #  diff_norm = calcNorm(eqn, results_diff)
  #  println("diff_norm = ", diff_norm)

   #= 
    # test entropy variables
    fname = "input_vals_vortexa.jl"
    println("\n\ntesting ", fname)
    mesh, sbp, eqn, opts = solvePDE(fname)

    @test  calcNorm(eqn, eqn.res_vec)  < 1e-9
    =#

    fname = "input_vals_vortex2.jl"
    println("\n\ntesting ", fname)
    mesh, sbp, eqn, opts = solvePDE(fname)
    @test  calcNorm(eqn, eqn.res_vec)  < 1e-9

  #=
    # test entropy variables
    fname = "input_vals_vortex2a.jl"
    mesh, sbp, eqn, opts = solvePDE(fname)
    @test  calcNorm(eqn, eqn.res_vec)  < 1e-9
  =#

    fname = "input_vals_vortex3dg.jl"
    println("\n\ntesting ", fname)
    mesh, sbp, eqn, opts = solvePDE(fname)
    @test  calcNorm(eqn, eqn.res_vec)  < 1e-9

  #=
    # test entropy variables
    fname = "input_vals_vortex3a.jl"
    mesh, sbp, eqn, opts = solvePDE(fname)
    @test  calcNorm(eqn, eqn.res_vec)  < 1e-9
  =#

    fname = "input_vals_vortex4.jl"
    println("\n\ntesting ", fname)
    mesh, sbp, eqn, opts = solvePDE(fname)
    @test  calcNorm(eqn, eqn.res_vec)  < 1e-9

  #=
    # test entropy variables
    fname = "input_vals_vortex4a.jl"
    mesh, sbp, eqn, opts = solvePDE(fname)
    @test  calcNorm(eqn, eqn.res_vec)  < 1e-9
  =#
  end

  return nothing
end

add_func1!(EulerTests, test_modes, [TAG_SHORTTEST])
