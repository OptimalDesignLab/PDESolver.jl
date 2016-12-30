
facts("--- Testing Sparse/Dense Jacobian ---") do

  resize!(ARGS, 1)
  ARGS[1] = "input_vals_vortex1.jl"
  println("\n\ntesting ", ARGS[1])
#  ARGS[1] = "input_vals_vortex5.jl"
  include("../src/solver/euler/startup.jl")

  @fact calcNorm(eqn, eqn.res_vec) --> less_than(1e-9)


  println("\n\n----- Testing jacobian vector product -----")
  # test jacobian vector product
  # load a new initial condition
  ICFunc = EulerEquationMod.ICDict["ICIsentropicVortex"]
  ICfunc(mesh, sbp, eqn, opts, eqn.q_vec)
  NonlinearSolvers.disassembleSolution(mesh, sbp, eqn, opts, eqn.q_vec)
  jac = SparseMatrixCSC(mesh.sparsity_bnds, typeof(real(eqn.res[1])))
  epsilon = 1e-20
  pert = complex(0, epsilon)
  newton_data = NonlinearSolvers.NewtonData(mesh, sbp, eqn, opts)
  NonlinearSolvers.calcJacobianSparse(newton_data, mesh, sbp, eqn, opts, EulerEquationMod.evalResidual, Array(Float64, 0,0,0), pert, jac)

  newton_data = NonlinearSolvers.NewtonData(mesh, sbp, eqn, opts)
  v = ones(mesh.numDof)  # vector to multiply jacobian against
  result1 = jac*v
  result2 = zeros(mesh.numDof)
  NonlinearSolvers.calcJacVecProd(newton_data, mesh, sbp, eqn, opts, pert, EulerEquationMod.evalResidual, v, result2)

  # check the two products are equal
  for i=1:mesh.numDof
    @fact result1[i] --> roughly(result2[i], atol=1e-13)
  end

#  results_all = hcat(result1, result2)
#  println("results_all = \n", results_all)

#  results_diff = result1 - result2
#  diff_norm = calcNorm(eqn, results_diff)
#  println("diff_norm = ", diff_norm)

  
  # test entropy variables
  ARGS[1] = "input_vals_vortexa.jl"
  println("\n\ntesting ", ARGS[1])
  mesh, sbp, eqn, opts = run_euler(ARGS[1])

  @fact calcNorm(eqn, eqn.res_vec) --> less_than(1e-9)
  

  resize!(ARGS, 1)
  ARGS[1] = "input_vals_vortex2.jl"
  println("\n\ntesting ", ARGS[1])
  mesh, sbp, eqn, opts = run_euler(ARGS[1])
  @fact calcNorm(eqn, eqn.res_vec) --> less_than(1e-9)

#=
  # test entropy variables
  ARGS[1] = "input_vals_vortex2a.jl"
  mesh, sbp, eqn, opts = run_euler(ARGS[1])
  @fact calcNorm(eqn, eqn.res_vec) --> less_than(1e-9)
=#

  resize!(ARGS, 1)
  ARGS[1] = "input_vals_vortex3.jl"
  println("\n\ntesting ", ARGS[1])
  mesh, sbp, eqn, opts = run_euler(ARGS[1])
  @fact calcNorm(eqn, eqn.res_vec) --> less_than(1e-9)

#=
  # test entropy variables
  ARGS[1] = "input_vals_vortex3a.jl"
  mesh, sbp, eqn, opts = run_euler(ARGS[1])
  @fact calcNorm(eqn, eqn.res_vec) --> less_than(1e-9)
=#

  resize!(ARGS, 1)
  ARGS[1] = "input_vals_vortex4.jl"
  println("\n\ntesting ", ARGS[1])
  mesh, sbp, eqn, opts = run_euler(ARGS[1])
  @fact calcNorm(eqn, eqn.res_vec) --> less_than(1e-9)

#=
  # test entropy variables
  ARGS[1] = "input_vals_vortex4a.jl"
  mesh, sbp, eqn, opts = run_euler(ARGS[1])
  @fact calcNorm(eqn, eqn.res_vec) --> less_than(1e-9)
=#
end

facts("-----  Testing rk4 -----") do

  include("test_rk4.jl")
end
