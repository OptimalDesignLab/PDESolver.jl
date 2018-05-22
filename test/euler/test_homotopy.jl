# tests for homotopy.jl

const test_homotopy_inputfile = "input_vals_channel.jl"
const test_homotopy_moddict = Dict{ASCIIString, Any}(
  "Flux_name" => "RoeFlux", 
  "use_DG" => true, 
  "IC_name" => "ICFreeStream",
  "BC1_name" => "FreeStreamBC",
  "jac_type" => 1,
  "jac_method" => 2,
  "new_fname" => "input_vals_channel_dg_homotopy")


function test_homotopy(mesh, sbp, eqn, opts)

  facts("----- Testing Homotopy operators ------") do
  # the initial condition is uniform flow, so the residual of the homotopy
  # should be zero

  res = zeros(eqn.res)
  fill!(eqn.res, 42)  # make sure eqn.res is not modified

  EulerEquationMod.calcHomotopyDiss(mesh, sbp, eqn, opts, res)

  for i=1:mesh.numEl
    @fact norm(res[:, :, i]) --> roughly(0.0, atol=1e-12)
  end

  for i=1:length(eqn.res)
    @fact eqn.res[i] --> 42
  end

  # make homotopy look like a physics
  function homotopy_physics_test(mesh, sbp, eqn, opts, t=0.0)
    evalHomotopy(mesh, sbp, eqn, opts, eqn.res)
  end

  ctx_residual = (homotopy_physics_test,)

  # test jacobian
  t = 0.0
  EulerEquationMod.ICExp(mesh, sbp, eqn, opts, eqn.q_vec)
  array1DTo3D(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)

  pc, lo = NonlinearSolvers.getHomotopyPCandLO(mesh, sbp, eqn, opts)
  ls_dense = StandardLinearSolver(pc, lo, eqn.comm, opts)
#  newton_data_dense, rhs_vec = NonlinearSolvers.setupNewton(mesh, mesh, sbp, eqn, opts, ls_dense)

  opts2 = copy(opts)
  opts2["jac_type"] = 2  # SparseMatrixCSC
  pc, lo = NonlinearSolvers.getHomotopyPCandLO(mesh, sbp, eqn, opts2)
  ls_sparse = StandardLinearSolver(pc, lo, eqn.comm, opts)

#  newton_data_sparse, rhs_vec_sparse = NonlinearSolvers.setupNewton(mesh, mesh, sbp, eqn, opts2, ls_sparse)

  lo2_dense = getBaseLO(ls_dense.lo)
  lo2_sparse = getBaseLO(ls_sparse.lo)
  @fact typeof(lo2_dense.A) <: Array --> true
  @fact typeof(lo2_sparse.A) <: SparseMatrixCSC -->  true

  calcLinearOperator(ls_dense.lo, mesh, sbp, eqn, opts, ctx_residual, t)
  calcLinearOperator(ls_sparse.lo, mesh, sbp, eqn, opts2, ctx_residual, t)
#  NonlinearSolvers.physicsJac(newton_data_dense, mesh, sbp, eqn, opts, jac, ctx_residual)
#  NonlinearSolvers.physicsJac(newton_data_sparse, mesh, sbp, eqn, opts2, jac_sparse, ctx_residual)

  jac_dense = lo2_dense.A
  jac_dense2 = full(lo2_sparse.A)


  for i=1:mesh.numDof
    for j=1:mesh.numDof
      @fact jac_dense2[j, i] --> roughly(jac_dense[j, i], atol=1e-12)
    end
  end

  #test matrix-free products
  opts2["jac_type"] = 4
  println("constructing mat-free linear operator")
  pc_free, lo_free = NonlinearSolvers.getHomotopyPCandLO(mesh, sbp, eqn, opts2)
  println("typeof(lo_free) = ", typeof(lo_free))
  println("calculating linear operaotr")
  calcLinearOperator(lo_free, mesh, sbp, eqn, opts2, ctx_residual, t)
  x = rand(mesh.numDof)
  b = zeros(x)
  b2 = zeros(x)

  println("applying linear operator")
  applyLinearOperator(ls_dense.lo, mesh, sbp, eqn, opts, ctx_residual, t, x, b)
  applyLinearOperator(lo_free, mesh, sbp, eqn, opts2, ctx_residual, t, x, b2)

  @fact norm(b - b2) --> roughly(0.0, atol=1e-13)


#=
  # explicit jacobian computation
  println("testing explicit jacobian calculation")
  opts2["jac_type"] = 2
  pc, lo = NonlinearSolvers.getHomotopyPCandLO(mesh, sbp, eqn, opts)
  pc, lo2 = NonlinearSolvers.getHomotopyPCandLO(mesh, sbp, eqn, opts)

  calcLinearOperator(lo, mesh, sbp, eqn, opts, ctx_residual, t)
  opts2["calc_jac_explicit"] = true
  calcLinearOperator(lo2, mesh, sbp, eqn, opts2, ctx_residual, t)

  A = getBaseLO(lo).A
  A2 = getBaseLO(lo2).A
  @fact norm(full(A) - full(A2)) --> roughly(0.0, atol=1e-13)
=#

 
end  # end do
  return nothing
end

add_func3!(EulerTests, test_homotopy, test_homotopy_inputfile, test_homotopy_moddict, [TAG_HOMOTOPY, TAG_SHORTTEST])

function test_homotopy_convergence()

  println("\ntesting homotopy jacobians")
  # run for 2 iterations, with explicit, coloring jacobians, make sure
  # residuals are the same
  opts_tmp = read_input_file("input_vals_homotopy.jl")
  opts_tmp["calc_jac_explicit"] = true
  opts_tmp["itermax"] = 2
  make_input(opts_tmp, "input_vals_homotopy_tmp.jl")
  opts_tmp["calc_jac_explicit"] = false
  make_input(opts_tmp, "input_vals_homotopy_tmp2.jl")

  mesh, sbp, eqn, opts = run_solver("input_vals_homotopy_tmp.jl")
  mesh2, sbp2, eqn2, opts2 = run_solver("input_vals_homotopy_tmp2.jl")

  @fact calcNorm(eqn, eqn.q) --> roughly(calcNorm(eqn2, eqn2.q), atol=1e-13)


  # now run a full case
  println("\nTesting homotopy")
  rmfile("convergence.dat")
  mesh, sbp, eqn, opts = run_solver("input_vals_homotopy.jl")

  @fact calcNorm(eqn, eqn.res_vec, strongres=true) --> less_than(opts["res_abstol"])
  data = readdlm("convergence.dat")

  @fact size(data, 1) --> less_than(30)  # 16 iterations to convergence

  return nothing
end

add_func1!(EulerTests, test_homotopy_convergence, [TAG_HOMOTOPY, TAG_LONGTEST])
