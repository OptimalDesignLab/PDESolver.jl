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
  EulerEquationMod.ICExp(mesh, sbp, eqn, opts, eqn.q_vec)
  disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)

  newton_data_dense, jac, rhs_vec = NonlinearSolvers.setupNewton(mesh, mesh, sbp, eqn, opts)

  opts2 = copy(opts)
  opts2["jac_type"] = 2  # SparseMatrixCSC
  newton_data_sparse, jac_sparse, rhs_vec_sparse = NonlinearSolvers.setupNewton(mesh, mesh, sbp, eqn, opts2)

  @fact typeof(jac) <: Array --> true
  @fact typeof(jac_sparse) <: SparseMatrixCSC -->  true

  NonlinearSolvers.physicsJac(newton_data_dense, mesh, sbp, eqn, opts, jac, ctx_residual)
  NonlinearSolvers.physicsJac(newton_data_sparse, mesh, sbp, eqn, opts2, jac_sparse, ctx_residual)

  jac_dense2 = full(jac_sparse)


  for i=1:mesh.numDof
    for j=1:mesh.numDof
      @fact jac_dense2[j, i] --> roughly(jac[j, i], atol=1e-12)
    end
  end


  return nothing
end

add_func3!(EulerTests, test_homotopy, test_homotopy_inputfile, test_homotopy_moddict, [TAG_HOMOTOPY, TAG_SHORTTEST])

function test_homotopy_convergence()

  rmfile("convergence.dat")
  mesh, sbp, eqn, opts = run_solver("input_vals_homotopy.jl")

  @fact calcNorm(eqn, eqn.res_vec, strongres=true) --> less_than(opts["res_abstol"])
  data = readdlm("convergence.dat")

  @fact size(data, 1) --> less_than(30)  # 16 iterations to convergence

  return nothing
end

add_func1!(EulerTests, test_homotopy_convergence, [TAG_HOMOTOPY, TAG_LONGTEST])
