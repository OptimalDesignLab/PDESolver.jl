# test jacobian calculation functions

"""
  Test the jacobian of individual terms
"""
function test_jac_terms()

  mesh, sbp, eqn, opts = run_solver("input_vals_jac2d.jl")
  mesh3, sbp3, eqn3, opts3 = run_solver("input_vals_jac3d.jl")

  test_pressure(eqn.params)
  test_pressure(eqn3.params)

  test_eulerflux(eqn.params)
  test_eulerflux(eqn3.params)

  nrm = [0.45, 0.55]
  nrm2 = -nrm

  println("testing all positive eigenvalues")
  q = Complex128[2.0, 3.0, 4.0, 7.0]
  qg = q + 1
  test_ad_inner(eqn.params, q, qg, nrm)

  println("testing all negative eigenvalues")
  q = Complex128[2.0, 3.0, 4.0, 7.0]
  qg = q + 1
  test_ad_inner(eqn.params, q, qg, nrm2)


  println("testing lambda1 entropy fix")
  q = Complex128[1.1, 0.47, 0.53, 2.2]
  qg = q + 1
  test_ad_inner(eqn.params, q, qg, nrm)
 
  println("testing lambda2 entropy fix")
  q = Complex128[1.1, -1.32, -1.34, 2.2] 
  qg = q + 1
  test_ad_inner(eqn.params, q, qg, nrm)
  test_ad_inner(eqn.params, q, qg, nrm2)

  println("testing lambda3 entropy fix")
  q = Complex128[1.1, -0.42, -0.45, 2.2]
  qg = q + 1
  test_ad_inner(eqn.params, q, qg, nrm)
  test_ad_inner(eqn.params, q, qg, nrm2)


  nrm = [0.45, 0.55, 0.65]
  nrm2 = -nrm

  println("testing all positive eigenvalues")
  q = Complex128[2.0, 3.0, 4.0, 5.0, 13.0]
  qg = q + 1
  test_ad_inner(eqn3.params, q, qg, nrm)

  println("testing all negative eigenvalues")
  q = Complex128[2.0, 3.0, 4.0, 5.0, 13.0]
  qg = q + 1
  test_ad_inner(eqn3.params, q, qg, nrm2)


  println("testing lambda1 entropy fix")
  # lambda1 entropy fix active
  q = Complex128[1.05, -1.1, -1.2, -1.3, 2.5] 
  qg = q + 1
  test_ad_inner(eqn3.params, q, qg, nrm)
  test_ad_inner(eqn3.params, q, qg, nrm2)

  println("testing lambda2 entropy fix")
  # lambda1 entropy fix active
  q = Complex128[1.05, 0.9, 1.2, -1.3, 2.5]
  qg = q + 1
  test_ad_inner(eqn3.params, q, qg, nrm)
  test_ad_inner(eqn3.params, q, qg, nrm2)

  println("testing lambda3 entropy fix")
  # lambda3 entropy fix active
  q = Complex128[1.05, -0.52, -0.47, -0.36, 8.5]
  qg = q + 1
  test_ad_inner(eqn3.params, q, qg, nrm)
  test_ad_inner(eqn3.params, q, qg, nrm2)

  test_jac_assembly(mesh, sbp, eqn, opts)

  println("testing 3d")
  test_jac_assembly(mesh3, sbp3, eqn3, opts3)

  return nothing
end


add_func1!(EulerTests, test_jac_terms, [TAG_SHORTTEST, TAG_JAC])

function test_pressure{Tdim}(params::AbstractParamType{Tdim})


  if Tdim == 2
    q = Complex128[2.0, 3.0, 4.0, 7.0]
  else
    q = Complex128[2.0, 3.0, 4.0, 5.0, 13.0]
  end

  numDofPerNode = length(q)

  p_dot = zeros(q)

  h = 1e-20
  pert = Complex128(0, h)
  for j=1:numDofPerNode
    q[j] += pert
    p = EulerEquationMod.calcPressure(params, q)
    p_dot[j] = imag(p)/h
    q[j] -= pert
  end

  p_dot2 = zeros(q)
  EulerEquationMod.calcPressure_diff(params, q, p_dot2)

  @fact maximum(abs(p_dot - p_dot2))  --> roughly(0.0, atol=1e-14)

  return nothing
end


function test_eulerflux{Tdim}(params::AbstractParamType{Tdim})

  if Tdim == 2
    nrm = [0.45, 0.55]
    q = Complex128[2.0, 3.0, 4.0, 7.0]
  else
    nrm = [0.45, 0.55, 0.65]
    q = Complex128[2.0, 3.0, 4.0, 5.0, 13.0]
  end

  numDofPerNode = length(q)

  aux_vars = Complex128[0.0]
  F = zeros(Complex128, numDofPerNode)

  res = zeros(Complex128, numDofPerNode, numDofPerNode)

  h =1e-20
  pert = Complex128(0, h)
  for j=1:numDofPerNode
    q[j] += pert
    EulerEquationMod.calcEulerFlux(params, q, aux_vars, nrm, F)

    for k=1:numDofPerNode
      res[k, j] = imag(F[k])/h
    end

    q[j] -= pert
  end


  res2 = zeros(res)
  EulerEquationMod.calcEulerFlux_diff(params, q, aux_vars, nrm, res2)

  @fact maximum(abs(res - res2)) --> roughly(0.0, atol=1e-14)
end

function test_ad_inner{Tdim}(params::AbstractParamType{Tdim}, qL, qR, nrm)

  # compute jacobian with complex step and AD, compare results

  numDofPerNode = length(qL)

  aux_vars = Complex128[0.0]
  F = zeros(Complex128, numDofPerNode)

  resL = zeros(Complex128, numDofPerNode, numDofPerNode)
  resR = zeros(Complex128, numDofPerNode, numDofPerNode)
  h = 1e-20
  pert = Complex128(0, h)
  for j=1:numDofPerNode
#    println("\nj = ", j)
    qL[j] += pert
    EulerEquationMod.RoeSolver(params, qL, qR, aux_vars, nrm, F)

    for k=1:numDofPerNode
      resL[k, j] = imag(F[k])/h
    end

    qL[j] -= pert
  end

  for j=1:numDofPerNode
#    println("\nj = ", j + 5)
    qR[j] += pert
    EulerEquationMod.RoeSolver(params, qL, qR, aux_vars, nrm, F)

    for k=1:numDofPerNode
      resR[k, j] = imag(F[k])/h
    end

    qR[j] -= pert
  end

  # AD version
#  println("\nTesting AD version")
  resL2 = zeros(resL)
  resR2 = zeros(resR)

  EulerEquationMod.RoeSolver_diff(params, qL, qR, aux_vars, nrm, resL2, resR2)
#=
  println("\n----- Comparing left results -----")
  println("resL = \n", resL)
  println("resL2 = \n", resL2)
  println("diff = \n", resL - resL2)
  println("diffnorm = \n", vecnorm(resL - resL2))

  println("\n----- Comparing right results -----")
  println("resR = \n", resR)
  println("resR2 = \n", resR2)
  println("diff = \n", resR - resR2)
  println("diffnorm = \n", vecnorm(resR - resR2))
=#

  @fact maximum(abs(resL - resL2)) --> roughly(0.0, atol=1e-14)
  @fact maximum(abs(resR - resR2)) --> roughly(0.0, atol=1e-14)

  return nothing
end

function test_jac_assembly(mesh, sbp, eqn, opts)

  # use a spatially varying solution
  icfunc = EulerEquationMod.ICDict["ICExp"]
  icfunc(mesh, sbp, eqn, opts, eqn.q_vec)
  disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)

  # get the correct differentiated flux function (this is needed because the
  # input file set calc_jac_explicit = false
  eqn.flux_func_diff = EulerEquationMod.FluxDict_diff[opts["Flux_name"]]

  # test volume terms
  opts["addBoundaryIntegrals"] = false
  opts["addFaceIntegrals"] = false
  
  jac1 = SparseMatrixCSC(mesh, Float64)
  jac2 = zeros(mesh.numDof, mesh.numDof)
  assembler = NonlinearSolvers._AssembleElementData(jac2, mesh, sbp, eqn, opts)

  # compute jacobian via coloring
  opts["calc_jac_explicit"] = false
  ctx_residual = (evalResidual,)
  NonlinearSolvers.physicsJac(mesh, sbp, eqn, opts, jac1, ctx_residual)

  # compute jacobian explicitly
  opts["calc_jac_explicit"] = true
  evalJacobian(mesh, sbp, eqn, opts, assembler)

  jac1d = full(jac1)

  @fact maximum(abs(jac1d - jac2)) --> roughly(0.0, atol=1e-14)


  # test face integrals
  opts["addVolumeIntegrals"] = false
  opts["addFaceIntegrals"] = true
  fill!(jac1, 0.0)
  fill!(jac2, 0.0)
  # compute jacobian via coloring
  opts["calc_jac_explicit"] = false
  ctx_residual = (evalResidual,)
  NonlinearSolvers.physicsJac(mesh, sbp, eqn, opts, jac1, ctx_residual)

  # compute jacobian explicitly
  opts["calc_jac_explicit"] = true
  evalJacobian(mesh, sbp, eqn, opts, assembler)

  jac1d = full(jac1)

  @fact maximum(abs(jac1d - jac2)) --> roughly(0.0, atol=1e-14)

  # test boundary integral
  # test face integrals
  opts["addVolumeIntegrals"] = false
  opts["addFaceIntegrals"] = false
  opts["addBoundaryIntegrals"] = true
  fill!(jac1, 0.0)
  fill!(jac2, 0.0)
  # compute jacobian via coloring
  opts["calc_jac_explicit"] = false
  ctx_residual = (evalResidual,)
  NonlinearSolvers.physicsJac(mesh, sbp, eqn, opts, jac1, ctx_residual)

  # compute jacobian explicitly
  opts["calc_jac_explicit"] = true
  evalJacobian(mesh, sbp, eqn, opts, assembler)

  jac1d = full(jac1)

  @fact maximum(abs(jac1d - jac2)) --> roughly(0.0, atol=1e-14)


  return nothing
end
