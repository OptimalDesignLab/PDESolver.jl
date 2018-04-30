# test jacobian calculation functions

"""
  Test the jacobian of individual terms.  This only tests the SBP Omega
  operators due to test time limits.
"""
function test_jac_terms()

  fname3 = "input_vals_jac3d.jl"
  mesh, sbp, eqn, opts = run_solver("input_vals_jac2d.jl")
  mesh3, sbp3, eqn3, opts3 = run_solver(fname3)
#=
  # SBPOmega, Petsc Mat
  fname4 = "input_vals_jac_tmp.jl"
  opts_tmp = read_input_file(fname3)
  opts_tmp["jac_type"] = 3
  opts_tmp["operator_type"] = "SBPOmega"
  make_input(opts_tmp, fname4)
  mesh4, sbp4, eqn4, opts4 = run_solver(fname4)
=#



  facts("----- Testing jacobian -----") do
    test_pressure(eqn.params)
    test_pressure(eqn3.params)

    test_eulerflux(eqn.params)
    test_eulerflux(eqn3.params)

    nrm = [0.45, 0.55]
    nrm2 = -nrm

    println("testing all positive eigenvalues")
    func = EulerEquationMod.RoeSolver
    func_diff = EulerEquationMod.RoeSolver_diff
    func2 = EulerEquationMod.calcLFFlux
    func2_diff = EulerEquationMod.calcLFFlux_diff
    q = Complex128[2.0, 3.0, 4.0, 7.0]
    qg = q + 1
    test_ad_inner(eqn.params, q, qg, nrm, func, func_diff)
    test_lambda(eqn.params, q, nrm)
    test_lambdasimple(eqn.params, q, qg, nrm)
    test_ad_inner(eqn.params, q, qg, nrm, func2, func2_diff)
    # make sure arrays are zerod out
    test_ad_inner(eqn.params, q, qg, nrm, func2, func2_diff)

    println("testing all negative eigenvalues")
    q = Complex128[2.0, 3.0, 4.0, 7.0]
    qg = q + 1
    test_ad_inner(eqn.params, q, qg, nrm2, func, func_diff)
    test_lambda(eqn.params, q, nrm2)
    test_lambdasimple(eqn.params, q, qg, nrm2)


    println("testing lambda1 entropy fix")
    q = Complex128[1.1, 0.47, 0.53, 2.2]
    qg = q + 1
    test_ad_inner(eqn.params, q, qg, nrm, func, func_diff)
   
    println("testing lambda2 entropy fix")
    q = Complex128[1.1, -1.32, -1.34, 2.2] 
    qg = q + 1
    test_ad_inner(eqn.params, q, qg, nrm, func, func_diff)
    test_ad_inner(eqn.params, q, qg, nrm2, func, func_diff)

    println("testing lambda3 entropy fix")
    q = Complex128[1.1, -0.42, -0.45, 2.2]
    qg = q + 1
    test_ad_inner(eqn.params, q, qg, nrm, func, func_diff)
    test_ad_inner(eqn.params, q, qg, nrm2, func, func_diff)


    nrm = [0.45, 0.55, 0.65]
    nrm2 = -nrm

    println("testing all positive eigenvalues")
    q = Complex128[2.0, 3.0, 4.0, 5.0, 13.0]
    qg = q + 1
    test_ad_inner(eqn3.params, q, qg, nrm, func, func_diff)
    test_lambda(eqn3.params, q, nrm)
    test_lambdasimple(eqn3.params, q, qg, nrm)

    println("testing all negative eigenvalues")
    q = Complex128[2.0, 3.0, 4.0, 5.0, 13.0]
    qg = q + 1
    test_ad_inner(eqn3.params, q, qg, nrm2, func, func_diff)
    test_lambda(eqn3.params, q, nrm2)
    test_lambdasimple(eqn3.params, q, qg, nrm2)


    println("testing lambda1 entropy fix")
    # lambda1 entropy fix active
    q = Complex128[1.05, -1.1, -1.2, -1.3, 2.5] 
    qg = q + 1
    test_ad_inner(eqn3.params, q, qg, nrm, func, func_diff)
    test_ad_inner(eqn3.params, q, qg, nrm2, func, func_diff)

    println("testing lambda2 entropy fix")
    # lambda1 entropy fix active
    q = Complex128[1.05, 0.9, 1.2, -1.3, 2.5]
    qg = q + 1
    test_ad_inner(eqn3.params, q, qg, nrm, func, func_diff)
    test_ad_inner(eqn3.params, q, qg, nrm2, func, func_diff)

    println("testing lambda3 entropy fix")
    # lambda3 entropy fix active
    q = Complex128[1.05, -0.52, -0.47, -0.36, 8.5]
    qg = q + 1
    test_ad_inner(eqn3.params, q, qg, nrm, func, func_diff)
    test_ad_inner(eqn3.params, q, qg, nrm2, func, func_diff)

    
    println("\ntesting jac assembly 2d")
    test_jac_assembly(mesh, sbp, eqn, opts)
    opts_tmp = copy(opts)
    test_jac_homotopy(mesh, sbp, eqn, opts_tmp)

    println("\ntesting jac assembly 3d")
    test_jac_assembly(mesh3, sbp3, eqn3, opts3)
    

  end
  return nothing
end


add_func1!(EulerTests, test_jac_terms, [TAG_SHORTTEST, TAG_JAC])


"""
  Tests assembling the jacobian of all the different operators
"""
function test_jac_terms_long()

  facts("----- Testing additional Jacobian calculation -----") do

    fname3 = "input_vals_jac3d.jl"
    # SBPGamma, Petsc Mat
    fname4 = "input_vals_jac_tmp.jl"
    opts_tmp = read_input_file(fname3)
    opts_tmp["jac_type"] = 3
    make_input(opts_tmp, fname4)
    mesh4, sbp4, eqn4, opts4 = run_solver(fname4)


    # SBPDiagonalE, Petsc Mat
    fname4 = "input_vals_jac_tmp.jl"
    opts_tmp = read_input_file(fname3)
    opts_tmp["jac_type"] = 3
    opts_tmp["operator_type"] = "SBPDiagonalE"
    opts_tmp["order"] = 2
    opts_tmp["write_dofs"] = true
    make_input(opts_tmp, fname4)
    mesh5, sbp5, eqn5, opts5 = run_solver(fname4)


    # SBPDiagonalE, SparseMatrixCSC
    fname4 = "input_vals_jac_tmp.jl"
    opts_tmp = read_input_file(fname3)
    opts_tmp["jac_type"] = 2
    opts_tmp["operator_type"] = "SBPDiagonalE"
    make_input(opts_tmp, fname4)
    mesh6, sbp6, eqn6, opts6 = run_solver(fname4)

    # SBPDiagonalE, Petsc Mat, use_Minv
    fname4 = "input_vals_jac_tmp.jl"
    opts_tmp = read_input_file(fname3)
    opts_tmp["jac_type"] = 3
    opts_tmp["operator_type"] = "SBPDiagonalE"
    opts_tmp["order"] = 2
    opts_tmp["use_Minv"] = true
    make_input(opts_tmp, fname4)
    mesh7, sbp7, eqn7, opts7 = run_solver(fname4)

    # SBPOmega, Petsc Mat
    fname4 = "input_vals_jac_tmp.jl"
    opts_tmp = read_input_file(fname3)
    opts_tmp["jac_type"] = 3
    opts_tmp["operator_type"] = "SBPOmega"
    make_input(opts_tmp, fname4)
    mesh8, sbp8, eqn8, opts8 = run_solver(fname4)

    # test various matrix and operator combinations
    println("testing mode 4")
    test_jac_general(mesh4, sbp4, eqn4, opts4)

    opts4["preallocate_jacobian_coloring"] = true
    test_jac_general(mesh4, sbp4, eqn4, opts4, is_prealloc_exact=false, set_prealloc=false)
    
    println("testing mode 5")
    test_jac_general(mesh5, sbp5, eqn5, opts5)
    # run the test twice to make sure the arrays are zeroed out properly
    println("testing mode 5 twice")
    test_jac_general(mesh5, sbp5, eqn5, opts5)

    println("testing mode 6")
    test_jac_general(mesh6, sbp6, eqn6, opts6)
 
    println("testing mode 7")
    test_jac_general(mesh7, sbp7, eqn7, opts7)
 
    println("testing mode 8")
    test_jac_general(mesh8, sbp8, eqn8, opts8)
  
    opts4["preallocate_jacobian_coloring"] = true
    test_jac_general(mesh8, sbp8, eqn8, opts8, is_prealloc_exact=true, set_prealloc=false)

  end

  return nothing
end

add_func1!(EulerTests, test_jac_terms_long, [TAG_LONGTEST, TAG_JAC])


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

function test_lambda{Tdim}(params::AbstractParamType{Tdim}, qL::AbstractVector,
                           nrm::AbstractVector)


  lambda_dot = zeros(qL)
  EulerEquationMod.getLambdaMax_diff(params, qL, nrm, lambda_dot)

  lambda_dot2 = zeros(qL)
  h=1e-20
  pert = Complex128(0, h)
  for i=1:length(lambda_dot)
    qL[i] += pert
    lambda_dot2[i] = imag(EulerEquationMod.getLambdaMax(params, qL, nrm))/h
    qL[i] -= pert
  end

  @fact norm(lambda_dot - lambda_dot2) --> roughly(0.0, atol=1e-13)


  return nothing
end

function test_lambdasimple{Tdim}(params::AbstractParamType{Tdim}, qL::AbstractVector,
                                 qR::AbstractVector,
                                 nrm::AbstractVector)


  lambda_dotL = zeros(qL)
  lambda_dotR = zeros(qL)
  EulerEquationMod.getLambdaMaxSimple_diff(params, qL, qR, nrm, lambda_dotL, lambda_dotR)

  lambda_dotL2 = zeros(qL)
  lambda_dotR2 = zeros(qL)
  h = 1e-20
  pert = Complex128(0, h)
  for i=1:length(lambda_dotL)
    qL[i] += pert
    lambda_dotL2[i] = imag(EulerEquationMod.getLambdaMaxSimple(params, qL, qR, nrm))/h
    qL[i] -= pert
  end

  for i=1:length(lambda_dotL)
    qR[i] += pert
    lambda_dotR2[i] = imag(EulerEquationMod.getLambdaMaxSimple(params, qL, qR, nrm))/h
    qR[i] -= pert
  end


  @fact norm(lambda_dotL - lambda_dotL2) --> roughly(0.0, atol=1e-13)
  @fact norm(lambda_dotR - lambda_dotR2) --> roughly(0.0, atol=1e-13)


  return nothing
end



"""
  Test a differentiated numerical flux function via complex step of the
  original function

  **Inputs**

   * params: a Params object
   * qL: left state
   * qR: right state
   * nrm: normal vector
   * func: the original function
   * func_diff: the differentiated version
"""
function test_ad_inner{Tdim}(params::AbstractParamType{Tdim}, qL, qR, nrm,
                             func, func_diff, output=false)

  # compute jacobian with complex step and AD, compare results

  numDofPerNode = length(qL)

  aux_vars = Complex128[0.0]
  F = zeros(Complex128, numDofPerNode)

  resL = zeros(Complex128, numDofPerNode, numDofPerNode)
  resR = zeros(Complex128, numDofPerNode, numDofPerNode)
  h = 1e-20
  pert = Complex128(0, h)
  for j=1:numDofPerNode
    qL[j] += pert
    func(params, qL, qR, aux_vars, nrm, F)

    for k=1:numDofPerNode
      resL[k, j] = imag(F[k])/h
    end

    qL[j] -= pert
  end

  for j=1:numDofPerNode
    qR[j] += pert
    func(params, qL, qR, aux_vars, nrm, F)

    for k=1:numDofPerNode
      resR[k, j] = imag(F[k])/h
    end

    qR[j] -= pert
  end

  # AD version
  resL2 = zeros(resL)
  resR2 = zeros(resR)

  func_diff(params, qL, qR, aux_vars, nrm, resL2, resR2)

  if output
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
  end


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
  
  jac1 = SparseMatrixCSC(mesh, Float64, COLORING, LinearSolvers.getFaceType(mesh.sbpface))
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


"""
  Test the entire jacobian assembly, for any type of jacobian matrix

  is_prealloc_exact: test that the jacobian preallocation is exact (Petsc only)
  set_prealloc: if true, set the preallocation of the jacobian to be tight
                for the explicitly computed jacobian
                if false, use the value currently in the dictionary
"""
function test_jac_general(mesh, sbp, eqn, opts; is_prealloc_exact=true, set_prealloc=true)

  # use a spatially varying solution
  icfunc = EulerEquationMod.ICDict["ICExp"]
  icfunc(mesh, sbp, eqn, opts, eqn.q_vec)
  disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)

  # get the correct differentiated flux function (this is needed because the
  # input file set calc_jac_explicit = false
  eqn.flux_func_diff = EulerEquationMod.FluxDict_diff[opts["Flux_name"]]


  startSolutionExchange(mesh, sbp, eqn, opts)

  opts["calc_jac_explicit"] = false
  pc1, lo1 = NonlinearSolvers.getNewtonPCandLO(mesh, sbp, eqn, opts)

  opts["calc_jac_explicit"] = true
  val_orig = opts["preallocate_jacobian_coloring"]
  if set_prealloc
    opts["preallocate_jacobian_coloring"] = false
  end
  pc2, lo2 = NonlinearSolvers.getNewtonPCandLO(mesh, sbp, eqn, opts)
  opts["preallocate_jacobian_coloring"] = val_orig

  jac1 = getBaseLO(lo1).A
  jac2 = getBaseLO(lo2).A

  assembler = NonlinearSolvers._AssembleElementData(getBaseLO(lo2).A, mesh, sbp, eqn, opts)

  opts["calc_jac_explicit"] = false
  println("calculating regular jacobian"); flush(STDOUT)
  println(STDERR, "calculating regular jacobian"); flush(STDERR)
  ctx_residual = (evalResidual,)
  NonlinearSolvers.physicsJac(mesh, sbp, eqn, opts, jac1, ctx_residual)

  # compute jacobian explicitly
  opts["calc_jac_explicit"] = true
  println("calculating explicit jacobian"); flush(STDOUT)
  println(STDERR, "calculating explicit jacobian"); flush(STDERR)

  evalJacobian(mesh, sbp, eqn, opts, assembler)

  assembly_begin(jac1, MAT_FINAL_ASSEMBLY)
  assembly_begin(jac2, MAT_FINAL_ASSEMBLY)

  # multiply against a random vector to make sure the jacobian is
  # the same
  for i=1:10
    x = rand(PetscScalar, mesh.numDof)
    b1 = zeros(PetscScalar, mesh.numDof)
    b2 = zeros(PetscScalar, mesh.numDof)

    t = 0.0
    applyLinearOperator(lo1, mesh, sbp, eqn, opts, ctx_residual, t, x, b1)
    applyLinearOperator(lo2, mesh, sbp, eqn, opts, ctx_residual, t, x, b2)

    @fact norm(b1 - b2) --> roughly(0.0, atol=1e-12)
  end

  A = getBaseLO(lo2).A
  if typeof(A) <: PetscMat
    matinfo = MatGetInfo(A, PETSc2.MAT_LOCAL)
    if is_prealloc_exact
      @fact matinfo.nz_unneeded --> 0
    else
      @fact matinfo.nz_unneeded --> greater_than(0)
    end
      

  end

  free(lo1)
  free(lo2)
  free(pc1)
  free(pc2)

  return nothing
end

function test_jac_homotopy(mesh, sbp, eqn, opts)

  println("\nTesting homotopy jacobian")

  # use a spatially varying solution
  icfunc = EulerEquationMod.ICDict["ICExp"]
  icfunc(mesh, sbp, eqn, opts, eqn.q_vec)
  disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)

  # get the correct differentiated flux function (this is needed because the
  # input file set calc_jac_explicit = false
  eqn.flux_func_diff = EulerEquationMod.FluxDict_diff["RoeFlux"]
  opts["homotopy_addBoundaryIntegrals"] = true
#=
  res1 = zeros(eqn.res)
  res2 = zeros(eqn.res)
  println("\ncomputing regular homotopy dissipation")
#  EulerEquationMod.calcHomotopyDiss(mesh, sbp, eqn, opts, res1)
  println("\ncomputing new homotopy dissipation")
  h = 1e-20
  pert = Complex128(0, h)
  eqn.q[1] += pert
  EulerEquationMod.calcHomotopyDiss(mesh, sbp, eqn, opts, res2)
  eqn.q[1] -= pert
=#
#=
  println("diffnorm = ", vecnorm(res1 - res2))
  println("res1 = \n", res1)
  println("res2 = \n", res2)
  println("diff = \n", res1 - res2)
  @assert vecnorm(res1 - res2) < 1e-13
=#
  startSolutionExchange(mesh, sbp, eqn, opts)

  println("constructing first operator")
  opts["calc_jac_explicit"] = false
  pc1, lo1 = NonlinearSolvers.getHomotopyPCandLO(mesh, sbp, eqn, opts)

  println("constructing second operator")
  opts["calc_jac_explicit"] = true
  pc2, lo2 = NonlinearSolvers.getHomotopyPCandLO(mesh, sbp, eqn, opts)

  jac1 = getBaseLO(lo1).A
  jac2 = getBaseLO(lo2).A

  assembler = NonlinearSolvers._AssembleElementData(getBaseLO(lo2).A, mesh, sbp, eqn, opts)

  function _evalHomotopy(mesh, sbp, eqn, opts, t)
    evalHomotopy(mesh, sbp, eqn, opts, eqn.res, t)
  end

  ctx_residual = (_evalHomotopy,)
  println("\nevaluating jacobians")

  opts["calc_jac_explicit"] = false
  println("calculating regular jacobian"); flush(STDOUT)
  println(STDERR, "calculating regular jacobian"); flush(STDERR)
  NonlinearSolvers.physicsJac(mesh, sbp, eqn, opts, jac1, ctx_residual)

  # compute jacobian explicitly
  opts["calc_jac_explicit"] = true
  println("calculating explicit jacobian"); flush(STDOUT)
  println(STDERR, "calculating explicit jacobian"); flush(STDERR)

  
  evalHomotopyJacobian(mesh, sbp, eqn, opts, assembler, lo2.lambda)

  assembly_begin(jac1, MAT_FINAL_ASSEMBLY)
  assembly_begin(jac2, MAT_FINAL_ASSEMBLY)

  # multiply against a random vector to make sure the jacobian is
  # the same
  for i=1:10
    x = rand(PetscScalar, mesh.numDof)
    b1 = zeros(PetscScalar, mesh.numDof)
    b2 = zeros(PetscScalar, mesh.numDof)

    t = 0.0
    applyLinearOperator(lo1, mesh, sbp, eqn, opts, ctx_residual, t, x, b1)
    applyLinearOperator(lo2, mesh, sbp, eqn, opts, ctx_residual, t, x, b2)

    @fact norm(b1 - b2) --> roughly(0.0, atol=1e-12)
  end

  free(lo1)
  free(lo2)
  free(pc1)
  free(pc2)

  println("finished testing Homotopy operators")
  return nothing
end
