"""
  Does basic testing of explicit jacobian calculation in parallel
"""
function test_jac_parallel()

  # SBPOmega
  fname = "input_vals_jac3dp.jl"
  fname2 = "input_vals_jac_tmp.jl"
  mesh, sbp, eqn, opts = run_solver(fname)

  MPI.Barrier(mesh.comm)
  test_jac_parallel_inner(mesh, sbp, eqn, opts)
  opts["preallocate_jacobian_coloring"] = true
  test_jac_parallel_inner(mesh, sbp, eqn, opts, is_prealloc_exact=true, set_prealloc=false)

  return nothing
end


add_func1!(EulerTests, test_jac_parallel, [TAG_SHORTTEST, TAG_JAC]) 


"""
  Does more thorough testing of jacobian calculation in parallel
"""
function test_jac_parallel_long()

  @testset "----- Testing jacobian assembly long -----" begin
    fname = "input_vals_jac3dp.jl"
    fname2 = "input_vals_jac_tmp.jl"

    myrank = MPI.Comm_rank(MPI.COMM_WORLD)

    f = open("testlog_$myrank.dat", "w")

    # SBPGamma
    if myrank == 0
      opts_tmp = read_input_file(fname)
      opts_tmp["operator_type"] = "SBPGamma"
      make_input(opts_tmp, fname2)
    end
    MPI.Barrier(MPI.COMM_WORLD)
    mesh4, sbp4, eqn4, opts4 = run_solver(fname2)


    # SBPDiagonalE
    if myrank == 0
      opts_tmp = read_input_file(fname)
      opts_tmp["operator_type"] = "SBPDiagonalE"
      make_input(opts_tmp, fname2)
    end
    MPI.Barrier(MPI.COMM_WORLD)
    mesh5, sbp5, eqn5, opts5 = run_solver(fname2)

    if myrank == 0
      opts_tmp = read_input_file(fname)
      opts_tmp["operator_type"] = "SBPDiagonalE"
      opts_tmp["use_Minv"] = true
      make_input(opts_tmp, fname2)
    end
    MPI.Barrier(MPI.COMM_WORLD)
    mesh6, sbp6, eqn6, opts6 = run_solver(fname2)

    # SBPGamma ES
    if myrank == 0
      opts_tmp = read_input_file(fname)
      opts_tmp["operator_type"] = "SBPGamma"
      opts_tmp["volume_integral_type"] = 2
      opts_tmp["Volume_flux_name"] = "IRFlux"
      opts_tmp["Flux_name"] = "IRFlux"
      opts_tmp["face_integral_type"] = 2
      opts_tmp["FaceElementIntegral_name"] = "ESLFFaceIntegral"
#      opts_tmp["addBoundaryIntegrals"] = false # DEBUGGING
#      opts_tmp["addVolumeIntegrals"] = false
#      opts_tmp["addFaceIntegrals"] = false
      make_input(opts_tmp, fname2)
    end
    MPI.Barrier(MPI.COMM_WORLD)
    mesh7, sbp7, eqn7, opts7 = run_solver(fname2)


    opts4_tmp = copy(opts4)
    test_jac_parallel_inner(mesh4, sbp4, eqn4, opts4)
    test_jac_homotopy(mesh4, sbp4, eqn4, opts4_tmp)
    test_revm_product(mesh4, sbp4, eqn4, opts4, f)
    test_revq_product(mesh4, sbp4, eqn4, opts4, f)


    test_jac_parallel_inner(mesh5, sbp5, eqn5, opts5)

    # run test twice to make sure arrays are zeroed out correctly
    test_jac_parallel_inner(mesh5, sbp5, eqn5, opts5)

    test_revm_product(mesh5, sbp5, eqn5, opts5, f)

    test_jac_parallel_inner(mesh6, sbp6, eqn6, opts6)

    test_jac_parallel_inner(mesh7, sbp7, eqn7, opts7)
    test_jac_parallel_inner(mesh7, sbp7, eqn7, opts7)

    println(f, "\nrunning test 7")
    test_revm_product(mesh7, sbp7, eqn7, opts7, f)
    test_revq_product(mesh7, sbp7, eqn7, opts7, f)
    close(f)
  end

  return nothing
end

add_func1!(EulerTests, test_jac_parallel_long, [TAG_LONGTEST, TAG_JAC, TAG_TMP]) 

#------------------------------------------------------------------------------
# functions that run individual tests


function rand_realpart(dims...)

  a = rand(Complex128, dims...)
  for i=1:length(a)
    a[i] = real(a[i])
  end

  return a
end



function test_jac_parallel_inner(mesh, sbp, eqn, opts; is_prealloc_exact=true, set_prealloc=true)

  # use a spatially varying solution
  icfunc = EulerEquationMod.ICDict["ICExp"]
  icfunc(mesh, sbp, eqn, opts, eqn.q_vec)
  array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)

  # get the correct differentiated flux function (this is needed because the
  # input file set calc_jac_explicit = false
  eqn.flux_func_diff = EulerEquationMod.FluxDict_diff[opts["Flux_name"]]
  eqn.volume_flux_func_diff = EulerEquationMod.FluxDict_diff[opts["Volume_flux_name"]]

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

  # compute jacobian via coloring
  opts["calc_jac_explicit"] = false
  ctx_residual = (evalResidual,)
  NonlinearSolvers.physicsJac(mesh, sbp, eqn, opts, jac1, ctx_residual)

  # compute jacobian explicitly
  opts["calc_jac_explicit"] = true
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

    @test isapprox( norm(b1 - b2), 0.0) atol=1e-12
  end

  A = getBaseLO(lo2).A
  if typeof(A) <: PetscMat
    matinfo = MatGetInfo(A, PETSc2.MAT_LOCAL)
    if is_prealloc_exact
      @test ( matinfo.nz_unneeded )== 0
    else
      @test  matinfo.nz_unneeded  > 0
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
  array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)

  # get the correct differentiated flux function (this is needed because the
  # input file set calc_jac_explicit = false
  eqn.flux_func_diff = EulerEquationMod.FluxDict_diff["RoeFlux"]

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
#=
  println("diffnorm = ", vecnorm(res1 - res2))
  println("res1 = \n", res1)
  println("res2 = \n", res2)
  println("diff = \n", res1 - res2)
  @assert vecnorm(res1 - res2) < 1e-13
=#
  startSolutionExchange(mesh, sbp, eqn, opts, wait=true)

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
  NonlinearSolvers.physicsJac(mesh, sbp, eqn, opts, jac1, ctx_residual)

  # compute jacobian explicitly
  opts["calc_jac_explicit"] = true

  
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

    @test isapprox( norm(b1 - b2), 0.0) atol=1e-12
  end

  free(lo1)
  free(lo2)
  free(pc1)
  free(pc2)

  return nothing
end


#DEBUGGING
function zeroSharedElements(mesh, sbp, eqn, opts, qvec)

  for i=1:mesh.npeers
    for el in mesh.local_element_lists[i]
      for j=1:mesh.numNodesPerElement
        for k=1:mesh.numDofPerNode
          qvec[mesh.dofs[k, j, el]] = 0
        end
      end
    end
  end

  return nothing
end




function test_revm_product(mesh, sbp, eqn, opts, f::IO)

  srand(1234)  # reproducable results

  h = 1e-20
  pert = Complex128(0, h)

  icfunc = EulerEquationMod.ICDict["ICExp"]
  icfunc(mesh, sbp, eqn, opts, eqn.q_vec)
  array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)
  eqn.q_vec .+= 0.01.*rand(size(eqn.q_vec))  # add a little noise, to make jump across
                                     # interfaces non-zero

  # fields: dxidx, jac, nrm_bndry, nrm_face, coords_bndry

  res_bar = rand_realpart(mesh.numDof)

  dxidx_dot       = rand_realpart(size(mesh.dxidx))
  jac_dot         = rand_realpart(size(mesh.jac))
  nrm_bndry_dot   = rand_realpart(size(mesh.nrm_bndry))
  nrm_face_dot    = rand_realpart(size(mesh.nrm_face_bar))
  coords_bndry_dot = rand_realpart(size(mesh.coords_bndry))
  nrm_sharedface_dot = Array{Array{Complex128, 3}}(mesh.npeers)
  for i=1:mesh.npeers
    nrm_sharedface_dot[i] = rand_realpart(size(mesh.nrm_sharedface[i]))
    #nrm_sharedface_dot[i] = zeros(Complex128, size(mesh.nrm_sharedface[i]))
  end
  nrm_sharedface_dot[1][1] = 1



  zeroBarArrays(mesh)

  mesh.dxidx        .+= pert*dxidx_dot
  mesh.jac          .+= pert*jac_dot
  mesh.nrm_bndry    .+= pert*nrm_bndry_dot
  mesh.nrm_face     .+= pert*nrm_face_dot
  mesh.coords_bndry .+= pert*coords_bndry_dot
  for i=1:mesh.npeers
    mesh.nrm_sharedface[i] .+= pert*nrm_sharedface_dot[i]
  end

  startSolutionExchange(mesh, sbp, eqn, opts)

  fill!(eqn.res, 0)
  evalResidual(mesh, sbp, eqn, opts)
  array3DTo1D(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
  val = sum(imag(eqn.res_vec)/h .* res_bar)

#  val = MPI.Allreduce(val, MPI.SUM, eqn.comm)

  mesh.dxidx        .-= pert*dxidx_dot
  mesh.jac          .-= pert*jac_dot
  mesh.nrm_bndry    .-= pert*nrm_bndry_dot
  mesh.nrm_face     .-= pert*nrm_face_dot
  mesh.coords_bndry .-= pert*coords_bndry_dot
  for i=1:mesh.npeers
    mesh.nrm_sharedface[i] .-= pert*nrm_sharedface_dot[i]
  end



  evalResidual_revm(mesh, sbp, eqn, opts, res_bar)
  val2 = 0
  val2 = sum(mesh.dxidx_bar .* dxidx_dot)              +
         sum(mesh.jac_bar .* jac_dot)                  +
         sum(mesh.nrm_bndry_bar .* nrm_bndry_dot)      +
         sum(mesh.nrm_face_bar .* nrm_face_dot)        +
         sum(mesh.coords_bndry_bar .* coords_bndry_dot)

  for i=1:mesh.npeers
    val2 += sum(mesh.nrm_sharedface_bar[i] .* nrm_sharedface_dot[i])
  end

#  val2 = MPI.Allreduce(val2, MPI.SUM, eqn.comm)

  println(f, "val = ", real(val))
  println(f, "val2 = ", real(val2))
  println(f, "diff = ", real(val - val2))
  println(f, "max dxidx_bar = ", maximum(abs.(mesh.dxidx_bar)))
  println(f, "max jac_bar = ", maximum(abs.(mesh.jac_bar)))
  println(f, "max nrm_bndry_bar = ", maximum(abs.(mesh.nrm_bndry_bar)))
  println(f, "max nrm_face_bar = ", maximum(abs.(mesh.nrm_face_bar)))
  println(f, "max coords_bndry_bar = ", maximum(abs.(mesh.coords_bndry_bar)))
  for i=1:mesh.npeers
    println(f, "sum nrm_sharedface_bar ", i, " = ", sum(mesh.nrm_sharedface_bar[i]))
    println(f, "sum(nrm_sharedface) ", i, " = ", sum(mesh.nrm_sharedface[i]))
    println(f, "sum(q_recv) ", i, " = ", sum(eqn.shared_data[i].q_recv))
  end
  @test abs(val - val2) < 1e-12


  return nothing
end


function test_revq_product(mesh, sbp, eqn, opts, f::IO)

  srand(1234)
  h = 1e-20
  pert = Complex128(0, h)

  icfunc = EulerEquationMod.ICDict["ICExp"]
  icfunc(mesh, sbp, eqn, opts, eqn.q_vec)
#  array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)
  eqn.q_vec .+= 0.01.*rand(size(eqn.q_vec))  # add a little noise, to make jump across
                                     # interfaces non-zero

  # fields: dxidx, jac, nrm_bndry, nrm_face, coords_bndry

  res_vec_bar = rand_realpart(mesh.numDof)
#  res_vec_bar = zeros(Complex128, mesh.numDof)
  q_vec_dot = rand_realpart(mesh.numDof)
#  q_vec_dot = zeros(Complex128, mesh.numDof)
  q_vec_bar = zeros(Complex128, mesh.numDof)

#  zeroSharedElements(mesh, sbp, eqn, opts, res_vec_bar)
#  zeroSharedElements(mesh, sbp, eqn, opts, q_vec_dot)


#=  
  if mesh.myrank == 0
    zeroSharedElements(mesh, sbp, eqn, opts, res_vec_bar)
    zeroSharedElements(mesh, sbp, eqn, opts, q_vec_dot)

    #for el in mesh.local_element_lists[1]
    el = 11
      for j=1:mesh.numNodesPerElement
        for k=1:mesh.numDofPerNode
          res_vec_bar[mesh.dofs[k, j, el]] = j + k
          q_vec_dot[mesh.dofs[k, j, el]] = j + k
        end
      end
#    end
  end


  if mesh.myrank == 1
    zeroSharedElements(mesh, sbp, eqn, opts, res_vec_bar)
  end

  if mesh.myrank == 2
    zeroSharedElements(mesh, sbp, eqn, opts, res_vec_bar)
  end

  if mesh.myrank == 3
    
    zeroSharedElements(mesh, sbp, eqn, opts, res_vec_bar)
    zeroSharedElements(mesh, sbp, eqn, opts, q_vec_dot)


    #for j=1:mesh.numNodesPerElement
    #  for k=1:mesh.numDofPerNode
    #    res_vec_bar[mesh.dofs[k, j, 2]] = j + k + mesh.myrank
    #    q_vec_dot[mesh.dofs[k, j, 2]] = j + k + mesh.myrank
    #  end
    #end
    
  end
=#  

  fill!(eqn.res, 0)
  eqn.q_vec .+= pert*q_vec_dot
  array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)

  if opts["parallel_type"] != 1
    println(f, "starting solution exchange")
    startSolutionExchange(mesh, sbp, eqn, opts)
  end

  evalResidual(mesh, sbp, eqn, opts)
  eqn.q_vec .-= pert*q_vec_dot
  array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)

  array3DTo1D(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
  val = sum(res_vec_bar .* imag(eqn.res_vec)/h)
  val = MPI.Allreduce(val, MPI.SUM, mesh.comm)


  evalResidual_revq(mesh, sbp, eqn, opts, res_vec_bar, q_vec_bar)
  val2 = sum(q_vec_bar .* q_vec_dot)
  val2 = MPI.Allreduce(val2, MPI.SUM, mesh.comm)

  println(f, "val = ", val)
  println(f, "val2 = ", val2)
  println(f, "diff = ", abs(val - val2))
  @test abs(val - val2) < 1e-12

  return nothing
end
