#using FactCheck
#using ODLCommonTools
#global const STARTUP_PATH = joinpath(Pkg.dir("PDESolver"), "src/solver/advection/startup_advection.jl")

facts("----- Testing Jacobian -----") do
   resize!(ARGS, 1)
   ARGS[1] = "input_vals_8el.jl"
   include(STARTUP_PATH)


  for el = 1:mesh.numEl
    println("----- Doing Finite Differences -----")
    ARGS[1] = "input_vals_8el.jl"
    include(STARTUP_PATH)

    jac_fd = zeros(Float64, 3,3)
    eps_fd = 1e-7
    # calculate jacobian of the first element

    fill!(eqn.res, 0.0)
    AdvectionEquationMod.evalAdvection(mesh, sbp, eqn, opts)
    res_0 = copy(reshape(eqn.res[1, :, el], 3))
    for i=1:3
      eqn.q[1, i, el] += eps_fd
      fill!(eqn.res, 0.0)
      AdvectionEquationMod.evalAdvection(mesh, sbp, eqn, opts)
      res_i = reshape(eqn.res[1, :, el], 3)
      for j=1:3
        jac_fd[j, i] = (res_i[j] - res_0[j])/eps_fd
      end

      #undo perturbation
      eqn.q[1, i, el] -= eps_fd
    end

    # now do complex step
    println("----- Doing Complex step -----")
    include(ARGS[1])
    arg_dict["run_type"] = 5
    f = open("input_vals_8elc.jl", "w")
    println(f, arg_dict)
    close(f)
    ARGS[1] = "input_vals_8elc.jl"
    include(STARTUP_PATH)


    jac_c = zeros(Float64, 3,3)
    eps_c = complex(0, 1e-20)
    for i=1:3
      eqn.q[1, i, el] += eps_c
      fill!(eqn.res, 0.0)
      AdvectionEquationMod.evalAdvection(mesh, sbp, eqn, opts)
      res_i = reshape(eqn.res[1, :, el], 3)
      for j=1:3
        jac_c[j, i] = imag(res_i[j])/abs(eps_c)
      end

      #undo perturbation
      eqn.q[1, i, el] -= eps_c
    end

    @fact jac_c --> roughly(jac_fd, atol=1e-6)
  end

  # back to finite differences
  println("----- Testing Finite Difference Jacobian -----")
  ARGS[1] = "input_vals_8el.jl"
  include(STARTUP_PATH)
  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)

  # needed for calls to NewtonData below
  Tsol = eltype(eqn.q)
  Tres = eltype(eqn.res)

  # now test full jacobian
  newton_data = NonlinearSolvers.NewtonData{Tsol, Tres}(mesh, sbp, eqn, opts)
  fill!(eqn.res, 0.0)
  AdvectionEquationMod.evalAdvection(mesh, sbp, eqn, opts)
  eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
  res_3d0 = copy(eqn.res)
  res_0 = copy(eqn.res_vec)
  jac = zeros(Float64, mesh.numDof, mesh.numDof)
  eps_fd = 1e-7
  fill!(eqn.res, 0.0)
  NonlinearSolvers.calcJacFD(newton_data, mesh, sbp, eqn, opts, AdvectionEquationMod.evalAdvection, res_0, eps_fd, jac)

#  jac_sparse = SparseMatrixCSC(mesh.sparsity_bounds, Float64)
  println("mesh.coloringDistance = ", mesh.coloringDistance)
  println("typeof(mesh) = ", typeof(mesh))
  println("mesh.pertNeighborEls = ", mesh.pertNeighborEls)
  println("mesh.dofs = ", mesh.dofs)
  jac_sparse = SparseMatrixCSC(mesh.sparsity_bnds, Float64)
  println("create jac_sparse")
  fill!(eqn.res, 0.0)
  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
  println("about to calculate jacobian")
  NonlinearSolvers.calcJacobianSparse(newton_data, mesh, sbp, eqn, opts, AdvectionEquationMod.evalAdvection, res_3d0, eps_fd, jac_sparse)
  println("finished calculating jacbian")

  jac_sparsefull = full(jac_sparse)
  jac_diff = jac - jac_sparsefull
  for i=1:mesh.numDof
    for j=1:mesh.numDof
      println("i, j ", i, ", ", j)
      @fact abs(jac_diff[j, i]) --> roughly(0.0, atol=1e-6)
    end
  end

  # back to complex step
  println("----- Testing Complex Step Jacobian -----")
  ARGS[1] = "input_vals_8elc.jl"
  arg_dict["run_type"] = 5  # something screwy is going on because this is necessary
  include(STARTUP_PATH)
  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)

  # now test full jacobian
  newton_data = NonlinearSolvers.NewtonData{Tsol, Tres}(mesh, sbp, eqn, opts)
  jac_c = zeros(Complex128, mesh.numDof, mesh.numDof)
  eps_c = complex(0, 1e-20)
  fill!(eqn.res, 0.0)
  NonlinearSolvers.calcJacobianComplex(newton_data, mesh, sbp, eqn, opts, AdvectionEquationMod.evalAdvection, eps_c, jac_c)

#  jac_csparse = SparseMatrixCSC(mesh.sparsity_bounds, Float64)
  jac_csparse = SparseMatrixCSC(mesh.sparsity_bnds, Float64)
  fill!(eqn.res, 0.0)
  res_3d0 = Array(Float64, 0, 0, 0)
  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
  NonlinearSolvers.calcJacobianSparse(newton_data, mesh, sbp, eqn, opts, AdvectionEquationMod.evalAdvection, res_3d0, eps_c, jac_csparse)


  jac_csparsefull = full(jac_csparse)
  jac_diff = jac_c - jac_csparsefull
  for i=1:mesh.numDof
    for j=1:mesh.numDof
      @fact abs(jac_diff[j, i]) --> roughly(0.0, atol=1e-12)
    end
  end


  # now check FD vs Complex step
  for i=1:mesh.numDof
    for j=1:mesh.numDof
      @fact abs(jac_c[i, j] - jac[i,j]) --> roughly(0.0, atol=1e-6)
    end
  end





end
