#using FactCheck
#global const STARTUP_PATH = joinpath(Pkg.dir("PDESolver"), "src/solver/advection/startup_advection.jl")

facts("----- Testing Jacobian -----") do
   resize!(ARGS, 1)
   ARGS[1] = "input_vals_8el.jl"
   include(STARTUP_PATH)


  for el = 1:mesh.numEl
    println("element ", el)
    println("----- Doing Finite Differences -----")
    ARGS[1] = "input_vals_8el.jl"
    include(STARTUP_PATH)

    jac_fd = zeros(Float64, 3,3)
    eps_fd = 1e-7
    # calculate jacobian of the first element

    fill!(eqn.res, 0.0)
    AdvectionEquationMod.evalAdvection(mesh, sbp, eqn, opts)
    res_0 = copy(reshape(eqn.res[1, :, el], 3))
    println("res_0 = ", res_0)
    for i=1:3
      eqn.q[1, i, el] += eps_fd
      fill!(eqn.res, 0.0)
      AdvectionEquationMod.evalAdvection(mesh, sbp, eqn, opts)
      res_i = reshape(eqn.res[1, :, el], 3)
      println("res_$i = ", res_i)
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
      println("res_$i = ", res_i)
      for j=1:3
        jac_c[j, i] = imag(res_i[j])/abs(eps_c)
      end

      #undo perturbation
      eqn.q[1, i, el] -= eps_c
    end

    println("jac_fd = \n", jac_fd)
    println("jac_c = \n", jac_c)
    println("diff = \n", jac_c - jac_fd)
    @fact jac_c => roughly(jac_fd, atol=1e-6)
  end

end
