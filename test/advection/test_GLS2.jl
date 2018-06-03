# test the second GLS implementation

"""
  Compute GLS in an inefficent way and compare to the efficient implementation
  in the code.  This is not a test function, but is called by test functions.
"""
function getGLS_test(mesh, sbp, eqn, opts)
#    println("----- Entered get GLS_test -----")
  Tmsh = eltype(mesh.dxidx)
  Tsol = eltype(eqn.q)
  Tres = eltype(eqn.res)

  for i=1:1  # DEBUGGING
    alpha_x = eqn.params.alpha_x
    alpha_y = eqn.params.alpha_y
    # do only element 1
    Dxi = diagm(1./sbp.w)*sbp.Q[:, :, 1]
    Deta = diagm(1./sbp.w)*sbp.Q[:, :, 2]

    dxidx_true = mesh.dxidx[:, :, 1, i]*mesh.jac[1,i]
    jac = mesh.jac[:, i]
    Dx = Dxi*dxidx_true[1, 1] + Deta*dxidx_true[2, 1]
    Dy = Dxi*dxidx_true[1, 2] + Deta*dxidx_true[2, 2]
    q = [eqn.q[1,1,i], eqn.q[1,2,i], eqn.q[1,3,i]]
    tau = zeros(Tres, 3,3)
    tau[1,1] = AdvectionEquationMod.getTau(alpha_x, alpha_y, jac[1], mesh.min_node_dist)
    tau[2,2] = AdvectionEquationMod.getTau(alpha_x, alpha_y, jac[2], mesh.min_node_dist)
    tau[3,3] = AdvectionEquationMod.getTau(alpha_x, alpha_y, jac[3], mesh.min_node_dist)


    tmp = alpha_x*Dx*q + alpha_y*Dy*q
    tmp2 = tau*diagm(sbp.w)*diagm(1./jac)*tmp
    gls_test = -(alpha_x*Dx.'*tmp2 + alpha_y*Dy.'*tmp2)

    fill!(eqn.res, 0.0)
    AdvectionEquationMod.applyGLS2(mesh, sbp, eqn, opts, eqn.src_func)

    gls_code = reshape(eqn.res[:, :, i], 3)
    @test isapprox( gls_code, gls_test) atol=1e-14
  end

#    println("----- finished get GLS_test -----")

  return nothing

end  # end function

"""
  Test GLS against reference implemntation, also tau for some cases
  where it is known exactly.
"""
function test_GLS2_term()
  @testset "----- Testing GLS2 -----" begin

    ARGS[1] = "input_vals_GLS2.jl"
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])

    eqn.params.alpha_x = 2.0
    eqn.params.alpha_y = 1.0
    dxidx_true = mesh.dxidx[:, :, 1, 2]*mesh.jac[1,2]
    println("dxidx_true = ", dxidx_true)
    tau_code = AdvectionEquationMod.getTau(2.0, 1.0, dxidx_true, 2)
    @test isapprox( tau_code, 1/sqrt(8)) atol=1e-14


    fill!(eqn.q, 1.0)
    fill!(eqn.res, 0.0)
    AdvectionEquationMod.applyGLS2(mesh, sbp, eqn, opts, eqn.src_func)
    vals_code = reshape(eqn.res[:, :, 1], 3)
    @test isapprox( vals_code, zeros(Float64, 3)) atol=1e-14


    
    
    
    # set q to something interesting
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        x = mesh.coords[1, j, i]
        y = mesh.coords[2, j, i]
        eqn.q[1, j, i] = sin(x)*cos(y)
      end
    end

    # set source term
    eqn.src_func = AdvectionEquationMod.SRCDict["SRC0"]
    getGLS_test(mesh, sbp, eqn, opts)

    eqn.q[1, :, 1] = [1., 0, 0]

    fill!(eqn.res, 0.0)
    AdvectionEquationMod.applyGLS2(mesh, sbp, eqn, opts, eqn.src_func)
    res_tmp = reshape(eqn.res[1, :, 1], 3)


    # now do complex step
    println("----- Doing Complex Step -----")
    ARGS[1] = "input_vals_GLS2.jl"
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])
    arg_dict["run_type"] = 5  # complex step run
    arg_dict["jac_method"] = 2  # complex step run
    f = open("input_vals_GLS2c.jl", "w")
    println(f, arg_dict)
    close(f)
    ARGS[1] = "input_vals_GLS2c.jl"
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])

    # set q to something interesting
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        x = mesh.coords[1, j, i]
        y = mesh.coords[2, j, i]
        eqn.q[1, j, i] = sin(x)*cos(y)
      end
    end

  end  # end facts block 

  return nothing
end  # end function

#test_GLS2_term()
add_func1!(AdvectionTests, test_GLS2_term, [TAG_SHORTTEST])

"""
  Finite difference and complex step test for GLS term
"""
function test_GLS2_jac()
  # finite difference checks
  @testset "----- Performing GLS2 Finite Difference Checks -----" begin

    ARGS[1] = "input_vals_GLS2.jl"
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])
    eqn.src_func = AdvectionEquationMod.SRCDict["SRC0"]
    eqn.params.alpha_x = 2.0
    eqn.params.alpha_y = 1.0
    # set q to something interesting
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        x = mesh.coords[1, j, i]
        y = mesh.coords[2, j, i]
        eqn.q[1, j, i] = sin(x)*cos(y)
      end
    end

    fill!(eqn.res, 0.0)

    jac_fd = zeros(Float64, 3, 3)

    eps_fd = 1e-7
    # do initial calculation
    fill!(eqn.res, 0.0)
    AdvectionEquationMod.applyGLS2(mesh, sbp, eqn, opts, eqn.src_func)
    res_0 = copy(reshape(eqn.res[1, :, 1], 3))
    for i=1:3
      eqn.q[1, i, 1] += eps_fd
      fill!(eqn.res, 0.0)
      getGLS_test(mesh, sbp, eqn, opts)
      res_i = reshape(eqn.res[1, :, 1], 3)
      for j=1:3
        jac_fd[j, i] = (res_i[j] - res_0[j])/eps_fd
      end

      # undo perturbation
      eqn.q[1, i, 1] -= eps_fd
    end

    eqn.q[1, :, 1] = [1., 0, 0]

    fill!(eqn.res, 0.0)
    AdvectionEquationMod.applyGLS2(mesh, sbp, eqn, opts, eqn.src_func)
    res_tmp = reshape(eqn.res[1, :, 1], 3)


    # now do complex step
    println("----- Doing Complex Step -----")
    ARGS[1] = "input_vals_GLS2.jl"
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])
    arg_dict["run_type"] = 5  # complex step run
    arg_dict["jac_method"] = 2
    f = open("input_vals_GLS2c.jl", "w")
    println(f, arg_dict)
    close(f)
    ARGS[1] = "input_vals_GLS2c.jl"
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])

    # set q to something interesting
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        x = mesh.coords[1, j, i]
        y = mesh.coords[2, j, i]
        eqn.q[1, j, i] = sin(x)*cos(y)
      end
    end

    eqn.params.alpha_x = 2.0
    eqn.params.alpha_y = 1.0

    # now we have complex equation object
    jac_c = zeros(Float64, 3,3)
    eps_complex = complex(0, 1e-20)
    for i=1:3 
      eqn.q[1, i, 1] += eps_complex
      fill!(eqn.res, 0.0)
      getGLS_test(mesh, sbp, eqn, opts)  # do checks on complex version
  #    AdvectionEquationMod.applyGLS2(mesh, sbp, eqn, opts)
      res_i = reshape(eqn.res[1, :, 1], 3)
      for j=1:3
        jac_c[j, i] = imag(res_i[j])/abs(eps_complex)
      end

      eqn.q[1, i, 1] -= eps_complex
    end

    @test isapprox( jac_c, jac_fd) atol=1e-6

  end  # end facts block

  return nothing
end  # end function

#test_GLS2_jac()
add_func1!(AdvectionTests, test_GLS2_jac, [TAG_SHORTTEST])
