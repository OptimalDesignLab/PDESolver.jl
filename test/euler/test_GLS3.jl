"""
  This function tests GLS for a given state loaded in the eqn object.
  This function is not a test function, but is called by test functions.
"""
function test_GLS{Tsol, Tres, Tmsh}(mesh::AbstractMesh{Tmsh}, sbp, eqn::AbstractSolutionData{Tsol, Tres}, opts)

  eqn.params.tau_type = 2
  Dxi = diagm(1./sbp.w)*sbp.Q[:, :, 1]
  Deta = diagm(1./sbp.w)*sbp.Q[:, :, 2]

  # create indices
  idx_range = Array(UnitRange{Int64}, mesh.numNodesPerElement)
  for i=1:mesh.numNodesPerElement
    start_idx = (i-1)*mesh.numDofPerNode + 1
    end_idx = i*mesh.numDofPerNode
    idx_range[i] = copy(start_idx:end_idx)
  end

#    println("idx_range = ", idx_range)

    size_block = mesh.numNodesPerElement*mesh.numDofPerNode
#    println("size_block = ", size_block)

    # get GLS from the code
    fill!(eqn.res, 0.0)
    EulerEquationMod.applyGLS3(mesh, sbp, eqn, opts)  

  # testing: only do one element
  for el =1:mesh.numEl
#    println("testing element ", el)

    # constant mapping elements only
    dxidx = zeros(2,2)
    dxidx_hat_el = sview(mesh.dxidx, :, :, 1, el)
    jac_el = sview(mesh.jac, :, el)
    q_el = reshape(copy(eqn.q[:, :, el]), size_block)

    for i=1:2
      for j=1:2
        dxidx[i,j] = dxidx_hat_el[i, j, 1]*jac_el[1]
      end
    end

    # calculate Dx, Dy
    Dx = dxidx[1,1]*Dxi + dxidx[2, 1]*Deta
    Dy = dxidx[1,2]*Dxi + dxidx[2, 2]*Deta

#    println("Dx = \n", Dx)
#    println("Dy = \n", Dy)

    # create block Dx, Dy
    Dx_tilde = zeros(Tmsh, size_block, size_block)
    Dy_tilde = zeros(Dx_tilde)

   
    for i=1:mesh.numNodesPerElement
      idx_i = idx_range[i]
#      println("idx_i = ", idx_i)
      for j=1:mesh.numNodesPerElement
        idx_j = idx_range[j]
#        println("  idx_j = ", idx_j)
        Dx_tilde[idx_i, idx_j] = Dx[i, j]*eye(mesh.numDofPerNode)
        Dy_tilde[idx_i, idx_j] = Dy[i, j]*eye(mesh.numDofPerNode)
      end
    end

#    println("Dx_tilde = \n", Dx_tilde)
#    println("Dy_tilde = \n", Dy_tilde)

    # create A1 tilde and A2 tilde
    A1_tilde = zeros(Tsol, size_block, size_block)
    A2_tilde = zeros(Tsol, size_block, size_block)

    for i=1:mesh.numNodesPerElement
      idx_i = idx_range[i]
      q_i = q_el[idx_i]

      A1 = sview(A1_tilde, idx_i, idx_i)
      EulerEquationMod.calcA1(eqn.params, q_i, A1)

      A2 = sview(A2_tilde, idx_i, idx_i)
      EulerEquationMod.calcA2(eqn.params, q_i, A2)
    end

#    println("A1 = \n", A1_tilde)
#    println("A2 = \n", A2_tilde)

    # create middle terms, including tau
    middle_tilde = zeros(Tres, size_block, size_block)

    for i=1:mesh.numNodesPerElement

      idx_i = idx_range[i]
      tau = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode)
      EulerEquationMod.getTau(eqn.params, jac_el[i], tau)  # tau number 2

      middle_tilde[idx_i, idx_i] = (sbp.w[i]/jac_el[i])*tau
    end


    # create the operator
    fancy_L  = A1_tilde*Dx_tilde + A2_tilde*Dy_tilde
    gls_operator = fancy_L.'*middle_tilde*fancy_L
    
    @fact isSymmetric(gls_operator, 1e-12) --> true
    #println("max asymmetry = ", maximum(abs(gls_operator - gls_operator.')))

    gls_test = -gls_operator*q_el

    gls_code = reshape(copy(eqn.res[:, :, el]), size_block)

    @fact gls_code --> roughly(gls_test, atol=1e-12)
#    println("gls_test = \n", gls_test)
#    println("gls_code = \n", gls_code)

    trial_term = fancy_L*q_el
    middle_term = middle_tilde*trial_term
    gls_term = -fancy_L.'*middle_term
#    println("element $el trial term = \n", trial_term)
#    println("elemetn $el middle_term = \n", middle_term)
#    println("element $el gls_term = \n", gls_term)

    @fact gls_term --> roughly(gls_test, atol=1e-12)
    # test matrix transpose
    tmp1 = A1_tilde.'*middle_term
    tmp2 = A2_tilde.'*middle_term
    gls_test2 = -(Dx_tilde.'*tmp1 + Dy_tilde.'*tmp2)
    @fact gls_test2 --> roughly(gls_test, atol=1e-12)


  end  # end loop over elements

  return nothing

end  # end function


# run the tests
"""
  Test GLS on a channel flow
"""
function test_gls_channel(mesh, sbp, eqn, opts)
  #=
  include("input_vals_channel.jl")
  arg_dict["solve"] = false
  arg_dict["variable_type"] = :entropy
  f = open("input_vals_channel_gls.jl", "w")
  println(f, "arg_dict = ")
  println(f, arg_dict)
  close(f)
  =#


  facts("----- Testing GLS3 channel -----") do
 #   resize!(ARGS, 1)
 #   ARGS[1] = "input_vals_channel_gls.jl"
 #   include(STARTUP_PATH)
    eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    test_GLS(mesh, sbp, eqn, opts)
  end  # end facts block

  return nothing
end

test_gls_channel_inputname = "input_vals_channel.jl"
test_gls_channel_moddict = Dict{ASCIIString, Any}("solve" => false, "variable_type" => :entropy, "new_fname" => "input_vals_channel_gls")
#test_gls_channel(mesh, sbp, eqn, opts)
add_func3!(EulerTests, test_gls_channel, test_gls_channel_inputname, test_gls_channel_moddict, [TAG_ENTROPYVARS])


"""
  Test GLS on the isentropic vortex for all degree operators
"""
function test_gls_vortex()
  for p = 1:4
    if true
      facts("----- Testing GLS3 p$p  on isentropic vortex -----") do
        # test on isentropic vortex
        include("input_vals_vortex3.jl")
        arg_dict["order"] = p
        arg_dict["solve"] = false
        arg_dict["variable_type"] = :entropy
        #TODO: use make_input function instead
        f = open("input_vals_vortex3_gls.jl", "w")
        print(f, "arg_dict = ")
        println(f, arg_dict)
        close(f)

        resize!(ARGS, 1)
        ARGS[1] = "input_vals_vortex3_gls.jl"
        include(STARTUP_PATH)
        eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
        test_GLS(mesh, sbp, eqn ,opts)
      end  # end facts block
    end
  end

  return nothing
end

#test_gls_vortex(mesh, sbp, eqn, opts)
add_func1!(EulerTests, test_gls_vortex, [TAG_ENTROPYVARS])

"""
  Test finite differencing and complex stepping of the GLS term
"""
function test_gls_fd()
  for p = 1:4
    if true
      facts("----- Performing GLS3 p$p finite difference checks -----") do
        ARGS[1] = "input_vals_vortex3_gls.jl"
        include(STARTUP_PATH)

        arg_dict["order"] = p
        f = open("input_vals_vortex3_gls.jl", "w")
        print(f, "arg_dict = ")
        println(f, arg_dict)
        close(f)
        ARGS[1] = "input_vals_vortex3_gls.jl"
        include(STARTUP_PATH)

        # rescale the problem
        for i = 1:length(eqn.q)
          eqn.q[i] = 1000*eqn.q[i]
        end

        # trick code into only doing the first element
        mesh.numEl = 1


        len = mesh.numDofPerNode*mesh.numNodesPerElement
        eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
        test_GLS(mesh, sbp, eqn, opts)
        jac_fd = zeros(len, len)
        eps_fd = 1e-8
        res0 = copy(reshape(eqn.res[:, :, 1], len))  # use res from previous run
        println("doing finite differences")
        for j=1:mesh.numNodesPerElement
          for i=1:mesh.numDofPerNode
            pos = (j-1)*mesh.numDofPerNode + i
            eqn.q[i, j, 1] += eps_fd
            test_GLS(mesh, sbp, eqn, opts)
            res_ij = copy(reshape(eqn.res[:, :, 1], len))
            jac_fd[:, pos] = (res_ij - res0)/eps_fd
            eqn.q[i, j, 1] -= eps_fd  # undo perturbation
          end
        end

        # do complex step
        println("doing complex step")
        arg_dict["run_type"] = 5
        arg_dict["jac_method"] = 2
        f = open("input_vals_vortex3c_gls.jl", "w")
        print(f, "arg_dict = ")
        println(f, arg_dict)
        close(f)
        ARGS[1] = "input_vals_vortex3c_gls.jl"
        include(STARTUP_PATH)

        # rescale the problem
        for i = 1:length(eqn.q)
          eqn.q[i] = 1000*eqn.q[i]
        end


        # trick code into only doing 1 element
        mesh.numEl = 1

        eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
        eqn.params.tau_type = p

        eps_c = 1e-20
        jac_c = zeros(len, len)
        for j=1:mesh.numNodesPerElement
          for i=1:mesh.numDofPerNode
            pos = (j-1)*mesh.numDofPerNode + i
            eqn.q[i, j, 1] += complex(0, eps_c)
            test_GLS(mesh, sbp, eqn, opts)
            res_ij = copy(reshape(eqn.res[:, :, 1], len))
            jac_c[:, pos] = imag(res_ij)/eps_c
            eqn.q[i, j, 1] -= complex(0, eps_c)  # undo perturbatino
          end
        end

        for j=1:len
          tol = 0.002
          for k = 1:len
            if abs(jac_fd[k,j]) > 1e-4
              @fact abs((jac_c[k, j] - jac_fd[k, j])/jac_c[k,j]) --> less_than(tol)
            end
          end
        end

    #=
        for j=1:12
          println("column $j")
          println("jac_c = \n", jac_c[:, j])
          println("jac_fd = \n", jac_fd[:, j])
          println("jac diff = \n", jac_c[:,j] - jac_fd[:, j])
        end    
    =#
      end  # end facts block
    end  # end if statement
  end  # end p=1:4

  return nothing

end  # end function

#test_gls_fd()
add_func1!(EulerTests, test_gls_fd, [TAG_COMPLEX])
