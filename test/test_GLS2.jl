push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/euler"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))

using PDESolver
#using Base.Test
using FactCheck
using ODLCommonTools
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
using EulerEquationMod
using ForwardDiff
using NonlinearSolvers   # non-linear solvers
using ArrayViews

global const STARTUP_PATH = joinpath(Pkg.dir("PDESolver"), "src/solver/euler/startup.jl")


facts ("----- Testing GLS2 -----") do

  include("input_vals_channel.jl")
  arg_dict["solve"] = false
  f = open("input_vals_channel_gls.jl", "w")
  println(f, arg_dict)
  close(f)

  resize!(ARGS, 1)
  ARGS[1] = "input_vals_channel_gls.jl"
  include(STARTUP_PATH)
  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)

  
  function test_gls(mesh, sbp, eqn, opts)
   
    Tsol = eltype(eqn.q)
    Tmsh = eltype(mesh.dxidx)
    Tres = eltype(eqn.res)
    # do GLS for first element only, compare the slow, explicit version 
    # with the fast index notation one
    q = zeros(Tsol, 12)
    q[1:4] = eqn.q[:, 1, 1]
    q[5:8] = eqn.q[:, 2, 1]
    q[9:12] = eqn.q[:, 3, 1]
    dxidx_hat = mesh.dxidx[:, :, :, 1]
    jac = mesh.jac[:, 1]
    dxidx = copy(dxidx_hat)
    for i=1:3
      dxidx[:, :, i] = dxidx[:, :, i]*jac[i]
    end

    function getTau_test(params, q, dxidx)
      A1 = zeros(Tsol, 4,4)
      A2 = zeros(Tsol, 4,4)
      A0inv = zeros(Tsol, 4,4)
      tau = zeros(Tres, 12, 12)
      range_idx = (1:4, 5:8, 9:12)
      for i=1:3
        idx_i = range_idx[i]
        tau_i = zeros(Tres, 4,4)
        q_i = q[idx_i]
        dxidx_i = dxidx[:, :, i] 
        EulerEquationMod.calcA1(params, q_i, A1)
        EulerEquationMod.calcA2(params, q_i, A2)

        tau_i[:, :] += (dxidx_i[1, 1]*dxidx_i[1, 1] + dxidx[2, 1]*dxidx[2, 1])*A1*A1
        
        tau_i[:, :] += (dxidx_i[1, 1]*dxidx_i[1, 2] + dxidx[2, 1]*dxidx[2, 2])*A1*A2

        tau_i[:, :] += (dxidx_i[1, 2]*dxidx_i[1, 1] + dxidx[2, 2]*dxidx[2, 1])*A2*A1
        tau_i[:, :] += (dxidx_i[1, 2]*dxidx_i[1, 2] + dxidx[2, 2]*dxidx[2, 2])*A2*A2

        EulerEquationMod.calcA0Inv(params, q_i, A0inv)

        D, V = eig(tau_i)
        D2 = diagm(D.^(-0.5))
        new_tau = V*D2*inv(V)
        new_tau2 = A0inv*new_tau
        tau[idx_i, idx_i] = real(new_tau2)

        # check agains the tau calculation in EulerEquationMod
        A_mat = zeros(Tsol, 4,4,2)
        A_mat[:, :, 1] = A1
        A_mat[:, :, 2] = A2


        tau_old = zeros(Tres, 4, 4)
        EulerEquationMod.getTau(params, A_mat, dxidx_i, tau_old)

        @fact tau_old => roughly(new_tau2, atol=1e-14)
      end

      return tau
    end

    tau_tilde = getTau_test(eqn.params, q, dxidx_hat)

    function getDtilde(sbp, dir::Integer)
      range_idx = (1:4, 5:8, 9:12)
      D = inv(diagm(sbp.w))*sbp.Q[:, :, dir]

      Dtilde = zeros(12, 12)
      for i=1:3
        idx_i = range_idx[i]
        for j=1:3
          idx_j = range_idx[j]
          Dtilde[idx_i, idx_j] = D[i,j]*eye(4,4)
        end
      end

      return Dtilde
    end

    D_tilde_xi = getDtilde(sbp, 1)
    D_tilde_eta = getDtilde(sbp, 2)

    # for steady channel, test Dtilde_x*q = 0?

    function getAtilde(params, q, dir)

      range_idx = (1:4, 5:8, 9:12)
      Atilde = zeros(Tsol, 12, 12)
      A = zeros(Tsol, 4,4)
      for i=1:3
        idx_i = range_idx[i]
        if dir == 1
          EulerEquationMod.calcA1(params, q[idx_i], A)
        else
          EulerEquationMod.calcA2(params, q[idx_i], A)
        end

        Atilde[idx_i, idx_i] = A
      end

      return Atilde
    end

    A1tilde = getAtilde(eqn.params, q, 1)
    A2tilde = getAtilde(eqn.params, q, 2)

    function getdxidx_tilde(dxidx, k, j)

      range_idx = (1:4, 5:8, 9:12)
      dxidx_tilde = zeros(Tmsh, 12, 12)

      for i=1:3
        idx_i = range_idx[i]
        dxidx_tilde[idx_i, idx_i] = dxidx[k, j, i]*eye(4,4)
      end

      return dxidx_tilde
    end

    # the versions not scaled 1/|J|, for the weighting space
    dxidxhat_tilde_11 = getdxidx_tilde(dxidx_hat, 1, 1)
    dxidxhat_tilde_12 = getdxidx_tilde(dxidx_hat, 1, 2)
    dxidxhat_tilde_21 = getdxidx_tilde(dxidx_hat, 2, 1)
    dxidxhat_tilde_22 = getdxidx_tilde(dxidx_hat, 2, 2)

    # the versions scaled by 1/|J|, for the trial space
    dxidx_tilde_11 = getdxidx_tilde(dxidx, 1, 1)
    dxidx_tilde_12 = getdxidx_tilde(dxidx, 1, 2)
    dxidx_tilde_21 = getdxidx_tilde(dxidx, 2, 1)
    dxidx_tilde_22 = getdxidx_tilde(dxidx, 2, 2)

    function getH(sbp)
      H = zeros(12, 12)
      range_idx = (1:4, 5:8, 9:12)

      for i=1:3
        idx_i = range_idx[i]
        H[idx_i, idx_i] = sbp.w[i]*eye(4)
      end

      return H
    end

    H_tilde = getH(sbp)

    # now compute the whole GLS term
    weighting_term = A1tilde*(dxidx_tilde_11*D_tilde_xi + dxidx_tilde_21*D_tilde_eta) +
                     A2tilde*(dxidx_tilde_12*D_tilde_xi + dxidx_tilde_22*D_tilde_eta)

    trial_term = A1tilde*(dxidxhat_tilde_11*D_tilde_xi*q + dxidxhat_tilde_21*D_tilde_eta*q) + 
                 A2tilde*(dxidxhat_tilde_12*D_tilde_xi*q + dxidxhat_tilde_22*D_tilde_eta*q)
    gls_test = weighting_term.'*H_tilde*tau_tilde*trial_term

    # now compute it in the code
    print("\n\n")
    fill!(eqn.res, 0.0)
    EulerEquationMod.applyGLS2(mesh, sbp, eqn, opts)
    gls_code = reshape(eqn.res[:, :, 1], 12)
#=
    println("gls_test = ", gls_test)
    println("gls_code = ", gls_code)
=#
    @fact gls_code => roughly(gls_test, atol=1e-14)
#    println("gls_code - gls_test = ", gls_code - gls_test)
#    print("\n\n")
#=
    println("\nPrinting intermediate quantities:\n")
    println("q = ", q)
    println("dxidx_hat = ", dxidx_hat)
    println("Dxi = ", D_tilde_xi)
    println("Deta = ", D_tilde_eta)
    println("A1tilde = ", A1tilde)
    println("A2tilde = ", A2tilde)
    println("dxidx_tilde_11 = ", dxidx_tilde_11)
    println("dxidx_tilde_12 = ", dxidx_tilde_12)
    println("dxidx_tilde_21 = ", dxidx_tilde_21)
    println("dxidx_tilde_22 = ", dxidx_tilde_22)
    println("qxi = ", D_tilde_xi*q)
    println("qeta = ", D_tilde_eta*q)
    println("dxidx = ", dxidx)
=#

    range_idx = (1:4, 5:8, 9:12)
    res_test = zeros(Tres, 12)
    for i=1:3
      idx_i = range_idx[i]
      for j=1:3
        idx_j = range_idx[j]
        # trial space
        tmp1 = D_tilde_xi[idx_j, :]*q
        tmp2 = D_tilde_eta[idx_j, :]*q
        tmp3 = dxidxhat_tilde_11[idx_j, idx_j]*tmp1 + dxidxhat_tilde_21[idx_j, idx_j]*tmp2
        tmp4 = dxidxhat_tilde_12[idx_j, idx_j]*tmp1 + dxidxhat_tilde_22[idx_j, idx_j]*tmp2
        tmp5 = A1tilde[idx_j, idx_j]*tmp3
        tmp6 = A2tilde[idx_j, idx_j]*tmp4
        tmp7 = tmp5 + tmp6
        tmp8 = tau_tilde[idx_j, idx_j]*tmp7
        tmp9 = H_tilde[idx_j, idx_j]*tmp8
        # weigthing space
        tmp10 = A1tilde[idx_j, idx_j].'*tmp9
        tmp11 = A2tilde[idx_j, idx_j].'*tmp9
        tmp12 = (dxidx_tilde_11[idx_i, idx_i]*D_tilde_xi[idx_j, idx_i] + dxidx_tilde_21[idx_i, idx_i]*D_tilde_eta[idx_j, idx_i])
        tmp13 = (dxidx_tilde_12[idx_i, idx_i]*D_tilde_xi[idx_j, idx_i] + dxidx_tilde_22[idx_i, idx_i]*D_tilde_eta[idx_j, idx_i])
        tmp14 = tmp12.'*tmp10
        tmp15 = tmp13.'*tmp11
        tmp16 = tmp14 + tmp15
         
        
        for k=1:4
          res_test[idx_i[k]] += tmp16[k]
        end
#=
        println("\n  qx = \n", tmp3) 
        println("\n  qy = \n", tmp4)
        println("\n  trial space x term = \n", tmp5)
        println("\n  trial space y term = \n", tmp6)
        println("\n  trial space term = \n", tmp7)
        println("\n  tau = \n", tau_tilde[idx_j, idx_j])
        println("\n  after multiplication by tau, reduction vector = \n", tmp8)
        println("\n  after multiplication by w, reduction vector = \n", tmp9)
        println("\n  after multiplication by Ax, reduction vector = \n", tmp10)
        println("\n  after multiplication by Ay, reduction vector = \n", tmp11)
        println("\n  dxidx_11 = \n", dxidx_tilde_11[idx_i, idx_i])
        println("\n  dxidx_21 = \n", dxidx_tilde_21[idx_i, idx_i])
        println("\n  dxidx_12 = \n", dxidx_tilde_12[idx_i, idx_i])
        println("\n  dxidx_22 = \n", dxidx_tilde_22[idx_i, idx_i])
        println("\n  D_xi = \n", D_tilde_xi[idx_j, idx_i])
        println("\n  D_eta = \n", D_tilde_eta[idx_j, idx_i])
        println("\n  x factor = \n", tmp12)
        println("\n  y factor = \n", tmp13)
        println("\n  after multiplication by Dx, reduction vector = \n", tmp14)
        println("\n  after multiplication by Dy, reduction vector = \n", tmp15)
        println("\n  final results = ", res_test)
=#
      end
    end

    @fact res_test => roughly(gls_test, atol=1e-14)
#    println("diff = ", gls_test - res_test)

    # another check of the block matrix approach
    tmp1 = D_tilde_xi*q
    tmp2 = D_tilde_eta*q
    tmp3 = dxidxhat_tilde_11*tmp1 + dxidxhat_tilde_21*tmp2
    tmp4 = dxidxhat_tilde_12*tmp1 + dxidxhat_tilde_22*tmp2
    tmp5 = A1tilde*tmp3
    tmp6 = A2tilde*tmp4
    tmp7 = tmp5 + tmp6
    tmp8 = tau_tilde*tmp7
    tmp9 = H_tilde*tmp8
    # weigthing space
    tmp10 = A1tilde.'*tmp9
    tmp11 = A2tilde.'*tmp9
    tmp12 = (dxidx_tilde_11*D_tilde_xi + dxidx_tilde_21*D_tilde_eta).'*tmp10
    tmp13 = (dxidx_tilde_12*D_tilde_xi + dxidx_tilde_22*D_tilde_eta).'*tmp11
    res_test2 = tmp12 + tmp13

    @fact res_test2 => roughly(gls_test, atol=1e-14)
#    println("block matrix diff = ", gls_test - res_test2)
#=
    # check the weighting term

    tmp8a = A1tilde.'
    tmp9a = A2tilde.'
    tmp10a = (dxidx_tilde_11*D_tilde_xi + dxidx_tilde_21*D_tilde_eta).'*tmp8a
    tmp11a = (dxidx_tilde_12*D_tilde_xi + dxidx_tilde_22*D_tilde_eta).'*tmp9a
    weighting_term2 = tmp10a + tmp11a
    println("weighting term 1-2 diff = ", weighting_term.' - weighting_term2)

    # another weighting term formulation
    tmp8a = A1tilde*(dxidx_tilde_11*D_tilde_xi + dxidx_tilde_21*D_tilde_eta)
    tmp9a = A2tilde*(dxidx_tilde_12*D_tilde_xi + dxidx_tilde_22*D_tilde_eta)
    weighting_term3 = (tmp8a + tmp9a).'

    println("weighting term 2-3 diff = ", weighting_term2 - weighting_term3)

    println("weighting term 1-3 diff = ", weighting_term.' - weighting_term3)

=#



    return gls_test, gls_code
  end  # end function test_gls

  
  println("----- Testing GLS2 on steady channel -----")
  # test on the steady channel case
  include("input_vals_channel.jl")
  arg_dict["solve"] = false
  f = open("input_vals_channel_gls.jl", "w")
  println(f, arg_dict)
  close(f)

  resize!(ARGS, 1)
  ARGS[1] = "input_vals_channel_gls.jl"
  include(STARTUP_PATH)
  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
  test_gls(mesh, sbp, eqn ,opts)

  println("----- Testing GLS2 on isentropic vortex -----")
  # test on isentropic vortex
  include("input_vals_vortex3.jl")
  arg_dict["solve"] = false
  f = open("input_vals_vortex3_gls.jl", "w")
  println(f, arg_dict)
  close(f)

  resize!(ARGS, 1)
  ARGS[1] = "input_vals_vortex3_gls.jl"
  include(STARTUP_PATH)
  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
  test_gls(mesh, sbp, eqn ,opts)

  println("----- Performing GLS2 finite difference checks -----")
  jac_fd = zeros(12, 12)
  eps_fd = 1e-6
  res0 = copy(reshape(eqn.res[:, :, 1], 12))  # use res from previous run
  println("doing finite differences")
  for j=1:3
    for i=1:4
      pos = (j-1)*4 + i
      println("pos = ", pos)
      eqn.q[i, j, 1] += eps_fd
      test_gls(mesh, sbp, eqn, opts)
      res_ij = copy(reshape(eqn.res[:, :, 1], 12))
      jac_fd[:, pos] = (res_ij - res0)/eps_fd
      eqn.q[i, j, 1] -= eps_fd  # undo perturbatino
    end
  end

  # now do complex step

  println("doing complex step")
  println("typeof(q_vec) = ", typeof(q_vec))
  arg_dict["run_type"] = 5
  f = open("input_vals_vortex3c_gls.jl", "w")
  println(f, arg_dict)
  close(f)
  ARGS[1] = "input_vals_vortex3c_gls.jl"
  include(STARTUP_PATH)

  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
  println("typeof(q_vec) = ", typeof(q_vec))

  eps_c = 1e-20
  jac_c = zeros(12, 12)
  println("typeof(eqn) = ", typeof(eqn))
  println("eqn.q[:, :, 1] = ", eqn.q[:, :, 1])
  for j=1:3
    for i=1:4
      pos = (j-1)*4 + i
      println("pos = ", pos)
      eqn.q[i, j, 1] += complex(0, eps_c)
      test_gls(mesh, sbp, eqn, opts)
      res_ij = copy(reshape(eqn.res[:, :, 1], 12))
      jac_c[:, pos] = imag(res_ij)/eps_c
      eqn.q[i, j, 1] -= complex(0, eps_c)  # undo perturbatino
    end
  end

  for j=1:12
    @fact jac_c[:, j] => roughly(jac_fd[:, j], atol=1e-5)
  end

  for j=1:12
    println("column $j")
    println("jac_c = \n", jac_c[:, j])
    println("jac_fd = \n", jac_fd[:, j])
    println("jac diff = \n", jac_c[:,j] - jac_fd[:, j])
  end    



end  # end facts do block
