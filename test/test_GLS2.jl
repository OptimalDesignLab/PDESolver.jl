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
    
    # do GLS for first element only, compare the slow, explicit version 
    # with the fast index notation one
    q = zeros(12)
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
      A1 = zeros(4,4)
      A2 = zeros(4,4)
      A0inv = zeros(4,4)
      tau = zeros(12, 12)
      range_idx = (1:4, 5:8, 9:12)
      for i=1:3
        println("i = ", i)
        idx_i = range_idx[i]
        tau_i = zeros(4,4)
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
        A_mat = zeros(4,4,2)
        A_mat[:, :, 1] = A1
        A_mat[:, :, 2] = A2


        tau_old = zeros(4, 4)
        EulerEquationMod.getTau(params, A_mat, dxidx_i, tau_old)

        println("tau_test = ", new_tau2)
        println("tau_code = ", tau_old)
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
      Atilde = zeros(12, 12)
      A = zeros(4,4)
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
      dxidx_tilde = zeros(12, 12)

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
    fill!(eqn.res, 0.0)
    EulerEquationMod.applyGLS2(mesh, sbp, eqn, opts)
    gls_code = reshape(eqn.res[:, :, 1], 12)

    println("gls_test = ", gls_test)
    println("gls_code = ", gls_code)

    @fact gls_code => roughly(gls_test, atol=1e-14)

    println("q = ", q)
    println("dxidx_hat = ", dxidx_hat)
    println("Dxi = ", D_tilde_xi)
    println("Deta = ", D_tilde_eta)
    println("A1tilde = ", A1tilde)
    println("A2tilde = ", A2tilde)
    println("qxi = ", D_tilde_xi*q)
    println("qeta = ", D_tilde_eta*q)
    println("dxidx = ", dxidx)

    tmp1 = dxidxhat_tilde_11*D_tilde_xi*q + dxidxhat_tilde_21*D_tilde_eta*q
    tmp2 = dxidxhat_tilde_12*D_tilde_xi*q + dxidxhat_tilde_22*D_tilde_eta*q
    tmp3 = A1tilde*tmp1
    tmp4 = A2tilde*tmp2
    tmp5 = tmp3 + tmp4
    tmp6 = tau_tilde*tmp5
    tmp7 = H_tilde*tmp6
    tmp8 = A1tilde.'*tmp7
    tmp9 = A2tilde.'*tmp7
    tmp10 = (dxidx_tilde_11*D_tilde_xi + dxidx_tilde_21*D_tilde_eta).'*tmp8
    tmp11 = (dxidx_tilde_12*D_tilde_xi + dxidx_tilde_22*D_tilde_eta).'*tmp9
    tmp12 = tmp10 + tmp11

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
    println("\nPrinting intermediate quantities:\n")
    println("qx = ", tmp1) 
    println("qy = ", tmp2)
    println("trial space x term = ", tmp3)
    println("trial space y term = ", tmp4)
    println("after multiplication by tau, reduction vector = \n", tmp6)
    println("after multiplication by w, reduction vector = \n", tmp7)
    println("after multiplication by Ax, reduction vector = \n", tmp8)
    println("after multiplication by Ay, reduction vector = \n", tmp9)
    println("after multiplication by Dx, reduction vector = \n", tmp10)
    println("after multiplication by Dy, reduction vector = \n", tmp11)
    println("final results = ", tmp12)

    @fact tmp12 => roughly(gls_test, atol=1e-14)
    println("diff = ", gls_test - tmp12)
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



end  # end facts do block
