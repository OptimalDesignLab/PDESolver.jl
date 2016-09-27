# test the entropy stable interface calculation

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
include( joinpath(Pkg.dir("PDESolver"), "src/solver/euler/complexify.jl"))
include( joinpath(Pkg.dir("PDESolver"), "src/input/make_input.jl"))
global const STARTUP_PATH = joinpath(Pkg.dir("PDESolver"), "src/solver/euler/startup.jl")
# insert a command line argument
resize!(ARGS, 1)


function calcESS{Tdim, Tsol, Tres, Tmsh}(params::AbstractParamType{Tdim}, sbpface, iface, 
                 qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol},
                 aux_vars::AbstractMatrix{Tres},
                 dxidx_face::Abstract3DArray{Tmsh}, functor::FluxType, 
                 resL::AbstractMatrix{Tres},
                 resR::AbstractMatrix{Tres})

  println("----- entered calcEss test -----")
  # form permutation matrices
  numDofPerNode, numNodesPerElement = size(qL)

  P_nu = zeros(numNodesPerElement, numNodesPerElement)
  P_kappa = zeros(P_nu)
  permvec_nu = sbpface.perm[:, iface.faceR]
  permvec_kappa = sbpface.perm[:, iface.faceL]
#  println("permvec_kappa = \n", permvec_kappa)
#  println("sbpface.perm = ", sbpface.perm)
#  println("sbpface.interp = ", sbpface.interp)
  for i=1:length(permvec_nu)
    P_nu[i, permvec_nu[i]] = 1
    P_kappa[i, permvec_kappa[i]] = 1
  end

#  println("P_nu = \n", P_nu)
#  println("P_kappa = \n", P_kappa)

  numFaceNodes = length(sbpface.wface)
#  println("numFaceNodes = ", numFaceNodes)
  Rprime = zeros(numFaceNodes, numNodesPerElement)
  Rprime[1:size(sbpface.interp,2), 1:size(sbpface.interp,1)] = sbpface.interp.'
#  println("Rprime = \n", Rprime)

  B = diagm(sbpface.wface)
#  println("sbpface.wface = ", sbpface.wface)
  Nx = zeros(Tmsh, numFaceNodes, numFaceNodes)
  facenormal = sbpface.normal[:, iface.faceL]

  println("dxidx_face = \n", dxidx_face)
  println("facenormal = ", facenormal)

  full_size = numDofPerNode*numNodesPerElement
  F = zeros(Tres, full_size, full_size)

  nrm = zeros(3)
  for dim=1:Tdim  # DEBUGGING loop 1:Tdim
    println("dim = ", dim)
    fill!(nrm, 0)
    nrm[dim] = 1
    # calculate F
    F_tmp = zeros(Tres, numDofPerNode)
    for i=1:numNodesPerElement
      qL_i = qL[:, i]
      aux_vars_i = aux_vars[:, i]
      idx_i = (numDofPerNode*(i-1) + 1):(numDofPerNode*i)
      for j=1:numNodesPerElement
        qR_j = qR[:, j]
        functor(params, qL_i, qR_j, aux_vars_i, nrm, F_tmp)
        idx_j = (numDofPerNode*(j-1) + 1):(numDofPerNode*j)

        # store into F
        F[idx_i, idx_j] = diagm(F_tmp)

      end
    end
#    println("F = \n", F)

    # calculate Nx
    for i=1:numFaceNodes
      nx = zero(Tmsh)
      for d=1:Tdim
        nx += facenormal[d]*dxidx_face[d, dim, i]
      end
      Nx[i, i] = nx
    end

#    println("Nx = ", Nx)
    println("Nx times wface = \n", Nx*B)

    # calculate the full operator

    E_expensiveL = P_kappa.'*Rprime.'*B*Nx*Rprime*P_nu
    println("E_expensiveL = \n", E_expensiveL)
    
    E_expensiveL_full = kron(E_expensiveL, eye(numDofPerNode, numDofPerNode))
    E_expensiveR_full = kron(E_expensiveL.', eye(numDofPerNode, numDofPerNode))

    termL = zeros(E_expensiveL_full)
    termR = zeros(E_expensiveR_full)
    @assert size(termL) == size(F)
    @assert size(E_expensiveL_full) == size(E_expensiveR_full)
    for i =1:length(E_expensiveL_full)
      termL[i] = E_expensiveL_full[i]*F[i]
      termR[i] = E_expensiveR_full[i]*F[i]
    end
    termL_reduced = sum(termL, 2)
    termR_reduced = sum(termR, 2)

    tmp = B*Nx*Rprime
    println("Nx times wface times Rprime = \n", tmp)

    tmp2 = Rprime.'*tmp
    println("Rprime transpoes times A = \n", tmp2)

    tmp3 = tmp2*P_nu
    println("postmulitply P_nu = \n", tmp3)

    tmp4 = P_kappa.'*tmp3
    println("P_gamma transpose = \n", P_kappa.')
    println("premultiply P_gamma = \n", tmp4)


#=
    # compute some terms for comparison
    I4 = eye(4,4)
    P_nu_full = kron(P_nu, I4)
    println("P_nu_full = \n", P_nu_full)
    permF = P_nu_full.*F
    println("permF = \n", permF)
    reduced_permF = sum(permF, 2)
    println("reduced_permF = \n", reduced_permF)

    tmp2 = kron(Rprime, I4)*reduced_permF
    println("after multiplication by Rprime, tmp2 = \n", tmp2)

    tmp2 = kron(B*Nx, I4)*tmp2
    println("after multiplication by B, Nx, tmp2 = \n", tmp2)

    tmp1 = kron(Rprime.', I4)*tmp2
    println("after multiplication by Rprime transpose, tmp1 = \n", tmp1)

    println("P_kappa transpose = \n", P_kappa.')
    println("P_kappa_big = \n", kron(P_kappa.', I4))
    tmp3 = kron(P_kappa.', I4)*tmp1
    println("after multiplication by Ptranspose, tmp3 = \n", tmp3)

    println("termL_reduced = \n", termL_reduced)

    println("---checking elementR ---")
    permF = kron(P_kappa, I4).*F
    println("permF = \n", permF)
    reduced_permF = sum(permF, 2)
    println("reduced_permF = \n", reduced_permF)

    tmp2 = kron(Rprime, I4)*reduced_permF
    println("after multiplication by Rprime, tmp2 = \n", tmp2)

    tmp2 = kron(B*Nx, I4)*tmp2
    println("after multiplication by B, Nx, tmp2 = \n", tmp2)

    tmp1 = kron(Rprime.', I4)*tmp2
    println("after multiplication by Rprime transpose, tmp1 = \n", tmp1)

    println("P_kappa transpose = \n", P_nu.')
    println("P_kappa_big = \n", kron(P_nu.', I4))
    tmp3 = kron(P_nu.', I4)*tmp1
    println("after multiplication by Ptranspose, tmp3 = \n", tmp3)

    println("termR_reduced = \n", termR_reduced)
=#


    # update res
    for i=1:numNodesPerElement
      for j=1:numDofPerNode
        idx = numDofPerNode*(i-1) + j
        resL[j, i] += termL_reduced[idx]
        resR[j, i] -= termR_reduced[idx]
      end
    end

  end  # end loop over Tdim

  return resL, resR
end

function runESSTest(mesh, sbp, eqn, opts)

  functor = EulerEquationMod.IRFlux()
  println("numInterfaces = ", mesh.numInterfaces)
  for i=1:mesh.numInterfaces  # DEBUGGING
    println("interface ", i)
    iface = mesh.interfaces[i]
    println("iface = \n", iface)
    qL = eqn.q[:, :, iface.elementL]
    qR = eqn.q[:, :, iface.elementR]
    println("qL = \n", qL)
    println("qR = \n", qR)
    aux_vars = eqn.aux_vars[:,:, iface.elementL]
    dxidx_face = mesh.dxidx_face[:, :, :, i]
    Tres = eltype(eqn.res)
    resL_code = zeros(Tres, size(qL))
    resR_code = zeros(Tres, size(qL))
    resL_test = zeros(Tres, size(qL))
    resR_test = zeros(Tres, size(qL))
    EulerEquationMod.calcESFaceIntegral(eqn.params, mesh.sbpface, iface, qL, qR, aux_vars, dxidx_face, functor, resL_code, resR_code)

    println("resL_code = ", resL_code)
    println("resR_code = ", resR_code)
    calcESS(eqn.params, mesh.sbpface, iface, qL, qR, aux_vars, dxidx_face, functor, resL_test, resR_test)
    println("resL_test = ", resL_test)
    println("resR_test = ", resR_test)

    for i=1:size(resL_code, 1)
      for j=1:size(resR_code, 2)
        @fact resL_code[i, j] --> roughly(resL_test[i, j], atol=1e-12)
        
        @fact resR_code[i, j] --> roughly(resR_test[i, j], atol=1e-12)
      end
    end

  end

end

facts("----- testing ESS -----") do
  ARGS[1] = "input_vals_channel_dg_large.jl"
  include(STARTUP_PATH)
  runESSTest(mesh, sbp, eqn, opts)

  ICFunc = EulerEquationMod.ICDict["ICExp"]
  ICFunc(mesh, sbp, eqn, opts, eqn.q_vec)
  disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
  runESSTest(mesh, sbp, eqn, opts)

end


