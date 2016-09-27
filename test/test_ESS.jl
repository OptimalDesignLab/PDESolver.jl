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

#  println("----- entered calcEss test -----")
  # form permutation matrices
  numDofPerNode, numNodesPerElement = size(qL)

  permvec_nu = sbpface.perm[:, iface.faceR]
  permvec_kappa = sbpface.perm[:, iface.faceL]
  P_nu = permMatrix(permvec_nu)
  P_kappa = permMatrix(permvec_kappa)


  numFaceNodes = length(sbpface.wface)
#  println("numFaceNodes = ", numFaceNodes)
  Rprime = zeros(numFaceNodes, numNodesPerElement)
  Rprime[1:size(sbpface.interp,2), 1:size(sbpface.interp,1)] = sbpface.interp.'
#  println("Rprime = \n", Rprime)

  B = diagm(sbpface.wface)
#  println("sbpface.wface = ", sbpface.wface)
  Nx = zeros(Tmsh, numFaceNodes, numFaceNodes)
  facenormal = sbpface.normal[:, iface.faceL]

  full_size = numDofPerNode*numNodesPerElement
  F = zeros(Tres, full_size, full_size)

  nrm = zeros(3)
  for dim=1:Tdim
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

    # calculate the full operator

    E_expensiveL = P_kappa.'*Rprime.'*B*Nx*Rprime*P_nu
    
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
  for i=1:mesh.numInterfaces
    iface = mesh.interfaces[i]
    qL = eqn.q[:, :, iface.elementL]
    qR = eqn.q[:, :, iface.elementR]
    aux_vars = eqn.aux_vars[:,:, iface.elementL]
    dxidx_face = mesh.dxidx_face[:, :, :, i]
    Tres = eltype(eqn.res)
    resL_code = zeros(Tres, size(qL))
    resR_code = zeros(Tres, size(qL))
    resL_test = zeros(Tres, size(qL))
    resR_test = zeros(Tres, size(qL))
    EulerEquationMod.calcESFaceIntegral(eqn.params, mesh.sbpface, iface, qL, qR, aux_vars, dxidx_face, functor, resL_code, resR_code)

    calcESS(eqn.params, mesh.sbpface, iface, qL, qR, aux_vars, dxidx_face, functor, resL_test, resR_test)

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


