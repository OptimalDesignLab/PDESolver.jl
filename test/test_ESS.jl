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

function interiorfacepenalty!{Tdim, Tsol, Tres, Tmsh}(
   params::AbstractParamType{Tdim}, sbpface::AbstractFace, iface::Interface,
   qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, 
   aux_vars::AbstractMatrix{Tres}, dxidx_face::Abstract3DArray{Tmsh},
   functor::FluxType, resL::AbstractMatrix{Tres}, resR::AbstractMatrix{Tres})

  fill!(resL, 0.0)
  fill!(resR, 0.0)
  resL2 = copy(resL)
  resR2 = copy(resR)
  println("----- entered interiorfacepenalty -----")
#  println("qL = \n", qL)
#  println("qR = \n", qR)

  Flux_tmp = params.flux_vals1
  Flux_tmp2 = params.flux_vals2
  nrm = params.nrm
  E_full = zeros(sbpface.stencilsize, sbpface.stencilsize, Tdim)
  F = zeros(4, sbpface.stencilsize, sbpface.stencilsize)
  for dim = 1:1  # DEBUGGING 1:TDIM
    fill!(nrm, 0.0)
    nrm[dim] = 1

    # loop over the nodes of "left" element that are in the stencil of interp
    for i = 1:sbpface.stencilsize
      p_i = sbpface.perm[i, iface.faceL]
      qi = qL[:, p_i]
      aux_vars_i = sview(aux_vars, :, p_i)  # !!!! why no aux_vars_j???

      # loop over the nodes of "right" element that are in the stencil of interp
      for j = 1:sbpface.stencilsize
        p_j = sbpface.perm[j, iface.faceR]
        qj = qR[:, p_j]
        # construct Eij
        Eij = zero(Tsol)  # should be Tres
        for k = 1:sbpface.numnodes
          # the computation of nrm_k could be moved outside i,j loops and saved
          # in an array of size [3, sbp.numnodes]
          nrm_k = zero(Tmsh)
          for d = 1:Tdim
            nrm_k += sbpface.normal[d, iface.faceL]*dxidx_face[d, dim, k]
          end
          kR = sbpface.nbrperm[k, iface.orient]
          Eij += sbpface.interp[i,k]*sbpface.interp[j,kR]*sbpface.wface[k]*nrm_k
        end  # end loop k
        E_full[p_i, p_j, dim] = Eij
        
        # compute flux and add contribution to left and right elements
        functor(params, qi, qj, aux_vars_i, nrm, Flux_tmp)
        F[:, p_i, p_j] = Flux_tmp[:]
        resL[:, p_i] -= Eij*Flux_tmp[:]
        resR[:, p_j] += Eij*Flux_tmp[:]

      end
    end

#    println("E_full = \n", E_full)
    for i=1:sbpface.stencilsize
      qi = sview(qL, :, i)
      aux_vars_i = sview(aux_vars, :, i)
      for j=1:sbpface.stencilsize
        qj = sview(qR, :, j)
        normdiff = norm(F[:, j, i] - F[:, i, j])

        functor(params, qi, qj, aux_vars_i, nrm, Flux_tmp)
        functor(params, qj, qi, aux_vars_i, nrm, Flux_tmp2)
        normdiff2 = norm(Flux_tmp - F[:, i, j])
        normdiff3 = norm(Flux_tmp - Flux_tmp2)
#        println("normdiff2 = ", normdiff2)
#        println("normdiff3 = ", normdiff3)
#        println("F ij, ji = \n", F[:, i, j], ", ", F[:, j, i])
        @assert (normdiff2 < 1e-12)
#        println("element L multiplying ", E_full[i, j], " with flux ", i, ", ", j, " and storing in res ", i)
#        println("element R multiplying ", E_full[j, i], " with flux ", i, ", ", j, " and storing to res ", i)
        resL2[:, i] -= E_full[i, j]*Flux_tmp[:]
        resR2[:, j] += E_full[j, i]*Flux_tmp[:]  # should be resR[:, j]

      end
    end

  end  # end loop Tdim


  return E_full
end



function calcESS{Tdim, Tsol, Tres, Tmsh}(params::AbstractParamType{Tdim}, sbpface, iface, 
                 qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol},
                 aux_vars::AbstractMatrix{Tres},
                 dxidx_face::Abstract3DArray{Tmsh}, functor::FluxType, 
                 resL::AbstractMatrix{Tres},
                 resR::AbstractMatrix{Tres})

  println("----- entered calcEss test -----")

  fill!(resL, 0.0)
  fill!(resR, 0.0)

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
  Rprime_nu = Rprime[sbpface.nbrperm[:, 1], :]
#  println("Rprime = \n", Rprime)

  B = diagm(sbpface.wface)
#  println("sbpface.wface = ", sbpface.wface)
  Nx = zeros(Tmsh, numFaceNodes, numFaceNodes)
  facenormal = sbpface.normal[:, iface.faceL]

  full_size = numDofPerNode*numNodesPerElement
  F = zeros(Tres, full_size, full_size)

  nrm = zeros(3)
  for dim=1:Tdim  #DEBUGGING 1:TDIM
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

    E_expensiveL = P_kappa.'*Rprime.'*B*Nx*Rprime_nu*P_nu
    E_expensiveL_full = kron(E_expensiveL, eye(numDofPerNode, numDofPerNode))
    E_expensiveR_full = kron(E_expensiveL.', eye(numDofPerNode, numDofPerNode))

    #  calculate terms for comparison
    interp_to_face = kron(Rprime_nu*P_nu, eye(4,4))

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
        # because this term is on the rhs, the signs are reversed
        resL[j, i] -= termL_reduced[idx]
        resR[j, i] += termR_reduced[idx]
      end
    end

  end  # end loop over Tdim

  return nothing
end

function psi_vec(params, q_vals)
  s = EulerEquationMod.calcEntropy(params, q_vals)
  rho = q_vals[1]
  u = q_vals[2]/rho
  v = q_vals[3]/rho
  psi = [rho*u, rho*v]
  return psi
end

function getPsi{Tsol}(params, qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, nrm::AbstractVector)

  numDofPerNode, numNodesPerElement = size(qL)
  psiL = zeros(Tsol, numNodesPerElement)
  psiR = zeros(psiL)
  dim = length(nrm)
  # calculate psi at volume nodes
  for i=1:numNodesPerElement
    qL_i = qL[:, i]
    qR_i = qR[:, i]

    psi_vecL = psi_vec(params, qL_i)
    psi_vecR = psi_vec(params, qR_i)
    for d=1:dim
      psiL[i] += nrm[d]*psi_vecL[d]
      psiR[i] += nrm[d]*psi_vecR[d]
    end
  end

  return psiL, psiR
end



function getEPsi{Tsol}(iface, sbpface, params,  qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, dir::Integer)
  # this computes the regular E * psi
  # not what is needed for lemma 3
  println("----- entered getEPsi -----")

  numDofPerNode, numNodesPerElement = size(qL)
  numFaceNodes = sbpface.numnodes
  
  psiL = zeros(Tsol, numNodesPerElement, 1)
  psiR = zeros(psiL)

  # calculate psi at volume nodes
  for i=1:numNodesPerElement
    qL_i = qL[:, i]
    qR_i = qR[:, i]

    psi = psi_vec(params, qL_i)
    psiL[i] = psi[dir]

    psi = psi_vec(params, qR_i)
    psiR[i] = psi[dir]
  end

  println("psiL = \n", psiL)
  println("psiR = \n", psiR)

  # interpolate to face
  psifaceL = zeros(Tsol, numFaceNodes, 1)
  psifaceR = zeros(psifaceL)

  bndryL = Boundary(1, iface.faceL)
  bndryR = Boundary(1, iface.faceR)

  boundaryinterpolate!(sbpface, [bndryL], psiL, psifaceL)
  boundaryinterpolate!(sbpface, [bndryR], psiR, psifaceR)

  println("psifaceL = \n", psifaceL)
  println("psifaceR = \n", psifaceR)

  # integrate and interpolate back to volume nodes
  resL = zeros(Tsol, numNodesPerElement, 1)
  resR = zeros(resL)

  boundaryintegrate!(sbpface, [bndryL], psifaceL, resL)
  boundaryintegrate!(sbpface, [bndryR], psifaceR, resR)

  println("resL = \n", resL)
  println("resR = \n", resR)

  return sum(resL), sum(resR)

end

function contractLHS{Tsol, Tres}(params, qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, resL::AbstractMatrix{Tres}, resR::AbstractMatrix{Tres})

#  println("----- entered contractLHS -----")

  val1 = contractLHS(params, qL, resL)
  val2 = contractLHS(params, qR, resR)
  return val1, val2
end

function contractLHS{Tsol, Tres}(params, qL::AbstractMatrix{Tsol}, resL::AbstractMatrix{Tres})

  numDofPerNode, numNodesPerElement = size(qL)
  resL_vec = reshape(resL, length(resL))

  wL = zeros(Tsol, numDofPerNode, numNodesPerElement)

  aux_vars = Array(Tsol, 1)
#  w_tmp = zeros(Tsol, numDofPerNode)
  for i=1:numNodesPerElement
    qL_i = qL[:, i]

    w_tmp = sview(wL, :, i)
    EulerEquationMod.convertToEntropy_(eqn.params, qL_i, w_tmp)
    scale!(w_tmp, 1./params.gamma_1)

  end

  println("wL = \n", wL)

#  println("wL = \n", wL)
#  println("wR = \n", wR)
#  println("resL = \n", resL)
#  println("resR = \n", resR)
  wL_vec = reshape(wL, length(wL))

  return dot(wL_vec, resL_vec)
end

function getEface(iface, sbpface, dxidx_face, dir::Integer)
  dim = size(dxidx_face, 1)
  EL = zeros(sbpface.stencilsize, sbpface.stencilsize)
  ER = zeros(EL)
  for i=1:sbpface.stencilsize
#    p_i = sbpface.perm[i, iface.faceL]
    p_iL = sbpface.perm[i, iface.faceL]
    p_iR = sbpface.perm[i, iface.faceR]
    for j=1:sbpface.stencilsize
#      p_j = sbpface.perm[j, iface.faceL]
      p_jL = sbpface.perm[j, iface.faceL]
      p_jR = sbpface.perm[j, iface.faceR]
      for k=1:sbpface.numnodes
        nrm_k = 0.0
        for d=1:dim
          nrm_k += sbpface.normal[d, iface.faceL]*dxidx_face[d, dir, k]
        end
        EL[i, p_j] += sbpface.interp[i, k]*sbpface.interp[j,k]*sbpface.wface[k]*nrm_k
        kR = sbpface.nbrperm[k, iface.orient]
        # need to consider nbrperm for p_i?
        ER[i, p_j] += -sbpface.interp[i, kR]*sbpface.interp[j, kR]*sbpface.wface[k]*nrm_k
      end
    end
  end

  return EL, ER
end

function reduceEface(iface, sbpface, dxidx_face, dir::Integer, psiL, psiR)

  dim = size(dxidx_face, 1)
  RHS1 = 0.0
  RHS2 = 0.0

  for i=1:sbpface.stencilsize
    for j=1:sbpface.stencilsize
      val_acc = 0.0
      psi_val = psiL[sbpface.perm[j, iface.faceL]]
      for k=1:sbpface.numnodes
        nrm_k = 0.0
        for d=1:dim
          nrm_k += sbpface.normal[d, iface.faceL]*dxidx_face[d, dir, k]
        end
        val = sbpface.interp[i,k]*sbpface.interp[j,k]*sbpface.wface[k]*nrm_k
        val_acc += val
        RHS1 += val*psi_val
        kR = sbpface.nbrperm[k, iface.orient]
        RHS2 -= sbpface.interp[i, kR]*sbpface.interp[j, kR]*sbpface.wface[k]*nrm_k*psiR[sbpface.perm[j, iface.faceR]]
      end
    end
  end

  return RHS1 + RHS2
end


  



function runESSTest(mesh, sbp, eqn, opts; test_boundaryintegrate=false)

  println("coords2 = \n", mesh.vert_coords[:, :, 2])
  println("coords7 = \n", mesh.vert_coords[:, :, 7])
  functor = EulerEquationMod.IRFlux()
  eqn.flux_func = functor
  fill!(eqn.res, 0.0)
  res_test = copy(eqn.res)
  res_test2 = copy(eqn.res)
  for i=1:mesh.numInterfaces
    iface = mesh.interfaces[i]
    println("iface ", i)
    println(iface)
    elL = iface.elementL
    elR = iface.elementR
    qL = eqn.q[:, :, iface.elementL]
    qR = eqn.q[:, :, iface.elementR]
    aux_vars = eqn.aux_vars[:,:, iface.elementL]
    dxidx_face = mesh.dxidx_face[:, :, :, i]
    Tres = eltype(eqn.res)
    resL_code = sview(eqn.res, :, :, elL)
    resR_code = sview(eqn.res, :, :, elR)
    resL_test = sview(res_test, :, :, elL)
    resR_test = sview(res_test, :, :, elR)
    resL_test2 = sview(res_test2, :, :, elL)
    resR_test2 = sview(res_test2, :, :, elR)

#    resL_code = zeros(Tres, size(qL))
#    resR_code = zeros(Tres, size(qL))
#    resL_test = zeros(Tres, size(qL))
#    resR_test = zeros(Tres, size(qL))

    EulerEquationMod.calcESFaceIntegral(eqn.params, mesh.sbpface, iface, qL, qR, aux_vars, dxidx_face, functor, resL_code, resR_code)

    calcESS(eqn.params, mesh.sbpface, iface, qL, qR, aux_vars, dxidx_face, functor, resL_test, resR_test)

    interiorfacepenalty!(eqn.params, mesh.sbpface, iface, qL, qR, aux_vars, dxidx_face, functor, resL_test2, resR_test2)

    for i=1:size(resL_code, 1)
      for j=1:size(resR_code, 2)
#        @fact resL_code[i, j] --> roughly(resL_test[i, j], atol=1e-12)
#        @fact resR_code[i, j] --> roughly(resR_test[i, j], atol=1e-12)

        @fact resL_code[i, j] --> roughly(resL_test2[i, j], atol=1e-12)
        @fact resR_code[i, j] --> roughly(resR_test2[i, j], atol=1e-12)
      end
    end

  end  # end loop over interfaces
#=

  # test the overall calling sequence
  fill!(eqn.res, 0.0)
  opts["face_integral_type"] = 2
  EulerEquationMod.evalFaceIntegrals(mesh, sbp, eqn, opts)

  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      for k=1:mesh.numDofPerNode
        @fact eqn.res[k, j, i] --> roughly(res_test[k, j, i], atol=1e-12)
      end
    end
  end
=#
  if test_boundaryintegrate
    println("testing boundaryintegrate")
    EulerEquationMod.interpolateFace(mesh, sbp, eqn, opts, eqn.q, eqn.q_face)
    EulerEquationMod.calcFaceFlux(mesh, sbp, eqn, eqn.flux_func, mesh.interfaces, eqn.flux_face)
    res2 = copy(eqn.res)
    fill!(res2, 0.0)
    interiorfaceintegrate!(mesh.sbpface, mesh.interfaces, eqn.flux_face, res2, SummationByParts.Subtract())

    # compare res to res2
    for i=1:mesh.numEl
      println("element ", i)
      for j=1:mesh.numNodesPerElement
        for k=1:mesh.numDofPerNode
          # was eqn.res on rhs
          @fact res2[k, j, i] --> roughly(eqn.res[k, j, i], atol=1e-12)
        end
      end
    end

    println("finished testing boundaryintegrate")
  end

  # verify lemma 3
  println("\nchecking lemma 3")
  for i=1:mesh.numInterfaces
    iface = mesh.interfaces[i]
    println("iface ", i)
    println(iface)
    elL = iface.elementL
    elR = iface.elementR
    qL = eqn.q[:, :, iface.elementL]
    qR = eqn.q[:, :, iface.elementR]
    println("qL = \n", qL)
    println("qR = \n", qR)

    aux_vars = eqn.aux_vars[:,:, iface.elementL]
    dxidx_face = mesh.dxidx_face[:, :, :, i]
    Tres = eltype(eqn.res)
    numDofPerNode = size(eqn.res, 1)
    numNodesPerElement = size(eqn.res, 2)
    resL = zeros(Tres, numDofPerNode, numNodesPerElement)
    resR = zeros(resL)

    println("dxidx_face = \n", dxidx_face)
    E_expensive = interiorfacepenalty!(eqn.params, mesh.sbpface, iface, qL, qR, aux_vars, dxidx_face, functor, resL, resR)
    println("E_expensive = \n", E_expensive)
    println("resL = \n", resL)
    println("resR = \n", resR)
    rhsL = zero(Tres)
    rhsR = zero(Tres)

    nrm = zeros(mesh.dim)
#    facenormal = sbpface.normal[:, iface.faceL]
#    calcBCNormal(eqn.params, dxidx_face[:, :, 1], facenormal, nrm) 

    rhs = 0.0
    for dir=1:1
      fill!(nrm, 0.0)
      nrm[dir] = 1.0
      psiL, psiR = getPsi(eqn.params, qL, qR, nrm)
      println("PsiL = \n", psiL)
      println("psiR = \n", psiR)
      rhs += reduceEface(iface, mesh.sbpface, dxidx_face, dir, psiL, psiR)
#      EL, ER = getEface(iface, mesh.sbpface, dxidx_face, dir)
#      rhsL += sum(EL*psiL)
#      rhsR += sum(ER*psiR)  # negative transpose
    end

    println("rhs = ", rhs)
#    println("rhsL = ", rhsL, ", rhsR = ", rhsR)

    lhsL, lhsR = contractLHS(eqn.params, qL, qR, resL, resR)
    println("lhs = ", lhsL + lhsR)
    @fact -(lhsL + lhsR) --> roughly(rhs, atol=1e-12)


#=
    # check first line of lemma 3 proof
    E_expensiveL = copy(E_expensive)
    E_expensiveR = copy(E_expensive).'
    FL = zeros(numDofPerNode, numNodesPerElement, numNodesPerElement)
    FR = zeros(FL)
    F_tmp = zeros(numDofPerNode)

    # calculate F
    for j=1:numNodesPerElement
      nrm = [1.0; 0.0]
      q_j = qL[:, j]
      aux_vars_j = aux_vars[:, j]
      for k=1:numNodesPerElement
        q_k = qR[:, k]
        functor(eqn.params, q_j, q_k, aux_vars_j, nrm, F_tmp)
        FL[:, j, k] = F_tmp

        aux_vars_j[1] = EulerEquationMod.calcPressure(eqn.params, q_k)
        functor(eqn.params, q_k, q_j, aux_vars_j, nrm, F_tmp)
        FR[:, k, j] = F_tmp
      end
    end

    lhs = zeros(numDofPerNode, numNodesPerElement, numNodesPerElement)
    rhs = zeros(lhs)
    for j=1:numNodesPerElement
      for k=1:numNodesPerElement
        lhs[:, j, k] = E_expensiveL[j, k]*FL[:, j, k]
        rhs[:, j, k] = E_expensiveR[j, k]*FR[:, j, k]
      end
    end

    diff = zeros(lhs)
    for j=1:numNodesPerElement
      for k=1:numNodesPerElement
        diff[:, j, k] = lhs[:, j, k] - rhs[:, j, k]
      end
    end

    sumL = sum(lhs)
    sumR = sum(rhs)
    rsumL = sum(resL)
    rsumR = sum(resR)
#    println("sumL = ", sumL, ", sumR = ", sumR, ", diff = ", sumL-sumR)
#    println("rsumL = ", rsumL, ", rsumR = ", rsumR, ", diff = ", rsumL+rsumR)


    @fact sumL --> roughly(sumR, atol=1e-12)
    @fact rsumL --> roughly(-rsumR, atol=1e-12)
=#
  end

  println("finished checking lemma 3")
#=
  println("checking symmetry of E_expensive")

  # reverse interfaces
  ifaces2 = copy(mesh.interfaces)
  for i=1:mesh.numInterfaces
    iface = ifaces2[i]
    ifaces2[i] = Interface(iface.elementR, iface.elementL, iface.faceR, iface.faceL, iface.orient)
  end
  ifaces_orig = mesh.interfaces
  mesh.interfaces = ifaces2

  dxidx_faceR, other_vals = PdePumiInterface.interpolateMapping(mesh)
  mesh.interfaces = ifaces_orig

  for i=1:mesh.numInterfaces
    iface = mesh.interfaces[i]
    println("iface ", i)
    println(iface)
    elL = iface.elementL
    elR = iface.elementR
    qL = eqn.q[:, :, iface.elementL]
    qR = eqn.q[:, :, iface.elementR]
    aux_vars = eqn.aux_vars[:,:, iface.elementL]
    dxidx_face = mesh.dxidx_face[:, :, :, i]
    Tres = eltype(eqn.res)
    numDofPerNode = size(eqn.res, 1)
    numNodesPerElement = size(eqn.res, 2)
    resL = zeros(Tres, numDofPerNode, numNodesPerElement)
    resR = zeros(resL)

    println("calculating elementL")
    E_expensive = EulerEquationMod.calcESFaceIntegral(eqn.params, mesh.sbpface, iface, qL, qR, aux_vars, dxidx_face, functor, resL, resR)
    E_exp

    # for non uniform flow, would have to recalculate aux_vars here too
    dxidx_face = dxidx_faceR[:, :, :, i]

    iface2 = ifaces2[i]
    println("iface2 = ")
    println(iface2)

    println("calculating elementR")
    E_expensive2
    E_expensive2 = EulerEquationMod.calcESFaceIntegral(eqn.params, mesh.sbpface, iface2, qR, qL, aux_vars, dxidx_face, functor, resR, resL)

    E_expensive3 = -E_expensive2.'

    @fact E_expensive3 --> roughly(E_expensive, atol=1e-12)


    print("\n")

  end

  println("finished checking symmetry of E_expensive")

=#

  # check lemma2 - volume terms
  println("checking lemma 2")
  fill!(eqn.res, 0.0)
  EulerEquationMod.calcVolumeIntegralsSplitForm(mesh, sbp, eqn, opts, functor)
  # check the lemma, then check that calcVolumeIntegralsSplitForm is 
  # computing the same thing

  S = zeros(mesh.numNodesPerElement, mesh.numNodesPerElement, mesh.dim)
  E = zeros(S)
  for dim=1:mesh.dim
    S[:, :, dim] = 0.5*(sbp.Q[:, :, dim] - sbp.Q[:, :, dim].')
    E[:, :, dim] = (sbp.Q[:, :, dim] + sbp.Q[:, :, dim].')
  end

  println("E = \n", E)

  F = zeros(mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numNodesPerElement, mesh.dim)
  res_el = zeros(mesh.numDofPerNode, mesh.numNodesPerElement)
  psi = zeros(mesh.numNodesPerElement, mesh.dim)
  nrm = zeros(2)
  for i=1:mesh.numEl
    println("element ", i)
    q_i = eqn.q[:, :, i]
    println("q_i = \n", q_i)
    aux_vars_i = eqn.aux_vars[:, :, i]
    dxidx = mesh.dxidx[:, :, :, i]
    # calculate F
    for j=1:mesh.numNodesPerElement
      q_j = q_i[:, j]
      aux_vars_j = aux_vars_i[:, j]
      for k=1:mesh.numNodesPerElement
        q_k = q_i[:, k]

        for d=1:mesh.dim
          for p=1:mesh.dim
            nrm[p] = dxidx[d, p, j]
          end
#          fill!(nrm, 0.0)
#          nrm[d] = 1.0

          F_tmp = sview(F, :, j, k, d)
          functor(eqn.params, q_j, q_k, aux_vars_j, nrm, F_tmp)
        end  # end loop d

      end  # end loop k
    end  # end loop j

    # calculate psi vector
    for j=1:mesh.numNodesPerElement
      tmp = psi_vec(eqn.params, q_i[:, j])
      for d=1:mesh.dim
        psi[j, d] = tmp[d]
      end
    end

    res_total = zeros(mesh.numDofPerNode, mesh.numNodesPerElement)
    for d=1:mesh.dim
      S_d = S[:, :, d]
      E_d = E[:, :, d]

      lhs_tmp1 = zeros(mesh.numDofPerNode, mesh.numNodesPerElement) 
      for j=1:mesh.numNodesPerElement
        for k=1:mesh.numNodesPerElement
          lhs_tmp1[:, j] += 2*S_d[j, k]*F[:, j, k, d]
        end
      end

      println("lhs_tmp1 = ", lhs_tmp1)

      lhs_reduced = contractLHS(eqn.params, q_i, lhs_tmp1)
      println("Psi = \n", psi)

      # now do rhs
      rhs_reduced = 0.0
      psi_nrm = dxidx[d, 1, 1]*psi[:, 1] + dxidx[d, 2, 1]*psi[:, 2]
      rhs_reduced = sum(E_d*psi_nrm)
#      println("E*psi = ", E_d*psi[:, d])

      println("lhs_reduced = ", lhs_reduced, ", rhs_reduced = ", rhs_reduced)

      @fact lhs_reduced --> roughly(-rhs_reduced, atol=1e-12)

      res_total[:, :] += lhs_tmp1
    end
    # check that this calculation is doing the same thing as the actual code
    @fact res_total --> roughly(-eqn.res[:, :, i], atol=1e-12)



  end  # end loop i

  println("finished checking lemma 2")





end

facts("----- testing ESS -----") do
  ARGS[1] = "input_vals_channel_dg_large.jl"
  include(STARTUP_PATH)

  # evaluate the residual to confirm it is zero
  EulerEquationMod.evalEuler(mesh, sbp, eqn, opts)
  println("eqn.res = \n", eqn.res)

  runESSTest(mesh, sbp, eqn, opts, test_boundaryintegrate=false)

  println("\nchecking ICExp")
  ICFunc = EulerEquationMod.ICDict["ICExp"]
  ICFunc(mesh, sbp, eqn, opts, eqn.q_vec)
  disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
  q = eqn.q
  q[:, :, 2] = q[:, :, 3]
  q[:, :, 7] = q[:, :, 4]
  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      eqn.aux_vars[1, j, i] = EulerEquationMod.calcPressure(eqn.params, eqn.q[:, j, i])
    end
  end
  runESSTest(mesh, sbp, eqn, opts)

end


