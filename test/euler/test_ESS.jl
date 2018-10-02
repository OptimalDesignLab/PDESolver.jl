# test the entropy stable interface calculation

"""
  A way of calculating the entropy stable face integral, used for comparison
  and debugging.
"""
function calcECFaceIntegralTest(params::AbstractParamType{Tdim},
        sbpface::AbstractFace,
        iface::Interface,
        qL::AbstractMatrix{Tsol},
        qR::AbstractMatrix{Tsol}, 
        aux_vars::AbstractMatrix{Tres},
        nrm_face::AbstractArray{Tmsh, 2},
        functor::FluxType,
        resL::AbstractMatrix{Tres}, 
        resR::AbstractMatrix{Tres}) where {Tdim, Tsol, Tres, Tmsh}

#  println("----- entered calcEss test -----")

  numDofPerNode, numNodesPerElement = size(qL)
  numFaceNodes = length(sbpface.wface)

  Rprime = params.Rprime

  F_tmp = zeros(Tres, numDofPerNode)
  workA = zeros(Tres, size(Rprime))
  workB = zeros(Tres, numNodesPerElement, numNodesPerElement)
  workC = zeros(Tres, numNodesPerElement, numNodesPerElement)
  nrm = zeros(Tmsh, size(nrm_face, 1))
  iperm = zeros(Int, size(sbpface.perm, 1))

  perm_nu = sview(sbpface.perm, :, iface.faceR)
  perm_gamma = sview(sbpface.perm, :, iface.faceL)

  # inverse (transpose) permutation vectors
  iperm_gamma = copy(iperm)
  inversePerm(perm_gamma, iperm_gamma)

  # get the face normals
  facenormal = sview(sbpface.normal, :, iface.faceL)
  for dim=1:Tdim
    fill!(nrm, 0.0)
    nrm[dim] = 1

    # Nx, wface times Rprime
    for i=1:numFaceNodes
      fac = sbpface.wface[i]*nrm_face[dim, i]
      # multiply by Rprime into A
      for j=1:numNodesPerElement
        # should nbrperm be after perm_nu?
        workA[i, j] = fac*Rprime[sbpface.nbrperm[i, iface.orient], j]
      end
    end

    # multiply by Rprime.'*A = B
    smallmatTmat!(Rprime, workA, workB)
    
    # post multiply by perm_nu
    applyPermColumn(perm_nu, workB, workC)

    # pre multiply by perm_gamma
    applyPermRow(iperm_gamma, workC, workB)

    # compute (B hadamard F)1
    for i=1:numNodesPerElement
      q_i = sview(qL, :, i)
      aux_vars_i = sview(aux_vars, :, i)
      for j=1:numNodesPerElement
        q_j = sview(qR, :, j)
        functor(params, q_i, q_j, aux_vars_i, nrm, F_tmp)
        B_val = workB[i, j]
        
        for p=1:numDofPerNode
          # because this term is on the rhs, the signs are reversed
          resL[p, i] -= B_val*F_tmp[p]
          resR[p, j] += B_val*F_tmp[p]
        end
      end
    end

  end  # end loop Tdim

  return copy(workB)
end
 
"""
  Calculate the entropy flux psi at a node
"""
function psi_vec(params::AbstractParamType{Tdim}, q_vals) where Tdim
  psi_vec = zeros(eltype(q_vals), Tdim)
  for i=1:Tdim
    psi_vec[i] = q_vals[i+1]
  end

  return psi_vec
end

"""
  Calculate the entropy flux psi for an element.  
"""
function getPsi(params, qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, nrm::AbstractVector) where Tsol

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


"""
  Calculate the boundary operator in the specified direction times
  the vector of entropy flux (psi) vals
"""
function getEPsi(iface, sbpface, params,  qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, dir::Integer) where Tsol
  # this computes the regular E * psi
  # not what is needed for lemma 3
#  println("----- entered getEPsi -----")

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

  # interpolate to face
  psifaceL = zeros(Tsol, numFaceNodes, 1)
  psifaceR = zeros(psifaceL)

  bndryL = Boundary(1, iface.faceL)
  bndryR = Boundary(1, iface.faceR)

  boundaryinterpolate!(sbpface, [bndryL], psiL, psifaceL)
  boundaryinterpolate!(sbpface, [bndryR], psiR, psifaceR)

  # integrate and interpolate back to volume nodes
  resL = zeros(Tsol, numNodesPerElement, 1)
  resR = zeros(resL)

  boundaryintegrate!(sbpface, [bndryL], psifaceL, resL)
  boundaryintegrate!(sbpface, [bndryR], psifaceR, resR)

  return sum(resL), sum(resR)

end

"""
  Contract resL and resR with entropy variables.
"""
function contractLHS(params, qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, resL::AbstractMatrix{Tres}, resR::AbstractMatrix{Tres}) where {Tsol, Tres}

#  println("----- entered contractLHS -----")

  val1 = contractLHS(params, qL, resL)
  val2 = contractLHS(params, qR, resR)
  return val1, val2
end

"""
  Contract a residual array for an element with the entropy variables
"""
function contractLHS(params, qL::AbstractMatrix{Tsol}, resL::AbstractMatrix{Tres}) where {Tsol, Tres}

  numDofPerNode, numNodesPerElement = size(qL)
  resL_vec = reshape(resL, length(resL))

  wL = zeros(Tsol, numDofPerNode, numNodesPerElement)

  aux_vars = Array{Tsol}(1)
#  w_tmp = zeros(Tsol, numDofPerNode)
  for i=1:numNodesPerElement
    qL_i = qL[:, i]

    w_tmp = sview(wL, :, i)
    EulerEquationMod.convertToEntropy_(params, qL_i, w_tmp)
    scale!(w_tmp, 1./params.gamma_1)

  end

  wL_vec = reshape(wL, length(wL))

  return dot(wL_vec, resL_vec)
end

function getEface(iface, sbpface, nrm_face, dir::Integer)
  # this is inconsistent with reduceEface, for reasons I don't understand
  error("getEface does not work correctly")

  dim = size(nrm_face, 1)
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
        nrm_k = nrm_face[dir, k]
        EL[i, p_j] += sbpface.interp[i, k]*sbpface.interp[j,k]*sbpface.wface[k]*nrm_k
        kR = sbpface.nbrperm[k, iface.orient]
        # need to consider nbrperm for p_i?
        ER[i, p_j] += -sbpface.interp[i, kR]*sbpface.interp[j, kR]*sbpface.wface[k]*nrm_k
      end
    end
  end

  return EL, ER
end

import EulerEquationMod.reduceEface
#=
"""
  Compute E for the current face (specified by dir) times psiL and psiR
"""
function reduceEface(iface, sbpface, nrm_face::AbstractMatrix, dir::Integer, psiL, psiR)

  dim = size(nrm_face, 1)
  RHS1 = 0.0
  RHS2 = 0.0

  for i=1:sbpface.stencilsize
    for j=1:sbpface.stencilsize
      val_acc = 0.0
      psi_val = psiL[sbpface.perm[j, iface.faceL]]
      for k=1:sbpface.numnodes
        nrm_k = nrm_face[dir, k]
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
=#

  
function entropyDissipativeRef(
              params::AbstractParamType{Tdim},
              sbpface::AbstractFace, iface::Interface,
              qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol},
              aux_vars::AbstractMatrix{Tres}, nrm_face::AbstractArray{Tmsh, 2},
              resL::AbstractMatrix{Tres}, resR::AbstractMatrix{Tres}) where {Tdim, Tsol, Tres, Tmsh}

#  println("----- entered entropyDissipativeRef -----")
  numDofPerNode = size(qL, 1)
  # convert to entropy variables
  wL = zeros(Tsol, numDofPerNode, sbpface.stencilsize)
  wR = zeros(wL)

  for i=1:sbpface.stencilsize
    qL_i = sview(qL, :, i)
    qR_i = sview(qR, :, i)
    wL_i = sview(wL, :, i)
    wR_i = sview(wR, :, i)
    EulerEquationMod.convertToIR(params, qL_i, wL_i)
    EulerEquationMod.convertToIR(params, qR_i, wR_i)
  end

  wLT = wL.'
  wRT = wR.'
  wLTP = zeros(wLT)  # permuted
  wRTP = zeros(wRT)

  applyPermRow(sbpface.perm[:, iface.faceL], wLT, wLTP)
  applyPermRow(sbpface.perm[:, iface.faceR], wRT, wRTP)

  wLTP_interp = sbpface.interp.'*wLTP
  wRTP_interp = sbpface.interp.'*wRTP

  # apply neigbor permutation to wR
  wRTNP_interp  = zeros(wRTP_interp)
  applyPermRow(vec(sbpface.nbrperm[:, iface.orient]), wRTP_interp, wRTNP_interp)

  # transpose back to numDofPerNode x numNodes
  wL_face = wLTP_interp.'
  wR_face = wRTNP_interp.'

  # get the middle term
  middle_term = zeros(Tres, numDofPerNode, numDofPerNode, sbpface.numnodes)
  A0 = zeros(Tsol, numDofPerNode, numDofPerNode)
  qL_i = zeros(Tsol, numDofPerNode)
  qR_i = zeros(Tsol, numDofPerNode)
  for i=1:sbpface.numnodes
    wL_i = wL_face[:, i]
    wR_i = wR_face[:, i]
    nrm = nrm_face[:, i]

    EulerEquationMod.convertToConservativeFromIR_(params, wL_i, qL_i)
    EulerEquationMod.convertToConservativeFromIR_(params, wR_i, qR_i)

    q_avg = 0.5*(qL_i + qR_i)
    EulerEquationMod.getIRA0(params, q_avg, A0)
    lambda_max = EulerEquationMod.getLambdaMaxSimple(params, qL_i, qR_i, nrm)

    middle_term[:, :, i] = lambda_max*sbpface.wface[i]*A0
  end

  # get the left term
  # this way is slower than the way the right term is computed, but more instructive
  PL = zeros(Int, sbpface.stencilsize, sbpface.stencilsize)
  PR = zeros(PL)
  permMatrix!(sbpface.perm[:, iface.faceL], PL)
  permMatrix!(sbpface.perm[:, iface.faceR], PR)
  Pnbr = permMatrix(vec(sbpface.nbrperm[:, iface.orient]))
  Rtranspose = sbpface.interp
  
  lhs_L = PL.'*Rtranspose
  lhs_R = PR.'*Rtranspose*Pnbr.'

  # construct w^T * lhs_L
  R = Rtranspose.'

  contractL3 = zeros(Tres, numDofPerNode, sbpface.numnodes)
  contractR3 = zeros(contractL3)
  for i=1:sbpface.numnodes
    for j=1:sbpface.stencilsize
      contractL3[:, i] += lhs_L[j, i]*wL[:, j]
      contractR3[:, i] += lhs_R[j, i]*wR[:, j]
    end
  end
  
  middle_sum = zeros(Tres, numDofPerNode, sbpface.numnodes)
  for i=1:sbpface.numnodes
    middle_sum[:, i] = middle_term[:, :, i]*(wL_face[:, i] - wR_face[:, i])
  end

  # update res
  delta_s = zero(Tres)
  for i=1:sbpface.numnodes
    delta_lhs = contractL3[:, i] - contractR3[:, i]
    delta_s += dot(delta_lhs, middle_sum[:, i])
  end

  for i=1:sbpface.numnodes
    for j=1:sbpface.stencilsize
      resL[:, j] -= lhs_L[j, i]*middle_term[:, :, i]*(wL_face[:, i] - wR_face[:, i])
      resR[:, j] += lhs_R[j, i]*middle_term[:, :, i]*(wL_face[:, i] - wR_face[:, i])
    end
  end

  return nothing
end


"""
  Test the entropy stable volume terms and face integrals against the
  analytical entropy flux
"""
function runECTest(mesh, sbp, eqn, opts, func_name="ECFaceIntegral"; test_ref=false)
# test_ref: whether or not to compare against the reference implementations above
# func_name: name of functor from FaceelementDict

  functor = EulerEquationMod.IRFlux()
  ec_integral = EulerEquationMod.FaceElementDict[func_name](mesh, eqn)
  eqn.flux_func = functor
  fill!(eqn.res, 0.0)
  res_test = copy(eqn.res)
  res_test2 = copy(eqn.res)
  for i=1:mesh.numInterfaces
    iface = mesh.interfaces[i]
    elL = iface.elementL
    elR = iface.elementR
    qL = eqn.q[:, :, iface.elementL]
    qR = eqn.q[:, :, iface.elementR]
    aux_vars = eqn.aux_vars[:,:, iface.elementL]
    nrm_face = mesh.nrm_face[:, :, i]
    Tres = eltype(eqn.res)
    resL_code = sview(eqn.res, :, :, elL)
    resR_code = sview(eqn.res, :, :, elR)
    resL_test = sview(res_test, :, :, elL)
    resR_test = sview(res_test, :, :, elR)
    resL_test2 = sview(res_test2, :, :, elL)
    resR_test2 = sview(res_test2, :, :, elR)


    ec_integral(eqn.params, mesh.sbpface, iface, qL, qR, aux_vars, nrm_face, functor, resL_code, resR_code)

    if test_ref
      calcECFaceIntegralTest(eqn.params, mesh.sbpface, iface, qL, qR, aux_vars, nrm_face, functor, resL_test2, resR_test2)

      @test isapprox( norm(vec(resL_code - resL_test2)), 0.0) atol=1e-12*length(resL_code)
      @test isapprox( norm(vec(resR_code - resR_test2)), 0.0) atol=1e-12*length(resR_code)

      #=
      for j=1:size(resL_code, 1)
        for k=1:size(resR_code, 2)

          @test isapprox( resL_code[j, k], resL_test2[j, k]) atol=1e-12
          @test isapprox( resR_code[j, k], resR_test2[j, k]) atol=1e-12
        end
      end
      =#

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
        @test isapprox( eqn.res[k, j, i], res_test[k, j, i]) atol=1e-12
      end
    end
  end
=#

  # verify lemma 3
  total_potentialflux = 0.0
  println("\nchecking lemma 3")
  for i=1:mesh.numInterfaces
    iface = mesh.interfaces[i]
    elL = iface.elementL
    elR = iface.elementR
    qL = eqn.q[:, :, iface.elementL]
    qR = eqn.q[:, :, iface.elementR]

    aux_vars = eqn.aux_vars[:,:, iface.elementL]
    nrm_face = mesh.nrm_face[:, :, i]
    Tres = eltype(eqn.res)
    numDofPerNode = size(eqn.res, 1)
    numNodesPerElement = size(eqn.res, 2)
    resL = zeros(Tres, numDofPerNode, numNodesPerElement)
    resR = zeros(resL)

    # calculate the integral of entropy flux from the residual
    ec_integral(eqn.params, mesh.sbpface, iface, qL, qR, aux_vars, nrm_face, functor, resL, resR)
    lhsL, lhsR = contractLHS(eqn.params, qL, qR, resL, resR)

    nrm = zeros(mesh.dim)

    # calculate the integral of entropy flux from q
    rhs = reduceEface(eqn.params, iface, mesh.sbpface, nrm_face, qL, qR)
#=
    rhs = 0.0
    for dir=1:mesh.dim
      fill!(nrm, 0.0)
      nrm[dir] = 1.0
      psiL, psiR = getPsi(eqn.params, qL, qR, nrm)
      rhs += reduceEface(iface, mesh.sbpface, nrm_face, dir, psiL, psiR)
    end
=#

    @test isapprox( -(lhsL + lhsR), rhs) atol=1e-12
    total_potentialflux -= lhsL + lhsR


  end

  println("finished checking lemma 3")

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

  F = zeros(mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numNodesPerElement, mesh.dim)
  res_el = zeros(mesh.numDofPerNode, mesh.numNodesPerElement)
  psi = zeros(mesh.numNodesPerElement, mesh.dim)
  nrm = zeros(mesh.dim)
  for i=1:mesh.numEl
    q_i = eqn.q[:, :, i]
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

    # calculate (S .* F)1
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


      # contract residual with entropy variables
      lhs_reduced = contractLHS(eqn.params, q_i, lhs_tmp1)

      # calculate expected potential flux from q variables
      rhs_reduced = 0.0
      psi_nrm = zeros(eltype(psi), size(psi, 1))
      for d2=1:mesh.dim
        psi_nrm[:] += dxidx[d, d2, 1]*psi[:, d2]# + dxidx[d, 2, 1]*psi[:, 2]
      end
      rhs_reduced = sum(E_d*psi_nrm)

      @test isapprox( lhs_reduced, -rhs_reduced) atol=1e-12

      res_total[:, :] += lhs_tmp1
    end
    # check that this calculation is doing the same thing as the actual code
    @test isapprox( res_total, -eqn.res[:, :, i]) atol=1e-12



  end  # end loop i

  println("finished checking lemma 2")


end

function runESTest(mesh, sbp, eqn, opts, penalty_name::String; test_ref=false, zero_penalty=false)
# run entropy stability tests
# test_ref: whether or not to compare against the reference implementations above
# zero_penalty: test whether the entropy stability penalty should be zero
# penalty name: name of penalty functor to call

  println("checking entropy dissipation")

  # the penalty functions don't need a flux, so pick an arbitrary one
  flux_func = EulerEquationMod.FluxDict["StandardFlux"]  
  # check calcLFEntropyPenaltyIntegral
  for i=1:mesh.numInterfaces
    iface = mesh.interfaces[i]
    elL = iface.elementL
    elR = iface.elementR
    qL = eqn.q[:, :, iface.elementL]
    qR = eqn.q[:, :, iface.elementR]
    aux_vars = eqn.aux_vars[:,:, iface.elementL]
    nrm_face = mesh.nrm_face[:, :, i]

    resL = zeros(mesh.numDofPerNode, mesh.numNodesPerElement)
    resR = zeros(resL)
    resL2 = zeros(resL)
    resR2 = zeros(resR)

#    resL3 = zeros(resL)
#    resR3 = zeros(resR)


    penalty_func = EulerEquationMod.FaceElementDict[penalty_name](mesh, eqn)
    lf_penalty_func = EulerEquationMod.FaceElementDict["ELFPenaltyFaceIntegral"](mesh, eqn)

    penalty_func(eqn.params, mesh.sbpface, iface, qL, qR, aux_vars, nrm_face, flux_func,  resL, resR)

    if test_ref
      entropyDissipativeRef(eqn.params, mesh.sbpface, iface, qL, qR, aux_vars, nrm_face, resL2, resR2)

      @test isapprox( norm(vec(resL - resL2)), 0.0) atol=1e-12
      @test isapprox( norm(vec(resR - resR2)), 0.0) atol=1e-12
 
      #=
      for j=1:size(resL, 1)
        for k=1:size(resL, 2)
          @test isapprox( abs(resL[j, k] - resL2[j, k]), 0.0) atol=1e-12
          @test isapprox( abs(resR[j, k] - resR2[j, k]), 0.0) atol=1e-12
        end
      end
      =#

    end

    # verify conservation
    resL_sum = sum(resL, 2)
    resR_sum = sum(resR, 2)

    @test isapprox( resL_sum, -resR_sum) atol=1e-13

    if zero_penalty
      @test isapprox( norm(vec(resL)), 0.0) atol=1e-13
      @test isapprox( norm(vec(resR)), 0.0) atol=1e-13
 
      #=
      for j=1:mesh.numNodesPerElement
        for p=1:mesh.numDofPerNode
          @test isapprox( resL[p, j], 0.0) atol=1e-13
          @test isapprox( resR[p, j], 0.0) atol=1e-13
        end
      end
      =#
    end
#=
    # verify equality with Lax-Friedrich
    for i=1:size(resL, 1)
      for j=1:size(resL, 2)
        @test isapprox( abs(resL[i, j] - resL3[i, j]), 0.0) atol=1e-12
        @test isapprox( abs(resR[i, j] - resR3[i, j]), 0.0) atol=1e-12
      end
    end
=#


    # contract with entropy variables, verify result is negative
    resL_reduced = zeros(mesh.numNodesPerElement)
    resR_reduced = zeros(mesh.numNodesPerElement)
    wL = zeros(mesh.numDofPerNode)
    wR = zeros(mesh.numDofPerNode)
    resL_reduced2 = zeros(mesh.numNodesPerElement)
    resR_reduced2 = zeros(mesh.numNodesPerElement)

    for j=1:mesh.numNodesPerElement
      EulerEquationMod.convertToEntropy_(eqn.params, qL[:, j], wL)
      EulerEquationMod.convertToEntropy_(eqn.params, qR[:, j], wR)
      scale!(wL, 1/eqn.params.gamma_1)
      scale!(wR, 1/eqn.params.gamma_1)

      resL_reduced[j] = dot(wL, resL[:, j])
      resR_reduced[j] = dot(wR, resR[:, j])
      
      resL_reduced2[j] = dot(wL, resL2[:, j])
      resR_reduced2[j] = dot(wR, resR2[:, j])

    end

#    println("resL_reduced = \n", resL_reduced)
#    println("resR_reduced = \n", resR_reduced)

#    println("sbp.w = ", sbp.w)
    delta_sL = sum(resL_reduced)
    delta_sR = sum(resR_reduced)
    delta_s = delta_sL + delta_sR
    @test  delta_s  < eps()

    if !(delta_s < eps())
      println("entropy growth at interface ", i)
      println("iface = ", iface)
      println("qL = \n", qL)
      println("qR = \n", qR)
    end
    if test_ref
      delta_sL2 = sum(resL_reduced2)
      delta_sR2 = sum(resR_reduced2)
      delta_s2 = delta_sL2 + delta_sR2
      @test  delta_s2  < eps()
    end

  end  # end loop over interfaces

  println("finished checking entropy dissipation")

end

"""
  Test the entropy-stable boundary condition
"""
function test_ESSBC(mesh, sbp, eqn, opts)

  func = EulerEquationMod.BCDict["noPenetrationESBC"](mesh, eqn)

  EulerEquationMod.interpolateBoundary(mesh, sbp, eqn, opts, eqn.q, eqn.q_bndry, eqn.aux_vars_bndry)

  bndry_flux = zeros(eltype(eqn.res), mesh.numDofPerNode)
  w = zeros(eltype(eqn.q), mesh.numDofPerNode)

  for i=1:mesh.numBoundaryFaces
    for j=1:mesh.numNodesPerFace
      q = ro_sview(eqn.q, :, j, i)
      aux_vars = ro_sview(eqn.aux_vars_bndry, :, j, i)
      coords = ro_sview(mesh.coords_bndry, :, j, i)
      nrm_xy = mesh.nrm_bndry[:, j, i]
      scale!(nrm_xy, 1/norm(nrm_xy))  # normalize

      func(eqn.params, q, aux_vars, coords, nrm_xy, bndry_flux)
      EulerEquationMod.convertToIR(eqn.params, q, w)

      # test: psi - w^T f(q) < 0
      psi_n = q[2]*nrm_xy[1] + q[3]*nrm_xy[2]
      @test  psi_n - dot(w, bndry_flux)  <= 1e-12
    end
  end

  return nothing
end


"""
  Used for debugging LW functions, not used for regular testing
"""
function testLW(mesh, sbp, eqn::EulerEquationMod.EulerData{Tsol, Tres, Tdim}, opts) where {Tsol, Tres, Tdim}
  # computes the Lax-Wendroff term using the maximum eigenvalue for all 
  # eigenvalue, turning it into Lax-Friedrich
  # compare agains LF
  params = eqn.params
  sbpface = mesh.sbpface

  Y = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode)
  A0 = copy(Y)
  S2 = zeros(Tsol, mesh.numDofPerNode)  # S is defined s.t. (YS)*(YS).' = A0
  Lambda = zeros(Tsol, mesh.numDofPerNode)
  nrm = sview(mesh.nrm_face, :, 1, 1)

  iface = mesh.interfaces[1]
  elL = iface.elementL
  elR = iface.elementR
  qL = eqn.q[:, 1, elL]
  qR = eqn.q[:, 1, elR]

  q_avg = 0.5*(qL + qR)

  # delta_w doesn't matter, so make an arbitrary vector
  delta_w = collect(Tsol, 1:length(qL))

  fill!(tmp3, 0.0)
  fill!(nrm, 0.0)

  # compute LW term
  lambda_net = 0.0

  nrm2 = nrm./norm(nrm)
  A1 = zeros(4,4)  # flux jacobian computed via sum of x, y directions
  A2 = zeros(4,4)  # flux jacobian in nrm direction computed directly
  for dim=1:Tdim
    # get the eigensystem in the current direction
    if dim == 1
      EulerEquationMod.calcEvecsx(params, q_avg, Y)
      EulerEquationMod.calcEvalsx(params, q_avg, Lambda)
      EulerEquationMod.calcEScalingx(params, q_avg, S2)
    elseif dim == 2
      EulerEquationMod.calcEvecsy(params, q_avg, Y)
      EulerEquationMod.calcEvalsy(params, q_avg, Lambda)
      EulerEquationMod.calcEScalingy(params, q_avg, S2)
    elseif dim == 3
      EulerEquationMod.calcEvecsz(params, q_avg, Y)
      EulerEquationMod.calcEvalsz(params, q_avg, Lambda)
      EulerEquationMod.calcEScalingz(params, q_avg, S2)
    end

    A1 += nrm[dim]*Y*diagm(Lambda)*inv(Y)
    # DEBUGGING: turn this into Lax-Friedrich
    lambda_max = maximum(absvalue(Lambda))
    lambda_net += absvalue(lambda_max*nrm[dim])
    fill!(Lambda, lambda_max)


    # compute the Lax-Wendroff term, returned in tmp2
    EulerEquationMod.applyEntropyLWUpdate(Y, Lambda, S2, delta_w, nrm[dim], tmp1, tmp2)
    # accumulate result
    for j=1:length(tmp3)
      tmp3[j] += tmp2[j]
    end
  end


  # compute LF term

  # this is not exactly the same eigenvalue calculation as lambda(q_avg)
#  lambda_max = EulerEquationMod.getLambdaMax(params, qL, qR, nrm)
  Un = (q_avg[2]*nrm[1] + q_avg[3]*nrm[2])/q_avg[1]
  p_avg = EulerEquationMod.calcPressure(params, q_avg)
  dA = sqrt(nrm[1]*nrm[1] + nrm[2]*nrm[2])
  a_avg = dA*sqrt(params.gamma*p_avg/q_avg[1])  # speed of sound
  lambda_max = absvalue(Un) + absvalue(a_avg)
  EulerEquationMod.getIRA0(params, q_avg, A0)

  # lambda_max * A0 * delta w
  smallmatvec!(A0, delta_w, tmp1)
  scale!(tmp1, lambda_max)

  # compute scaling parameter alpha
  alpha = lambda_max/lambda_net
  scale!(tmp3, alpha)

  println("LW = \n", tmp3)
  println("LF = \n", tmp1)
  println("diff = \n", tmp3 - tmp1)

  println("lambda_net = ", lambda_net)
  println("lambda_max = ", lambda_max)
  println("alpha = ", alpha)
  println("nrm = ", nrm)

  q2 = convert(Array{Complex128, 1}, q_avg)
  h = 1e-20
  pert = Complex128(0, h)
  F = zeros(Complex128, length(q_avg))
  aux_vars = zeros(Complex128, 1)
  for i=1:length(q_avg)
    q2[i] += pert
    EulerEquationMod.calcEulerFlux(params, q2, aux_vars, nrm, F)
    A2[:, i] = imag(F)/h
    q2[i] -= pert
  end

  println("A1 = \n", A1)
  println("A2 = \n", A2)
  println("ratio = \n", A1./A2)
  println("diff = \n", A1 - A2)


  return nothing
end

function applyPoly(mesh, sbp, eqn::EulerData{Tsol, Tres}, opts, p) where {Tsol, Tres}
# set the solution to be a polynomial of degree p of the entropy variables

  v_vals = zeros(Tsol, mesh.numDofPerNode)
  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      x = mesh.coords[1, j, i]
      y = mesh.coords[2, j, i]

      v_vals[1] = 4*(x^p) + 4*(y^p) - 100
      v_vals[2] = 2*(x^p) + 2*(y^p) + 1
      v_vals[3] = 3*(x^p) + 3*(y^p) + 2
      v_vals[4] = 4*(x^p) + 4*(y^p) - 10
      if mesh.dim == 3
        v_vals[5] = 4*(x^p) + 4*(y^p) - 10
      end

#      println("v_vals = \n", v_vals)
      q_i = sview(eqn.q, :, j, i)
      EulerEquationMod.convertToConservativeFromIR_(eqn.params, v_vals, q_i)
      eqn.aux_vars[1, j, i] = EulerEquationMod.calcPressure(eqn.params, q_i)
    end
  end
        

  return nothing
end

function factRes0(mesh, sbp, eqn, opts)

  @test isapprox( norm(vec(eqn.res)), 0.0) atol=1e-13
  #=
  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      for k=1:mesh.numDofPerNode
        @test isapprox( eqn.res[k, j, i], 0.0) atol=1e-13
      end
    end
  end
  =#

  return nothing
end


function test_ESS()
  @testset "----- testing ESS -----" begin
    fname = "input_vals_channel_dg_large.jl"
    mesh, sbp, eqn, opts = solvePDE(fname)
    # evaluate the residual to confirm it is zero
    EulerEquationMod.evalResidual(mesh, sbp, eqn, opts)
    penalty_lf = "ELFPenaltyFaceIntegral"
    penalty_lw2 = "ELW2PenaltyFaceIntegral"

    for dim =2:3
      println("\n\n\ndim = ", dim)
      println("testing dimensions ", dim)
      if dim == 2
        meshname = "SRCMESHES/tri8l.smb"
      else
        meshname = "SRCMESHES/tet8cubep.smb"
      end

      for p=1:4
        println("\n\np = ", p)
        println("test p = ", p)
        if dim == 3 && p > 2  # 3D Omega operators for p > 2 don't exists yet
          continue
        end
        println("testing p = ", p)
        opts["dimensions"] = dim
        opts["smb_name"] = meshname
        opts["operator_type"] = "SBPOmega"
        opts["order"] = p
        opts["numBC"] = 0
        fname = "input_vals_ESS_test2"
        make_input(opts, fname)
        fname = fname*".jl"
        mesh, sbp, eqn, opts = solvePDE(fname)
       

        runECTest(mesh, sbp, eqn, opts, test_ref=true)

        ICFunc = EulerEquationMod.ICDict["ICExp"]
        ICFunc(mesh, sbp, eqn, opts, eqn.q_vec)
        scale!(eqn.q_vec, 0.01)
        array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)
        for i=1:mesh.numEl
          for j=1:mesh.numNodesPerElement
            eqn.aux_vars[1, j, i] = EulerEquationMod.calcPressure(eqn.params, eqn.q[:, j, i])
          end
        end
        runECTest(mesh, sbp, eqn, opts, test_ref=true)
        println("testing LF dissipation")
        runESTest(mesh, sbp, eqn, opts, penalty_lf, test_ref=true)
        println("testing LW2 dissipation")
        runESTest(mesh, sbp, eqn, opts, penalty_lw2, test_ref=false)
        println("finished testing LW dissipation")
        if dim == 2
          println("testing ESBC")
          test_ESSBC(mesh, sbp, eqn, opts)
        end


    #    testLW(mesh, sbp, eqn, opts)
        # check polynomial
        println("testing polynomial")
        applyPoly(mesh, sbp, eqn, opts, p)
        runECTest(mesh, sbp, eqn, opts, test_ref=true)
        runESTest(mesh, sbp, eqn, opts, penalty_lf, test_ref=true, zero_penalty=true)
        if dim == 2
          println("testing ESBC")
          test_ESSBC(mesh, sbp, eqn, opts)
        end


        # check full calling sequence
        fill!(eqn.res, 0.0)

       penalty_functor = EulerEquationMod.FaceElementDict["ELFPenaltyFaceIntegral"](mesh, eqn)
        println("about to check full calling sequence")
        println("typeof(mesh) = ", typeof(mesh))
        println("typeof(penalty_functor) = ", typeof(penalty_functor))
        println("size(A0) = ", size(penalty_functor.kernel.A0))
        EulerEquationMod.getFaceElementIntegral(mesh, sbp, eqn, penalty_functor, eqn.flux_func, mesh.sbpface, mesh.interfaces)

        factRes0(mesh, sbp, eqn, opts)


        # check gamma operators
        opts["operator_type"] = "SBPGamma"
        fname = "input_vals_ESS_test2"
        make_input(opts, fname)
        fname = fname*".jl"
        mesh, sbp, eqn, opts = solvePDE(fname)
        
        ICFunc = EulerEquationMod.ICDict["ICExp"]
        ICFunc(mesh, sbp, eqn, opts, eqn.q_vec)
        scale!(eqn.q_vec, 0.01)
        array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)
        for i=1:mesh.numEl
          for j=1:mesh.numNodesPerElement
            eqn.aux_vars[1, j, i] = EulerEquationMod.calcPressure(eqn.params, eqn.q[:, j, i])
          end
        end
        runECTest(mesh, sbp, eqn, opts, test_ref=false)
        runESTest(mesh, sbp, eqn, opts, penalty_lf, test_ref=false)

        # check polynomial
        applyPoly(mesh, sbp, eqn, opts, p)
        runECTest(mesh, sbp, eqn, opts, test_ref=false)
        runESTest(mesh, sbp, eqn, opts, penalty_lf, test_ref=false, zero_penalty=true)

        # check full calling sequence
        fill!(eqn.res, 0.0)
        penalty_functor = EulerEquationMod.FaceElementDict["ELFPenaltyFaceIntegral"](mesh, eqn)
        EulerEquationMod.getFaceElementIntegral(mesh, sbp, eqn, penalty_functor, eqn.flux_func, mesh.sbpface, mesh.interfaces)

        factRes0(mesh, sbp, eqn, opts)

      end  # end loop over p
    end  # end loop over dim

    println("\nTesting diagonal E")
    # test a diagonal E operator
    dim = 3
    p = 2
    meshname = "SRCMESHES/tet8cubep.smb"
    println("testing p = ", p)
    opts["dimensions"] = dim
    opts["smb_name"] = meshname
    opts["operator_type"] = "SBPDiagonalE"
    opts["order"] = p
    opts["numBC"] = 0
    fname = "input_vals_ESS_test2"
    make_input(opts, fname)
    fname = fname*".jl"
    mesh, sbp, eqn, opts = solvePDE(fname)
   

    runECTest(mesh, sbp, eqn, opts, test_ref=false)

    # because there is no interpolation, add a small amount to the solution
    # on one side of the interface (so delta_w will be non-zero)
    for iface in mesh.interfaces
      for j in mesh.sbpface.perm
        for i=1:mesh.numDofPerNode
          eqn.q[i, j, iface.elementL] += 0.001
        end
      end
    end
    #TODO: test entropy dissipation
    println("testing LF dissipation")
    runESTest(mesh, sbp, eqn, opts, penalty_lf, test_ref=false)
    println("testing LW2 dissipation")
    runESTest(mesh, sbp, eqn, opts, penalty_lw2, test_ref=false)
    println("finished testing LW dissipation")



  end  # end facts block

  return nothing
end  # end function

#test_ESS()
add_func1!(EulerTests, test_ESS, [TAG_FLUX, TAG_ESS, TAG_SHORTTEST])
