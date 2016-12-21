# test the entropy stable interface calculation

"""
  A way of calculating the entropy stable face integral, used for comparison
  and debugging.
"""
function calcESFaceIntegralTest{Tdim, Tsol, Tres, Tmsh}(
                                params::AbstractParamType{Tdim},
                                sbpface::AbstractFace,
                                iface::Interface,
                                qL::AbstractMatrix{Tsol},
                                qR::AbstractMatrix{Tsol}, 
                                aux_vars::AbstractMatrix{Tres},
                                dxidx_face::Abstract3DArray{Tmsh},
                                functor::FluxType,
                                resL::AbstractMatrix{Tres}, 
                                resR::AbstractMatrix{Tres})

#  println("----- entered calcEss test -----")

  numDofPerNode, numNodesPerElement = size(qL)
  numFaceNodes = length(sbpface.wface)

  Rprime = params.Rprime

  F_tmp = params.flux_vals1

  workA = params.A
  workB = sview(params.B, :, :, 1)
  workC = sview(params.B, :, :, 2)
  nrm = params.nrm

  perm_nu = sview(sbpface.perm, :, iface.faceR)
  perm_gamma = sview(sbpface.perm, :, iface.faceL)

  # inverse (transpose) permutation vectors
  iperm_gamma = params.iperm 
  inversePerm(perm_gamma, iperm_gamma)

  # get the face normals
  facenormal = sview(sbpface.normal, :, iface.faceL)
  for dim=1:Tdim
    fill!(nrm, 0.0)
    nrm[dim] = 1

    # Nx, wface times Rprime
    for i=1:numFaceNodes
      nrm_i = zero(Tmsh)
      for d=1:Tdim
        nrm_i += facenormal[d]*dxidx_face[d, dim, i]
      end
      fac = sbpface.wface[i]*nrm_i
      # multiply by Rprime into A
      for j=1:numNodesPerElement
        # should nbrperm be after perm_nu?
        workA[i, j] = fac*Rprime[sbpface.nbrperm[i], j]
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
function psi_vec(params, q_vals)
  s = EulerEquationMod.calcEntropy(params, q_vals)
  rho = q_vals[1]
  u = q_vals[2]/rho
  v = q_vals[3]/rho
  psi = [rho*u, rho*v]
  return psi
end

"""
  Calculate the entropy flux psi for an element.  
"""
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


"""
  Calculate the boundary operator in the specified direction times
  the vector of entropy flux (psi) vals
"""
function getEPsi{Tsol}(iface, sbpface, params,  qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, dir::Integer)
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
function contractLHS{Tsol, Tres}(params, qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, resL::AbstractMatrix{Tres}, resR::AbstractMatrix{Tres})

#  println("----- entered contractLHS -----")

  val1 = contractLHS(params, qL, resL)
  val2 = contractLHS(params, qR, resR)
  return val1, val2
end

"""
  Contract a residual array for an element with the entropy variables
"""
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

  wL_vec = reshape(wL, length(wL))

  return dot(wL_vec, resL_vec)
end

function getEface(iface, sbpface, dxidx_face, dir::Integer)
  # this is inconsistent with reduceEface, for reasons I don't understand
  error("getEface does not work correctly")

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

"""
  Compute E for the current face (specified by dir) times psiL and psiR
"""
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


  


"""
  Test the entropy stable volume terms and face integrals against the
  analytical entropy flux
"""
function runESSTest(mesh, sbp, eqn, opts; test_boundaryintegrate=false)

  functor = EulerEquationMod.IRFlux()
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
    dxidx_face = mesh.dxidx_face[:, :, :, i]
    Tres = eltype(eqn.res)
    resL_code = sview(eqn.res, :, :, elL)
    resR_code = sview(eqn.res, :, :, elR)
    resL_test = sview(res_test, :, :, elL)
    resR_test = sview(res_test, :, :, elR)
    resL_test2 = sview(res_test2, :, :, elL)
    resR_test2 = sview(res_test2, :, :, elR)


    EulerEquationMod.calcESFaceIntegral(eqn.params, mesh.sbpface, iface, qL, qR, aux_vars, dxidx_face, functor, resL_code, resR_code)

    calcESFaceIntegralTest(eqn.params, mesh.sbpface, iface, qL, qR, aux_vars, dxidx_face, functor, resL_test2, resR_test2)

    for i=1:size(resL_code, 1)
      for j=1:size(resR_code, 2)

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
    dxidx_face = mesh.dxidx_face[:, :, :, i]
    Tres = eltype(eqn.res)
    numDofPerNode = size(eqn.res, 1)
    numNodesPerElement = size(eqn.res, 2)
    resL = zeros(Tres, numDofPerNode, numNodesPerElement)
    resR = zeros(resL)

    # calculate the integral of entropy flux from the residual
    E_expensive = calcESFaceIntegralTest(eqn.params, mesh.sbpface, iface, qL, qR, aux_vars, dxidx_face, functor, resL, resR)
    lhsL, lhsR = contractLHS(eqn.params, qL, qR, resL, resR)

    nrm = zeros(mesh.dim)

    # calculate the integral of entropy flux from q
    rhs = 0.0
    for dir=1:mesh.dim
      fill!(nrm, 0.0)
      nrm[dir] = 1.0
      psiL, psiR = getPsi(eqn.params, qL, qR, nrm)
      rhs += reduceEface(iface, mesh.sbpface, dxidx_face, dir, psiL, psiR)
    end


    @fact -(lhsL + lhsR) --> roughly(rhs, atol=1e-12)
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
  nrm = zeros(2)
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
      psi_nrm = dxidx[d, 1, 1]*psi[:, 1] + dxidx[d, 2, 1]*psi[:, 2]
      rhs_reduced = sum(E_d*psi_nrm)

      @fact lhs_reduced --> roughly(-rhs_reduced, atol=1e-12)

      res_total[:, :] += lhs_tmp1
    end
    # check that this calculation is doing the same thing as the actual code
    @fact res_total --> roughly(-eqn.res[:, :, i], atol=1e-12)



  end  # end loop i

  println("finished checking lemma 2")





end

"""
  Running the entropy stable tests.
"""
function test_ESS()
  facts("----- testing ESS -----") do
    ARGS[1] = "input_vals_channel_dg_large.jl"
    include(STARTUP_PATH)
    # evaluate the residual to confirm it is zero
    EulerEquationMod.evalEuler(mesh, sbp, eqn, opts)

    println("checking channel flow")
    runESSTest(mesh, sbp, eqn, opts, test_boundaryintegrate=false)

    println("\nchecking ICExp")
    ICFunc = EulerEquationMod.ICDict["ICExp"]
    ICFunc(mesh, sbp, eqn, opts, eqn.q_vec)
    scale!(eqn.q_vec, 0.01)
    disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        eqn.aux_vars[1, j, i] = EulerEquationMod.calcPressure(eqn.params, eqn.q[:, j, i])
      end
    end
    runESSTest(mesh, sbp, eqn, opts)

  end  # end facts block

  return nothing
end  # end functions

#test_ESS()
add_func1!(EulerTests, test_ESS, [TAG_FLUX])
