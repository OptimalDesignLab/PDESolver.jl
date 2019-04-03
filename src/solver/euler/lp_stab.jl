# Local Projection (LP) stabilization

import SummationByParts: UnaryFunctor, Add, Subtract

#=
"""
  Computes the LPS matrix P = I - L*L^T*H
"""
function getLPSMatrix(sbp::AbstractOperator)

  #TODO: do this correctly
  P = zeros(sbp.numnodes, sbp.numnodes)
  for i=1:length(P)
    P[i] = 0.01*i
  end

  return P
end
=#
function getLPSMatrix(sbp::AbstractOperator)

  diss = zeros(sbp.numnodes, sbp.numnodes)

  x = calcnodes(sbp)
  V = SummationByParts.calcvandermondproriol(x.', sbp.degree)

  # Minv = diagm(1./diag(V.'*diagm(sbp.w)*V))
  # Minv is the idenity, since the V's are orthogonal in L2
  diss[:,:] = eye(sbp.numnodes) - (V*V.')*diagm(sbp.w)

  return diss
end

"""
  Main function for applying LPS stabilization

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
"""
function applyLPStab(mesh, sbp, eqn, opts)

  lps_data = eqn.params.lps_data
  entropy_vars = lps_data.entropy_vars

  _applyLPStab(mesh, sbp, eqn, opts, lps_data, entropy_vars)

  return nothing
end

@noinline function _applyLPStab(mesh, sbp, eqn, opts, data::LPSData{Tsol},
                       entropy_vars::AbstractVariables) where {Tsol}

  w_i = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerElement)
  for i=1:mesh.numEl
    q_i = ro_sview(eqn.q, :, :, i)
    dxidx = ro_sview(mesh.dxidx, :, :, :, i)
    res_i = sview(eqn.res, :, :, i)

    # convert to entropy variables
    for j=1:mesh.numNodesPerElement
      q_j = ro_sview(q_i, :, j)
      w_j = sview(w_i, :, j)
      convertToEntropy(entropy_vars, eqn.params, q_j, w_j)
    end

    applyLPSKernel(eqn.params, data, sbp, entropy_vars, q_i, w_i, dxidx,res_i, i)
  end

  return nothing
end


"""
  Compute res += op( P^T * H * A * P * u_in)  for each component of the system.

  **Inputs**

   * params: ParamType
   * data: [`LPSData`](@ref)
   * sbp: SBP operator
   * q: conservative variables for the element, `numDofPerNode` x
        `numNodesPerElement`
   * u_in: values to multiply against, same size as `q`
   * op: UnaryFunctor specifying whether to add or subtract from `res`.  Default
         subtract.

  **Inputs/Outputs**

   * res: array to update with result, same size as `q`

  **Aliasing**

  Note that `q` and `u_in` can alias
"""
function applyLPSKernel(params::ParamType{Tdim}, data::LPSData, sbp::AbstractOperator,
                        entropy_vars::AbstractVariables,
                        q::AbstractMatrix{Tsol}, u_in::AbstractMatrix,
                        dxidx::Abstract3DArray{Tmsh},
                        res::AbstractMatrix{Tres}, elnum, op::UnaryFunctor=Subtract()
                        ) where {Tsol, Tres, Tmsh, Tdim}
  # For scalar equations, the operator is applied -epsilon * P^T M P * u
  # For vector equations, P needs to be applied to all equations as once:
  # utmp^T = P*u^T
  # Instead, work with utmp = (P*u^T)^T = u*P^T

  numDofPerNode, numNodesPerElement = size(q)

  t1 = zeros(Tsol, numDofPerNode, numNodesPerElement)
  t2 = zeros(Tsol, numDofPerNode, numNodesPerElement)
  A0 = zeros(Tsol, numDofPerNode, numDofPerNode)
  dir = zeros(Tmsh, Tdim)

  # apply P
  smallmatmatT!(u_in, data.P, t1)

  if elnum == 1
    println("\nEntered applyLPSKernel")
    println("t1_dot = \n", imag(t1[1,1])./1e-20)
  end

  # apply mass matrix and scaling
  @simd for i=1:numNodesPerElement
    q_i = ro_sview(q, :, i)
    getA0(entropy_vars, params, q_i, A0)
    t1_i = ro_sview(t1, :, i)
    t2_i = sview(t2, :, i)

    # scaling by wave speed
    lambda = zero(Tres)
    for d1=1:Tdim
      for d2=1:Tdim
        dir[d2] = dxidx[d1, d2, i]
      end
      lambda += getLambdaMax(params, q_i, dir)
    end
    fac = lambda*sbp.w[i]/Tdim

    # multiply by A0 and fac at the same time
    smallmatvec_kernel!(A0, t1_i, t2_i, fac, 0)
  end

  if elnum == 1
    println("t2_dot = \n", imag(t2[1,1])./1e-20)
  end

  # apply P^T
  smallmatmat!(t2, data.P, t1)
  @simd for i=1:numNodesPerElement
    @simd for j=1:numDofPerNode
      res[j, i] += op(data.alpha*t1[j, i])
    end
  end

  if elnum == 1
    println("t3_dot = \n", imag(res[1,1])./1e-20)
  end

  return nothing
end


#------------------------------------------------------------------------------
# Differentiated version

"""
  Differentiated version of `applyLPStab`
"""
function applyLPStab_diff(mesh, sbp, eqn, opts, assem::AssembleElementData)

  lps_data = eqn.params.lps_data
  entropy_vars = lps_data.entropy_vars
  
  _applyLPStab_diff(mesh, sbp, eqn, opts, lps_data, entropy_vars, assem)

  return nothing
end


@noinline function _applyLPStab_diff(mesh::AbstractMesh, sbp::AbstractOperator,
                                 eqn::EulerData, opts,
                                 data::LPSData{Tsol, Tres},
                                 entropy_vars::AbstractVariables,
                                 assem::AssembleElementData) where {Tsol, Tres}

  res_jac = eqn.params.calc_volume_integrals_data.res_jac

  w_in = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerElement)
  w_jac = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode,
                      mesh.numNodesPerElement)
  for i=1:mesh.numEl
    q_i = ro_sview(eqn.q, :, :, i)
    dxidx = ro_sview(mesh.dxidx, :, :, :, i)
    # convert to entropy variables
    for j=1:mesh.numNodesPerElement
      q_j = ro_sview(q_i, :, j)
      w_j = sview(w_in, :, j)
      w_jac_j = sview(w_jac, :, :, j)
      convertToEntropy(entropy_vars, eqn.params, q_j, w_j)
      getA0inv(entropy_vars, eqn.params, q_j, w_jac_j)
    end

    fill!(res_jac, 0)
    applyLPSKernel_diff(eqn.params, data, sbp, entropy_vars, q_i, w_in, w_jac, dxidx,
                        res_jac, i)

    assembleElement(assem, mesh, i, res_jac)
  end

  return nothing
end



"""
  Differentiated version of `applyEntropyKernel`

  **Inputs**

   * params
   * data: [`LPSData`](@ref)
   * sbp
   * q
   * u_in_jac: 3D jacobian of `u_in` wrt `q`, `numDofPerNode` x `numDofPerNode`
               x `numNodesPerElement`
   * op

  **Inputs/Outputs**

   * res_jac: 4D jacobian to be updated with the result.
"""
function applyLPSKernel_diff(params::ParamType{Tdim}, data::LPSData,
                        sbp::AbstractOperator, entropy_vars::AbstractVariables,
                        q_el::AbstractMatrix{Tsol}, u_in::AbstractMatrix,
                        u_in_jac::Abstract3DArray,
                        dxidx::Abstract3DArray{Tmsh},
                        res_jac::Abstract4DArray{Tres}, elnum, op::UnaryFunctor=Subtract()
                        ) where {Tsol, Tres, Tmsh, Tdim}


  numDofPerNode, numNodesPerElement = size(q_el)
#=
  # compute derivative contribution from v
  @simd for p=1:numNodesPerElement
    @simd for q=1:numNodesPerElement
      # calculate the A[p, q]
      Apq = zero(Tsol)
      @simd for k=1:numNodesPerElement
        Apq += data.P[k, p]*sbp.w[k]*data.P[k, q]
      end
      @simd for j=1:numDofPerNode
        @simd for i=1:numDofPerNode
          res_jac[i, j, p, q] += op(data.alpha*Apq*u_in_jac[i, j, q])  # * ee[1, p]
        end
      end
    end
  end
=#
  t1 = zeros(Tsol, numDofPerNode, numNodesPerElement)
  t1_jac = zeros(Tsol, numDofPerNode, numDofPerNode, numNodesPerElement,
                       numNodesPerElement)
  t2 = zeros(Tsol, numDofPerNode, numNodesPerElement)
  t2_jac = zeros(Tsol, numDofPerNode, numDofPerNode, numNodesPerElement,
                       numNodesPerElement)
  t3_jac = zeros(Tsol, numDofPerNode, numDofPerNode, numNodesPerElement,
                       numNodesPerElement)
  A0 = zeros(Tsol, numDofPerNode, numDofPerNode)
  A0_diff = zeros(Tsol, numDofPerNode, numDofPerNode, numDofPerNode)
  q_dot = zeros(numDofPerNode, numDofPerNode)
  for i=1:numDofPerNode
    q_dot[i, i] = 1
  end

  lambdas = zeros(Tres, numNodesPerElement)
  lambdas_jac = zeros(Tres, numDofPerNode, numNodesPerElement)
  dir = zeros(Tmsh, Tdim)

  smallmatmatT!(u_in, data.P, t1)
  for q=1:numNodesPerElement
    for p=1:numNodesPerElement
      for j=1:numDofPerNode
        for i=1:numDofPerNode
          t1_jac[i, j, p, q] = data.P[p, q]*u_in_jac[i, j, q]
        end
      end
    end
  end

  if elnum == 1
    println("\nEntered applyLPSKernel_diff")
    println("t1_dot - \n", real(t1_jac[1, 1, 1, 1]))
  end
  # TODO: try to fuse p and q loops with above

  # multiply by du/dw
  for q=1:numNodesPerElement
    for p=1:numNodesPerElement
      q_p = ro_sview(q_el, :, p)
      getA0(entropy_vars, params, q_p, A0)
      for j=1:numDofPerNode
        for i=1:numDofPerNode
          for k=1:numDofPerNode
            # this should be A0[i, k], but A0 is symmetric
            t2_jac[i, j, p, q] += A0[k, i]*t1_jac[k, j, p, q]
          end
        end
      end
    end  # end p

    # contribution of d A0/dq
    q_q = ro_sview(q_el, :, q)
    getA0_diff(entropy_vars, params, q_q, q_dot, A0, A0_diff)
    for j=1:numDofPerNode
      for i=1:numDofPerNode
        for k=1:numDofPerNode
          # A0 is symmetric, so d A0[i, k]/dq = d A0[k, i]/dq
          t2_jac[i, j, q, q] += A0_diff[k, i, j]*t1[k, q]
        end
      end
    end

  end  # end q


  # compute t2
  # apply mass matrix and scaling
  @simd for i=1:numNodesPerElement
    q_i = ro_sview(q_el, :, i)
    getA0(entropy_vars, params, q_i, A0)
    t1_i = ro_sview(t1, :, i)
    t2_i = sview(t2, :, i)

    # scaling by wave speed
    lambda = zero(Tres)
    for d1=1:Tdim
      for d2=1:Tdim
        dir[d2] = dxidx[d1, d2, i]
      end
      lambda += getLambdaMax(params, q_i, dir)
    end
    fac = lambda*sbp.w[i]/Tdim

    # multiply by A0 and fac at the same time
    smallmatvec_kernel!(A0, t1_i, t2_i, fac, 0)
  end



  # compute lambda and its derivative at all volume nodes
  lambda_dot_p = zeros(Tres, numDofPerNode)
  for p=1:numNodesPerElement
    q_p = ro_sview(q_el, :, p)
    for d1=1:Tdim
      for d2=1:Tdim
        dir[d2] = dxidx[d1, d2, p]
      end
      lambdas[p] += getLambdaMax_diff(params, q_p, dir, lambda_dot_p)
      for i=1:numDofPerNode
        lambdas_jac[i, p] += lambda_dot_p[i]
      end
    end
    lambdas[p] /= Tdim
    for i=1:numDofPerNode
      lambdas_jac[i, p] /= Tdim
    end
  end

  # compute t2
  @simd for i=1:numNodesPerElement
    q_i = ro_sview(q_el, :, i)
    getA0(entropy_vars, params, q_i, A0)
    t1_i = ro_sview(t1, :, i)
    t2_i = sview(t2, :, i)
    fac = 1

    # multiply by A0 and fac at the same time
    smallmatvec_kernel!(A0, t1_i, t2_i, fac, 0)
  end


  # multipy by lambda and H
  for q=1:numNodesPerElement
    for p=1:numNodesPerElement
      fac = lambdas[p]*sbp.w[p]
      for j=1:numDofPerNode
        for i=1:numDofPerNode
          t3_jac[i, j, p, q] = fac*t2_jac[i, j, p, q]
        end
      end
    end  # end p

    for i=1:numDofPerNode
      for j=1:numDofPerNode
        t3_jac[i, j, q, q] += sbp.w[q]*lambdas_jac[j, q]*t2[i, q]
      end
    end

  end  # end q

  if elnum == 1
    println("t2_dot - \n", real(t2_jac[1, 1, 1, 1]))
  end

  # multiply by P^T
  for q=1:numNodesPerElement
    for p=1:numNodesPerElement
      for k=1:numNodesPerElement
        for j=1:numDofPerNode
          for i=1:numDofPerNode
            res_jac[i, j, p, q] += op(data.alpha*data.P[k, p]*t3_jac[i, j, k, q])
          end
        end
      end
    end
  end

  if elnum == 1
    println("t4_dot - \n", real(res_jac[1, 1, 1, 1]))
  end

  return nothing
end
