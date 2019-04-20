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
    @simd for j=1:mesh.numNodesPerElement
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

  @unpack data t1 t2 A0 dir
  numDofPerNode, numNodesPerElement = size(q)

  # apply P
  smallmatmatT!(u_in, data.P, t1)

  # apply mass matrix and scaling
  @simd for i=1:numNodesPerElement
    q_i = ro_sview(q, :, i)
    getA0(entropy_vars, params, q_i, A0)
    t1_i = ro_sview(t1, :, i)
    t2_i = sview(t2, :, i)

    # scaling by wave speed
    lambda = zero(Tres)
    @simd for d1=1:Tdim
      @simd for d2=1:Tdim
        dir[d2] = dxidx[d1, d2, i]
      end
      lambda += getLambdaMax(params, q_i, dir)
    end
    fac = lambda*sbp.w[i]/Tdim

    # multiply by A0 and fac at the same time
    smallmatvec_kernel!(A0, t1_i, t2_i, fac, 0)
  end

  # apply P^T
  smallmatmat!(t2, data.P, t1)
  @simd for i=1:numNodesPerElement
    @simd for j=1:numDofPerNode
      res[j, i] += op(data.alpha*t1[j, i])
    end
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
    @simd for j=1:mesh.numNodesPerElement
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

  @unpack data t1 t1_jac t2 t2_jac t3_jac A0 A0_diff q_dot lambdas lambdas_jac
  @unpack data dir lambda_dot_p

  # compute required primal quantities (t1 and t2)
  smallmatmatT!(u_in, data.P, t1)

  @simd for i=1:numNodesPerElement
    q_i = ro_sview(q_el, :, i)
    getA0(entropy_vars, params, q_i, A0)
    t1_i = ro_sview(t1, :, i)
    t2_i = sview(t2, :, i)
    fac = 1

    # multiply by A0 and fac at the same time
    smallmatvec_kernel!(A0, t1_i, t2_i, fac, 0)
  end

  # compute lambda and its derivative at all volume nodes
  fill!(lambdas, 0); fill!(lambdas_jac, 0)
  @simd for p=1:numNodesPerElement
    q_p = ro_sview(q_el, :, p)
    @simd for d1=1:Tdim
      @simd for d2=1:Tdim
        dir[d2] = dxidx[d1, d2, p]
      end
      lambdas[p] += getLambdaMax_diff(params, q_p, dir, lambda_dot_p)
      @simd for i=1:numDofPerNode
        lambdas_jac[i, p] += lambda_dot_p[i]
      end
    end
    lambdas[p] /= Tdim
    @simd for i=1:numDofPerNode
      lambdas_jac[i, p] /= Tdim
    end
  end


  #--------------------------
  # compute derivative

  @simd for q=1:numNodesPerElement

    # t1
    @simd for p=1:numNodesPerElement
      @simd for j=1:numDofPerNode
        @simd for i=1:numDofPerNode
          t1_jac[i, j, p] = data.P[p, q]*u_in_jac[i, j, q]
        end
      end
    end

    # A0 * t1_dot
    @simd for p=1:numNodesPerElement
      q_p = ro_sview(q_el, :, p)
      getA0(entropy_vars, params, q_p, A0)
      @simd for j=1:numDofPerNode
        @simd for i=1:numDofPerNode
          t2_jac[i, j, p] = 0
          @simd for k=1:numDofPerNode
            # this should be A0[i, k], but A0 is symmetric
            t2_jac[i, j, p] += A0[k, i]*t1_jac[k, j, p]
          end
        end
      end
    end  # end p

    # contribution of d A0/dq * t1
    q_q = ro_sview(q_el, :, q)
    getA0_diff(entropy_vars, params, q_q, q_dot, A0, A0_diff)
    @simd for j=1:numDofPerNode
      @simd for i=1:numDofPerNode
        @simd for k=1:numDofPerNode
          # A0 is symmetric, so d A0[i, k]/dq = d A0[k, i]/dq
          t2_jac[i, j, q] += A0_diff[k, i, j]*t1[k, q]
        end
      end
    end

    # multiply by lambda and H
    @simd for p=1:numNodesPerElement
      fac = lambdas[p]*sbp.w[p]
      @simd for j=1:numDofPerNode
        @simd for i=1:numDofPerNode
          t3_jac[i, j, p] = fac*t2_jac[i, j, p]
        end
      end
    end  # end p


    @simd for i=1:numDofPerNode
      @simd for j=1:numDofPerNode
        t3_jac[i, j, q] += sbp.w[q]*lambdas_jac[j, q]*t2[i, q]
      end
    end

    # multiply by P^T
    @simd for p=1:numNodesPerElement
      @simd for k=1:numNodesPerElement
        @simd for j=1:numDofPerNode
          @simd for i=1:numDofPerNode
            res_jac[i, j, p, q] += op(data.alpha*data.P[k, p]*t3_jac[i, j, k])
          end
        end
      end
    end

  end  # end q

  return nothing
end


#------------------------------------------------------------------------------
# revq

"""
  Main function for applying LPS stabilization revq

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
"""
function applyLPStab_revq(mesh, sbp, eqn, opts)

  lps_data = eqn.params.lps_data
  entropy_vars = lps_data.entropy_vars

  _applyLPStab_revq(mesh, sbp, eqn, opts, lps_data, entropy_vars)

  return nothing
end

@noinline function _applyLPStab_revq(mesh, sbp, eqn, opts, data::LPSData{Tsol, Tres},
                       entropy_vars::AbstractVariables) where {Tsol, Tres}

  w_i = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerElement)
  w_i_bar = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)
  A0inv = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode)
  for i=1:mesh.numEl
    q_i = ro_sview(eqn.q, :, :, i); q_i_bar = sview(eqn.q_bar, :, :, i)
    dxidx = ro_sview(mesh.dxidx, :, :, :, i)
    res_i = sview(eqn.res, :, :, i); res_i_bar = sview(eqn.res_bar, :, :, i)

    # convert to entropy variables
    @simd for j=1:mesh.numNodesPerElement
      q_j = ro_sview(q_i, :, j)
      w_j = sview(w_i, :, j)
      convertToEntropy(entropy_vars, eqn.params, q_j, w_j)
    end

    fill!(w_i_bar, 0)
    applyLPSKernel_revq(eqn.params, data, sbp, entropy_vars, q_i, q_i_bar,
                        w_i, w_i_bar, dxidx, res_i, res_i_bar, i)

    for j=1:mesh.numNodesPerElement
      q_j = ro_sview(q_i, :, j)
      getA0inv(entropy_vars, eqn.params, q_j, A0inv)
      smallmatvec_kernel!(A0inv, sview(w_i_bar, :, j), sview(q_i_bar, :, j), 1, 1)
    end
  end

  return nothing
end




function applyLPSKernel_revq(params::ParamType{Tdim}, data::LPSData,
                        sbp::AbstractOperator,
                        entropy_vars::AbstractVariables,
                        q::AbstractMatrix{Tsol}, q_bar::AbstractMatrix,
                        u_in::AbstractMatrix, u_in_bar::AbstractMatrix,
                        dxidx::Abstract3DArray{Tmsh},
                        res::AbstractMatrix{Tres}, res_bar::AbstractMatrix,
                        elnum, op::UnaryFunctor=Subtract()
                        ) where {Tsol, Tres, Tmsh, Tdim}
  # For scalar equations, the operator is applied -epsilon * P^T M P * u
  # For vector equations, P needs to be applied to all equations as once:
  # utmp^T = P*u^T
  # Instead, work with utmp = (P*u^T)^T = u*P^T

  @unpack data t1 t2 A0 dir
  numDofPerNode, numNodesPerElement = size(q)

  # apply P
  smallmatmatT!(u_in, data.P, t1)
#=
  # apply mass matrix and scaling
  @simd for i=1:numNodesPerElement
    q_i = ro_sview(q, :, i)
    getA0(entropy_vars, params, q_i, A0)
    t1_i = ro_sview(t1, :, i)
    t2_i = sview(t2, :, i)

    # scaling by wave speed
    lambda = zero(Tres)
    @simd for d1=1:Tdim
      @simd for d2=1:Tdim
        dir[d2] = dxidx[d1, d2, i]
      end
      lambda += getLambdaMax(params, q_i, dir)
    end
    fac = lambda*sbp.w[i]/Tdim

    # multiply by A0 and fac at the same time
    smallmatvec_kernel!(A0, t1_i, t2_i, fac, 0)
  end
=#
#=
  # apply P^T
  smallmatmat!(t2, data.P, t1)
  @simd for i=1:numNodesPerElement
    @simd for j=1:numDofPerNode
      res[j, i] += op(data.alpha*t1[j, i])
    end
  end
=#
  #--------------------------
  # reverse sweep
  @unpack data t1_bar t2_bar tmp tmp_bar A0_bar

  for i=1:numNodesPerElement
    for j=1:numDofPerNode
      t1_bar[j, i] = op(data.alpha*res_bar[j, i])
    end
  end

  smallmatmatT!(t1_bar, data.P, t2_bar)

  fill!(t1_bar, 0)

  for i=1:numNodesPerElement

    # need fac, tmp, A0
    q_i = ro_sview(q, :, i); q_i_bar = sview(q_bar, :, i)
    t1_i = ro_sview(t1, :, i)

    getA0(entropy_vars, params, q_i, A0)

    lambda = zero(Tres)
    @simd for d1=1:Tdim
      @simd for d2=1:Tdim
        dir[d2] = dxidx[d1, d2, i]
      end
      lambda += getLambdaMax(params, q_i, dir)
    end
    fac = lambda*sbp.w[i]/Tdim

    smallmatvec_kernel!(A0, t1_i, tmp, 1, 0)

    # reverse mode
    fill!(tmp_bar, 0)
    fac_bar = zero(Tres)
    for j=1:numDofPerNode
      tmp_bar[j] = fac*t2_bar[j, i]
      fac_bar   += tmp[j]*t2_bar[j, i]
    end

    t1_i_bar = sview(t1_bar, :, i)
    smallmatvec!(A0, tmp_bar, t1_i_bar)
    for j=1:numDofPerNode
      for k=1:numDofPerNode
        A0_bar[k, j] = tmp_bar[k]*t1[j, i]
      end
    end

    lambda_bar = sbp.w[i]*fac_bar/Tdim
    @simd for d1=1:Tdim
      @simd for d2=1:Tdim
        dir[d2] = dxidx[d1, d2, i]
      end
      getLambdaMax_revq(params, q_i, q_i_bar, dir, lambda_bar)
    end

    getA0_revq(entropy_vars, params, q_i, q_i_bar, A0, A0_bar)
  end

  smallmatmat!(t1_bar, data.P, u_in_bar)


  return nothing
end


#------------------------------------------------------------------------------
# revm

"""
  Main function for applying LPS stabilization revm

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
"""
function applyLPStab_revm(mesh, sbp, eqn, opts)

  lps_data = eqn.params.lps_data
  entropy_vars = lps_data.entropy_vars

  _applyLPStab_revm(mesh, sbp, eqn, opts, lps_data, entropy_vars)

  return nothing
end

@noinline function _applyLPStab_revm(mesh, sbp, eqn, opts, data::LPSData{Tsol},
                       entropy_vars::AbstractVariables) where {Tsol}

  w_i = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerElement)
  for i=1:mesh.numEl
    q_i = ro_sview(eqn.q, :, :, i)
    dxidx = ro_sview(mesh.dxidx, :, :, :, i);
    dxidx_bar = sview(mesh.dxidx_bar, :, :, :, i)
    res_i = sview(eqn.res, :, :, i); res_i_bar = sview(eqn.res_bar, :, :, i)

    # convert to entropy variables
    @simd for j=1:mesh.numNodesPerElement
      q_j = ro_sview(q_i, :, j)
      w_j = sview(w_i, :, j)
      convertToEntropy(entropy_vars, eqn.params, q_j, w_j)
    end

    applyLPSKernel_revm(eqn.params, data, sbp, entropy_vars, q_i, w_i,
                        dxidx, dxidx_bar, res_i, res_i_bar, i)
  end

  return nothing
end




function applyLPSKernel_revm(params::ParamType{Tdim}, data::LPSData,
                      sbp::AbstractOperator,
                      entropy_vars::AbstractVariables,
                      q::AbstractMatrix{Tsol}, u_in::AbstractMatrix,
                      dxidx::Abstract3DArray{Tmsh}, dxidx_bar::Abstract3DArray,
                      res::AbstractMatrix{Tres}, res_bar::AbstractMatrix,
                      elnum, op::UnaryFunctor=Subtract()
                        ) where {Tsol, Tres, Tmsh, Tdim}
  # For scalar equations, the operator is applied -epsilon * P^T M P * u
  # For vector equations, P needs to be applied to all equations as once:
  # utmp^T = P*u^T
  # Instead, work with utmp = (P*u^T)^T = u*P^T

  @unpack data t1 t2 A0 dir
  numDofPerNode, numNodesPerElement = size(q)

  # apply P
  smallmatmatT!(u_in, data.P, t1)
#=
  # apply mass matrix and scaling
  @simd for i=1:numNodesPerElement
    q_i = ro_sview(q, :, i)
    getA0(entropy_vars, params, q_i, A0)
    t1_i = ro_sview(t1, :, i)
    t2_i = sview(t2, :, i)

    # scaling by wave speed
    lambda = zero(Tres)
    @simd for d1=1:Tdim
      @simd for d2=1:Tdim
        dir[d2] = dxidx[d1, d2, i]
      end
      lambda += getLambdaMax(params, q_i, dir)
    end
    fac = lambda*sbp.w[i]/Tdim

    # multiply by A0 and fac at the same time
    smallmatvec_kernel!(A0, t1_i, t2_i, fac, 0)
  end
=#
#=
  # apply P^T
  smallmatmat!(t2, data.P, t1)
  @simd for i=1:numNodesPerElement
    @simd for j=1:numDofPerNode
      res[j, i] += op(data.alpha*t1[j, i])
    end
  end
=#
  #--------------------------
  # reverse sweep
  @unpack data t1_bar t2_bar tmp dir_bar

  for i=1:numNodesPerElement
    for j=1:numDofPerNode
      t1_bar[j, i] = op(data.alpha*res_bar[j, i])
    end
  end

  smallmatmatT!(t1_bar, data.P, t2_bar)

  for i=1:numNodesPerElement

    # need fac, tmp, A0
    q_i = ro_sview(q, :, i)
    t1_i = ro_sview(t1, :, i)

    getA0(entropy_vars, params, q_i, A0)
    smallmatvec_kernel!(A0, t1_i, tmp, 1, 0)

    # reverse mode
    fac_bar = zero(Tres)
    for j=1:numDofPerNode
      fac_bar   += tmp[j]*t2_bar[j, i]
    end

    lambda_bar = sbp.w[i]*fac_bar/Tdim
    @simd for d1=1:Tdim
      @simd for d2=1:Tdim
        dir[d2] = dxidx[d1, d2, i]
        dir_bar[d2] = 0
      end
      getLambdaMax_revm(params, q_i, dir, dir_bar, lambda_bar)

      @simd for d2=1:Tdim
        dxidx_bar[d1, d2, i] += dir_bar[d2]
      end
    end
  end

  return nothing
end


