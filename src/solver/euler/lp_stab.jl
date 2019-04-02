# Local Projection (LP) stabilization

import SummationByParts: UnaryFunctor, Add, Subtract

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
    res_i = sview(eqn.res, :, :, i)

    # convert to entropy variables
    for j=1:mesh.numNodesPerElement
      q_j = ro_sview(q_i, :, j)
      w_j = sview(w_i, :, j)
      convertToEntropy(entropy_vars, eqn.params, q_j, w_j)
    end

    applyLPSKernel(eqn.params, data, sbp, q_i, w_i, res_i)

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
function applyLPSKernel(params::ParamType, data::LPSData, sbp::AbstractOperator,
                        q::AbstractMatrix{Tsol}, u_in::AbstractMatrix,
                        res::AbstractMatrix, op::UnaryFunctor=Subtract()
                        ) where {Tsol}
  # For scalar equations, the operator is applied -epsilon * P^T M P * u
  # For vector equations, P needs to be applied to all equations as once:
  # utmp^T = P*u^T
  # Instead, work with utmp = (P*u^T)^T = u*P^T

  numDofPerNode, numNodesPerElement = size(q)

  t1 = zeros(Tsol, numDofPerNode, numNodesPerElement)
  t2 = zeros(Tsol, numDofPerNode, numNodesPerElement)

  # apply P
  smallmatmatT!(u_in, data.P, t1)

  # apply mass matrix
  @simd for i=1:numNodesPerElement
    fac = sbp.w[i]  # TODO: should have 1/|J| factor?  The paper says no, but...
    @simd for j=1:numDofPerNode
      t1[j, i] *= fac
    end
  end

  # apply P^T
  smallmatmat!(t1, data.P, t2)
  @simd for i=1:numNodesPerElement
    @simd for j=1:numDofPerNode
      res[j, i] += op(data.alpha*t2[j, i])
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

  w_jac = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode,
                      mesh.numNodesPerElement)
  for i=1:mesh.numEl
    q_i = ro_sview(eqn.q, :, :, i)
    # convert to entropy variables
    for j=1:mesh.numNodesPerElement
      q_j = ro_sview(q_i, :, j)
      w_jac_j = sview(w_jac, :, :, j)
      getA0inv(entropy_vars, eqn.params, q_j, w_jac_j)
    end

    fill!(res_jac, 0)
    applyLPSKernel_diff(eqn.params, data, sbp, q_i, w_jac, res_jac)

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
function applyLPSKernel_diff(params::ParamType, data::LPSData,
                        sbp::AbstractOperator,
                        q::AbstractMatrix{Tsol}, u_in_jac::Abstract3DArray,
                        res_jac::Abstract4DArray, op::UnaryFunctor=Subtract()
                        ) where {Tsol}


  numDofPerNode, numNodesPerElement = size(q)

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

  return nothing
end
