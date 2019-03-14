# the elementwise-constant diffusion model used for shock capturing


#------------------------------------------------------------------------------
# Required functions

"""
  Implements an element-wise constant diffusion model.
"""
function applyDiffusionTensor(obj::ShockDiffusion, sbp::AbstractOperator,
                    params::ParamType, q::AbstractMatrix,
                    w::AbstractMatrix,
                    coords::AbstractMatrix, dxidx::Abstract3DArray,
                    jac::AbstractVector,
                    i::Integer, nodes::AbstractVector,
                    dx::Abstract3DArray, flux::Abstract3DArray{Tres}
                   ) where {Tres}

  # For shock capturing it is much easier because the viscoscity is diagonal
  # for C[:, :, j, j] and zero for the cross terms

  numDofPerNode, numNodesPerElement, dim = size(flux)

  Se = obj.Se
  ee = obj.ee
  getShockSensor(params, sbp, obj.sensor, q, coords, dxidx, jac, Se, ee)

  @simd for d=1:dim
    @simd for j in nodes
      ee_j = ee[d, j]
      @simd for k=1:numDofPerNode
        flux[k, j, d] = ee_j*dx[k, j, d]
      end
    end
  end

  return nothing
end

#------------------------------------------------------------------------------
# Functions that do individual computations for the differentiated code

"""
  This function applies Lambda do a jacobian, ie.

    dt2/dq = Lambda * dt1/dq.

  Note that Lambda is also a function of q, and this function does *not* compute
  the dLambda/dq term

  **Inputs**

   * obj: ShockDiffusion
   * sbp: SBP Operator
   * params: ParamType
   * q: conservative variables for this element
   * w: entropy variables for this element
   * coords: xyz coordinates of nodes
   * dxidx: scaled mapping jacobian
   * jac: mapping jacobian determinant
   * elnum: element number
   * nodes: vector of nodes to compute the output for.
   * t1_dot: a 5D Jacobian
   
  **Inputs/Outputs**

   * t2_dot: a  5D Jacobian to be overwritten
"""
function applyLambda_diff(obj::ShockDiffusion, sbp::AbstractOperator,
                          params::ParamType, q::AbstractMatrix,
                          w::AbstractMatrix,
                          coords::AbstractMatrix, dxidx::Abstract3DArray,
                          jac::AbstractVector,
                          elnum::Integer, nodes::AbstractVector,
                          t1_dot::Abstract5DArray, t2_dot::Abstract5DArray{Tres}
                         ) where {Tres}

  numDofPerNode = size(t1_dot, 1)
  dim = size(t1_dot, 3)
  numNodesPerElement = size(t1_dot, 4)

  Se = obj.Se
  ee = obj.ee
  getShockSensor(params, sbp, obj.sensor, q, coords, dxidx, jac, Se, ee)

  @simd for q=1:numNodesPerElement
    @simd for p in nodes
      @simd for d=1:dim
        ee_d = ee[d, p]
        @simd for j=1:numDofPerNode
          @simd for i=1:numDofPerNode
            t2_dot[i, j, d, p, q] = ee_d*t1_dot[i, j, d, p, q]
          end
        end
      end
    end
  end

  return nothing
end


"""
  This function compute the chain rule term for Lambda, ie.

  dt2/dq += dLambda/dq * t1

  **Inputs**

   * obj: ShockDiffusion
   * sbp: SBP operator
   * params: ParamType
   * q: conservative variables for the element
   * w: entropy variables for the element
   * coords: xyz coordinates of nodes
   * dxidx: scaled mapping jacobian
   * jac: mapping jacobian determinant
   * elnum: element number
   * nodes: vector of nodes to compute the output for
   * t1: array that Lambda is multiplied against, `numDofPerNode` x
        `numNodesPerElement` x `dim`

  **Inputs/Outputs**

   * t2_dot: 5D Jacobian to update
"""
function applyLambdaDot_diff(obj::ShockDiffusion, sbp::AbstractOperator,
                             params::ParamType,
                             q_el::AbstractMatrix, w_el::AbstractMatrix,
                             coords::AbstractMatrix, dxidx::Abstract3DArray,
                             jac::AbstractVector,
                             elnum::Integer, nodes::AbstractVector,
                             t1::Abstract3DArray, t2_dot::Abstract5DArray{Tres}) where {Tres}
  
  numDofPerNode = size(t1, 1)
  dim = size(t2_dot, 3)
  numNodesPerElement = size(t1, 2)

  Se_dot = obj.Se_dot
  ee_dot = obj.ee_dot
  is_constant = getShockSensor_diff(params, sbp, obj.sensor, q_el,
                                    coords, dxidx, jac, Se_dot, ee_dot)

  if !is_constant
    @simd for q=1:numNodesPerElement
      @simd for p in nodes
        @simd for d=1:dim
          @simd for j=1:numDofPerNode
            ee_dot_j = ee_dot[d, j, p, q]
            @simd for i=1:numDofPerNode
              t2_dot[i, j, d, p, q] += ee_dot_j*t1[i, p, d]
            end
          end
        end
      end
    end
  end

  return nothing
end


#------------------------------------------------------------------------------
# Required functions

"""
  Differentiated version of [`ApplyDiffusionTensor`](@ref)  for `ShockDiffusion`
"""
function applyDiffusionTensor_diff(obj::ShockDiffusion, sbp::AbstractOperator,
                    params::ParamType, 
                    q_el::AbstractMatrix, w_el::AbstractMatrix,
                    coords::AbstractMatrix, dxidx::Abstract3DArray,
                    jac::AbstractVector,
                    i::Integer, nodes::AbstractVector, t1::Abstract3DArray,
                    t1_dot::Abstract5DArray, t2_dot::Abstract5DArray)

  # this is simpler than a physically accurate diffusion because the diffusion
  # value is constant for the entire element.
  # It is complicated slightly that the diffusion value might be a function
  # of the solution at all nodes of the element

  # do Lambda * t1_dot
  applyLambda_diff(obj, sbp, params, q_el, w_el, coords, dxidx, jac, i, nodes,
                   t1_dot, t2_dot)

  # do dLambda/dq * t1
  applyLambdaDot_diff(obj, sbp, params, q_el, w_el, coords, dxidx, jac,
                      i, nodes, t1, t2_dot)

  return nothing
end

"""
  Method for the case where t1 is a function of 2 states, qL and qR
"""
function applyDiffusionTensor_diff(obj::ShockDiffusion, sbp::AbstractOperator,
                    params::ParamType,
                    q_el::AbstractMatrix, w_el::AbstractMatrix,
                    coords::AbstractMatrix, dxidx::Abstract3DArray,
                    jac::AbstractVector,
                    i::Integer, nodes::AbstractVector, t1::Abstract3DArray,
                    t1_dotL::Abstract5DArray, t1_dotR::Abstract5DArray,
                    t2_dotL::Abstract5DArray, t2_dotR::Abstract5DArray)

  # if isleft:
  # compute dt2/dqL = Lambda(qL) * dt1/dqL + d Lambda/dqL * t2
  #     and dt2/dqR = Lambda(qL) * dt1/dqR
  # else: add the d Lambda/dqR term to t2_dotR
  applyLambda_diff(obj, sbp, params, q_el, w_el, coords, dxidx, jac, i, nodes,
                   t1_dotL, t2_dotL)
  applyLambda_diff(obj, sbp, params, q_el, w_el, coords, dxidx, jac, i, nodes,
                   t1_dotR, t2_dotR)

  # Lambda is only a function of qL (we assume) so only add its contribution
  # to t2_dotL
  applyLambdaDot_diff(obj, sbp, params, q_el, w_el, coords, dxidx, jac, i,
                      nodes, t1, t2_dotL)

  return nothing
end
