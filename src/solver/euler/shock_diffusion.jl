# the elementwise-constant diffusion model used for shock capturing

"""
  Applies the diffusion tensor to a given array, for an element 

  More specifically,

  ```
  for i=1:numNodesPerElement
    for d=1:mesh.dim
      for d2=1:mesh.dim
        flux[:, i, d] = C[:, :, d1, d2]*dx[:, i, d2]
      end
    end
  end
  ```
  where C[:, :, :, :] is the diffusion tensor computed at state w[:, i]

  **Inputs**

   * obj: the [`AbstractDiffusion`](@ref) object
   * w: the `numDofPerNode` x `numNodesPerElement` array of entropy variables
        for the element
   * i: the element number.  This is useful for diffusions where the coeffients
        are precomputed and stored in an array
   * dx: the values to multiply against, `numDofPerNode` x `numNodesPerElement`
         x `dim`

  **Inputs/Outputs**

   * flux: array to overwrite with the result, same size as `dx`

"""
function applyDiffusionTensor(obj::ShockDiffusion, w::AbstractMatrix,
                    i::Integer, dx::Abstract3DArray, flux::Abstract3DArray)

  # For shock capturing it is much easier because the viscoscity is diagonal
  # for C[:, :, j, j] and zero for the cross terms

  numDofPerNode, numNodesPerElement, dim = size(flux)
  ee = obj.ee[i]
  @simd for d=1:dim
    @simd for j=1:numNodesPerElement
      @simd for k=1:numDofPerNode
        flux[k, j, d] = ee*dx[k, j, d]
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
   * w: entropy variables for this element
   * elnum: element number
   * t1_dot: a 5D Jacobian
   
  **Inputs/Outputs**

   * t2_dot: a  5D Jacobian to be overwritten
"""
function applyLambda_diff(obj::ShockDiffusion, w::AbstractMatrix, elnum::Integer,
                          t1_dot::Abstract5DArray, t2_dot::Abstract5DArray)

  ee = obj.ee[elnum]
  @simd for i=1:length(t2_dot)
    t2_dot[i] = ee*t1_dot[i]
  end

  return nothing
end


"""
  This function compute the chain rule term for Lambda, ie.

  dt2/dq += dLambda/dq * t1

  **Inputs**

   * obj: ShockDiffusion
   * w: entropy variables for the element
   * elnum: element number
   * t1: array that Lambda is multiplied against, `numDofPerNode` x
        `numNodesPerElement` x `dim`

  **Inputs/Outputs**

   * t2_dot: 5D Jacobian to update
"""
function applyLambdaDot_diff(obj::ShockDiffusion, w::AbstractMatrix,
                             elnum::Integer,
                             t1::Abstract3DArray, t2_dot::Abstract5DArray)
  
  numDofPerNode = size(t1, 1)
  dim = size(t2_dot, 3)
  numNodesPerElement = size(t1, 2)

  if obj.is_nonlinear[elnum]
    @simd for q=1:numNodesPerElement
      @simd for p=1:numNodesPerElement
        @simd for d=1:dim
          @simd for i=1:numDofPerNode
            @simd for j=1:numDofPerNode
              t2_dot[i, j, d, p, q] += obj.ee_dot[j, q, elnum]*t1[i, p, d]
            end
          end
        end
      end
    end
  end

  return nothing
end


#------------------------------------------------------------------------------
# Externally visible functions

"""
  Differentiated version of [`ApplyDiffusionTensor`](@ref)

  **Inputs**

   * obj: the [`AbstractDiffusion`](@ref) object
   * w: the `numDofPerNode` x `numNodesPerElement` array of entropy variables
        for the element
   * i: the element number.  This is useful for diffusions where the coeffients
        are precomputed and stored in an array
   * t1: the values to multiply against, `numDofPerNode` x `numNodesPerElement`
         x `dim`
   * t1_dot: the derivative of t1, `numDofPerNode` x `numDofPerNode` x `dim`
            x `numNodesPerElement` x `numNodesPerElement`

  **Inputs/Outputs**

   * t2_dot: the derivative of t1, `numDofPerNode` x `numDofPerNode` x `dim`
            x `numNodesPerElement` x `numNodesPerElement`, overwritten


"""
function applyDiffusionTensor_diff(obj::ShockDiffusion, w::AbstractMatrix,
                    i::Integer, t1::Abstract3DArray, t1_dot::Abstract5DArray,
                    t2_dot::Abstract5DArray)
#TODO: consider an inplace version?

  # this is simpler than a physically accurate diffusion because the diffusion
  # value is constant for the entire element.
  # It is complicated slightly that the diffusion value might be a function
  # of the solution at all nodes of the element

  # do Lambda * t1_dot
  applyLambda_diff(obj, w, i, t1_dot, t2_dot)

  # do dLambda/dq * t1
  applyLambdaDot_diff(obj, w, i, t1, t2_dot)

  return nothing
end

"""
  Method for the case where t1 is a function of 2 states, qL and qR, ie.

  t2_dotL = Lambda * t1_dotL + dLambda/dq * t1
  t2_dotR = Lambda * t1_dotR

  Note that the chain rule term dLambda/dq is added to t2_dotL, so if using
  this function for elementR of a given interface, the argument `t1_dotL`
  should really be the derivative of `t1` with respect to qR (and similarly
  `t2_dotL` and `t2_dotR` should be reversed.

  **Inputs**

   * obj: [`ShockDiffusion`](@ref)
   * w: entropy variables for the element, `numDofPerNode` x `numNodesPerElement`
   * i: element number
   * t1: the values Lambda is multiplied against for the primal operation,
         `numDofPerNode` x `numNodesPerElement` x `dim`
   * t1_dotL: 5D Jacobian containing the derivative of t1 wrt the left element
   * t1_dotR: 5D Jacobian containing the derivative of t1 wrt the right element

  **Inputs/Outputs**

   * t2_dotL: 5D Jacobian overwritten with result for elementL
   * t2_dotR: 5D Jacobian overwritten with result for elementR
"""
function applyDiffusionTensor_diff(obj::ShockDiffusion, w::AbstractMatrix,
                    i::Integer, t1::Abstract3DArray, t1_dotL::Abstract5DArray,
                    t1_dotR::Abstract5DArray, t2_dotL::Abstract5DArray,
                    t2_dotR::Abstract5DArray)

  # if isleft:
  # compute dt2/dqL = Lambda(qL) * dt1/dqL + d Lambda/dqL * t2
  #     and dt2/dqR = Lambda(qL) * dt1/dqR
  # else: add the d Lambda/dqR term to t2_dotR
  applyLambda_diff(obj, w, i, t1_dotL, t2_dotL)
  applyLambda_diff(obj, w, i, t1_dotR, t2_dotR)

  # Lambda is only a function of qL (we assume) so only add its contribution
  # to t2_dotL
  applyLambdaDot_diff(obj, w, i, t1, t2_dotL)

  return nothing
end
