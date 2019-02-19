# Interface for any diffusion model

"""
  Abstract type for diffusion tensors.  Used mostly for shock capturing.
  This type specifies the diffusion tensor.  It must implement the functions
  below (see the documentation for these specific function signatures
  for more detail:

  ```
    applyDiffusionTensor(obj::ShockDiffusion, w::AbstractMatrix,
                    i::Integer, nodes::AbstractVector,
                    dx::Abstract3DArray, flux::Abstract3DArray)
  ```
"""
abstract type AbstractDiffusion end


#------------------------------------------------------------------------------
# Functions each diffusion model must implement
# Users should not generally call these functions, use the API functions below
"""
  Applies the diffusion tensor to a given array, for an element 

  More specifically,

  ```
  for i in nodes
    for d=1:mesh.dim
      for d2=1:mesh.dim
        flux[:, i, d] = C[:, :, d1, d2]*dx[:, i, d2]
      end
    end
  end
  ```
  where C[:, :, :, :] is the diffusion tensor computed at state w[:, i].
  Note that only the nodes of the element specified by the `nodes` vector
  should be computed.  The values of `dx` not in `nodes` are undefined.

  **Inputs**

   * obj: the [`AbstractDiffusion`](@ref) object
   * w: the `numDofPerNode` x `numNodesPerElement` array of entropy variables
        for the element
   * i: the element number.  This is useful for diffusions where the coeffients
        are precomputed and stored in an array
   * nodes: a vector of integers describing which nodes of `flux` need to be
            populated.
   * dx: the values to multiply against, `numDofPerNode` x `numNodesPerElement`
         x `dim`

  **Inputs/Outputs**

   * flux: array to overwrite with the result, same size as `dx`

"""
function applyDiffusionTensor(obj::AbstractDiffusion, w::AbstractMatrix,
                    i::Integer, nodes::AbstractVector,
                    dx::Abstract3DArray, flux::Abstract3DArray)

  error("abstract fallback for applyDiffusionTensor() reached.  Did you forget to extend applyDiffusionTensor with a new method?")
end


"""
  Differentiated version of [`ApplyDiffusionTensor`](@ref)

  **Inputs**

   * obj: the [`AbstractDiffusion`](@ref) object
   * w: the `numDofPerNode` x `numNodesPerElement` array of entropy variables
        for the element
   * i: the element number.  This is useful for diffusions where the coeffients
        are precomputed and stored in an array
   * nodes: vector of integers identifying  nodes at which to compute the
            output (4th dimension of `t1_dot` and `t2_dot`)
   * t1: the values to multiply against, `numDofPerNode` x `numNodesPerElement`
         x `dim`
   * t1_dot: the derivative of t1, `numDofPerNode` x `numDofPerNode` x `dim`
            x `numNodesPerElement` x `numNodesPerElement`

  **Inputs/Outputs**

   * t2_dot: the derivative of t1, `numDofPerNode` x `numDofPerNode` x `dim`
            x `numNodesPerElement` x `numNodesPerElement`, overwritten


"""
function applyDiffusionTensor_diff(obj::AbstractDiffusion, w::AbstractMatrix,
                    i::Integer, nodes::AbstractVector, t1::Abstract3DArray,
                    t1_dot::Abstract5DArray, t2_dot::Abstract5DArray)

  error("abstract fallback for applyDiffusionTensor_diff() reached.  Did you forget to extend applyDiffusionTensor_diff() with a new method?")

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
   * nodes: vector of integers identifying  nodes at which to compute the
            output (4th dimension of `t1_dot` and `t2_dot`)
   * t1: the values Lambda is multiplied against for the primal operation,
         `numDofPerNode` x `numNodesPerElement` x `dim`
   * t1_dotL: 5D Jacobian containing the derivative of t1 wrt the left element
   * t1_dotR: 5D Jacobian containing the derivative of t1 wrt the right element

  **Inputs/Outputs**

   * t2_dotL: 5D Jacobian overwritten with result for elementL
   * t2_dotR: 5D Jacobian overwritten with result for elementR
"""
function applyDiffusionTensor_diff(obj::AbstractDiffusion, w::AbstractMatrix,
                    i::Integer, nodes::AbstractVector, t1::Abstract3DArray,
                    t1_dotL::Abstract5DArray,
                    t1_dotR::Abstract5DArray, t2_dotL::Abstract5DArray,
                    t2_dotR::Abstract5DArray)


  error("abstract fallback for applyDiffusionTensor_diff() reached.  Did you forget to extend applyDiffusionTensor_diff() with a new method?")
end


#------------------------------------------------------------------------------
# API functions, common to all diffusion models (implemented in terms of the
# above)

"""
  This function applies the diffusion tensor to only the nodes in the stencil of
  R for a given face.  This is useful when `flux` will be interpolated to the
  face in the next step, and can be significantly faster than applying the
  diffusion tensor to all nodes of the element.

  **Inputs**

   * obj: the [`AbstractDiffusion`](@ref) object
   * w: the `numDofPerNode` x `numNodesPerElement` array of entropy variables
        for the element
   * i: the element number.  This is useful for diffusions where the coeffients
        are precomputed and stored in an array
   * sbpface: an `AbstractFace` object that describes the stencil of R
   * face: integer describing which face of the element will be used.
   * dx: the values to multiply against, `numDofPerNode` x `numNodesPerElement`
         x `dim`

  **Inputs/Outputs**

   * flux: array to overwrite with the result, same size as `dx`


"""
function applyDiffusionTensor(obj::AbstractDiffusion, w::AbstractMatrix,
                    i::Integer, sbpface::AbstractFace, face::Integer,
                    dx::Abstract3DArray, flux::Abstract3DArray)

  nodes = sview(sbpface.perm, :, face)
  applyDiffusionTensor(obj, w, i, nodes, dx, flux)

  return nothing
end

"""
  This method applies the diffusion tensor to all nodes of the element.

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
function applyDiffusionTensor(obj::AbstractDiffusion, w::AbstractMatrix,
                    i::Integer,
                    dx::Abstract3DArray, flux::Abstract3DArray)

  numNodesPerElement = size(dx, 2)
  nodes = 1:numNodesPerElement
  applyDiffusionTensor(obj, w, i, nodes, dx, flux)

  return nothing
end

"""
  This method computes the Jacobian of applying the diffusion tensor to
  only the nodes on the stencil of R.

  **Inputs**

   * obj: the [`AbstractDiffusion`](@ref) object
   * w: the `numDofPerNode` x `numNodesPerElement` array of entropy variables
        for the element
   * i: the element number.  This is useful for diffusions where the coeffients
        are precomputed and stored in an array
   * sbpface: the `AbstractFace` that identifies the stencil of R
   * face: the local face number
   * t1: the values to multiply against, `numDofPerNode` x `numNodesPerElement`
         x `dim`
   * t1_dot: the derivative of t1, `numDofPerNode` x `numDofPerNode` x `dim`
            x `numNodesPerElement` x `numNodesPerElement`

  **Inputs/Outputs**

   * t2_dot: the derivative of t1, `numDofPerNode` x `numDofPerNode` x `dim`
            x `numNodesPerElement` x `numNodesPerElement`, overwritten
"""
function applyDiffusionTensor_diff(obj::AbstractDiffusion, w::AbstractMatrix,
                    i::Integer, sbpface::AbstractFace, face::Integer,
                    t1::Abstract3DArray,
                    t1_dot::Abstract5DArray, t2_dot::Abstract5DArray)

  nodes = sview(sbpface.perm, :, face)
  applyDiffusionTensor_diff(obj, w, i, nodes, t1, t1_dot, t2_dot)

  return nothing
end


"""
  This method computes the Jacobian of applying the diffusion tensor to
  all nodes of the element

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
function applyDiffusionTensor_diff(obj::AbstractDiffusion, w::AbstractMatrix,
                    i::Integer,
                    t1::Abstract3DArray,
                    t1_dot::Abstract5DArray, t2_dot::Abstract5DArray)

  nodes = 1:size(t1, 2)
  applyDiffusionTensor_diff(obj, w, i, nodes, t1, t1_dot, t2_dot)

  return nothing
end

"""
  Method for computing the Jacobian of applying the diffusion tensor to
  only the nodes in the stencil of R, for the case where `t1` depends
  on both `qL` and `qR`.

  **Inputs**

   * obj: [`ShockDiffusion`](@ref)
   * w: entropy variables for the element, `numDofPerNode` x `numNodesPerElement`
   * i: element number
   * sbpface: the `AbstractFace` that describes the stencil of R
   * face: the local face number
   * t1: the values Lambda is multiplied against for the primal operation,
         `numDofPerNode` x `numNodesPerElement` x `dim`
   * t1_dotL: 5D Jacobian containing the derivative of t1 wrt the left element
   * t1_dotR: 5D Jacobian containing the derivative of t1 wrt the right element

  **Inputs/Outputs**

   * t2_dotL: 5D Jacobian overwritten with result for elementL
   * t2_dotR: 5D Jacobian overwritten with result for elementR
"""
function applyDiffusionTensor_diff(obj::AbstractDiffusion, w::AbstractMatrix,
                    i::Integer, sbpface::AbstractFace, face::Integer,
                    t1::Abstract3DArray,
                    t1_dotL::Abstract5DArray, t1_dotR::Abstract5DArray,
                    t2_dotL::Abstract5DArray, t2_dotR::Abstract5DArray)

  nodes = sview(sbpface.perm, :, face)
  applyDiffusionTensor_diff(obj, w, i, nodes, t1, t1_dotL, t1_dotR,
                            t2_dotL, t2_dotR)

  return nothing
end

"""
  Method for computing the Jacobian of applying the diffusion tensor to
  all the nodes of the element, for the case where `t1` depends
  on both `qL` and `qR`.

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
function applyDiffusionTensor_diff(obj::AbstractDiffusion, w::AbstractMatrix,
                    i::Integer,
                    t1::Abstract3DArray,
                    t1_dotL::Abstract5DArray, t1_dotR::Abstract5DArray,
                    t2_dotL::Abstract5DArray, t2_dotR::Abstract5DArray)

  nodes = 1:size(t1_dotL, 4)
  applyDiffusionTensor_diff(obj, w, i, nodes, t1, t1_dotL, t1_dotR,
                            t2_dotL, t2_dotR)

  return nothing
end



