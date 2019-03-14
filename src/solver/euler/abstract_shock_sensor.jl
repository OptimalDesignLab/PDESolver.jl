# abstract type and required methods for a shock sensor

"""
  Abstract type for all shock sensors.  The purpose of shock sensors it to tell
  whether or not there is a shock, and if so, what the viscoscity coefficient
  should be.  The type of artificial viscoscity to add is determined by the
  [`AbstractShockCapturing`](@ref).  Each type must implement:

  
  ```
  getShockSensor
  getShockSensor_diff
  ```

  Each type may optionally implement

  ```
  isShockElement
  ```
"""
abstract type AbstractShockSensor end


"""
  Returns true if the element has non-zero diffusion in it, false otherwise.

  This function has a fallback defined in terms of `getShockSensor`, but types
  may be able to implement this more efficiently.

  **Inputs**

   * params: a `AbstractParamType`
   * sbp: SBP operator
   * sensor: the `AbstractShockSensor` object
   * q: solution on a particular element, `numDofPerNode` x `numNodesPerElement`
   * coords: `dim` x `numNodesPerElement` array containing the coordinates
             of each node of the element
   * dxidx: `dim` x `dim` x `numNodesPerElement` array containing the (scaled)
            mapping jacobian dxi/dx
   * jac: mapping jacobian determinant for each node of the element, length
          `numNodesPerElement`

  **Outputs**

   * is_shock: bool
"""
function isShockElement(params::AbstractParamType{Tdim}, sbp::AbstractOperator,
                        sensor::AbstractShockSensor,
                        q::AbstractMatrix{Tsol}, coords::AbstractMatrix,
                        dxidx::Abstract3DArray, jac::AbstractVector{Tmsh},
                        ) where {Tsol, Tmsh, Tdim}

  Tres = promote_type(Tsol, Tmsh)
  numDofPerNode, numNodesPerElement = size(q)
  dim = size(dxidx, 1)
  Se = EmptyArray{Tres}(dim, numNodesPerElement)
  ee = EmptyArray{Tres}(dim, numNodesPerElement)

  return getShockSensor(params, sbp, sensor, q, coords, dxidx, jac, Se, ee)
end

"""
  Computes the shock sensor value and the associated viscoscity.  The
  viscoscity can be different in each Cartesian and at each node of the
  element,  A full diffusion tensor is not supported.

  **Inputs**

   * params: a `AbstractParamType`
   * sbp: SBP operator
   * sensor: the `AbstractShockSensor` object
   * q: solution on a particular element, `numDofPerNode` x `numNodesPerElement`
   * coords: `dim` x `numNodesPerElement` array containing the coordinates
             of each node of the element
   * dxidx: `dim` x `dim` x `numNodesPerElement` array containing the (scaled)
            mapping jacobian dxi/dx
   * jac: mapping jacobian determinant for each node of the element, length
          `numNodesPerElement`

  **Inputs/Outputs**

   * Se: `dim` x `numNodesPerElement` containing the shock sensor in each
         direction at each node.  
   * ee: the viscoscity coefficient, same size as `Se`

  **Outputs**

   * is_nonzero: Bool, true if ee is zero at all nodes, false otherwise

"""
function getShockSensor(params::AbstractParamType{Tdim}, sbp::AbstractOperator,
                          sensor::AbstractShockSensor,
                          q::AbstractMatrix{Tsol}, coords::AbstractMatrix,
                          dxidx::Abstract3DArray, jac::AbstractVector{Tmsh},
                          Se::AbstractMatrix, ee::AbstractMatrix
                         ) where {Tsol, Tmsh, Tdim}

  error("abstract method for getShockSensor() called: did you forget to extend it with a new method for shock sensor $(typeof(sensor))")

end


"""
  Differentiated version of [`getShockSensor`](@ref).  Most of the inputs
  and outputs are the same as that function, see its documentation.

  **Inputs**

   * params
   * sbp
   * q
   * coords
   * dxidx
   * jac

  **Inputs/Outputs**

   * Se_jac: `dim` x `numDofPerNode` x `numNodesPerElement` x 
             `numNodesPerElement` array to be overwritten with jacobian of
             `Se` wrt `q`
   * ee_jac: array to be overwritten with jacobian of `ee` wrt `q`, same size as
             `Se_jac`

  **Outputs**
  
   * is_const: true if the derivative of ee wrt q is zero.  This is useful
               for skipping terms in the jacobian calculation

"""
function getShockSensor_diff(params::AbstractParamType{Tdim}, sbp::AbstractOperator,
                      sensor::AbstractShockSensor,
                      q::AbstractMatrix{Tsol},
                      coords::AbstractMatrix, dxidx::Abstract3DArray,
                      jac::AbstractVector{Tmsh},
                      Se_jac::Abstract4DArray{Tres},
                      ee_jac::Abstract4DArray{Tres}
                     ) where {Tsol, Tmsh, Tres, Tdim}

  error("abstract method for getShockSensor_diff() called: did you forget to extend it with a new method for shock sensor $(typeof(sensor))")

end
