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
"""
abstract type AbstractShockSensor end


"""
  Computes the shock sensor value and the associated viscoscity

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

   * Se: the shock sensor value
   * ee: the viscoscity coefficient (constant for the entire element)

"""
function getShockSensor(params::AbstractParamType{Tdim}, sbp::AbstractOperator,
                          sensor::AbstractShockSensor,
                          q::AbstractMatrix{Tsol}, coords::AbstractMatrix,
                          dxidx::Abstract3DArray, jac::AbstractVector{Tmsh},
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

   * Se_jac: `numDofPerNode` x `numNodesPerElement` array to be overwritten
            with jacobian of `Se` wrt `q`
   * ee_jac: `numDofPerNode` x `numNodesPerElement` array to be overwritten
             with jacobian of `ee` wrt `q`

  **Outputs**
  
   * Se
   * ee
   * is_const: true if the derivative of ee wrt q is zero.  This is useful
               for skipping terms in the jacobian calculation

"""
function getShockSensor_diff(params::AbstractParamType{Tdim}, sbp::AbstractOperator,
                      sensor::AbstractShockSensor,
                      q::AbstractMatrix{Tsol},
                      coords::AbstractMatrix, dxidx::Abstract3DArray,
                      jac::AbstractVector{Tmsh},
                      Se_jac::AbstractMatrix{Tres},
                      ee_jac::AbstractMatrix{Tres}
                     ) where {Tsol, Tmsh, Tres, Tdim}

  error("abstract method for getShockSensor_diff() called: did you forget to extend it with a new method for shock sensor $(typeof(sensor))")

end
