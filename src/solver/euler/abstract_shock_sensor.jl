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
  updateMetrics
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
   * elnum: element number
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
                        q::AbstractMatrix{Tsol}, elnum::Integer,
                        coords::AbstractMatrix,
                        dxidx::Abstract3DArray, jac::AbstractVector{Tmsh},
                        ) where {Tsol, Tmsh, Tdim}

  Tres = promote_type(Tsol, Tmsh)
  numDofPerNode, numNodesPerElement = size(q)
  dim = size(dxidx, 1)
  Se = EmptyArray{Tres}(dim, numNodesPerElement)
  ee = EmptyArray{Tres}(dim, numNodesPerElement)

  return getShockSensor(params, sbp, sensor, q, elnum, coords, dxidx, jac, Se,
                        ee)
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
   * elnum: element number
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
                          q::AbstractMatrix{Tsol}, elnum::Integer,
                          coords::AbstractMatrix,
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
   * elnum
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
                      q::AbstractMatrix{Tsol}, elnum::Integer,
                      coords::AbstractMatrix, dxidx::Abstract3DArray,
                      jac::AbstractVector{Tmsh},
                      Se_jac::Abstract4DArray{Tres},
                      ee_jac::Abstract4DArray{Tres}
                     ) where {Tsol, Tmsh, Tres, Tdim}

  error("abstract method for getShockSensor_diff() called: did you forget to extend it with a new method for shock sensor $(typeof(sensor))")

end


"""
  Function to be called after the mesh metrics are updated.  This gives the
 eshock sensor the opportunity to recalculate mesh-wide quantities.

  **Inputs**

   * mesh
   * sbp
   * opts
   * sensor: an [`AbstractShockSensor`](@ref)
"""
function updateMetrics(mesh, sbp, opts, sensor::AbstractShockSensor)

end

"""
  Reverse mode with respect to `q` of the shock sensor

  **Inputs**

   * params
   * sbp
   * sensor: `AbstractShockSensor`
   * q
   * elnum
   * coords
   * dxidx
   * jac
   * ee_bar: adjoint part of `ee`, same size

  **Inputs/Outputs**
   
   * ee
   * q_bar: to be updated with result, same size as `q`

  **Outputs**

   * is_nonzero: Bool, true if ee is zero at all nodes, false otherwise
"""
function getShockSensor_revq(params::AbstractParamType, sbp::AbstractOperator,
                        sensor::AbstractShockSensor,
                        q::AbstractMatrix, q_bar::AbstractMatrix,
                        elnum::Integer,
                        coords::AbstractMatrix,
                        dxidx::Abstract3DArray, jac::AbstractVector,
                        ee_mat::AbstractMatrix, ee_bar::AbstractMatrix
                        )

  error("abstract method for getShockSensor_revq() called: did you forget to extend it with a new method for shock sensor $(typeof(sensor))")

end


"""
  Reverse mode wrt metrics of the shock sensor

  **Inputs**

   * params
   * sbp
   * sensor
   * q
   * elnum
   * coords
   * dxidx
   * jac
   * ee_bar: seed matrix for back-propigation

  **Inputs/Outputs**

   * coords_bar
   * dxidx_bar
   * jac_bar
"""
function getShockSensor_revm(params::AbstractParamType, sbp::AbstractOperator,
                        sensor::AbstractShockSensor,
                        q::AbstractMatrix,
                        elnum::Integer,
                        coords::AbstractMatrix, coords_bar::AbstractMatrix,
                        dxidx::Abstract3DArray, dxidx_bar::Abstract3DArray,
                        jac::AbstractVector, jac_bar::AbstractVector,
                        ee_bar::AbstractMatrix
                       )

  error("abstract method for getShockSensor_revm() called: did you forget to extend it with a new method for shock sensor $(typeof(sensor))")

end


"""
  Performs any initialization required by the shock sensor at the beginning of
  a reverse mode wrt metrics calculation.  This is mainly used for zeroing out
  arrays inside the shock sensor that will be += during the computation.

  This function has an (empty) default function, so shock sensors only need to
  implement this function if they have to do setup work.

  **Inputs**

   * sensor: an `AbstractShockSensor`
"""
function initForRevm(sensor::AbstractShockSensor)

  return nothing
end

"""
  Finish back propigation with respect to the metrics by back propgating any
  array in the `sensor` to the `mesh`.

  This function has an (empty) default function, so shock sensors only need to
  implement this function if they have work to do.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * sensor: `AbstractShockSensor`
"""
function finishRevm(mesh::AbstractMesh, sbp::AbstractOperator, eqn::EulerData,
                    opts, sensor::AbstractShockSensor)

  return nothing
end
