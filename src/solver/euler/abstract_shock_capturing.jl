
"""
  Abstract type for all shock capturing methods.  This type determines what
  kind of dissipation to add when there is a shock present.  Ideally, any
  shock capturing scheme should be able to be paired with any shock sensor.

  New shock capturing scheme should generally not subtype this type directly,
  instead they should subtype either [`AbstractVolumeShockCapturing`](@ref)
  or [`AbstractFaceShockCapturing`](@ref).

  Optional methods:

  ```
    getShockSensor
    getShockSensor
    getDiffusion
  ```

  `getShockSensor` and `setShockSensor only need to be extended if the
  `AbstractShockCapturing` does not store the sensor in the `.sensor` field.
  
  `getDiffusion` only needs to be extended if the `AbstractShockCapturing` does
  not store the sensor in the `.diffusion` field.

"""
abstract type AbstractShockCapturing end


"""
  Returns the shock sensor.  This function is not required to be type-stable

  **Inputs**

   * obj: an `AbstractShockCapturing`

  **Outputs**

   * sensor: an `AbstractShockSensor`
"""
function getShockSensor(obj::AbstractShockCapturing)

  return obj.sensor
end

"""
  Returns the `AbstractDiffusion` object.  This function is not required to
  be type stable.

  **Inputs**

   * obj: an `AbstractShockCapturing`

  **Outputs**

   * diffusion: an `AbstractDiffusion`
"""
function getDiffusion(obj::AbstractShockCapturing)

  return obj.diffusion
end

"""
  Counterpart of `getShockSensor`, sets a new shock sensor to be used with
  the existing `AbstractShockCapturing`.

  **Inputs**

   * obj: `AbstractShockCapturing`
   * sensor: the new `AbstractShockSensor`
"""
function setShockSensor(obj::AbstractShockCapturing, sensor::AbstractShockSensor)
  obj.sensor = sensor
end

#------------------------------------------------------------------------------
# Volume shock capturing

"""
  Abstract type for shock capturing methods that only have volume terms.
  These scheme are substantially simpler and do not require constructing
  a reduced mesh.  The interface for this type has not been finalized

  This type must implement functions:

  ```
    calcShockCapturing
    calcShockCapturing_diff
  ```
"""
abstract type AbstractVolumeShockCapturing <: AbstractShockCapturing end


"""
  This function should update the residual with the shock capturing term

  **Inputs**

   * mesh: mesh object
   * sbp: SBP operator
   * eqn: EulerData
   * opts: options dictonary
   * sensor: [`AbstractShockSensor`](@ref)
   * capture: the shock capturing scheme
"""
function calcShockCapturing(mesh::AbstractMesh, sbp::AbstractOperator,
                            eqn::EulerData, opts,
                            sensor::AbstractShockSensor,
                            capture::AbstractVolumeShockCapturing)

  error("abstract fallback for calcShockCapturing() called. Did you forget to extend it with a new method for your AbstractFaceShockCapturing?")

end


"""
  This function should compute the Jacobian of the shock capturing term

  **Inputs**

   * mesh: mesh object
   * sbp: SBP operator
   * eqn: EulerData
   * opts: options dictionary
   * sensor: an [`AbstractShockSensor`](@ref)
   * capture: the shock capturing object
   * assem: an [`AssembleElementData`](@ref)

  Note that many `capture` types require the same sensor to be passed into
  their constructor as to `calcShockCapturing`.

"""
function calcShockCapturing_diff(mesh::AbstractMesh, sbp::AbstractOperator,
                                 eqn::EulerData, opts,
                                 sensor::AbstractShockSensor,
                                 capture::AbstractVolumeShockCapturing,
                                 assem::AssembleElementData)

  error("abstract fallback for calcShockCapturing_diff() called. Did you forget to extend it with a new method for your AbstractFaceShockCapturing?")
end


#------------------------------------------------------------------------------
# AbstractFaceShockCapturing

"""
  Abstract type for shock capturing methods that do face integrals (these
  scheme may do volume integrals as well).  These schemes require constructing
  a reduced mesh for the elements that have shocks in them.

  These types must implement these functions:
  ```
    allocateArrays
    calcShockCapturing
    calcShockCapturing_diff
  ```

"""
abstract type AbstractFaceShockCapturing <: AbstractShockCapturing end


"""
  This function is called every time the `shockmesh` is updated, and performs
  any setup work that depends on the `shockmesh`


  **Inputs**

   * capture: [`LDGShockCapturing`](@ref)
   * mesh
   * shockmesh: a `ShockedElements` object, fully initialized

"""
function allocateArrays(capture::AbstractFaceShockCapturing,
                        mesh::AbstractMesh, shockmesh::ShockedElements)

  error("abstract fallback for allocateArrays() called. Did you forget to extend it with a new method for your AbstractFaceShockCapturing?")
end


"""
  Update `eqn.res` with the shock capturing term

  **Inputs**

   * mesh
   * sbp
   * eqn: `eqn.res` should be updated
   * opts
   * capture: the shock capturing scheme to be used (this argument must
              be specialized)
   * shockmesh

"""
function calcShockCapturing(mesh::AbstractMesh, sbp::AbstractOperator,
                           eqn::EulerData, opts,
                           capture::AbstractFaceShockCapturing,
                           shockmesh::ShockedElements)

  error("abstract fallback for calcShockCapturing() called. Did you forget to extend it with a new method for your AbstractFaceShockCapturing?")
end

"""
  Assemble the Jacobian of the shock capturing term into the Jacobian

  **Inputs**

   * mesh
   * sbp
   * eqn: 
   * opts
   * capture: the shock capturing scheme to be used (this argument must
              be specialized)
   * shockmesh: a `ShockedElements`
   * assem:  an `AssembleElementData`

"""
function calcShockCapturing_diff(mesh::AbstractMesh, sbp::AbstractOperator,
                             eqn::EulerData, opts,
                             capture::AbstractFaceShockCapturing,
                             shockmesh::ShockedElements,
                             assem::AssembleElementData)

  error("abstract fallback for calcShockCapturing_diff() called. Did you forget to extend it with a new method for your AbstractFaceShockCapturing?")
end
