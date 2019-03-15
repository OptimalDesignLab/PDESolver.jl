# differentiated shock capturing functions

# main entry point
"""
  Main function for assembling the shock capturing terms into the Jacobian.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * capture: an [`AbstractShockCapturing`](@ref)
   * assem: an [`AssembleElementData`](@ref)
"""
function applyShockCapturing_diff(mesh::AbstractMesh, sbp::AbstractOperator,
                             eqn::EulerData, opts,
                             capture::AbstractVolumeShockCapturing,
                             assem::AssembleElementData)

  # nothing to do here, call the implementation
  sensor = getShockSensor(capture)
  calcShockCapturing_diff(mesh, sbp, eqn, opts, sensor, capture, assem)

  return nothing
end


function applyShockCapturing_diff(mesh::AbstractMesh, sbp::AbstractOperator,
                             eqn::EulerData{Tsol, Tres}, opts,
                             capture::AbstractFaceShockCapturing,
                             assem::AssembleElementData) where {Tsol, Tres}

  sensor = getShockSensor(capture)
  _applyShockCapturing_diff(mesh, sbp, eqn, opts, sensor, capture, assem)

  return nothing
end

function _applyShockCapturing_diff(mesh::AbstractMesh, sbp::AbstractOperator,
                             eqn::EulerData{Tsol, Tres}, opts,
                             sensor::AbstractShockSensor,
                             capture::AbstractFaceShockCapturing,
                             assem::AssembleElementData) where {Tsol, Tres}

  if mesh.commsize > 1 && 
    getParallelData(eqn.shared_data) != PARALLEL_DATA_ELEMENT

    error("shock capturing requires PARALLEL_DATA_ELEMENT")
  end

  assertReceivesWaited(eqn.shared_data)  # parallel communication for the regular
                                     # face integrals should already have
                                     # finished parallel communication


  # re-using the shockmesh from one iteration to the next makes allows
  # the algorithm to re-use existing arrays that depend on the number of
  # elements with the shock in it.
  shockmesh = eqn.params.shockmesh
  reset(shockmesh)
  for i=1:mesh.numEl
    q_i = ro_sview(eqn.q, :, :, i)
    coords_i = ro_sview(mesh.coords, :, :, i)
    dxidx_i = ro_sview(mesh.dxidx, :, :, :, i)
    jac_i = ro_sview(mesh.jac, :, i)

    val = isShockElement(eqn.params, sbp, sensor, q_i, coords_i, dxidx_i,
                            jac_i)
    if val
      # push to shockmesh
      push!(shockmesh, i)
    end
  end

  completeShockElements(mesh, shockmesh)

  allocateArrays(capture, mesh, shockmesh)
  
  # call shock capturing scheme
  calcShockCapturing_diff(mesh, sbp, eqn, opts, capture, shockmesh, assem)

  return nothing
end



#------------------------------------------------------------------------------




