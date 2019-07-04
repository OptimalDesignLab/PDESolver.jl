# shock capturning functions


"""
  Main entry point for shock capturing.  Applies the specified shock capturing
  scheme to every element

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * data: an [`AbstractShockCaputring`](@ref) object
"""
function applyShockCapturing(mesh::AbstractMesh, sbp::AbstractOperator,
                             eqn::EulerData, opts,
                             capture::AbstractVolumeShockCapturing)

  # unlike AbstractFaceShockCapturing, nothing to do here

  sensor = getShockSensor(capture)
  calcShockCapturing(mesh, sbp, eqn, opts, sensor, capture)

  return nothing
end

function applyShockCapturing(mesh::AbstractMesh, sbp::AbstractOperator,
                             eqn::EulerData, opts,
                             capture::AbstractFaceShockCapturing)

  sensor = getShockSensor(capture)
  shockmesh = eqn.params.shockmesh
  setupShockCapturing(mesh, sbp, eqn, opts, sensor, capture, shockmesh)

  # call shock capturing scheme
  calcShockCapturing(mesh, sbp, eqn, opts, capture, shockmesh)


  return nothing
end

"""
  This function sets up the `ShockedElements` object and calls `allocateArrays`
  on the `AbstractFaceShockCapturing` object afterwards.  This completes
  the preparations required before the shock capturing term can be calculated

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * sensor
   * capture
"""
@noinline function setupShockCapturing(mesh::AbstractMesh, sbp::AbstractOperator,
                             eqn::EulerData, opts,
                             sensor::AbstractShockSensor,
                             capture::AbstractFaceShockCapturing,
                             shockmesh::ShockedElements)

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
  reset(shockmesh)
  numEl_shock = 0
  for i=1:mesh.numEl
    q_i = ro_sview(eqn.q, :, :, i)
    coords_i = ro_sview(mesh.coords, :, :, i)
    dxidx_i = ro_sview(mesh.dxidx, :, :, :, i)
    jac_i = ro_sview(mesh.jac, :, i)

    val = isShockElement(eqn.params, sbp, sensor, q_i, i, coords_i, dxidx_i,
                            jac_i)
    if val
      # push to shockmesh
      #println("element ", i, " viscoscity = ", ee)
      numEl_shock += 1
      push!(shockmesh, i)
    end
  end

  #println("counted ", numEl_shock, " elements with nonzero viscoscity")
  completeShockElements(mesh, shockmesh)
  @assert shockmesh.numShock == numEl_shock

  allocateArrays(capture, mesh, shockmesh)

  return nothing
end



#------------------------------------------------------------------------------
# Debugging/testing

using PumiInterface
using apf

function writeShockSensorField(mesh, sbp, eqn::EulerData{Tsol, Tres}, opts,
                               sensor::AbstractShockSensor) where {Tsol, Tres}

  fname = "shock_sensor"

  f_ptr = apf.findField(mesh.m_ptr, fname)
  if f_ptr == C_NULL
    fshape = apf.getConstantShapePtr(mesh.dim)
    f_ptr = apf.createPackedField(mesh, fname, 2, fshape)
  end

  vals = zeros(Float64, 2)
  numEl_shock = 0

  Se = zeros(Tres, mesh.dim, mesh.numNodesPerElement)
  ee = zeros(Tres, mesh.dim, mesh.numNodesPerElement)
  for i=1:mesh.numEl
    q_i = sview(eqn.q, :, :, i)
    coords_i = ro_sview(mesh.coords, :, :, i)
    dxidx_i = ro_sview(mesh.dxidx, :, :, :, i)
    jac_i = sview(mesh.jac, :, i)

    is_nonzero = getShockSensor(eqn.params, sbp, sensor, q_i, i, coords_i, dxidx_i, 
                            jac_i, Se, ee)


    #val1 = computeL2Norm(eqn.params, sbp, jac_i, Se)
    #val2 = computeL2Norm(eqn.params, sbp, jac_i, ee)
    val1 = zero(Tres); val2 = zero(Tres)
    for j=1:mesh.numNodesPerElement
      for d=1:mesh.dim
        val1 += Se[d, j]
        val2 += ee[d, j]
      end
    end

    val1 /= mesh.dim*mesh.numNodesPerElement
    val2 /= mesh.dim*mesh.numNodesPerElement

    vals[1] = val1
    vals[2] = val2
    apf.setComponents(f_ptr, mesh.elements[i], 0, vals)
    if is_nonzero
      numEl_shock += 1
    end
  end

  numEl_shock = MPI.Allreduce(numEl_shock, MPI.SUM, eqn.comm)

  if mesh.myrank == 0
    println("wrote shock field with ", numEl_shock, " elements with non-zero viscoscity")
  end

  return nothing
end



