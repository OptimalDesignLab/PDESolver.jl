# shock capturning functions


"""
  Main entry point for shock capturing.  Applies the specified shock capturing
  scheme to every element

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * sensor: an [`AbstractShockSensor`](@ref) object
   * data: an [`AbstractShockCaputring`](@ref) object
"""
function applyShockCapturing(mesh::AbstractMesh, sbp::AbstractOperator,
                             eqn::EulerData, opts,
                             sensor::AbstractShockSensor,
                             capture::AbstractVolumeShockCapturing)

  for i=1:mesh.numEl
    q_i = ro_sview(eqn.q, :, :, i)
    jac_i = ro_sview(mesh.jac, :, i)
    res_i = sview(eqn.res, :, :, i)

    calcShockCapturing(eqn.params, sbp, sensor, capture, q_i, jac_i, res_i)
  end

  return nothing
end


function applyShockCapturing(mesh::AbstractMesh, sbp::AbstractOperator,
                             eqn::EulerData, opts,
                             sensor::AbstractShockSensor,
                             capture::AbstractFaceShockCapturing)

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
  numEl_shock = 0
  for i=1:mesh.numEl
    q_i = ro_sview(eqn.q, :, :, i)
    jac_i = ro_sview(mesh.jac, :, i)

    Se, ee = getShockSensor(eqn.params, sbp, sensor, q_i, jac_i)
    if ee > 0
      # push to shockmesh
      println("element ", i, " viscoscity = ", ee)
      numEl_shock += 1
      push!(shockmesh, i, ee)
    end
  end

  println("counted ", numEl_shock, " elements with nonzero viscoscity")
  completeShockElements(mesh, shockmesh)
  @assert shockmesh.numShock == numEl_shock

  # compute the shock viscoscity for shared elements
  for peer=1:shockmesh.npeers
    peer_full = shockmesh.peer_indices[peer]
    data = eqn.shared_data[peer_full]
    metrics = mesh.remote_metrics[peer_full]

    for i in shockmesh.shared_els[peer]
      i_full = getSharedElementIndex(shockmesh, mesh, peer, i)
      q_i = ro_sview(data.q_recv, :, :, i_full)
      jac_i = ro_sview(metrics.jac, :, i_full)

      Se, ee = getShockSensor(eqn.params, sbp, sensor, q_i, jac_i)
      setViscoscity(shockmesh, i, ee)
    end
  end


  allocateArrays(capture, mesh, shockmesh)

  # call shock capturing scheme
  calcShockCapturing(mesh, sbp, eqn, opts, capture, shockmesh)

  return nothing
end


#------------------------------------------------------------------------------
# Shock capturing using Prof. Hicken's (as yet unnamed) projection method,
# using the PP shock sensor and viscoscity coefficient

"""
Returns a nxn matrix that projects out modes lower than sbp.degree
"""
function getFilterOperator!(sbp::TriSBP{T}, diss::AbstractArray{T,2}) where {T}
  x = calcnodes(sbp)
  # loop over ortho polys up to degree d
  #println("TEMP: filter operator set to sbp.degree-1!!!!!")
  d = 0
#  d = sbp.degree
  V = zeros(T, (sbp.numnodes, convert(Int, (d+1)*(d+2)/2)) )
  ptr = 0
  for r = 0:d
    for j = 0:r
      i = r-j
      V[:,ptr+1] = SummationByParts.OrthoPoly.proriolpoly(vec(x[1,:]), vec(x[2,:]), i, j)
      ptr += 1
    end
  end
  #F = eye(sbp.numnodes) - (pinv(V.')*V.').'
  #Minv = diagm(1./diag(V.'*diagm(sbp.w)*V))
  # Minv is the idenity, since the V's are orthogonal in L2
  diss[:,:] = eye(sbp.numnodes) - (V*V.')*diagm(sbp.w)
end

#TODO: this doesn't work yet
function getFilterOperator!(sbp::TetSBP{T}, diss::AbstractArray{T,2}) where {T}

  fill!(diss, 0)
end


function calcShockCapturing(params::ParamType, sbp::AbstractOperator,
                             sensor::AbstractShockSensor,
                             capture::ProjectionShockCapturing,
                             q::AbstractMatrix, jac::AbstractVector,
                             res::AbstractMatrix)

  numDofPerNode, numNodesPerElement = size(q)
  @unpack capture t1 t2 w

  Se, ee = getShockSensor(params, sbp, sensor, q, jac)
  if ee > 0  # if there is a shock
    # For scalar equations, the operator is applied -epsilon * P^T M P * u
    # For vector equations, P needs to be applied to all equations as once:
    # utmp^T = P*u^T
    # Instead, work with utmp = (P*u^T)^T = u*P^T

    # convert to entropy variables to make this term entropy-stable
    #w = zeros(eltype(q), numDofPerNode, numNodesPerElement)
    for i=1:numNodesPerElement
      w_i = sview(w, :, i)
      q_i = sview(q, :, i)
      convertToIR(params, q_i, w_i)
    end

    # apply P
    smallmatmatT!(w, capture.filt, t1)


    # apply mass matrix
    @simd for i=1:numNodesPerElement
      fac = sbp.w[i]/jac[i]
      @simd for j=1:numDofPerNode
        t1[j, i] *= fac
      end
    end

    # apply P^T
    smallmatmat!(t1, capture.filt, t2)
    @simd for i=1:numNodesPerElement
      @simd for j=1:numDofPerNode
        res[j, i] -= ee*t2[j, i]
      end
    end

  end  # end if

  return nothing
end


#------------------------------------------------------------------------------
# Debugging/testing

using PumiInterface
using apf

function writeShockSensorField(mesh, sbp, eqn, opts, sensor::AbstractShockSensor)

  fname = "shock_sensor"

  f_ptr = apf.findField(mesh.m_ptr, fname)
  if f_ptr == C_NULL
    fshape = apf.getConstantShapePtr(mesh.dim)
    f_ptr = apf.createPackedField(mesh.m_ptr, fname, 2, fshape)
  end

  vals = zeros(Float64, 2)
  numEl_shock = 0
  for i=1:mesh.numEl
    q_i = sview(eqn.q, :, :, i)
    jac_i = sview(mesh.jac, :, i)

    Se, ee = getShockSensor(eqn.params, sbp, sensor, q_i, jac_i)
    vals[1] = Se
    vals[2] = ee
    apf.setComponents(f_ptr, mesh.elements[i], 0, vals)
    if ee > 0
      numEl_shock += 1
    end
  end

  println("wrote shock field with ", numEl_shock, " elements with non-zero viscoscity")

  return nothing
end



