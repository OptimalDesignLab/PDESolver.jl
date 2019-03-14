# Shock capturing using Prof. Hicken's (as yet unnamed) projection method,
#
function calcShockCapturing(mesh::AbstractMesh, sbp::AbstractOperator,
                             eqn::EulerData, opts,
                             sensor::AbstractShockSensor,
                             capture::ProjectionShockCapturing)

  for i=1:mesh.numEl
    q_i = ro_sview(eqn.q, :, :, i)
    coords_i = ro_sview(mesh.coords, :, :, i)
    dxidx_i = ro_sview(mesh.dxidx, :, :, :, i)
    jac_i = ro_sview(mesh.jac, :, i)
    res_i = sview(eqn.res, :, :, i)

    calcShockCapturing(eqn.params, sbp, sensor, capture, q_i, coords_i,
                       dxidx_i, jac_i, res_i)
  end

  return nothing
end

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
                             q::AbstractMatrix, coords::AbstractMatrix,
                             dxidx::Abstract3DArray, jac::AbstractVector,
                             res::AbstractMatrix{Tres}) where {Tres}

  numDofPerNode, numNodesPerElement = size(q)
  dim = size(dxidx, 1)
  @unpack capture t1 t2 w

  #TODO: stash these somewhere
  Se = zeros(Tres, dim, numNodesPerElement)
  ee = zeros(Tres, dim, numNodesPerElement)
  is_nonconst = getShockSensor(params, sbp, sensor, q, coords, dxidx, jac, Se,
                               ee)
  if is_nonconst # if there is a shock
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
        # this uses only the x direction viscoscity
        res[j, i] -= ee[1, i]*t2[j, i]
      end
    end

  end  # end if

  return nothing
end


