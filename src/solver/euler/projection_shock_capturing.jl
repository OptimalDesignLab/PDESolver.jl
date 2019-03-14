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

    projectionShockCapturing(eqn.params, sbp, sensor, capture, q_i, coords_i,
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


function projectionShockCapturing(params::ParamType, sbp::AbstractOperator,
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


#------------------------------------------------------------------------------
# Differentiated version

function calcShockCapturing_diff(mesh::AbstractMesh, sbp::AbstractOperator,
                                 eqn::EulerData, opts,
                                 sensor::AbstractShockSensor,
                                 capture::ProjectionShockCapturing,
                                 assem::AssembleElementData)


  data = eqn.params.calc_volume_integrals_data
  res_jac = data.res_jac
  fill!(res_jac, 0)

  for i=1:mesh.numEl
    q_i = ro_sview(eqn.q, :, :, i)
    coords_i = ro_sview(mesh.coords, :, :, i)
    dxidx_i = ro_sview(mesh.dxidx, :, :, :, i)
    jac_i = ro_sview(mesh.jac, :, i)

    nonzero_jac = projectionShockCapturing_diff(eqn.params, sbp, sensor,
                              capture, q_i, coords_i, dxidx_i, jac_i, res_jac)
  
    # assembling into a sparse matrix is non-trivially expensive, don't do
    # it unless this element has shock capturing active
    if nonzero_jac
      if eqn.params.use_Minv == 1
        applyMinvElement(jac_i, sbp.w, res_jac)
      end

      # assemble element level jacobian into the residual
      assembleElement(assem, mesh, i, res_jac)
      fill!(res_jac, 0)
    end  # if nonzero_jac
  
  end  # end i

  return nothing
end

"""
  Differentiated version of `projectionShockCapturing` for
  [`ProjectionShockCapturing`](@ref).
"""
function projectionShockCapturing_diff(params::ParamType, sbp::AbstractOperator,
                                  sensor::AbstractShockSensor,
                                  capture::ProjectionShockCapturing,
                                  u::AbstractMatrix, 
                                  coords::AbstractMatrix,
                                  dxidx::Abstract3DArray,
                                  jac::AbstractVector{Tmsh},
                                  res_jac::AbstractArray{Tres, 4}) where {Tmsh, Tres}

  numDofPerNode, numNodesPerElement = size(u)
  dim = size(coords, 1)
  @unpack capture t1 t2 w Se_jac ee_jac A0inv

  #TODO: stash these somewhere
  Se = zeros(Tres, dim, numNodesPerElement)
  ee = zeros(Tres, dim, numNodesPerElement)
  #TODO: make shock capturing and shock sensing independent choices
  is_shock = isShockElement(params, sbp, sensor, u, coords, dxidx, jac)

  if is_shock
    fill!(Se_jac, 0); fill!(ee_jac, 0)
    # only compute the derivative if there is a shock
    getShockSensor(params, sbp, sensor, u, coords, dxidx, jac, Se, ee)
    ee_constant = getShockSensor_diff(params, sbp, sensor, u, coords,
                                              dxidx, jac, Se_jac, ee_jac)

    # the operator (for a scalar equation) is A = P^T * M * P * v, so
    # dR[p]/v[q] = (P^T * M * P)[p, q].  It then needs to be converted back
    # to conservative variables

    # compute derivative contribution from v
    @simd for p=1:numNodesPerElement
      @simd for q=1:numNodesPerElement
        getIRA0inv(params, sview(u, :, q), A0inv)
        # calculate the A[p, q]
        Apq = zero(Tres)
        @simd for k=1:numNodesPerElement
          Apq += capture.filt[k, p]*(sbp.w[k]/jac[k])*capture.filt[k, q]
        end
        @simd for j=1:numDofPerNode
          @simd for i=1:numDofPerNode
            # this uses only the x direction viscoscity
            res_jac[i, j, p, q] = -ee[1, p]*Apq*A0inv[i, j]
          end
        end
      end
    end

    # compute derivative contribution from ee
    if !ee_constant

      @simd for i=1:numNodesPerElement
        w_i = sview(w, :, i)
        q_i = sview(u, :, i)
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

      @simd for p=1:numNodesPerElement
        @simd for q=1:numNodesPerElement
          @simd for j=1:numDofPerNode
            @simd for i=1:numDofPerNode
              res_jac[i, j, p, q] -= ee_jac[1, j, p, q]*t2[i, p]
            end
          end
        end
      end

    end  # if ee_constant


  end  # end if

  return is_shock
end


