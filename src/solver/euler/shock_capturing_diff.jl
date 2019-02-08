# differentiated shock capturing functions

# main entry point
"""
  Main function for assembling the shock capturing terms into the Jacobian.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * sensor: an [`AbstractShockSensor`](@ref)
   * capture: an [`AbstractShockCapturing`](@ref)
   * assem: an [`AssembleElementData`](@ref)
"""
function applyShockCapturing_diff(mesh::AbstractMesh, sbp::AbstractOperator,
                             eqn::EulerData, opts,
                             sensor::AbstractShockSensor,
                             capture::AbstractVolumeShockCapturing,
                             assem::AssembleElementData)


  data = eqn.params.calc_volume_integrals_data
  res_jac = data.res_jac
  fill!(res_jac, 0)

  for i=1:mesh.numEl
    q_i = sview(eqn.q, :, :, i)
    jac_i = sview(mesh.jac, :, i)

    nonzero_jac = calcShockCapturing_diff(eqn.params, sbp, sensor, capture,
                                            q_i, jac_i, res_jac)
  
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



function applyShockCapturing_diff(mesh::AbstractMesh, sbp::AbstractOperator,
                             eqn::EulerData, opts,
                             sensor::AbstractShockSensor,
                             capture::AbstractFaceShockCapturing,
                             assem::AssembleElementData)

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
    jac_i = ro_sview(mesh.jac, :, i)

    Se, ee = getShockSensor(eqn.params, sbp, sensor, q_i, jac_i)
    if ee > 0
      # push to shockmesh
      push!(shockmesh, i, ee)
    end
  end

  completeShockElements(mesh, shockmesh)

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
  
  # compute the derivative of the shock sensor for the shocked elements
  ee_dot = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, shockmesh.numEl)
  Se_dot_i = zeros(Tres, mesh.numDofPerNOde, mesh.numNodesPerElement)
  is_nonlinear = BitArray(shockmesh.numEl)
  for i=1:shockmesh.numEl
    i_full = shockmesh.elnums_all[i]
    q_i = ro_sview(eqn.q, :, :, i_full)
    jac_i = ro_sview(mesh.jac, :, :, i_full)
    ee_dot_i = sivew(ee_dot, :, :, i)

    Se, ee, is_constant = getShockSensor_diff(eqn.params, sbp, sensor, q_i,
                                              jac_i, Se_dot_i, ee_dot_i)
    is_nonlinear[i] = !is_constant
  end

  #TODO: perhaps remove ee from shockmesh?
  setDiffusionArray_diff(capture.diffusion, shockmesh.ee, ee_dot, is_nonlinear)

  # call shock capturing scheme
  calcShockCapturing_diff(mesh, sbp, eqn, opts, capture, shockmesh)

  return nothing
end



#------------------------------------------------------------------------------


"""
  Differentiated version of `calcShockCapturing` for
  [`ProjectionShockCapturing`](@ref).
"""
function calcShockCapturing_diff(params::ParamType, sbp::AbstractOperator,
                                  sensor::AbstractShockSensor,
                                  capture::ProjectionShockCapturing,
                                  u::AbstractMatrix, jac::AbstractVector{Tmsh},
                                  res_jac::AbstractArray{Tres, 4}) where {Tmsh, Tres}

  numDofPerNode, numNodesPerElement = size(u)
  @unpack capture t1 t2 w Se_jac ee_jac A0inv

  #TODO: make shock capturing and shock sensing independent choices
  Se, ee = getShockSensor(params, sbp, sensor, u, jac)

  if ee > 0
    fill!(Se_jac, 0); fill!(ee_jac, 0)
    # only compute the derivative if there is a shock
    Se, ee, ee_constant = getShockSensor_diff(params, sbp, sensor, u, jac,
                                              Se_jac, ee_jac)

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
            res_jac[i, j, p, q] = -ee*Apq*A0inv[i, j]
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
              res_jac[i, j, p, q] -= ee_jac[j, q]*t2[i, p]
            end
          end
        end
      end

    end  # if ee_constant


  end  # end if

  return ee > 0
end



