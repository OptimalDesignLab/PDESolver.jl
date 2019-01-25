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

    applyShockCapturing(eqn.params, sbp, sensor, capture, q_i, jac_i, res_i)
  end

  return nothing
end


function applyShockCapturing(mesh::AbstractMesh, sbp::AbstractOperator,
                             eqn::EulerData, opts,
                             sensor::AbstractShockSensor,
                             capture::AbstractFaceShockCapturing)

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

  # call LDG
  applyShockCapturing(mesh, sbp, eqn, opts, capture, shockmesh)

  return nothing
end


#------------------------------------------------------------------------------
# Method from: Persson and Peraire, "Sub-Cell Shock Capturing for DG Methods"
#              AIAA 2006
# Some additional details are from Barter's Ph.D Thesis, "Shock Capturing
# with PDE-Based Artificial Viscosity for an Adaptive, Higher-Order DG
# Finite Element Method", MIT 2008.


"""
  Computes the shock sensor and the numerical viscoscity for the Persson
  and Perairi paper

  **Inputs**

   * params
   * sbp
   * q: solution on a particular element, `numDofPerNode` x `numNodesPerElement`
   * jac: mapping jacobian determinant for each node of the element, length
          `numNodesPerElement`

  **Outputs**

   * Se: the shock sensor value
   * ee: the viscoscity coefficient (constant for the entire element)
"""
function getShockSensor(params::ParamType, sbp::AbstractOperator,
                          sensor::ShockSensorPP,
                          q::AbstractMatrix{Tsol}, jac::AbstractVector{Tmsh},
                         ) where {Tsol, Tmsh}

  Tres = promote_type(Tsol, Tmsh)
  numNodesPerElement = size(q, 2)

  @unpack sensor up up_tilde up1_tilde s0 kappa e0
  #fill!(up_tilde, 0); fill!(up1_tilde, 0)

  # use density as the key variable
  for i=1:numNodesPerElement
    up[i] = q[1, i]
  end

  getFilteredSolution(params, sensor.Vp, up, up_tilde)
  getFilteredSolution(params, sensor.Vp1, up, up1_tilde)

  # compute the inner product
  num = zero(Tres)
  den = zero(Tres)
  for i=1:numNodesPerElement
    fac = sbp.w[i]/jac[i]
    delta_u = up_tilde[i] - up1_tilde[i]

    num += delta_u*fac*delta_u
    # use the filtered variables for (u, u).  This is a bit different than
    # finite element methods, where the original solution has a basis, and the
    # norm in any basis should be the same.  Here we use the filtered u rather
    # than the original because it is probably smoother.
    den += up_tilde[i]*fac*up_tilde[i]
  end


  Se = num/den
  se = log10(Se)
  # should this be a separate function from computing Se?
  if se < s0 - kappa
    ee = zero(Tres)
  elseif se > s0 - kappa && se < s0 + kappa
    ee = 0.5*e0*(1 + sinpi( (se - s0)/(2*kappa)))
  else
    ee = Tres(e0)
  end

  return Se, ee
end


"""
  Filter the solution for the shock sensor

  **Inputs**

   * params
   * vand: the [`VandermondeData`](@ref) to use
   * u_nodal: vector, length `numNodesPerElement` containing the solution at the
              nodes of the element

  **Inputs/Outputs**

   * u_filt: vector to be overwritten with the filtered solution, same length
             as `u_nodal`
"""
function getFilteredSolution(params::ParamType, vand::VandermondeData,
                             u_nodal::AbstractVector, u_filt::AbstractVector)

  # u_modal = Vpinv * u_nodal
  # u_filt = Vp * u_modal, thus this can be done in one step by:
  # u_filt = (Vp*Vpinv)*u_nodal
  smallmatvec!(vand.filt, u_nodal, u_filt)
        
  return nothing
end


#------------------------------------------------------------------------------
# ShockSensorNone

function getShockSensor(params::ParamType, sbp::AbstractOperator,
                          sensor::ShockSensorNone,
                          q::AbstractMatrix{Tsol}, jac::AbstractVector{Tmsh},
                         ) where {Tsol, Tmsh}

  error("getShockSensor called for ShockSensorNone: did you forget to specify the shock capturing scheme?")
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


function applyShockCapturing(params::ParamType, sbp::AbstractOperator,
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
  for i=1:mesh.numEl
    q_i = sview(eqn.q, :, :, i)
    jac_i = sview(mesh.jac, :, i)

    Se, ee = getShockSensor(eqn.params, sbp, sensor, q_i, jac_i)
    vals[1] = Se
    vals[2] = ee
    apf.setComponents(f_ptr, mesh.elements[i], 0, vals)
  end

  return nothing
end



