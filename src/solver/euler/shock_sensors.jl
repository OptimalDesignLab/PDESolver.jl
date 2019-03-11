#TODO: move shock sensor implementation to own file

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
function getShockSensor(params::ParamType{Tdim}, sbp::AbstractOperator,
                          sensor::ShockSensorPP,
                          q::AbstractMatrix{Tsol}, coords::AbstractMatrix,
                          dxidx::Abstract3DArray, jac::AbstractVector{Tmsh},
                         ) where {Tsol, Tmsh, Tdim}

  Tres = promote_type(Tsol, Tmsh)
  numNodesPerElement = size(q, 2)

  @unpack sensor up up_tilde up1_tilde s0 kappa e0
  #fill!(up_tilde, 0); fill!(up1_tilde, 0)

  # use density as the key variable
  for i=1:numNodesPerElement
    up[i] = q[1, i]
  end

  getFilteredSolutions(params, sensor.Vp, up, up_tilde, up1_tilde)
  #getFilteredSolution(params, sensor.Vp, up, up_tilde)
  #getFilteredSolution(params, sensor.Vp1, up, up1_tilde)

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

    # see if this makes the sensor less sensitive
    #den += up[i]*fac*up[i]
  end

  Se = num/den
  se = log10(Se)
  
  # should this be a separate function from computing Se?
  if se < s0 - kappa
    ee = Tres(0.0)
  elseif se > s0 - kappa && se < s0 + kappa
    ee = 0.5*e0*(1 + sinpi( (se - s0)/(2*kappa)))
  else
    ee = Tres(e0)
  end

  # scale by lambda_max * h/p to get subcell resolution
  lambda_max = zero(Tsol)
  h_avg = zero(Tmsh)
  for i=1:numNodesPerElement
    q_i = sview(q, :, i)
    lambda_max += getLambdaMax(params, q_i)
    h_avg += sbp.w[i]/jac[i]
  end

  lambda_max /= numNodesPerElement
  h_avg = h_avg^(1/Tdim)

  ee *= lambda_max*h_avg/sbp.degree
  
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

"""
  Computed the two filtered solutions required by the shock sensor (well, only
  one of them is required, the other one may or may not be useful)

  **Inputs**

   * params
   * vand: `VandermondeData`
   * u_nodal: `numDofPerNode` x `numNodesPerElement` nodal solution

  **Inputs/Outputs**

   * u_filt1: The result of `u = Vp*pinv(Vp)*u_nodal`  This projects the
              solution onto the orthogonal basis and then projects it back.
              This may or may not be useful
   * u_filt2: The result of projecting `u` onto the orthgonal basis and then
              projecting all but the highest mode(s) back to a nodal back to
              the nodal basis
"""
function getFilteredSolutions(params::ParamType, vand::VandermondeData,
                              u_nodal, u_filt1, u_filt2)


  smallmatvec!(vand.filt, u_nodal, u_filt1)
  smallmatvec!(vand.filt1, u_nodal, u_filt2)

  return nothing
end

#------------------------------------------------------------------------------
# ShockSensorNone

function getShockSensor(params::ParamType, sbp::AbstractOperator,
                          sensor::ShockSensorNone,
                          q::AbstractMatrix{Tsol}, coords::AbstractMatrix,
                          dxidx::Abstract3DArray, jac::AbstractVector{Tmsh},
                         ) where {Tsol, Tmsh}

  error("getShockSensor called for ShockSensorNone: did you forget to specify the shock capturing scheme?")
end

#------------------------------------------------------------------------------
# ShockSensorEverywhere

function getShockSensor(params::ParamType, sbp::AbstractOperator,
                        sensor::ShockSensorEverywhere{Tsol, Tres},
                        q::AbstractMatrix, coords::AbstractMatrix,
                        dxidx::Abstract3DArray, jac::AbstractVector{Tmsh},
                        ) where {Tsol, Tres, Tmsh}

  return Tsol(1.0), Tres(1.0)
end

#------------------------------------------------------------------------------
# ShockSensorHIso

function getShockSensor(params::ParamType{Tdim}, sbp::AbstractOperator,
                        sensor::ShockSensorHIso{Tsol, Tres},
                        q::AbstractMatrix, coords::AbstractMatrix,
                        dxidx::Abstract3DArray, jac::AbstractVector{Tmsh},
                         ) where {Tsol, Tres, Tmsh, Tdim}

  # compute | div(F) |
  # Do this in Cartesian coordinates because it makes the differentiation easier
  numDofPerNode, numNodesPerElement = size(q)
  res = sensor.res

  computeStrongResidual(params, sbp, sensor.strongdata, q, dxidx, jac, res)
  val = computeL2Norm(params, sbp, jac, res)
  h_avg = computeElementVolume(params, sbp, jac)

  h_fac = h_avg^((2 - sensor.beta)/Tdim)

  ee = Tres(sensor.C_eps*h_fac*val)
#=
  if ee < 1e-12  # make this zero so the shockmesh machinery knows to exclude
                 # this element, saving cost
    ee = Tres(0)
  end
=#


  # calling | div(f) | as the shock sensor, which is somewhat arbitrary
  return val, ee
end

"""
  Computes del * F.  where del is the diverage in the xyz directions,
  for a single element.

  **Inputs**

   * params
   * sbp
   * data: `StrongFormData` to hold temporary arrays
   * q: solution on the entire element, `numDofPerNode` x `numNodesPerElement`
   * dxidx: scaled mapping jacobian for entire element, `dim` x `dim` x
            `numNodesPerElement`
   * jac: mapping jacobian determinant for entire element, `numNodesPerElement`

  **Inputs/Outputs**

   * res: array to overwrite with residual, same size as `q`.
"""
function computeStrongResidual(params::ParamType{Tdim}, sbp::AbstractOperator,
                               data::StrongFormData{Tsol, Tres},
                               q::AbstractMatrix,
                               dxidx::Abstract3DArray,
                               jac::AbstractVector,
                               res::AbstractMatrix) where {Tdim, Tsol, Tres}

  @unpack data flux nrm aux_vars work
  fill!(res, 0); fill!(nrm, 0)

  numDofPerNode, numNodesPerElement = size(q)

  for i=1:numNodesPerElement
    q_i = sview(q, :, i)
    for d=1:Tdim
      nrm[d] = 1
      flux_d = sview(flux, :, i, d)

      calcEulerFlux(params, q_i, aux_vars, nrm, flux_d)

      nrm[d] = 0
    end
  end

  applyDx(sbp, flux, dxidx, jac, work, res)

  return nothing
end

"""
  Computes the integral L2 norm for a single element

  **Inputs**

   * params
   * sbp
   * jac: mapping jacobian determinant for the element, `numNodesPerElement`
   * res: `numDofPerNode` x `numNodesPerElement` array to compute the norm
          of.  
"""
function computeL2Norm(params::ParamType{Tdim}, sbp::AbstractOperator,
                       jac::AbstractVector, res::AbstractMatrix{Tres}
                      ) where {Tres, Tdim}

  numDofPerNode, numNodesPerElement = size(res)

  # compute norm
  val = zero(Tres)
  for i=1:numNodesPerElement
    fac = sbp.w[i]/jac[i]
    for j=1:numDofPerNode
      val += res[j, i]*fac*res[j, i]
    end
  end

  val = sqrt(val)

  return val
end


"""
  Computes volume of element by integrating the mapping jacobian determinant
"""
function computeElementVolume(params::ParamType{Tdim}, sbp::AbstractOperator,
                         jac::AbstractVector{Tmsh}) where {Tdim, Tmsh}

  # compute norm and mesh size h
  vol = zero(Tmsh)
  for i=1:length(jac)
    vol += sbp.w[i]/jac[i]
  end

  return vol
end




#------------------------------------------------------------------------------
# ShockSensorBO

function getShockSensor(params::ParamType{Tdim}, sbp::AbstractOperator,
                        sensor::ShockSensorBO{Tsol, Tres},
                        q::AbstractMatrix, coords::AbstractMatrix,
                        dxidx::Abstract3DArray, jac::AbstractVector{Tmsh},
                         ) where {Tsol, Tres, Tmsh, Tdim}

  numDofPerNode, numNodesPerElement = size(q)

  lambda_max = zero(Tsol)
  h_avg = zero(Tmsh)
  for i=1:numNodesPerElement
    q_i = sview(q, :, i)
    lambda_max += getLambdaMax(params, q_i)
    h_avg += sbp.w[i]/jac[i]
  end

  lambda_max /= numNodesPerElement
  h_avg = h_avg^(1/Tdim)

  return lambda_max, sensor.alpha*lambda_max*h_avg/sbp.degree
end

#------------------------------------------------------------------------------
# ShockSensorHHO



function getShockSensor(params::ParamType{Tdim}, sbp::AbstractOperator,
                        sensor::ShockSensorHHO{Tsol, Tres},
                        q::AbstractMatrix, coords::AbstractMatrix,
                        dxidx::Abstract3DArray, jac::AbstractVector{Tmsh},
                         ) where {Tsol, Tres, Tmsh, Tdim}
#NOTE: the real sensor produces spatially varying viscoscity.  Here we take
#      the norm so it fits into the elementwise-constant viscoscity model

  numDofPerNode, numNodesPerElement = size(q)

  @unpack sensor p_dot press_el press_dx work res
#  p_dot = zeros(Tsol, numDofPerNode)
#  press_el = zeros(Tsol, 1, numNodesPerElement)
#  press_dx = zeros(Tres, 1, numNodesPerElement, Tdim)
#  work = zeros(Tres, 1, numNodesPerElement, Tdim)
#  res = zeros(Tres, numDofPerNode, numNodesPerElement)

  fill!(press_dx, 0); fill!(res, 0)

  #TODO: in the real sensor, this should include an anisotropy factor
  h_avg = computeElementVolume(params, sbp, jac)
  h_fac = (h_avg^(1/Tdim))/(sbp.degree + 1)

  # compute Rm
  computeStrongResidual(params, sbp, sensor.strongdata, q, dxidx, jac, res)

  # compute Rp
  for i=1:numNodesPerElement
    q_i = ro_sview(q, :, i)
    press = calcPressure_diff(params, q_i, p_dot)
    press_el[1, i] = press

    for j=1:numDofPerNode
      res[j, i] *= p_dot[j]/press
    end
  end

  # compute pressure switch fp
  applyDx(sbp, press_el, dxidx, jac, work, press_dx)
  for i=1:numNodesPerElement
    val = zero(Tres)
    for d=1:Tdim
      val += press_dx[1, i, d]*press_dx[1, i, d]
    end
    val = sqrt(val)*h_fac/(press_el[1, i] + 1e-12)

    # lump fp in with |Rm| for now
    for j=1:numDofPerNode
      res[j, i] *= val
    end
  end

  # compute norm to make this elementwise constant
  val = computeL2Norm(params, sbp, jac, res)

  ee = val*h_fac*h_fac*sensor.C_eps
  #ee = val*h_fac*sensor.C_eps

  return val, ee
end
