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
# ShockSensorEverywhere

function getShockSensor(params::ParamType, sbp::AbstractOperator,
                        sensor::ShockSensorEverywhere{Tsol, Tres},
                        q::AbstractMatrix, jac::AbstractVector,
                        ) where {Tsol, Tres}

  return Tsol(1.0), Tres(1.0)
end



