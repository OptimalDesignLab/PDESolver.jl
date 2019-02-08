# differentiated version of shock sensors

#------------------------------------------------------------------------------
# ShockSensorPP

"""
  Differentiated versionn of `getShockSensor` for [`ShockSensorPP`](@ref)

  **Inputs**

   * params
   * sbp
   * q
   * jac

  **Inputs/Outputs**

   * Se_jac: `numDofPerNode` x `numNodesPerElement` array to be overwritten
            with jacobian of `Se` wrt `q`
   * ee_jac: `numDofPerNode` x `numNodesPerElement` array to be overwritten
             with jacobian of `ee` wrt `q`

  **Outputs**
  
   * Se
   * ee
   * is_const: true if the derivative of ee wrt q is zero.  This is useful
               for skipping terms in the jacobian calculation
"""
function getShockSensor_diff(params::ParamType, sbp::AbstractOperator,
                      sensor::ShockSensorPP,
                      q::AbstractMatrix{Tsol},
                      jac::AbstractVector{Tmsh},
                      Se_jac::AbstractMatrix{Tres},
                      ee_jac::AbstractMatrix{Tres}) where {Tsol, Tmsh, Tres}
# computes the Jacobian of Se and ee wrt q. Does not take in q_dot because
# that would be way more expensive
# Se_jac and ee_jac are overwritten
# The third output argument tells if ee_jac is all zeros (hopefully a common
# case).

  numDofPerNode, numNodesPerElement = size(q)
  @unpack sensor up up_tilde up1_tilde s0 kappa e0 num_dot den_dot
  fill!(num_dot, 0); fill!(den_dot, 0)

  @simd for i=1:numNodesPerElement
    up[i] = q[1, i]
  end

  # for getFiltered solution, the matrix itself is the Jacobian
  getFilteredSolution(params, sensor.Vp, up, up_tilde)
  getFilteredSolution(params, sensor.Vp1, up, up1_tilde)
  up_tilde_dotT = sensor.Vp.filtT  # use transposed because of memory order
  up1_tilde_dotT = sensor.Vp1.filtT

  # compute the inner product
  num = zero(Tres)
  den = zero(Tres)

  @simd for i=1:numNodesPerElement
    fac = sbp.w[i]/jac[i]
    delta_u = up_tilde[i] - up1_tilde[i]


    num += delta_u*fac*delta_u
    @simd for j=1:numNodesPerElement
      delta_u_dot = up_tilde_dotT[j, i] - up1_tilde_dotT[j, i]
      num_dot[j] += 2*fac*delta_u*delta_u_dot
    end

    # use the filtered variables for (u, u).  This is a bit different than
    # finite element methods, where the original solution has a basis, and the
    # norm in any basis should be the same.  Here we use the filtered u rather
    # than the original because it is probably smoother.
    den += up_tilde[i]*fac*up_tilde[i]
    @simd for j=1:numNodesPerElement
      den_dot[j] += 2*fac*up_tilde[i]*up_tilde_dotT[j, i]
    end
  end

  Se = num/den
  fac2 = 1/(den*den)
  @simd for i=1:numNodesPerElement
    Se_jac[1, i] = (num_dot[i]*den - den_dot[i]*num)*fac2
  end


  se = log10(Se)
  eejac_zero = true
  
  if se < s0 - kappa
    ee = zero(Tres)
    fill!(ee_jac, 0)
  elseif se > s0 - kappa && se < s0 + kappa
    ee = 0.5*e0*(1 + sinpi( (se - s0)/(2*kappa)))

    # derivative of ee wrt Se (not se)
    fac3 = 0.5*e0*cospi( (se - s0)/(2*kappa) ) * (Float64(pi)/(2*kappa*log(10)*Se))
    fill!(ee_jac, 0)
    @simd for i=1:numNodesPerElement
      ee_jac[1, i] = fac3*Se_jac[1, i]
    end
    eejac_zero = false
  else
    ee = Tres(e0)
    fill!(ee_jac, 0)
  end

  return Se, ee, eejac_zero
end


#------------------------------------------------------------------------------
# ShockSensorNone

function getShockSensor_diff(params::ParamType, sbp::AbstractOperator,
                      sensor::ShockSensorNone,
                      q::AbstractMatrix{Tsol},
                      jac::AbstractVector{Tmsh},
                      Se_jac::AbstractMatrix{Tres},
                      ee_jac::AbstractMatrix{Tres}) where {Tsol, Tmsh, Tres}

  error("getShockSensori_diff called for ShockSensorNone: did you forget to specify the shock capturing scheme?")
end

#------------------------------------------------------------------------------
# ShockSensorEverywhere

function getShockSensor_diff(params::ParamType, sbp::AbstractOperator,
                      sensor::ShockSensorEverywhere,
                      q::AbstractMatrix{Tsol},
                      jac::AbstractVector{Tmsh},
                      Se_jac::AbstractMatrix{Tres},
                      ee_jac::AbstractMatrix{Tres}) where {Tsol, Tmsh, Tres}

  fill!(Se_jac, 0)
  fill!(ee_jac, 0)

  return Tres(1.0), Tres(1.0), true
end
