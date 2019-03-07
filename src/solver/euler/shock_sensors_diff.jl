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
function getShockSensor_diff(params::ParamType{Tdim}, sbp::AbstractOperator,
                      sensor::ShockSensorPP,
                      q::AbstractMatrix{Tsol},
                      coords::AbstractMatrix, dxidx::Abstract3DArray,
                      jac::AbstractVector{Tmsh},
                      Se_jac::AbstractMatrix{Tres},
                      ee_jac::AbstractMatrix{Tres}
                     ) where {Tsol, Tmsh, Tres, Tdim}
# computes the Jacobian of Se and ee wrt q. Does not take in q_dot because
# that would be way more expensive
# Se_jac and ee_jac are overwritten
# The third output argument tells if ee_jac is all zeros (hopefully a common
# case).

  numDofPerNode, numNodesPerElement = size(q)
  @unpack sensor up up_tilde up1_tilde s0 kappa e0 num_dot den_dot ee_dot lambda_max_dot
  fill!(num_dot, 0); fill!(den_dot, 0)

  @simd for i=1:numNodesPerElement
    up[i] = q[1, i]
  end

  # for getFiltered solution, the matrix itself is the Jacobian
  getFilteredSolutions(params, sensor.Vp, up, up_tilde, up1_tilde)
#  getFilteredSolution(params, sensor.Vp, up, up_tilde)
#  getFilteredSolution(params, sensor.Vp1, up, up1_tilde)
  up_tilde_dotT = sensor.Vp.filtT  # use transposed because of memory order
  #up1_tilde_dotT = sensor.Vp1.filtT
  up1_tilde_dotT = sensor.Vp1.filt1T

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

    #den += up[i]*fac*up[i]
    #den_dot[i] += 2*fac*up[i]
  end

  Se = num/den
  fac2 = 1/(den*den)
  @simd for i=1:numNodesPerElement
    Se_jac[1, i] = (num_dot[i]*den - den_dot[i]*num)*fac2
  end


  se = log10(Se)
  eejac_zero = true
  ee_dot = zeros(Tres, numNodesPerElement)
  
  if se < s0 - kappa
    #ee = zero(Tres)
    ee = Tres(0.0)
    fill!(ee_dot, 0)
  elseif se > s0 - kappa && se < s0 + kappa
    ee = 0.5*e0*(1 + sinpi( (se - s0)/(2*kappa)))

    # derivative of ee wrt Se (not se)
    fac3 = 0.5*e0*cospi( (se - s0)/(2*kappa) ) * (Float64(pi)/(2*kappa*log(10)*Se))
    fill!(ee_dot, 0)
    @simd for i=1:numNodesPerElement
      ee_dot[i] = fac3*Se_jac[1, i]
    end
    eejac_zero = false
  else
    ee = Tres(e0)
    fill!(ee_dot, 0)
    eejac_zero = false
  end

  # multiply by lambda_max * h/p to get subcell resolution
  lambda_max = zero(Tsol)
  lambda_max_dot = zeros(Tres, numDofPerNode, numNodesPerElement)  #TODO: use preallocated version
  h_avg = zero(Tmsh)
  for i=1:numNodesPerElement
    q_i = sview(q, :, i)
    lambda_max_dot_i = sview(lambda_max_dot, :, i)
    lambda_max += getLambdaMax_diff(params, q_i, lambda_max_dot_i)
    h_avg += sbp.w[i]/jac[i]
  end

  lambda_max /= numNodesPerElement
  scale!(lambda_max_dot, 1/numNodesPerElement)
  h_avg = h_avg^(1/Tdim)

  fill!(ee_jac, 0)
  for i=1:numNodesPerElement
    ee_jac[1, i] = lambda_max*h_avg*ee_dot[i]/sbp.degree
    for j=1:numDofPerNode
      ee_jac[j, i] += lambda_max_dot[j, i]*ee*h_avg/sbp.degree
    end
  end
  ee *= lambda_max*h_avg/sbp.degree

  return Se, ee, eejac_zero
end


#------------------------------------------------------------------------------
# ShockSensorNone

function getShockSensor_diff(params::ParamType, sbp::AbstractOperator,
                      sensor::ShockSensorNone,
                      q::AbstractMatrix{Tsol},
                      coords::AbstractMatrix, dxidx::Abstract3DArray,
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
                      coords::AbstractMatrix, dxidx::Abstract3DArray,
                      jac::AbstractVector{Tmsh},
                      Se_jac::AbstractMatrix{Tres},
                      ee_jac::AbstractMatrix{Tres}) where {Tsol, Tmsh, Tres}

  fill!(Se_jac, 0)
  fill!(ee_jac, 0)

  return Tres(1.0), Tres(1.0), true
end


#------------------------------------------------------------------------------
# ShockSensorHIso

function getShockSensor_diff(params::ParamType{Tdim}, sbp::AbstractOperator,
                      sensor::ShockSensorHIso,
                      q::AbstractMatrix{Tsol},
                      coords::AbstractMatrix, dxidx::Abstract3DArray,
                      jac::AbstractVector{Tmsh},
                      Se_jac::AbstractMatrix{Tres},
                      ee_jac::AbstractMatrix{Tres}) where {Tsol, Tmsh, Tres, Tdim}

  # compute | div(F) |
  # Do this in Cartesian coordinates because it makes the differentiation easier
  numDofPerNode, numNodesPerElement = size(q)

  @unpack sensor flux nrm aux_vars work res
  fill!(res, 0); fill!(nrm, 0)

  # this only needs to be 4D, but applyOperatorJac doesn't have a method
  # for that
  flux_jac = zeros(Tres, numDofPerNode, numDofPerNode, Tdim,
                   numNodesPerElement, numNodesPerElement)

  Dx = zeros(numNodesPerElement, numNodesPerElement, Tdim)
  res_jac = zeros(Tres, numDofPerNode, numDofPerNode, numNodesPerElement,
                        numNodesPerElement)
  for i=1:numNodesPerElement
    q_i = sview(q, :, i)
    for d=1:Tdim
      nrm[d] = 1
      flux_d = sview(flux, :, i, d)
      fluxjac_d = sview(flux_jac, :, :, d, i, i)

      calcEulerFlux(params, q_i, aux_vars, nrm, flux_d)
      calcEulerFlux_diff(params, q_i, aux_vars, nrm, fluxjac_d)

      nrm[d] = 0
    end
  end

  applyDx(sbp, flux, dxidx, jac, work, res)

  calcDx(sbp, dxidx, jac, Dx)
  applyOperatorJac(Dx, flux_jac, res_jac)

  # compute norm and mesh size h
  val = zero(Tres)
  fill!(Se_jac, 0)
  h_avg = zero(Tmsh)
  for i=1:numNodesPerElement
    fac = sbp.w[i]/jac[i]  # not sure about this factor of 1/|J|, because
                           # this is the strong form residual
    for j=1:numDofPerNode
      val += res[j, i]*fac*res[j, i]
    end
    h_avg += fac
  end

  # jacobian of val
  for q=1:numNodesPerElement
    for p=1:numNodesPerElement
      fac = sbp.w[p]/jac[p]
      for j=1:numDofPerNode
        for i=1:numDofPerNode
          Se_jac[j, q] += 2*res[i, p]*fac*res_jac[i, j, p, q]
        end
      end
    end
  end

  # was 2 - beta
  h_fac = h_avg^((1 - sensor.beta)/Tdim)

  ee = Tres(sensor.C_eps*h_fac*val)
  for q=1:numNodesPerElement
    for j=1:numDofPerNode
      ee_jac[j, q] = sensor.C_eps*h_fac*Se_jac[j, q]
    end
  end


  # calling | div(f) | as the shock sensor, which is somewhat arbitrary
  return val, ee, false
end

#------------------------------------------------------------------------------
# ShockSensorBO

function getShockSensor_diff(params::ParamType{Tdim}, sbp::AbstractOperator,
                      sensor::ShockSensorBO,
                      q::AbstractMatrix{Tsol},
                      coords::AbstractMatrix, dxidx::Abstract3DArray,
                      jac::AbstractVector{Tmsh},
                      Se_jac::AbstractMatrix{Tres},
                      ee_jac::AbstractMatrix{Tres}) where {Tsol, Tmsh, Tres, Tdim}

  numDofPerNode, numNodesPerElement = size(q)

  lambda_max = zero(Tsol)
  #lambda_max_dot = zeros(Tres, numDofPerNode, numNodesPerElement)
  h_avg = zero(Tmsh)
  for i=1:numNodesPerElement
    q_i = sview(q, :, i)
    lambda_max_dot_i = sview(Se_jac, :, i)
    lambda_max += getLambdaMax_diff(params, q_i, lambda_max_dot_i)
    h_avg += sbp.w[i]/jac[i]
  end

  lambda_max /= numNodesPerElement
  scale!(Se_jac, 1/numNodesPerElement)
  h_avg = h_avg^(1/Tdim)

  for i=1:numNodesPerElement
    for j=1:numDofPerNode
      ee_jac[j, i] = sensor.alpha*h_avg*Se_jac[j, i]/sbp.degree
    end
  end

  return lambda_max, sensor.alpha*h_avg*lambda_max/sbp.degree, false
end
