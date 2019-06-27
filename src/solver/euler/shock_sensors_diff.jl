# differentiated version of shock sensors

#------------------------------------------------------------------------------
# ShockSensorPP

"""
  Differentiated versionn of `getShockSensor` for [`ShockSensorPP`](@ref)
"""
function getShockSensor_diff(params::ParamType{Tdim}, sbp::AbstractOperator,
                      sensor::ShockSensorPP,
                      q::AbstractMatrix{Tsol}, elnum::Integer,
                      coords::AbstractMatrix, dxidx::Abstract3DArray,
                      jac::AbstractVector{Tmsh},
                      Se_jac::Abstract4DArray{Tres},
                      ee_jac::Abstract4DArray{Tres}
                     ) where {Tsol, Tmsh, Tres, Tdim}
# computes the Jacobian of Se and ee wrt q. Does not take in q_dot because
# that would be way more expensive
# Se_jac and ee_jac are overwritten
# The third output argument tells if ee_jac is all zeros (hopefully a common
# case).

  numDofPerNode, numNodesPerElement = size(q)
  dim = size(coords, 1)

  @unpack sensor up up_tilde up1_tilde s0 kappa e0 num_dot den_dot ee_dot lambda_max_dot
  fill!(num_dot, 0); fill!(den_dot, 0); fill!(Se_jac, 0); fill!(ee_jac, 0)

  @simd for i=1:numNodesPerElement
    up[i] = q[1, i]
  end

  # for getFiltered solution, the matrix itself is the Jacobian
  getFilteredSolutions(params, sensor.Vp, up, up_tilde, up1_tilde)
#  getFilteredSolution(params, sensor.Vp, up, up_tilde)
#  getFilteredSolution(params, sensor.Vp1, up, up1_tilde)
  up_tilde_dotT = sensor.Vp.filtT  # use transposed because of memory order
  up1_tilde_dotT = sensor.Vp.filt1T
  #up1_tilde_dotT = sensor.Vp1.filtT
  #up1_tilde_dotT = sensor.Vp1.filt1T

  # compute the inner product
  num = zero(Tres)
  den = zero(Tres)

  if sensor.use_filtered
    @simd for i=1:numNodesPerElement
      fac = sbp.w[i]/jac[i]
      delta_u = up_tilde[i] - up1_tilde[i]

      num += delta_u*fac*delta_u
      @simd for j=1:numNodesPerElement
        delta_u_dot = up_tilde_dotT[j, i] - up1_tilde_dotT[j, i]
        num_dot[j] += 2*fac*delta_u*delta_u_dot
      end

      den += up_tilde[i]*fac*up_tilde[i]
      @simd for j=1:numNodesPerElement
        den_dot[j] += 2*fac*up_tilde[i]*up_tilde_dotT[j, i]
      end
    end
  else
    @simd for i=1:numNodesPerElement
      fac = sbp.w[i]/jac[i]
      delta_u = up[i] - up1_tilde[i]

      num += delta_u*fac*delta_u

      @simd for j=1:numNodesPerElement
        delta_u_dot = -up1_tilde_dotT[j, i]
        if j == i
          delta_u_dot += 1
        end

        num_dot[j] += 2*fac*delta_u*delta_u_dot
      end

      den += up[i]*fac*up[i]
      den_dot[i] += 2*fac*up[i]
    end
  end

  Se = num/den
  fac2 = 1/(den*den)
#  fac3 = 1/(2*Se)
  @simd for i=1:numNodesPerElement
    Se_jac[1, 1, 1, i] = (num_dot[i]*den - den_dot[i]*num)*fac2#*fac3
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
      ee_dot[i] = fac3*Se_jac[1, 1, 1, i]
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
  @simd for i=1:numNodesPerElement
    q_i = sview(q, :, i)
    lambda_max_dot_i = sview(lambda_max_dot, :, i)
    lambda_max += getLambdaMax_diff(params, q_i, lambda_max_dot_i)
    h_avg += sbp.w[i]/jac[i]
  end

  lambda_max /= numNodesPerElement
  scale!(lambda_max_dot, 1/numNodesPerElement)
  h_avg = h_avg^(1/Tdim)

  fill!(ee_jac, 0)
  @simd for i=1:numNodesPerElement
    ee_jac[1, 1, 1, i] = lambda_max*h_avg*ee_dot[i]/sbp.degree
    @simd for j=1:numDofPerNode
      ee_jac[1, j, 1, i] += lambda_max_dot[j, i]*ee*h_avg/sbp.degree
    end
  end
  ee *= lambda_max*h_avg/sbp.degree

  # fill in the rest of the jac arrays
  @simd for q=1:numNodesPerElement
    @simd for p=1:numNodesPerElement
      @simd for j=1:numDofPerNode
        @simd for i=1:dim
          Se_jac[i, j, p, q] = Se_jac[1, j, 1, q]
          ee_jac[i, j, p, q] = ee_jac[1, j, 1, q]
        end
      end
    end
  end

  return eejac_zero
end


function getShockSensor_revq(params::ParamType{Tdim}, sbp::AbstractOperator,
                        sensor::ShockSensorPP{Tsol, Tres},
                        q::AbstractMatrix, q_bar::AbstractMatrix,
                        elnum::Integer,
                        coords::AbstractMatrix,
                        dxidx::Abstract3DArray, jac::AbstractVector{Tmsh},
                        ee_mat::AbstractMatrix, ee_mat_bar::AbstractMatrix
                        ) where {Tsol, Tres, Tmsh, Tdim}

  numNodesPerElement = size(q, 2)

  @unpack sensor up up_tilde up1_tilde s0 kappa e0 up_tilde_bar up1_tilde_bar
  @unpack sensor up_bar
  fill!(up_tilde_bar, 0); fill!(up1_tilde_bar, 0); fill!(up_bar, 0)


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
  if sensor.use_filtered
    for i=1:numNodesPerElement
      fac = sbp.w[i]/jac[i]
      delta_u = up_tilde[i] - up1_tilde[i]

      num += delta_u*fac*delta_u
      den += up_tilde[i]*fac*up_tilde[i]
    end
  else
    for i=1:numNodesPerElement
      fac = sbp.w[i]/jac[i]
      delta_u = up[i] - up1_tilde[i]

      num += delta_u*fac*delta_u
      den += up[i]*fac*up[i]
    end
  end

  Se = num/den
  se = log10(Se)
 
  # should this be a separate function from computing Se?
  if se < s0 - kappa
    ee = Tres(0.0)
    is_nonzero = false
  elseif se > s0 - kappa && se < s0 + kappa
    ee = 0.5*e0*(1 + sinpi( (se - s0)/(2*kappa)))
    is_nonzero = true
  else
    ee = Tres(e0)
    is_nonzero = true
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

  ee_final = ee*lambda_max*h_avg/sbp.degree


  # use elementwise constant
  fill!(ee_mat, ee)

  #---------------------
  # reverse sweep
  ee_final_bar = zero(Tres)
  for i=1:numNodesPerElement
    for d=1:Tdim
      ee_final_bar += ee_mat_bar[d, i]
    end
  end

  ee_bar = zero(Tres); lambda_max_bar = zero(Tres)
  ee_bar         += ee_final_bar*lambda_max*h_avg/sbp.degree
  lambda_max_bar += ee_final_bar*ee*h_avg/sbp.degree

  lambda_max_bar /= numNodesPerElement
  for i=1:numNodesPerElement
    q_i = ro_sview(q, :, i)
    q_bar_i = sview(q_bar, :, i)
    getLambdaMax_revq(params, q_i, q_bar_i, lambda_max_bar)
  end

  # if s < s0 - kappa or > s0 + kappa, then ee is either 1 or 0, so there is
  # no q dependence and no need for further computation
  if se > s0 - kappa && se < s0 + kappa

    # derivative of ee wrt Se (not se)
    Se_bar = ee_bar*0.5*e0*cospi( (se - s0)/(2*kappa) ) * (Float64(pi)/(2*kappa*log(10)*Se))

    num_bar = Se_bar/den
    den_bar = -num*Se_bar/(den*den)

#    up_tilde_bar = zeros(Tsol, numNodesPerElement)
#    up1_tilde_bar = zeros(Tsol, numNodesPerElement)a

    if sensor.use_filtered
      for i=1:numNodesPerElement
        fac = sbp.w[i]/jac[i]
        delta_u = up_tilde[i] - up1_tilde[i]

        delta_u_bar       = 2*delta_u*fac*num_bar
        up_tilde_bar[i]  += 2*up_tilde[i]*fac*den_bar
        up_tilde_bar[i]  += delta_u_bar
        up1_tilde_bar[i] -= delta_u_bar
      end
    else
      for i=1:numNodesPerElement
        fac = sbp.w[i]/jac[i]
        delta_u = up[i] - up1_tilde[i]

        delta_u_bar       = 2*delta_u*fac*num_bar
        up_bar[i]        += 2*up[i]*fac*den_bar
        up_bar[i]        += delta_u_bar
        up1_tilde_bar[i] -= delta_u_bar
      end
    end

    getFilteredSolutions_rev(params, sensor.Vp, up_bar, up_tilde_bar, up1_tilde_bar)

    for i=1:numNodesPerElement
      q_bar[1, i] += up_bar[i]
    end
  end

  return is_nonzero
end


function getShockSensor_revm(params::ParamType{Tdim}, sbp::AbstractOperator,
                        sensor::ShockSensorPP{Tsol, Tres},
                        q::AbstractMatrix,
                        elnum::Integer,
                        coords::AbstractMatrix, coords_bar::AbstractMatrix,
                        dxidx::Abstract3DArray, dxidx_bar::Abstract3DArray,
                        jac::AbstractVector{Tmsh}, jac_bar::AbstractVector,
                        ee_mat_bar::AbstractMatrix
                       ) where {Tsol, Tres, Tmsh, Tdim}

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
  if sensor.use_filtered
    for i=1:numNodesPerElement
      fac = sbp.w[i]/jac[i]
      delta_u = up_tilde[i] - up1_tilde[i]

      num += delta_u*fac*delta_u
      den += up_tilde[i]*fac*up_tilde[i]

      # see if this makes the sensor less sensitive
      #den += up[i]*fac*up[i]
    end
  else
    for i=1:numNodesPerElement
      fac = sbp.w[i]/jac[i]
      delta_u = up[i] - up1_tilde[i]

      num += delta_u*fac*delta_u
      den += up[i]*fac*up[i]
    end

  end

  Se = num/den
  se = log10(Se)
 
  # should this be a separate function from computing Se?
  if se < s0 - kappa
    ee = Tres(0.0)
    is_nonzero = false
  elseif se > s0 - kappa && se < s0 + kappa
    ee = 0.5*e0*(1 + sinpi( (se - s0)/(2*kappa)))
    is_nonzero = true
  else
    ee = Tres(e0)
    is_nonzero = true
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
  h_avg2 = h_avg^(1/Tdim)

  ee_final = ee*lambda_max*h_avg2/sbp.degree


  # use elementwise constant
  #fill!(ee_mat, ee)

  #---------------------
  # reverse sweep
  ee_final_bar = zero(Tres)
  for i=1:numNodesPerElement
    for d=1:Tdim
      ee_final_bar += ee_mat_bar[d, i]
    end
  end

  ee_bar = zero(Tres); lambda_max_bar = zero(Tres); h_avg2_bar = zero(Tres)
  ee_bar         += ee_final_bar*lambda_max*h_avg2/sbp.degree
  h_avg2_bar += ee*lambda_max*ee_final_bar/sbp.degree

  h_avg_bar = (1/Tdim)*h_avg2_bar*h_avg^((1/Tdim) - 1)

  for i=1:numNodesPerElement
    jac_bar[i] += -sbp.w[i]*h_avg_bar/(jac[i]*jac[i])
  end
#=
  lambda_max_bar /= numNodesPerElement
  for i=1:numNodesPerElement
    q_i = ro_sview(q, :, i)
    q_bar_i = sview(q_bar, :, i)
    getLambdaMax_revq(params, q_i, q_bar_i, lambda_max_bar)
  end
=#
  # if s < s0 - kappa or > s0 + kappa, then ee is either 1 or 0, so there is
  # no metric dependence and no need for further computation
  if se > s0 - kappa && se < s0 + kappa

    # derivative of ee wrt Se (not se)
    Se_bar = ee_bar*0.5*e0*cospi( (se - s0)/(2*kappa) ) * (Float64(pi)/(2*kappa*log(10)*Se))

    num_bar = Se_bar/den
    den_bar = -num*Se_bar/(den*den)

    #TODO: why didn't this error?
    if sensor.use_filtered
      for i=1:numNodesPerElement
        fac = sbp.w[i]/jac[i]
        delta_u = up_tilde[i] - up1_tilde[i]

        fac_bar     = delta_u*delta_u*num_bar
        fac_bar    += up_tilde[i]*up_tilde[i]*den_bar
        jac_bar[i] += -sbp.w[i]*fac_bar/(jac[i]*jac[i])
      end
    else
      for i=1:numNodesPerElement
        fac = sbp.w[i]/jac[i]
        delta_u = up[i] - up1_tilde[i]

        fac_bar     = delta_u*delta_u*num_bar
        fac_bar    += up[i]*up[i]*den_bar
        jac_bar[i] += -sbp.w[i]*fac_bar/(jac[i]*jac[i])
      end

    end
  end

  return is_nonzero
end


"""
  Computed the two filtered solutions required by the shock sensor (well, only
  one of them is required, the other one may or may not be useful)

  **Inputs**

   * params
   * vand: `VandermondeData`
   * u_nodal_bar:

  **Inputs/Outputs**

   * u_filt1_bar:
   * u_filt2_bar
"""
function getFilteredSolutions_rev(params::ParamType, vand::VandermondeData,
                              u_nodal, u_filt1, u_filt2)


  smallmatTvec_kernel!(vand.filt, u_filt1, u_nodal, 1, 1)
  smallmatTvec_kernel!(vand.filt1, u_filt2,  u_nodal, 1, 1)

  return nothing
end




#------------------------------------------------------------------------------
# ShockSensorNone

function getShockSensor_diff(params::ParamType, sbp::AbstractOperator,
                      sensor::ShockSensorNone,
                      q::AbstractMatrix{Tsol}, elnum::Integer,
                      coords::AbstractMatrix, dxidx::Abstract3DArray,
                      jac::AbstractVector{Tmsh},
                      Se_jac::Abstract4DArray{Tres},
                      ee_jac::Abstract4DArray{Tres}) where {Tsol, Tmsh, Tres}

  error("getShockSensor_diff called for ShockSensorNone: did you forget to specify the shock capturing scheme?")
end

#------------------------------------------------------------------------------
# ShockSensorEverywhere

function getShockSensor_diff(params::ParamType, sbp::AbstractOperator,
                      sensor::ShockSensorEverywhere,
                      q::AbstractMatrix{Tsol}, elnum::Integer,
                      coords::AbstractMatrix, dxidx::Abstract3DArray,
                      jac::AbstractVector{Tmsh},
                      Se_jac::Abstract4DArray{Tres},
                      ee_jac::Abstract4DArray{Tres}) where {Tsol, Tmsh, Tres}

  fill!(Se_jac, 0)
  fill!(ee_jac, 0)

  return true
end

function getShockSensor_revq(params::ParamType, sbp::AbstractOperator,
                        sensor::ShockSensorEverywhere{Tsol, Tres},
                        q::AbstractMatrix, q_bar::AbstractMatrix,
                        elnum::Integer,
                        coords::AbstractMatrix,
                        dxidx::Abstract3DArray, jac::AbstractVector{Tmsh},
                        ee_mat::AbstractMatrix, ee_bar::AbstractMatrix
                        ) where {Tsol, Tres, Tmsh}

  # ee is constant, so q_bar += 0*ee_bar

  return true
end


#------------------------------------------------------------------------------
# ShockSensorVelocity

function getShockSensor_diff(params::ParamType, sbp::AbstractOperator,
                      sensor::ShockSensorVelocity,
                      q::AbstractMatrix{Tsol}, elnum::Integer,
                      coords::AbstractMatrix, dxidx::Abstract3DArray,
                      jac::AbstractVector{Tmsh},
                      Se_jac::Abstract4DArray{Tres},
                      ee_jac::Abstract4DArray{Tres}) where {Tsol, Tmsh, Tres}

  numNodesPerElement = size(q, 2)
  dim = size(coords, 1)

  fill!(Se_jac, 0)
  fill!(ee_jac, 0)

  # derivative with momentum in direction d = d, all others zero
  alpha = sensor.alpha
  @simd for p=1:numNodesPerElement
    @simd for d=1:dim
      Se_jac[d, d+1, p, p] = d*alpha
      ee_jac[d, d+1, p, p] = d*alpha
    end
  end

  return false
end


function getShockSensor_revq(params::ParamType, sbp::AbstractOperator,
                        sensor::ShockSensorVelocity{Tsol, Tres},
                        q::AbstractMatrix, q_bar::AbstractMatrix,
                        elnum::Integer,
                        coords::AbstractMatrix,
                        dxidx::Abstract3DArray, jac::AbstractVector{Tmsh},
                        ee_mat::AbstractMatrix, ee_bar::AbstractMatrix
                        ) where {Tsol, Tres, Tmsh}

  # ee is constant, so q_bar += 0*ee_bar
  numNodesPerElement = size(q, 2)
  dim = size(coords, 1)

  numNodesPerElement = size(q, 2)
  dim = size(coords, 1)

  @simd for i=1:numNodesPerElement
    @simd for d=1:dim
      ee_mat[d, i] = d*q[d+1, i]
      q_bar[d+1, i] += ee_bar[d, i]*d
    end
  end

  return true
end


function getShockSensor_revm(params::ParamType{Tdim}, sbp::AbstractOperator,
                        sensor::ShockSensorVelocity{Tsol, Tres},
                        q::AbstractMatrix,
                        elnum::Integer,
                        coords::AbstractMatrix, coords_bar::AbstractMatrix,
                        dxidx::Abstract3DArray, dxidx_bar::Abstract3DArray,
                        jac::AbstractVector{Tmsh}, jac_bar::AbstractVector,
                        ee_bar::AbstractMatrix
                       ) where {Tsol, Tres, Tmsh, Tdim}


  # nothing to do here

  return nothing
end

#------------------------------------------------------------------------------
# ShockSensorHIso

function getShockSensor_diff(params::ParamType{Tdim}, sbp::AbstractOperator,
                      sensor::ShockSensorHIso,
                      q::AbstractMatrix{Tsol}, elnum::Integer,
                      coords::AbstractMatrix, dxidx::Abstract3DArray,
                      jac::AbstractVector{Tmsh},
                      Se_jac::Abstract4DArray{Tres},
                      ee_jac::Abstract4DArray{Tres}) where {Tsol, Tmsh, Tres, Tdim}

  # compute | div(F) |
  # Do this in Cartesian coordinates because it makes the differentiation easier
  numDofPerNode, numNodesPerElement = size(q)
  dim = size(coords, 1)
  @unpack sensor res res_jac

  _Se_jac = zeros(Tres, numDofPerNode, numNodesPerElement)
  computeStrongResidual_diff(params, sbp, sensor.strongdata,
                             q, coords, dxidx, jac, res, res_jac)
  val = computeL2Norm_diff(params, sbp, jac, res, res_jac, _Se_jac)

  h_avg = computeElementVolume(params, sbp, jac)
  h_fac = h_avg^((2 - sensor.beta)/Tdim)

  ee = Tres(sensor.C_eps*h_fac*val)
  @simd for q=1:numNodesPerElement
    @simd for j=1:numDofPerNode
      ee_jac[1, j, 1, q] = sensor.C_eps*h_fac*_Se_jac[j, q]
    end
  end

  # fill in the rest of the jac arrays
  @simd for q=1:numNodesPerElement
    @simd for p=1:numNodesPerElement
      @simd for j=1:numDofPerNode
        @simd for i=1:dim
          Se_jac[i, j, p, q] = _Se_jac[j, q]
          ee_jac[i, j, p, q] = ee_jac[1, j, 1, q]
        end
      end
    end
  end

  # calling | div(f) | as the shock sensor, which is somewhat arbitrary
  return false
end

"""
  Differentiated version

  **Inputs**

   * params
   * sbp
   * data: `StrongFormData` to hold temporary arrays
   * q: solution on the entire element, `numDofPerNode` x `numNodesPerElement`
   * dxidx: scaled mapping jacobian for entire element, `dim` x `dim` x
            `numNodesPerElement`
   * jac: mapping jacobian determinant for entire element, `numNodesPerElement`

  **Inputs/Outputs**

   * val_jac: jacobian of `val` wrt `q`, same size as `q`

  **Outputs**

   * val
"""
function computeStrongResidual_diff(params::ParamType{Tdim},
                      sbp::AbstractOperator,
                      data::StrongFormData{Tsol, Tres},
                      q::AbstractMatrix,
                      coords::AbstractMatrix, dxidx::Abstract3DArray,
                      jac::AbstractVector{Tmsh},
                      res::AbstractMatrix, res_jac::Abstract4DArray) where {Tsol, Tres, Tmsh, Tdim}


  numDofPerNode, numNodesPerElement = size(q)
  @unpack data flux nrm aux_vars work flux_jac Dx
  fill!(res, 0); fill!(nrm, 0); fill!(flux_jac, 0); fill!(res_jac, 0)

  @simd for i=1:numNodesPerElement
    q_i = sview(q, :, i)
    @simd for d=1:Tdim
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

  return nothing
end


function computeStrongResidual_revq(params::ParamType{Tdim},
                               sbp::AbstractOperator,
                               data::StrongFormData{Tsol, Tres},
                               q::AbstractMatrix, q_bar::AbstractMatrix,
                               dxidx::Abstract3DArray,
                               jac::AbstractVector,
                               res_bar::AbstractMatrix
                              ) where {Tdim, Tsol, Tres}

  @unpack data flux nrm aux_vars work flux_bar
  fill!(nrm, 0); fill!(flux_bar, 0)

  numDofPerNode, numNodesPerElement = size(q)
#=
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
=#
  # reverse sweep
  applyDx_revq(sbp, flux_bar, dxidx, jac, work, res_bar)

  @simd for i=1:numNodesPerElement
    q_i = sview(q, :, i)
    q_bar_i = sview(q_bar, :, i)
    @simd for d=1:Tdim
      nrm[d] = 1
      flux_bar_d = ro_sview(flux_bar, :, i, d)

      calcEulerFlux_revq(params, q_i, q_bar_i, aux_vars, nrm, flux_bar_d)

      nrm[d] = 0
    end
  end

  return nothing
end


function computeStrongResidual_revm(params::ParamType{Tdim},
                             sbp::AbstractOperator,
                             data::StrongFormData{Tsol, Tres},
                             q::AbstractMatrix,
                             dxidx::Abstract3DArray, dxidx_bar::Abstract3DArray,
                             jac::AbstractVector, jac_bar::AbstractVector,
                             res_bar::AbstractMatrix
                            ) where {Tdim, Tsol, Tres}

  @unpack data flux nrm aux_vars work flux_bar
  fill!(nrm, 0); fill!(flux_bar, 0)

  numDofPerNode, numNodesPerElement = size(q)

  @simd for i=1:numNodesPerElement
    q_i = sview(q, :, i)
    @simd for d=1:Tdim
      nrm[d] = 1
      flux_d = sview(flux, :, i, d)

      calcEulerFlux(params, q_i, aux_vars, nrm, flux_d)

      nrm[d] = 0
    end
  end

  #applyDx(sbp, flux, dxidx, jac, work, res)

  # reverse sweep
  applyDx_revm(sbp, flux, flux_bar, dxidx, dxidx_bar, jac, jac_bar, work, res_bar)

  return nothing
end


"""
  **Inputs**

   * params
   * sbp
   * jac: mapping jacobian determinant for the element, `numNodesPerElement`
   * res: `numDofPerNode` x `numNodesPerElement` array to compute the norm
          of.
   * res_jac: 4D jacobian of `res` wrt `q`

  **Inputs/Outputs**

   * val_jac: array to be overwritten with jacobian of `val` wrt `q`

  **Outputs**

   * val: L2 norm
"""
function computeL2Norm_diff(params::ParamType{Tdim}, sbp::AbstractOperator,
                            jac::AbstractVector,
                            res::AbstractMatrix{Tres}, res_jac::Abstract4DArray,
                            val_jac::AbstractMatrix) where {Tres, Tdim}

  numDofPerNode, numNodesPerElement = size(res)

  # compute norm
  val = zero(Tres)
  fill!(val_jac, 0)
  @simd for i=1:numNodesPerElement
    fac = sbp.w[i]/jac[i]
    @simd for j=1:numDofPerNode
      val += res[j, i]*fac*res[j, i]
    end
  end

  val = sqrt(val)

  # jacobian of val (including the sqrt)
  fac_val = 1/(2*val)
  @simd for q=1:numNodesPerElement
    @simd for p=1:numNodesPerElement
      fac = sbp.w[p]/jac[p]
      @simd for j=1:size(res_jac, 2)
        @simd for i=1:size(res_jac, 1)
          val_jac[j, q] += fac_val*2*res[i, p]*fac*res_jac[i, j, p, q]
        end
      end
    end
  end

  return val
end



#------------------------------------------------------------------------------
# ShockSensorBO

function getShockSensor_diff(params::ParamType{Tdim}, sbp::AbstractOperator,
                      sensor::ShockSensorBO,
                      q::AbstractMatrix{Tsol}, elnum::Integer,
                      coords::AbstractMatrix, dxidx::Abstract3DArray,
                      jac::AbstractVector{Tmsh},
                      Se_jac::Abstract4DArray{Tres},
                      ee_jac::Abstract4DArray{Tres}) where {Tsol, Tmsh, Tres, Tdim}

  numDofPerNode, numNodesPerElement = size(q)
  dim = size(coords, 1)

  lambda_max = zero(Tsol)
  _Se_jac = zeros(Tres, numDofPerNode, numNodesPerElement)
  #lambda_max_dot = zeros(Tres, numDofPerNode, numNodesPerElement)
  h_avg = zero(Tmsh)
  @simd for i=1:numNodesPerElement
    q_i = sview(q, :, i)
    lambda_max_dot_i = sview(_Se_jac, :, i)
    lambda_max += getLambdaMax_diff(params, q_i, lambda_max_dot_i)
    h_avg += sbp.w[i]/jac[i]
  end

  lambda_max /= numNodesPerElement
  scale!(_Se_jac, 1/numNodesPerElement)
  h_avg = h_avg^(1/Tdim)

  @simd for i=1:numNodesPerElement
    @simd for j=1:numDofPerNode
      ee_jac[1, j, 1, i] = sensor.alpha*h_avg*_Se_jac[j, i]/sbp.degree
    end
  end

  # fill in the rest of the jac arrays
  @simd for q=1:numNodesPerElement
    @simd for p=1:numNodesPerElement
      @simd for j=1:numDofPerNode
        @simd for i=1:dim
          Se_jac[i, j, p, q] = _Se_jac[j, q]
          ee_jac[i, j, p, q] = ee_jac[1, j, 1, q]
        end
      end
    end
  end

  return false
end


#------------------------------------------------------------------------------
# ShockSensorHHO


function getShockSensor_diff(params::ParamType{Tdim}, sbp::AbstractOperator,
                      sensor::ShockSensorHHO,
                      q_el::AbstractMatrix{Tsol}, elnum::Integer,
                      coords::AbstractMatrix, dxidx::Abstract3DArray,
                      jac::AbstractVector{Tmsh},
                      Se_jac::Abstract4DArray{Tres},
                      ee_jac::Abstract4DArray{Tres}) where {Tsol, Tmsh, Tres, Tdim}


  numDofPerNode, numNodesPerElement = size(q_el)
  dim = size(coords, 1)

  @unpack sensor p_jac press_el press_dx work res res_jac p_hess Dx px_jac val_dot Rp fp Rp_jac fp_jac

  fill!(Rp, 0)
  fill!(Rp_jac, 0); fill!(fp_jac, 0)
  fill!(press_dx, 0); fill!(res, 0); fill!(px_jac, 0)

  computeStrongResidual_diff(params, sbp, sensor.strongdata,
                             q_el, coords, dxidx, jac, res, res_jac)

  # compute pressure and pressure_jac @simd for entire element
  @simd for p=1:numNodesPerElement
    q_p = ro_sview(q_el, :, p)
    p_dot = sview(p_jac, 1, :, p)
    press_el[1, p] = calcPressure_diff(params, q_p, p_dot)
  end

  # compute Rp
  # jacobian first: p_dot/p * dR/dq + R*d/dq (p_dot/p)
  # first term
  
  @simd for q=1:numNodesPerElement
    @simd for p=1:numNodesPerElement
      press = press_el[1, p]
      p_dot = sview(p_jac, 1, :, p)
      @simd for j=1:numDofPerNode
        @simd for i=1:numDofPerNode
          #res_jac[i, j, p, q] *= p_dot[i]/press
          Rp_jac[p, j, q] += res_jac[i, j, p, q]*p_dot[i]/press
        end
      end
    end
  end


  # second term
  @simd for p=1:numNodesPerElement
    q_p = sview(q_el, :, p)
    press = press_el[1, p]
    p_dot = sview(p_jac, 1, :, p)
    calcPressure_hess(params, q_p, p_hess)
    @simd for j=1:numDofPerNode
      @simd for i=1:numDofPerNode
        Rp_jac[p, j, p] += res[i, p]*(-p_dot[i]*p_dot[j]/(press*press) + p_hess[i, j]/press)
      end
    end
  end

  # residual itself
  @simd for i=1:numNodesPerElement
    q_i = ro_sview(q_el, :, i)
    press = press_el[1, i]
    p_dot = sview(p_jac, 1, :, i)

    @simd for j=1:numDofPerNode
      Rp[i] += p_dot[j]*res[j, i]/press
    end
  end

  # compute pressure switch fp
  calcDx(sbp, dxidx, jac, Dx)
  applyOperatorJac(Dx, p_jac, px_jac)


  applyDx(sbp, press_el, dxidx, jac, work, press_dx)
  @simd for p=1:numNodesPerElement
    val = zero(Tres)
    fill!(val_dot, 0)
    @simd for d=1:Tdim
      val += press_dx[1, p, d]*press_dx[1, p, d]
      @simd for q=1:numNodesPerElement
        @simd for j=1:numDofPerNode
          val_dot[j, q] += 2*press_dx[1, p, d]*px_jac[1, j, d, p, q]
        end
      end
    end  # end d
    press = press_el[1, p]

    val2 = sqrt(val)
    scale!(val_dot, 1/(2*val2))
  
    t1 = (press + 1e-12)
    val3 = val2/t1

    # contribution of derivative of val2
    @simd for q=1:numNodesPerElement
      @simd for j=1:numDofPerNode
        val_dot[j, q] = val_dot[j, q]/t1
      end
    end

    # contribution of derivative of p
    @simd for j=1:numDofPerNode
      val_dot[j, p] += val2*(-p_jac[1, j, p]/(t1*t1))
    end

    fp[p] = val3
    @simd for q=1:numNodesPerElement
      @simd for j=1:numDofPerNode
        fp_jac[j, p, q] = val_dot[j, q]
      end
    end
  end

  @simd for q=1:numNodesPerElement
    @simd for p=1:numNodesPerElement
      @simd for j=1:numDofPerNode
        @simd for d=1:Tdim
          h_fac = sensor.h_k_tilde[d, elnum]
          Se_jac[d, j, p, q] = h_fac*(fp[p]*sign_c(Rp[p])*Rp_jac[p, j, q] + 
                                      fp_jac[j, p, q]*absvalue(Rp[p]))
          ee_jac[d, j, p, q] = Se_jac[d, j, p, q]*h_fac*h_fac*sensor.C_eps
        end
      end
    end
  end

  return false
end


function getShockSensor_revq(params::ParamType{Tdim}, sbp::AbstractOperator,
                        sensor::ShockSensorHHO{Tsol, Tres},
                        q::AbstractMatrix, q_bar::AbstractMatrix,
                        elnum::Integer,
                        coords::AbstractMatrix,
                        dxidx::Abstract3DArray, jac::AbstractVector{Tmsh},
                        ee_mat::AbstractMatrix, ee_bar::AbstractMatrix
                       ) where {Tsol, Tres, Tmsh, Tdim}

  numDofPerNode, numNodesPerElement = size(q)

  @unpack sensor p_dot press_el press_dx work res Rp fp
  fill!(press_dx, 0); fill!(res, 0); fill!(Rp, 0)

  # compute Rm
  computeStrongResidual(params, sbp, sensor.strongdata, q, dxidx, jac, res)

  # compute Rp
  @simd for i=1:numNodesPerElement
    q_i = ro_sview(q, :, i)
    press = calcPressure_diff(params, q_i, p_dot)
    press_el[1, i] = press

    @simd for j=1:numDofPerNode
      Rp[i] += p_dot[j]*res[j, i]/press
    end
  end

  # compute pressure switch fp (excluding the h_k factor)
  applyDx(sbp, press_el, dxidx, jac, work, press_dx)
  @simd for i=1:numNodesPerElement
    val = zero(Tres)
    @simd for d=1:Tdim
      val += press_dx[1, i, d]*press_dx[1, i, d]
    end
    fp[i] = sqrt(val)/(press_el[1, i] + 1e-12)
  end

  @simd for i=1:numNodesPerElement
    @simd for d=1:Tdim
      h_fac = sensor.h_k_tilde[d, elnum]
      Se = h_fac*fp[i]*absvalue(Rp[i])
      ee_mat[d, i] = Se*h_fac*h_fac*sensor.C_eps
    end
  end


  # reverse sweep
  @unpack sensor fp_bar Rp_bar press_el_bar press_dx_bar p_dot_bar res_bar
  fill!(fp_bar, 0); fill!(Rp_bar, 0); fill!(press_el_bar, 0)
  fill!(press_dx_bar, 0); fill!(res_bar, 0)

  @simd for i=1:numNodesPerElement
    @simd for d=1:Tdim
      h_fac = sensor.h_k_tilde[d, elnum]
      Se_bar = h_fac*h_fac*sensor.C_eps*ee_bar[d, i]
      fp_bar[i] += h_fac*absvalue(Rp[i])*Se_bar
      Rp_bar[i] += h_fac*fp[i]*sign_c(Rp[i])*Se_bar
    end
  end

  # fp_bar
  @simd for i=1:numNodesPerElement
    val = zero(Tres)
    val_bar = zero(Tres)
    @simd for d=1:Tdim
      val += press_dx[1, i, d]*press_dx[1, i, d]
    end
    val_bar += fp_bar[i]/(2*sqrt(val)*(press_el[1, i] + 1e-12))
    press_el_bar[1, i] += -sqrt(val)*fp_bar[i]/( 
                          (press_el[1, i] + 1e-12)*(press_el[1, i] + 1e-12) )
    @simd for d=1:Tdim
      press_dx_bar[1, i, d] += 2*press_dx[1, i, d]*val_bar
    end
  end
  applyDx_revq(sbp, press_el_bar, dxidx, jac, work, press_dx_bar)


  # Rp_bar
   @simd for i=1:numNodesPerElement
    q_i = ro_sview(q, :, i)
    q_bar_i = sview(q_bar, :, i)
    calcPressure_diff(params, q_i, p_dot)
    press = press_el[1, i]
    fill!(p_dot_bar, 0)

    @simd for j=1:numDofPerNode
      #Rp[i] += p_dot[j]*res[j, i]/press
      p_dot_bar[j]       += res[j, i]*Rp_bar[i]/press
      res_bar[j, i]      += p_dot[j]*Rp_bar[i]/press
      press_el_bar[1, i] += -p_dot[j]*res[j, i]*Rp_bar[i]/(press*press)
    end

    calcPressure_diff_revq(params, q_i, q_bar_i, p_dot_bar, press_el_bar[1, i])
  end

  computeStrongResidual_revq(params, sbp, sensor.strongdata, q, q_bar, dxidx,
                             jac, res_bar)

  return true
end


function getShockSensor_revm(params::ParamType{Tdim}, sbp::AbstractOperator,
                        sensor::ShockSensorHHO{Tsol, Tres},
                        q::AbstractMatrix,
                        elnum::Integer,
                        coords::AbstractMatrix, coords_bar::AbstractMatrix,
                        dxidx::Abstract3DArray, dxidx_bar::Abstract3DArray,
                        jac::AbstractVector{Tmsh}, jac_bar::AbstractVector,
                        ee_bar::AbstractMatrix
                       ) where {Tsol, Tres, Tmsh, Tdim}

  numDofPerNode, numNodesPerElement = size(q)

  @unpack sensor p_dot press_el press_dx work res Rp fp
  fill!(press_dx, 0); fill!(res, 0); fill!(Rp, 0)

  # compute Rm
  computeStrongResidual(params, sbp, sensor.strongdata, q, dxidx, jac, res)

  # compute Rp
  @simd for i=1:numNodesPerElement
    q_i = ro_sview(q, :, i)
    press = calcPressure_diff(params, q_i, p_dot)
    press_el[1, i] = press

    @simd for j=1:numDofPerNode
      Rp[i] += p_dot[j]*res[j, i]/press
    end
  end

  # compute pressure switch fp (excluding the h_k factor)
  applyDx(sbp, press_el, dxidx, jac, work, press_dx)
  @simd for i=1:numNodesPerElement
    val = zero(Tres)
    @simd for d=1:Tdim
      val += press_dx[1, i, d]*press_dx[1, i, d]
    end
    fp[i] = sqrt(val)/(press_el[1, i] + 1e-12)
  end
#=
  @simd for i=1:numNodesPerElement
    @simd for d=1:Tdim
      h_fac = sensor.h_k_tilde[d, elnum]
      Se = h_fac*fp[i]*absvalue(Rp[i])
      ee_mat[d, i] = Se*h_fac*h_fac*sensor.C_eps
    end
  end
=#

  # reverse sweep
  @unpack sensor fp_bar Rp_bar press_dx_bar p_dot_bar press_el_bar res_bar
  fill!(fp_bar, 0); fill!(Rp_bar, 0); fill!(press_el_bar, 0)
  fill!(press_dx_bar, 0); fill!(res_bar, 0)

  @simd for i=1:numNodesPerElement
    @simd for d=1:Tdim
      h_fac = sensor.h_k_tilde[d, elnum]
      Se_bar = h_fac*h_fac*sensor.C_eps*ee_bar[d, i]
      fp_bar[i] += h_fac*absvalue(Rp[i])*Se_bar
      Rp_bar[i] += h_fac*fp[i]*sign_c(Rp[i])*Se_bar
      sensor.h_k_tilde_bar[d, elnum] += 3*fp[i]*absvalue(Rp[i])*Se_bar
    end
  end

  # fp_bar
  @simd for i=1:numNodesPerElement
    val = zero(Tres)
    val_bar = zero(Tres)
    @simd for d=1:Tdim
      val += press_dx[1, i, d]*press_dx[1, i, d]
    end
    val_bar += fp_bar[i]/(2*sqrt(val)*(press_el[1, i] + 1e-12))
    #press_el_bar[1, i] += -sqrt(val)*fp_bar[i]/( 
    #                      (press_el[1, i] + 1e-12)*(press_el[1, i] + 1e-12) )
    @simd for d=1:Tdim
      press_dx_bar[1, i, d] += 2*press_dx[1, i, d]*val_bar
    end
  end
  # press_el_bar unneeded, but required by applyDx_revm
  applyDx_revm(sbp, press_el, press_el_bar, dxidx, dxidx_bar, jac, jac_bar,
               work, press_dx_bar)


  # Rp_bar
   @simd for i=1:numNodesPerElement
    q_i = ro_sview(q, :, i)
    calcPressure_diff(params, q_i, p_dot)
    press = press_el[1, i]
    fill!(p_dot_bar, 0)

    @simd for j=1:numDofPerNode
      #Rp[i] += p_dot[j]*res[j, i]/press
      #p_dot_bar[j]       += res[j, i]*Rp_bar[i]/press
      res_bar[j, i]      += p_dot[j]*Rp_bar[i]/press
      #press_el_bar[1, i] += -p_dot[j]*res[j, i]*Rp_bar[i]/(press*press)
    end

    #calcPressure_diff_revq(params, q_i, q_bar_i, p_dot_bar, press_el_bar[1, i])
  end

  computeStrongResidual_revm(params, sbp, sensor.strongdata, q, dxidx,
                             dxidx_bar, jac, jac_bar, res_bar)

  return true
end


#------------------------------------------------------------------------------

function getShockSensor_diff(params::ParamType{Tdim}, sbp::AbstractOperator,
                      sensor::ShockSensorHHOConst,
                      q_el::AbstractMatrix{Tsol}, elnum::Integer,
                      coords::AbstractMatrix, dxidx::Abstract3DArray,
                      jac::AbstractVector{Tmsh},
                      Se_jac::Abstract4DArray{Tres},
                      ee_jac::Abstract4DArray{Tres}) where {Tsol, Tmsh, Tres, Tdim}


  numDofPerNode, numNodesPerElement = size(q_el)
  dim = size(coords, 1)

  @unpack sensor p_jac press_el press_dx work res res_jac p_hess Dx px_jac val_dot Rp fp  epsilon Rp_jac fp_jac epsilon_dot val2_dot

  fill!(Rp, 0)
  fill!(Rp_jac, 0); fill!(fp_jac, 0)
  fill!(press_dx, 0); fill!(res, 0); fill!(px_jac, 0)

  computeStrongResidual_diff(params, sbp, sensor.strongdata,
                             q_el, coords, dxidx, jac, res, res_jac)

  # compute pressure and pressure_jac @simd for entire element
  @simd for p=1:numNodesPerElement
    q_p = ro_sview(q_el, :, p)
    p_dot = sview(p_jac, 1, :, p)
    press_el[1, p] = calcPressure_diff(params, q_p, p_dot)
  end

  # compute Rp
  # jacobian first: p_dot/p * dR/dq + R*d/dq (p_dot/p)
  # first term
  
  @simd for q=1:numNodesPerElement
    @simd for p=1:numNodesPerElement
      press = press_el[1, p]
      p_dot = sview(p_jac, 1, :, p)
      @simd for j=1:numDofPerNode
        @simd for i=1:numDofPerNode
          #res_jac[i, j, p, q] *= p_dot[i]/press
          Rp_jac[p, j, q] += res_jac[i, j, p, q]*p_dot[i]/press
        end
      end
    end
  end


  # second term
  @simd for p=1:numNodesPerElement
    q_p = sview(q_el, :, p)
    press = press_el[1, p]
    p_dot = sview(p_jac, 1, :, p)
    calcPressure_hess(params, q_p, p_hess)
    @simd for j=1:numDofPerNode
      @simd for i=1:numDofPerNode
        Rp_jac[p, j, p] += res[i, p]*(-p_dot[i]*p_dot[j]/(press*press) + p_hess[i, j]/press)
      end
    end
  end

  # residual itself
  @simd for i=1:numNodesPerElement
    q_i = ro_sview(q_el, :, i)
    press = press_el[1, i]
    p_dot = sview(p_jac, 1, :, i)

    @simd for j=1:numDofPerNode
      Rp[i] += p_dot[j]*res[j, i]/press
    end
  end

  # compute pressure switch fp
  calcDx(sbp, dxidx, jac, Dx)
  applyOperatorJac(Dx, p_jac, px_jac)


  applyDx(sbp, press_el, dxidx, jac, work, press_dx)
  @simd for p=1:numNodesPerElement
    val = zero(Tres)
    fill!(val_dot, 0)
    @simd for d=1:Tdim
      val += press_dx[1, p, d]*press_dx[1, p, d]
      @simd for q=1:numNodesPerElement
        @simd for j=1:numDofPerNode
          val_dot[j, q] += 2*press_dx[1, p, d]*px_jac[1, j, d, p, q]
        end
      end
    end  # end d
    press = press_el[1, p]

    val2 = sqrt(val)
    scale!(val_dot, 1/(2*val2))
  
    t1 = (press + 1e-12)
    val3 = val2/t1

    # contribution of derivative of val2
    @simd for q=1:numNodesPerElement
      @simd for j=1:numDofPerNode
        val_dot[j, q] = val_dot[j, q]/t1
      end
    end

    # contribution of derivative of p
    @simd for j=1:numDofPerNode
      val_dot[j, p] += val2*(-p_jac[1, j, p]/(t1*t1))
    end

    fp[p] = val3
    @simd for q=1:numNodesPerElement
      @simd for j=1:numDofPerNode
        fp_jac[j, p, q] = val_dot[j, q]
      end
    end
  end

  @simd for i=1:numNodesPerElement
    epsilon[1, i] = fp[i]*Rp[i]
  end

  @simd for q=1:numNodesPerElement
    @simd for p=1:numNodesPerElement
      @simd for j=1:numDofPerNode
        epsilon_dot[1, j, p, q] = fp[p]*Rp_jac[p, j, q] + Rp[p]*fp_jac[j, p, q]
      end
    end
  end

  #val = computeL2Norm_diff(params, sbp, jac, epsilon, epsilon_dot, val2_dot)
  fill!(val2_dot, 0)
  @simd for q=1:numNodesPerElement
    @simd for p=1:numNodesPerElement
      @simd for j=1:numDofPerNode
        val2_dot[j, q] += absvalue_deriv(epsilon[1, p])*sbp.w[p]*epsilon_dot[1, j, p, q]
      end
    end
  end

  @simd for q=1:numNodesPerElement
    @simd for p=1:numNodesPerElement
      @simd for j=1:numDofPerNode
        @simd for d=1:Tdim
          h_fac = sensor.h_k_tilde[d, elnum]
          Se_jac[d, j, p, q] = h_fac*val2_dot[j, q]
          ee_jac[d, j, p, q] = Se_jac[d, j, p, q]*h_fac*h_fac*sensor.C_eps
        end
      end
    end
  end

  return false
end



function calcAnisoFactors_revm(mesh::AbstractMesh, sbp, opts,
                          hk_all::AbstractMatrix{T}, hk_all_bar) where {T}

  # compute h_k tilde = h_k/(p + 1), where h_k takes into account the
  # anisotropy of the element in directory k
  # compute the p_k first, then compute h_k, where p_k is the integral of the
  # face normal
  @assert size(hk_all) == (mesh.dim, mesh.numEl)

  fill!(hk_all, 0)

  # do all faces (interior, boundary, shared)
  @simd for i=1:mesh.numInterfaces
    iface_i = mesh.interfaces[i]
    nrm_i = ro_sview(mesh.nrm_face, :, :, i)

    @simd for j=1:mesh.numNodesPerFace
      @simd for k=1:mesh.dim
        fac = mesh.sbpface.wface[j]*absvalue2(nrm_i[k, j])
        hk_all[k, iface_i.elementL] += fac
        hk_all[k, iface_i.elementR] += fac
      end
    end
  end

  @simd for i=1:mesh.numBoundaryFaces
    bndry_i = mesh.bndryfaces[i]
    nrm_i = ro_sview(mesh.nrm_bndry, :, :, i)

    @simd for j=1:mesh.numNodesPerFace
      @simd for k=1:mesh.dim
        fac = mesh.sbpface.wface[j]*absvalue2(nrm_i[k, j])
        hk_all[k, bndry_i.element] += fac
      end
    end
  end

  @simd for peer=1:mesh.npeers
    @simd for i=1:length(mesh.shared_interfaces[peer])
      iface_i = mesh.shared_interfaces[peer]
      nrm_i = ro_sview(mesh.nrm_sharedface[peer], :, :, i)

      @simd for j=1:mesh.numNodesPerFace
        @simd for k=1:mesh.dim
          fac = mesh.sbpface.wface[j]*absvalue2(nrm_i[k, j])
          hk_all[k, iface_i.elementL] += fac
        end
      end
    end
  end

# don't do this now because the reverse sweep needs hk_all to have the
# same values as during the same time during the forward calculation
#=
  # hk_all now contains the p_i in each direction
  # Now compute h_k tilde from p_i
  for i=1:mesh.numEl

    jac_i = ro_sview(mesh.jac, :, i)
    vol = zero(T)
    for j=1:mesh.numNodesPerElement
      vol += sbp.w[j]/jac_i[j]
    end

    h_k = vol^(1/mesh.dim)
    pk_prod = one(T)
    for d=1:mesh.dim
      pk_prod *= hk_all[d, i]
    end

    for d=1:mesh.dim
      hk_all[d, i] = h_k*pk_prod/(hk_all[d, i]*(sbp.degree + 1))
    end
  end
=#

  # reverse sweep
  @simd for i=1:mesh.numEl

    jac_i = ro_sview(mesh.jac, :, i)
    vol = zero(T)
    @simd for j=1:mesh.numNodesPerElement
      vol += sbp.w[j]/jac_i[j]
    end

    h_k = vol^(1/mesh.dim)
    pk_prod = one(T)
    @simd for d=1:mesh.dim
      pk_prod *= hk_all[d, i]
    end
    pk_prod_d = pk_prod^(1/mesh.dim)

    pk_prod_d_bar = zero(T)
    h_k_bar = zero(T)
    @simd for d=1:mesh.dim
      #hk_all[d, i] = h_k*pk_prod_d/(hk_all[d, i]*(sbp.degree + 1))
      tmp = hk_all[d, i]*(sbp.degree + 1)
      h_k_bar += pk_prod_d*hk_all_bar[d, i]/tmp
      pk_prod_d_bar += h_k*hk_all_bar[d, i]/tmp
      hk_all_bar[d, i] = -h_k*pk_prod_d*hk_all_bar[d, i]/(hk_all[d, i]*tmp)
    end


    pk_prod_bar = (pk_prod_d_bar/mesh.dim)*(pk_prod)^((1/mesh.dim) - 1)
    @simd for d=mesh.dim:-1:1
      pk_prod /= hk_all[d, i]  # reset primal value
      hk_all_bar[d, i] += pk_prod*pk_prod_bar
      pk_prod_bar = hk_all[d, i]*pk_prod_bar
    end

    jac_bar_i = sview(mesh.jac_bar, :, i)
    vol_bar = (h_k_bar/mesh.dim)*vol^((1/mesh.dim) - 1)
    @simd for j=1:mesh.numNodesPerElement
      jac_bar_i[j] += -sbp.w[j]*vol_bar/(jac_i[j]*jac_i[j])
    end

  end  # end i

  @simd for i=1:mesh.numInterfaces
    iface_i = mesh.interfaces[i]
    nrm_i = ro_sview(mesh.nrm_face, :, :, i)
    nrm_bar_i = sview(mesh.nrm_face_bar, :, :, i)

    @simd for j=1:mesh.numNodesPerFace
      @simd for k=1:mesh.dim
        fac = absvalue2_deriv(nrm_i[k, j])*mesh.sbpface.wface[j]
        nrm_bar_i[k, j] += fac*hk_all_bar[k, iface_i.elementL]
        nrm_bar_i[k, j] += fac*hk_all_bar[k, iface_i.elementR]
      end
    end
  end

  @simd for i=1:mesh.numBoundaryFaces
    bndry_i = mesh.bndryfaces[i]
    nrm_i = ro_sview(mesh.nrm_bndry, :, :, i)
    nrm_bar_i = sview(mesh.nrm_bndry_bar, :, :, i)

    @simd for j=1:mesh.numNodesPerFace
      @simd for k=1:mesh.dim
        fac = mesh.sbpface.wface[j]*absvalue2_deriv(nrm_i[k, j])
        nrm_bar_i[k, j] += fac*hk_all_bar[k, bndry_i.element]
      end
    end
  end

  @simd for peer=1:mesh.npeers
    @simd for i=1:length(mesh.shared_interfaces[peer])
      iface_i = mesh.shared_interfaces[peer]
      nrm_i = ro_sview(mesh.nrm_sharedface[peer], :, :, i)
      nrm_bar_i = sview(mesh.nrm_sharedface_bar[peer], :, :, i)

      @simd for j=1:mesh.numNodesPerFace
        @simd for k=1:mesh.dim
          fac = mesh.sbpface.wface[j]*absvalue2_deriv(nrm_i[k, j])
          nrm_bar_i[k, j] += fac*hk_all_bar[k, iface_i.elementL]
        end
      end
    end
  end

  # finish the forward computation, so hk_all has the right values in it on
  # exit.
  # hk_all now contains the p_i in each direction
  # Now compute h_k tilde from p_i
  @simd for i=1:mesh.numEl

    jac_i = ro_sview(mesh.jac, :, i)
    vol = zero(T)
    @simd for j=1:mesh.numNodesPerElement
      vol += sbp.w[j]/jac_i[j]
    end

    h_k = vol^(1/mesh.dim)
    pk_prod = one(T)
    @simd for d=1:mesh.dim
      pk_prod *= hk_all[d, i]
    end
    pk_prod = pk_prod^(1/mesh.dim)

    @simd for d=1:mesh.dim
      hk_all[d, i] = h_k*pk_prod/(hk_all[d, i]*(sbp.degree + 1))
    end
  end

  return nothing
end


#------------------------------------------------------------------------------
# SensorHApprox

# this shock sensor is not a function of q, so the _diff and _revq methods
# are not required

# even though this shared the name with all the other getShockSensor methods,
# it signature and outputs are different.  This is sort-of ok because this
# is a special shock sensor and should never be used as the main shock sensor.
function getShockSensor_revm(params::ParamType{Tdim}, sbp::AbstractOperator,
                        sensor::ShockSensorHApprox{Tsol, Tres},
                        elnum::Integer, jac::AbstractVector{Tmsh},
                        jac_bar::AbstractVector, epsL_bar::Number
                       ) where {Tsol, Tres, Tmsh, Tdim}


  numNodesPerElement = length(jac)

  shared_idx = elnum - first(sensor.shared_els) + 1  # index of shared element
  if (elnum in sensor.local_els) ||
     (elnum in sensor.shared_els && sensor.shared_isShocked[shared_idx])

    # compute element size h for both elements
    h = zero(Tmsh)
    @simd for i=1:numNodesPerElement
      h += sbp.w[i]/jac[i]
    end

    # compute estimates viscosity
    lambda_max = sensor.lambda_max
    fac = sensor.alpha*lambda_max/sbp.degree
    #eps_L = fac*(h^(1/Tdim))

    #----------------------------------------
    # reverse sweep
    
    h_bar = fac*epsL_bar*(h^(1/Tdim - 1))/Tdim

    @simd for i=1:numNodesPerElement
      jac_bar[i] += -sbp.w[i]*h_bar/(jac[i]*jac[i])
    end
  end

end


