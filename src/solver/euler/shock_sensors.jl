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
"""
function getShockSensor(params::ParamType{Tdim}, sbp::AbstractOperator,
                          sensor::ShockSensorPP,
                          q::AbstractMatrix{Tsol}, elnum::Integer,
                          coords::AbstractMatrix,
                          dxidx::Abstract3DArray, jac::AbstractVector{Tmsh},
                          Se_mat::AbstractMatrix, ee_mat::AbstractMatrix
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

  ee *= lambda_max*h_avg/sbp.degree


  # use elementwise constant
  fill!(Se_mat, Se)
  fill!(ee_mat, ee)
  
  return is_nonzero
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
                          q::AbstractMatrix{Tsol}, elnum::Integer,
                          coords::AbstractMatrix,
                          dxidx::Abstract3DArray, jac::AbstractVector{Tmsh},
                          Se_mat::AbstractMatrix, ee_mat::AbstractMatrix
                         ) where {Tsol, Tmsh}

  error("getShockSensor called for ShockSensorNone: did you forget to specify the shock capturing scheme?")
end


function isShockElement(params::ParamType, sbp::AbstractOperator,
                          sensor::ShockSensorNone,
                          q::AbstractMatrix{Tsol}, elnum::Integer,
                          coords::AbstractMatrix,
                          dxidx::Abstract3DArray, jac::AbstractVector{Tmsh},
                         ) where {Tsol, Tmsh}

  error("isShockElement called for ShockSensorNone: did you forget to specify the shock capturing scheme?")
end


#------------------------------------------------------------------------------
# ShockSensorEverywhere

function getShockSensor(params::ParamType, sbp::AbstractOperator,
                        sensor::ShockSensorEverywhere{Tsol, Tres},
                        q::AbstractMatrix, elnum::Integer,
                        coords::AbstractMatrix,
                        dxidx::Abstract3DArray, jac::AbstractVector{Tmsh},
                        Se_mat::AbstractMatrix, ee_mat::AbstractMatrix
                        ) where {Tsol, Tres, Tmsh}

  fill!(Se_mat, 1.0)
  fill!(ee_mat, 1.0)

  return true
end

function isShockElement(params::ParamType, sbp::AbstractOperator,
                        sensor::ShockSensorEverywhere{Tsol, Tres},
                        q::AbstractMatrix, elnum::Integer,
                        coords::AbstractMatrix,
                        dxidx::Abstract3DArray, jac::AbstractVector{Tmsh},
                        ) where {Tsol, Tres, Tmsh}


  return true
end

#------------------------------------------------------------------------------
# ShockSensorVelocity

function getShockSensor(params::ParamType, sbp::AbstractOperator,
                        sensor::ShockSensorVelocity{Tsol, Tres},
                        q::AbstractMatrix, elnum::Integer,
                        coords::AbstractMatrix,
                        dxidx::Abstract3DArray, jac::AbstractVector{Tmsh},
                        Se_mat::AbstractMatrix, ee_mat::AbstractMatrix
                        ) where {Tsol, Tres, Tmsh}


  numNodesPerElement = size(q, 2)
  dim = size(coords, 1)

  for i=1:numNodesPerElement
    for d=1:dim
      Se_mat[d, i] = d*q[d+1, i]
      ee_mat[d, i] = d*q[d+1, i]
    end
  end

  return true
end

function isShockElement(params::ParamType, sbp::AbstractOperator,
                        sensor::ShockSensorVelocity{Tsol, Tres},
                        q::AbstractMatrix, elnum::Integer,
                        coords::AbstractMatrix,
                        dxidx::Abstract3DArray, jac::AbstractVector{Tmsh},
                        ) where {Tsol, Tres, Tmsh}

  return true
end

#------------------------------------------------------------------------------
# ShockSensorHIso

function getShockSensor(params::ParamType{Tdim}, sbp::AbstractOperator,
                        sensor::ShockSensorHIso{Tsol, Tres},
                        q::AbstractMatrix, elnum::Integer,
                        coords::AbstractMatrix,
                        dxidx::Abstract3DArray, jac::AbstractVector{Tmsh},
                        Se_mat::AbstractMatrix, ee_mat::AbstractMatrix
                         ) where {Tsol, Tres, Tmsh, Tdim}

  # compute | div(F) |
  # Do this in Cartesian coordinates because it makes the differentiation easier
  numDofPerNode, numNodesPerElement = size(q)
  res = sensor.res

  computeStrongResidual(params, sbp, sensor.strongdata, q, dxidx, jac, res)
  # Note: Hartmann's paper referrs to Jaffre's paper which says to use the
  #       maximum value of res over the element, but that is not differentiable,
  #       so take the L2 norm instead
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

  fill!(Se_mat, val)
  fill!(ee_mat, ee)

  # calling | div(f) | as the shock sensor, which is somewhat arbitrary
  return true
end

function isShockElement(params::ParamType{Tdim}, sbp::AbstractOperator,
                        sensor::ShockSensorHIso{Tsol, Tres},
                        q::AbstractMatrix, elnum::Integer,
                        coords::AbstractMatrix,
                        dxidx::Abstract3DArray, jac::AbstractVector{Tmsh},
                       ) where {Tsol, Tres, Tmsh, Tdim}

  return true
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
                        q::AbstractMatrix, elnum::Integer,
                        coords::AbstractMatrix,
                        dxidx::Abstract3DArray, jac::AbstractVector{Tmsh},
                        Se_mat::AbstractMatrix, ee_mat::AbstractMatrix
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

  fill!(Se_mat, lambda_max)
  fill!(ee_mat, sensor.alpha*lambda_max*h_avg/sbp.degree)

  return true
end

function isShockElement(params::ParamType{Tdim}, sbp::AbstractOperator,
                        sensor::ShockSensorBO{Tsol, Tres},
                        q::AbstractMatrix, elnum::Integer,
                        coords::AbstractMatrix,
                        dxidx::Abstract3DArray, jac::AbstractVector{Tmsh},
                       ) where {Tsol, Tres, Tmsh, Tdim}

  return true
end

#------------------------------------------------------------------------------
# ShockSensorHHO



function getShockSensor(params::ParamType{Tdim}, sbp::AbstractOperator,
                        sensor::ShockSensorHHO{Tsol, Tres},
                        q::AbstractMatrix, elnum::Integer,
                        coords::AbstractMatrix,
                        dxidx::Abstract3DArray, jac::AbstractVector{Tmsh},
                        Se_mat::AbstractMatrix, ee_mat::AbstractMatrix
                       ) where {Tsol, Tres, Tmsh, Tdim}

  numDofPerNode, numNodesPerElement = size(q)

  @unpack sensor p_dot press_el press_dx work res Rp fp
  fill!(press_dx, 0); fill!(res, 0); fill!(Rp, 0)

  # compute Rm
  computeStrongResidual(params, sbp, sensor.strongdata, q, dxidx, jac, res)

  # compute Rp
  for i=1:numNodesPerElement
    q_i = ro_sview(q, :, i)
    press = calcPressure_diff(params, q_i, p_dot)
    press_el[1, i] = press

    for j=1:numDofPerNode
      Rp[i] += p_dot[j]*res[j, i]/press
    end
  end

  # compute pressure switch fp (excluding the h_k factor)
  applyDx(sbp, press_el, dxidx, jac, work, press_dx)
  for i=1:numNodesPerElement
    val = zero(Tres)
    for d=1:Tdim
      val += press_dx[1, i, d]*press_dx[1, i, d]
    end
    fp[i] = sqrt(val)/(press_el[1, i] + 1e-12)
  end

  for i=1:numNodesPerElement
    for d=1:Tdim
      h_fac = sensor.h_k_tilde[d, elnum]
      Se_mat[d, i] = h_fac*fp[i]*absvalue(Rp[i])
      ee_mat[d, i] = Se_mat[d, i]*h_fac*h_fac*sensor.C_eps
    end
  end

  return true
end


#------------------------------------------------------------------------------
# ShockSensorHHOConst

function getShockSensor(params::ParamType{Tdim}, sbp::AbstractOperator,
                        sensor::ShockSensorHHOConst{Tsol, Tres},
                        q::AbstractMatrix, elnum::Integer,
                        coords::AbstractMatrix,
                        dxidx::Abstract3DArray, jac::AbstractVector{Tmsh},
                        Se_mat::AbstractMatrix, ee_mat::AbstractMatrix
                       ) where {Tsol, Tres, Tmsh, Tdim}

  numDofPerNode, numNodesPerElement = size(q)

  @unpack sensor p_dot press_el press_dx work res Rp fp epsilon
  fill!(press_dx, 0); fill!(res, 0); fill!(Rp, 0)

  # compute Rm
  computeStrongResidual(params, sbp, sensor.strongdata, q, dxidx, jac, res)

  # compute Rp
  for i=1:numNodesPerElement
    q_i = ro_sview(q, :, i)
    press = calcPressure_diff(params, q_i, p_dot)
    press_el[1, i] = press

    for j=1:numDofPerNode
      Rp[i] += p_dot[j]*res[j, i]/press
    end
  end

  # compute pressure switch fp (excluding the h_k factor)
  applyDx(sbp, press_el, dxidx, jac, work, press_dx)
  for i=1:numNodesPerElement
    val = zero(Tres)
    for d=1:Tdim
      val += press_dx[1, i, d]*press_dx[1, i, d]
    end
    fp[i] = sqrt(val)/(press_el[1, i] + 1e-12)
  end

  # compute L2 norm of fp * Rp
  for i=1:numNodesPerElement
    epsilon[1, i] = fp[i]*Rp[i]
  end

  val = computeL2Norm(params, sbp, jac, epsilon)

  for i=1:numNodesPerElement
    for d=1:Tdim
      h_fac = sensor.h_k_tilde[d, elnum]
      Se_mat[d, i] = h_fac*val
      ee_mat[d, i] = Se_mat[d, i]*h_fac*h_fac*sensor.C_eps
    end
  end

  return true
end



function isShockElement(params::ParamType{Tdim}, sbp::AbstractOperator,
                        sensor::HHOSensors{Tsol, Tres},
                        q::AbstractMatrix, elnum::Integer,
                        coords::AbstractMatrix,
                        dxidx::Abstract3DArray, jac::AbstractVector{Tmsh},
                       ) where {Tsol, Tres, Tmsh, Tdim}
  return true
end


"""
  Compute the h_k tilde factors from the paper
"""
function calcAnisoFactors(mesh::AbstractMesh, sbp, opts,
                          hk_all::AbstractMatrix{T}) where {T}

  # compute h_k tilde = h_k/(p + 1), where h_k takes into account the
  # anisotropy of the element in directory k
  # compute the p_k first, then compute h_k, where p_k is the integral of the
  # face normal

  @assert size(hk_all) == (mesh.dim, mesh.numEl)

  fill!(hk_all, 0)

  # do all faces (interior, boundary, shared)
  for i=1:mesh.numInterfaces
    iface_i = mesh.interfaces[i]
    nrm_i = ro_sview(mesh.nrm_face, :, :, i)

    for j=1:mesh.numNodesPerFace
      for k=1:mesh.dim
        fac = mesh.sbpface.wface[j]*absvalue2(nrm_i[k, j])
        hk_all[k, iface_i.elementL] += fac
        hk_all[k, iface_i.elementR] += fac
      end
    end
  end


  for i=1:mesh.numBoundaryFaces
    bndry_i = mesh.bndryfaces[i]
    nrm_i = ro_sview(mesh.nrm_bndry, :, :, i)

    for j=1:mesh.numNodesPerFace
      for k=1:mesh.dim
        fac = mesh.sbpface.wface[j]*absvalue2(nrm_i[k, j])
        hk_all[k, bndry_i.element] += fac
      end
    end
  end

  for peer=1:mesh.npeers
    for i=1:length(mesh.shared_interfaces[peer])
      iface_i = mesh.shared_interfaces[peer]
      nrm_i = ro_sview(mesh.nrm_sharedface[peer], :, :, i)

      for j=1:mesh.numNodesPerFace
        for k=1:mesh.dim
          fac = mesh.sbpface.wface[j]*absvalue2(nrm_i[k, j])
          hk_all[k, iface_i.elementL] += fac
        end
      end
    end
  end

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
    pk_prod = pk_prod^(1/mesh.dim)

    for d=1:mesh.dim
      hk_all[d, i] = h_k*pk_prod/(hk_all[d, i]*(sbp.degree + 1))
    end
  end

  return nothing
end






