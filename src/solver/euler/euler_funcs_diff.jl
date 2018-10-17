# differentiated version of functions in euler_funcs.jl

"""
  Applies the element level inverse mass matrix to an element-level jacobian
  contribution

  **Inputs**

   * jac_el: vector, length `numNodesPerElement`, containing the determinant
              of the mapping jacobian at each node of the element
   * weights: the SBP integration weights for each node of the element,
               length `numNodesPerElement

  **Inputs/Outputs**

   * res_jac: array containing element level Jacobian, to be multiplied by
              the inverse mass matrix, size `numDofPerNode` x `numDofPerNode`
              x `numNodesPerElement` x `numNodesPerElement`
"""
function applyMinvElement(jac_el::AbstractVector, weights::AbstractVector,
                          res_jac::AbstractArray{T, 4}) where {T}

  # multiply by Minv if needed
  for q=1:size(res_jac, 4)
    for p=1:size(res_jac, 3)
      val = jac_el[p]/weights[p]  # entry in Minv
      @simd for m=1:size(res_jac, 2)
        @simd for n=1:size(res_jac, 1)
          res_jac[n, m, p, q] *= val
        end
      end
    end
  end

  return nothing
end


"""
  Computes the derivataive of the volume term of the Roe scheme with respect
  to `q`, ie

  d/dq (Q^T * f)

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * assembler: used to assemble the contribution of each element into the
                Jacobian
"""
function calcVolumeIntegrals_nopre_diff(
                                   mesh::AbstractMesh{Tmsh},
                                   sbp::AbstractSBP,
                                   eqn::EulerData{Tsol, Tres, Tdim},
                                   opts,
                                   assembler::AssembleElementData) where {Tmsh, Tsol, Tres, Tdim}



  data = eqn.params.calc_volume_integrals_data
  @unpack data flux_jac res_jac nrm

  for i=1:mesh.numEl
    fill!(flux_jac, 0)
    fill!(res_jac, 0)
    for j=1:mesh.numNodesPerElement
      q_j = ro_sview(eqn.q, :, j, i)
      aux_vars_j = ro_sview(eqn.aux_vars, :, j, i)

      # compute dF/dq
      for k=1:Tdim
        fluxjac_k = sview(flux_jac, :, :, j, k)

        # get the direction vector
        for p=1:Tdim
          nrm[p] = mesh.dxidx[k, p, j, i]
        end
        calcEulerFlux_diff(eqn.params, q_j, aux_vars_j, nrm, fluxjac_k)
      end  # end loop k
    end  # end loop j

    # compute dR/dq
    for k=1:Tdim
      weakDifferentiateElement_jac!(sbp, k, sview(flux_jac, :, :, :, k), res_jac, SummationByParts.Add(), true)
    end

    
    if eqn.params.use_Minv == 1
      jac_el = sview(mesh.jac, :, i)
      applyMinvElement(jac_el, sbp.w, res_jac)
    end

    # assemble element level jacobian into the residual
    assembleElement(assembler, mesh, i, res_jac)
#    fill!(res_jac, 0.0)
    # flux_jac gets overwritten, so no need to zero it 

  end  # end loop i

  return nothing
end  # end function

"""
  Computes the derivative of the strong form volume terms with
  respect to `q`, ie.

  d/dq (-Q * f)

  but only the mesh.numDofPerNode x mesh.numDofPerNode diagonal block

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * assembler: used to assemble the contribution of each element into the
                Jacobian
  
"""
function calcVolumeIntegralsStrong_nopre_diff(
                                   mesh::AbstractMesh{Tmsh},
                                   sbp::AbstractSBP,
                                   eqn::EulerData{Tsol, Tres, Tdim},
                                   opts,
                                   assembler::AssembleElementData) where {Tmsh, Tsol, Tres, Tdim}


  @assert eqn.params.use_Minv != 1  # use_Minv not supported (TODO: why?)

  data = eqn.params.calc_volume_integrals_data
  @unpack data flux_jac res_jac nrm

  for i=1:mesh.numEl
    fill!(flux_jac, 0)
    fill!(res_jac, 0)
    for j=1:mesh.numNodesPerElement
      q_j = ro_sview(eqn.q, :, j, i)
      aux_vars_j = ro_sview(eqn.aux_vars, :, j, i)

      # compute dF/dq
      for k=1:Tdim
        fluxjac_k = sview(flux_jac, :, :, j, k)

        # get the direction vector
        for p=1:Tdim
          nrm[p] = mesh.dxidx[k, p, j, i]
        end
        calcEulerFlux_diff(eqn.params, q_j, aux_vars_j, nrm, fluxjac_k)
      end  # end loop k
    end  # end loop j

    # compute dR/dq
    for k=1:Tdim
      weakDifferentiateElement_jac!(sbp, k, sview(flux_jac, :, :, :, k), res_jac, SummationByParts.Subtract(), false)
    end

    
    if eqn.params.use_Minv == 1
      jac_el = sview(mesh.jac, :, i)
      applyMinvElement(jac_el, sbp.w, res_jac)
    end

    # assemble element level jacobian into the residual
    assembleElement(assembler, mesh, i, res_jac)

  end  # end loop i

  return nothing
end  # end function

"""
  Differentiated version of [`calcVolumeIntegralsSplitForm`](@ref)

  **Inputs**
   * mesh
   * sbp
   * eqn
   * opts
   * functor: the numerical flux function F, of type FluxType_diff
"""
function calcVolumeIntegralsSplitForm_diff(
                mesh::AbstractMesh{Tmsh},
                sbp::AbstractSBP,
                eqn::EulerData{Tsol, Tres, Tdim}, opts,
                functor::FluxType_diff,
                assembler::AssembleElementData) where {Tmsh, Tsol, Tres, Tdim}

  if opts["use_staggered_grid"]
    error("Jacobian of staggered grid volume integrals not supported")
  else
    # do the curvilinear version even if msh is linear
    calcVolumeIntegralsSplitFormCurvilinear_diff(mesh, sbp, eqn, opts, functor,
                                                 assembler)
  end

  return nothing
end



"""
  Computes the Jacobian contribution from [`calcVolumeIntegralsSplitFormCurvilinear`](@ref).
  
  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * functor: [`FluxType_diff`](@ref). 
   * assembler: used to assemble the contribution of each element into the
                Jacobian
 
"""
function calcVolumeIntegralsSplitFormCurvilinear_diff(
                mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP,
                eqn::EulerData{Tsol, Tres, Tdim}, opts,
                functor::FluxType_diff,
                assembler::AssembleElementData) where {Tmsh, Tsol, Tres, Tdim}

  dxidx = mesh.dxidx
  res = eqn.res
  aux_vars = eqn.aux_vars
  params = eqn.params

  data = params.calc_volume_integrals_data
  @unpack data nrmD F_d Sx res_jac

  flux_jacL = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode, Tdim)
  flux_jacR = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode, Tdim)

  # S is calculated in x-y-z, so the normal vectors should be the unit normals
  fill!(nrmD, 0.0)
  for d=1:Tdim
    nrmD[d, d] = 1
  end

  for i=1:mesh.numEl
    # get S for this element
    dxidx_i = ro_sview(dxidx, :, :, :, i)
    calcSCurvilinear(sbp, dxidx_i, Sx)

    fill!(res_jac, 0.0)
    for j=1:mesh.numNodesPerElement
      q_j = ro_sview(eqn.q, :, j, i)
      aux_vars_j = ro_sview(aux_vars, :, j, i)
      for k=1:(j-1)  # loop over lower triangle of S
        q_k = ro_sview(eqn.q, :, k, i)

        # calculate the numerical flux functions in all Tdim
        # directions at once
        fill!(flux_jacL, 0.0)
        fill!(flux_jacR, 0.0)
        functor(params, q_j, q_k, aux_vars_j, nrmD, flux_jacL, flux_jacR)

        @simd for d=1:Tdim
          @simd for p=1:mesh.numDofPerNode
            @simd for q=1:mesh.numDofPerNode
              res_jac[q, p, j, j] -= 2*Sx[j, k, d]*flux_jacL[q, p, d]
              res_jac[q, p, j, k] -= 2*Sx[j, k, d]*flux_jacR[q, p, d]

              res_jac[q, p, k, j] += 2*Sx[j, k, d]*flux_jacL[q, p, d]
              res_jac[q, p, k, k] += 2*Sx[j, k, d]*flux_jacR[q, p, d]
            end
          end
        end

      end  # end k loop
    end  # end j loop

    if eqn.params.use_Minv == 1
      jac_el = sview(mesh.jac, :, i)
      applyMinvElement(jac_el, sbp.w, res_jac)
    end

    # assemble element level jacobian into the residual
    assembleElement(assembler, mesh, i, res_jac)
  end  # end i loop

  return nothing
end





"""
  Computes the jacobian of [`calcEulerFlux`](@ref) with respect to `q`.
  Methods are available for 2D and 3D

  The caller must zero out the output array (if required)

  **Inputs**

   * params: ParamType, conservative variables only
   * q: vector of conservative variables at node
   * aux_vars: auxiliary variables at the node
   * dir: direction vector (possibly scaled) to compute the flux jacobian in

  **Inputs/Outputs**

   * Fjac: flux jacobian, numDofPerNode x numDofPerNode, summed into

"""
function calcEulerFlux_diff(params::ParamType{2, :conservative},
                      q::AbstractArray{Tsol,1},
                      aux_vars::AbstractArray{Tres, 1},
                      dir::AbstractArray{Tmsh},  Fjac::AbstractArray{Tsol,2}) where {Tmsh, Tsol, Tres}
# calculates the Euler flux in a particular direction at a point
# eqn is the equation type
# q is the vector (of length 4), of the conservative variables at the point
# aux_vars is the vector of auxiliary variables at the point
# dir is a vector of length 2 that specifies the direction
# F is populated with the flux Jacobian
# 2D  only


  p_dot = params.eulerfluxdata.p_dot
  press = calcPressure_diff(params, q, p_dot)
#  press = getPressure(aux_vars)
#  press = @getPressure(aux_vars)
  fac = 1/q[1]
  U = (q[2]*dir[1] + q[3]*dir[2])*fac
  U_dot1 = -(q[2]*dir[1] + q[3]*dir[2])*fac*fac
  U_dot2 = dir[1]*fac
  U_dot3 = dir[2]*fac

  # F[1] = q[1]*U
  # F[2] = q[2]*U + dir[1]*press
  # F[3] = q[3]*U + dir[2]*press
  # F[4] = (q[4] + press)*U
  Fjac[1, 1] += U + q[1]*U_dot1
  Fjac[2, 1] +=     q[2]*U_dot1 + dir[1]*p_dot[1]
  Fjac[3, 1] +=     q[3]*U_dot1 + dir[2]*p_dot[1]
  Fjac[4, 1] +=     q[4]*U_dot1 + press*U_dot1 + U*p_dot[1]

  Fjac[1, 2] +=     q[1]*U_dot2
  Fjac[2, 2] += U + q[2]*U_dot2 + dir[1]*p_dot[2]
  Fjac[3, 2] +=     q[3]*U_dot2 + dir[2]*p_dot[2]
  Fjac[4, 2] +=     q[4]*U_dot2 + press*U_dot2 + U*p_dot[2]

  Fjac[1, 3] +=     q[1]*U_dot3
  Fjac[2, 3] +=     q[2]*U_dot3 + dir[1]*p_dot[3]
  Fjac[3, 3] += U + q[3]*U_dot3 + dir[2]*p_dot[3]
  Fjac[4, 3] +=     q[4]*U_dot3 + press*U_dot3 + U*p_dot[3]

  Fjac[1, 4] += 0
  Fjac[2, 4] += dir[1]*p_dot[4]
  Fjac[3, 4] += dir[2]*p_dot[4]
  Fjac[4, 4] += U + U*p_dot[4]

  return nothing

end


function calcEulerFlux_diff(params::ParamType{3},
                      q::AbstractArray{Tsol,1},
                      aux_vars::AbstractArray{Tres, 1},
                      dir::AbstractArray{Tmsh},  Fjac::AbstractArray{Tsol,2}) where {Tmsh, Tsol, Tres}
# calculates the Euler flux in a particular direction at a point
# eqn is the equation type
# q is the vector (of length 4), of the conservative variables at the point
# aux_vars is the vector of auxiliary variables at the point
# dir is a vector of length 2 that specifies the direction
# F is populated with the flux Jacobian
# 2D  only


  p_dot = params.eulerfluxdata.p_dot
  press = calcPressure_diff(params, q, p_dot)
#  press = getPressure(aux_vars)
#  press = @getPressure(aux_vars)
  fac = 1/q[1]
  U = (q[2]*dir[1] + q[3]*dir[2] + q[4]*dir[3])*fac
  U_dot1 = -(q[2]*dir[1] + q[3]*dir[2] + q[4]*dir[3])*fac*fac
  U_dot2 = dir[1]*fac
  U_dot3 = dir[2]*fac
  U_dot4 = dir[3]*fac

  # F[1] = q[1]*U
  # F[2] = q[2]*U + dir[1]*press
  # F[3] = q[3]*U + dir[2]*press
  # F[4] = q[4]*U + dir[3]*press
  # F[4] = (q[5] + press)*U
  Fjac[1, 1] += U + q[1]*U_dot1
  Fjac[2, 1] +=     q[2]*U_dot1 + dir[1]*p_dot[1]
  Fjac[3, 1] +=     q[3]*U_dot1 + dir[2]*p_dot[1]
  Fjac[4, 1] +=     q[4]*U_dot1 + dir[3]*p_dot[1]
  Fjac[5, 1] +=     q[5]*U_dot1 + press*U_dot1 + U*p_dot[1]

  Fjac[1, 2] +=     q[1]*U_dot2
  Fjac[2, 2] += U + q[2]*U_dot2 + dir[1]*p_dot[2]
  Fjac[3, 2] +=     q[3]*U_dot2 + dir[2]*p_dot[2]
  Fjac[4, 2] +=     q[4]*U_dot2 + dir[3]*p_dot[2]
  Fjac[5, 2] +=     q[5]*U_dot2 + press*U_dot2 + U*p_dot[2]

  Fjac[1, 3] +=     q[1]*U_dot3
  Fjac[2, 3] +=     q[2]*U_dot3 + dir[1]*p_dot[3]
  Fjac[3, 3] += U + q[3]*U_dot3 + dir[2]*p_dot[3]
  Fjac[4, 3] +=     q[4]*U_dot3 + dir[3]*p_dot[3]
  Fjac[5, 3] +=     q[5]*U_dot3 + press*U_dot3 + U*p_dot[3]

  Fjac[1, 4] +=     q[1]*U_dot4
  Fjac[2, 4] +=     q[2]*U_dot4 + dir[1]*p_dot[4]
  Fjac[3, 4] +=     q[3]*U_dot4 + dir[2]*p_dot[4]
  Fjac[4, 4] += U + q[4]*U_dot4 + dir[3]*p_dot[4]
  Fjac[5, 4] +=     q[5]*U_dot4 + press*U_dot4 + U*p_dot[4]

  Fjac[1, 5] += 0
  Fjac[2, 5] += dir[1]*p_dot[5]
  Fjac[3, 5] += dir[2]*p_dot[5]
  Fjac[4, 5] += dir[3]*p_dot[5]
  Fjac[5, 5] += U + U*p_dot[5]

  return nothing

end

"""
  Computes the gradient of pressure with respect to `q` at a node.
  Methods are available in 2D and 3D

  **Inputs**

   * params: ParamType, conservative variables only
   * q: vector of conservative variables at the node

  **Inputs/Outputs**

   * pdot: vector of length numDofPerNode, overwritten with derivative of `p` 
           wrt `q` (overwritten)
"""
function calcPressure_diff(params::ParamType{2, :conservative},
                      q::AbstractArray{Tsol,1}, p_dot::AbstractVector{Tsol} ) where Tsol
  # calculate pressure for a node
  # q is a vector of length 4 of the conservative variables

  t1 = 1/(q[1]*q[1])
  t2 = q[2]*q[2]
  t3 = q[3]*q[3]

  p_dot[1] = (params.gamma_1)*( 0.5*(t2*t1 + t3*t1))
  p_dot[2] = -(params.gamma_1)*(q[2]/q[1])
  p_dot[3] = -(params.gamma_1)*(q[3]/q[1])
  p_dot[4] = params.gamma_1

  return  (params.gamma_1)*(q[4] - 0.5*(t2 + t3)/q[1])
end



function calcPressure_diff(params::ParamType{3, :conservative},
                      q::AbstractArray{Tsol,1}, p_dot::AbstractVector{Tsol} ) where Tsol
  # calculate pressure for a node

  t1 = 1/(q[1]*q[1])
  t2 = q[2]*q[2]
  t3 = q[3]*q[3]
  t4 = q[4]*q[4]

  p_dot[1] =  params.gamma_1*( 0.5*(t2 + t3 + t4)*t1)
  p_dot[2] = -params.gamma_1*(q[2]/q[1])
  p_dot[3] = -params.gamma_1*(q[3]/q[1])
  p_dot[4] = -params.gamma_1*(q[4]/q[1])
  p_dot[5] =  params.gamma_1

  return (params.gamma_1)*(q[5] - 0.5*(t2 + t3 + t4)/q[1])
end


"""
  Reverse mode of the pressure calculation

  **Inputs**

   * params: ParamType object
   * q: vector of solution variables at the node (length numDofPerNode)
   * p_bar: seed value for pressure

  **Inputs/Outputs**

   * q_bar: vector to be updated (not overwritten) with the result
"""
function calcPressure_revq(params::ParamType{2, :conservative},
                           q::AbstractArray{Tsol, 1}, q_bar::AbstractArray{Tsol, 1},
                           p_bar::Number) where {Tsol}

  q_bar[1] +=  params.gamma_1*0.5*(q[2]*q[2] + q[3]*q[3])/(q[1]*q[1])*p_bar
  q_bar[2] += -params.gamma_1*q[2]*p_bar/q[1]
  q_bar[3] += -params.gamma_1*q[3]*p_bar/q[1]
  q_bar[4] +=  params.gamma_1*p_bar

  return nothing
end


function calcPressure_revq(params::ParamType{3, :conservative},
                           q::AbstractArray{Tsol, 1}, q_bar::AbstractArray{Tsol, 1},
                           p_bar::Number) where {Tsol}

  q_bar[1] +=  params.gamma_1*0.5*(q[2]*q[2] + q[3]*q[3] + q[4]*q[4])/(q[1]*q[1])*p_bar
  q_bar[2] += -params.gamma_1*q[2]*p_bar/q[1]
  q_bar[3] += -params.gamma_1*q[3]*p_bar/q[1]
  q_bar[4] += -params.gamma_1*q[4]*p_bar/q[1]
  q_bar[5] +=  params.gamma_1*p_bar

  return nothing
end


"""
  Differentiated version of [`getLambdaMax`](@ref)

  **Inputs**

   * params: ParamType
   * qL: vector of conservative variables at a node
   * dir: direction vector (can be scaled)

  **Inputs/Outputs**

   * qL_dot: derivative of lambda max wrt qL (overwritten)

  **Outputs**

   * lambda_max: maximum eigenvalue
"""
function getLambdaMax_diff(params::ParamType{2},
                      qL::AbstractVector{Tsol},
                      dir::AbstractVector{Tmsh},
                      lambda_dot::AbstractVector{Tres}) where {Tsol, Tres, Tmsh}

  gamma = params.gamma
  Un = zero(Tres)
  dA = zero(Tmsh)
  rhoLinv = 1/qL[1]
  rhoLinv_dotL1 = -rhoLinv*rhoLinv

  p_dot = params.get_lambda_max_data.p_dot
  pL = calcPressure_diff(params, qL, p_dot)
  aL = sqrt(gamma*pL*rhoLinv)  # speed of sound
  t1 = gamma*rhoLinv/(2*aL)
  t2 = gamma*pL/(2*aL)
  aL_dotL1 = t1*p_dot[1] + t2*rhoLinv_dotL1
  aL_dotL2 = t1*p_dot[2]
  aL_dotL3 = t1*p_dot[3]
  aL_dotL4 = t1*p_dot[4]


  Un_dotL1 = dir[1]*qL[2]*rhoLinv_dotL1
  Un_dotL2 = dir[1]*rhoLinv
  Un += dir[1]*qL[2]*rhoLinv

  Un_dotL1 += dir[2]*qL[3]*rhoLinv_dotL1
  Un_dotL3 = dir[2]*rhoLinv
  Un += dir[2]*qL[3]*rhoLinv

  for i=1:2
    dA += dir[i]*dir[i]
  end

  dA = sqrt(dA)

  lambda_max = absvalue(Un) + dA*aL
  lambda_dot[1] = dA*aL_dotL1
  lambda_dot[2] = dA*aL_dotL2
  lambda_dot[3] = dA*aL_dotL3
  lambda_dot[4] = dA*aL_dotL4

  if Un > 0
    lambda_dot[1] += Un_dotL1
    lambda_dot[2] += Un_dotL2
    lambda_dot[3] += Un_dotL3
  else
    lambda_dot[1] -= Un_dotL1
    lambda_dot[2] -= Un_dotL2
    lambda_dot[3] -= Un_dotL3
  end


  return lambda_max
end



function getLambdaMax_diff(params::ParamType{3},
                      qL::AbstractVector{Tsol},
                      dir::AbstractVector{Tmsh},
                      lambda_dot::AbstractVector{Tres}) where {Tsol, Tres, Tmsh}

  gamma = params.gamma
  Un = zero(Tres)
  dA = zero(Tmsh)
  rhoLinv = 1/qL[1]
  rhoLinv_dotL1 = -rhoLinv*rhoLinv

  p_dot = params.get_lambda_max_data.p_dot
  pL = calcPressure_diff(params, qL, p_dot)
  aL = sqrt(gamma*pL*rhoLinv)  # speed of sound
  t1 = gamma*rhoLinv/(2*aL)
  t2 = gamma*pL/(2*aL)
  aL_dotL1 = t1*p_dot[1] + t2*rhoLinv_dotL1
  aL_dotL2 = t1*p_dot[2]
  aL_dotL3 = t1*p_dot[3]
  aL_dotL4 = t1*p_dot[4]
  aL_dotL5 = t1*p_dot[5]


  Un_dotL1 = dir[1]*qL[2]*rhoLinv_dotL1
  Un_dotL2 = dir[1]*rhoLinv
  Un += dir[1]*qL[2]*rhoLinv

  Un_dotL1 += dir[2]*qL[3]*rhoLinv_dotL1
  Un_dotL3 = dir[2]*rhoLinv
  Un += dir[2]*qL[3]*rhoLinv

  Un_dotL1 += dir[3]*qL[4]*rhoLinv_dotL1
  Un_dotL4 = dir[3]*rhoLinv
  Un += dir[3]*qL[4]*rhoLinv


  for i=1:3
    dA += dir[i]*dir[i]
  end

  dA = sqrt(dA)

  lambda_max = absvalue(Un) + dA*aL
  lambda_dot[1] = dA*aL_dotL1
  lambda_dot[2] = dA*aL_dotL2
  lambda_dot[3] = dA*aL_dotL3
  lambda_dot[4] = dA*aL_dotL4
  lambda_dot[5] = dA*aL_dotL5

  if Un > 0
    lambda_dot[1] += Un_dotL1
    lambda_dot[2] += Un_dotL2
    lambda_dot[3] += Un_dotL3
    lambda_dot[4] += Un_dotL4
  else
    lambda_dot[1] -= Un_dotL1
    lambda_dot[2] -= Un_dotL2
    lambda_dot[3] -= Un_dotL3
    lambda_dot[4] -= Un_dotL4
  end


  return lambda_max
end


"""
  Reverse mode of [`getLambdaMax`](@ref) with respect to the metrics

  **Inputs**

   * params
   * qL
   * dir
   * lambda_bar: seed value for reverse mode

  **Inputs/Outputs**

   * nrm_bar: vector to be updated (not overwritten) with the result
"""
function getLambdaMax_revm(params::ParamType{Tdim}, 
                      qL::AbstractVector{Tsol}, 
                      dir::AbstractVector{Tmsh},
                      dir_bar::AbstractVector{Tmsh}, lambda_bar::Number) where {Tsol, Tmsh, Tdim}

  Tres = promote_type(Tsol, Tmsh)
  gamma = params.gamma
  Un = zero(Tres)
  dA = zero(Tmsh)
  rhoLinv = 1/qL[1]

  
  pL = calcPressure(params, qL)
  aL = sqrt(gamma*pL*rhoLinv)  # speed of sound

  for i=1:Tdim
    Un += dir[i]*qL[i+1]*rhoLinv
    dA += dir[i]*dir[i]
  end

  dA2 = sqrt(dA)

  lambda_max = absvalue(Un) + dA2*aL

  # reverse sweep
  fac = Un > 0 ? 1 : -1
  Un_bar = fac*lambda_bar
  dA2_bar = aL*lambda_bar
#  aL_bar = dA2*lambda_bar

  dA_bar = dA2_bar/(2*dA2)

  for i=1:Tdim
    # Un
    dir_bar[i] += qL[i+1]*rhoLinv*Un_bar
    # dA
    dir_bar[i] += 2*dir[i]*dA_bar
  end


  return lambda_max
end


"""
  Reverse mode wrt q of [`getLambdaMax`](@ref)

  **Inputs**

   * params
   * qL
   * dir
   * lambda_bar: reverse mode seed value

  **Inputs/Outputs**

   * qL_bar: adjoint part of qL, will be updated (not overwritten)
"""
function getLambdaMax_revq(params::ParamType{Tdim}, 
                           qL::AbstractVector{Tsol}, 
                           qL_bar::AbstractVector{Tsol},
                           dir::AbstractVector{Tmsh},
                           lambda_bar::Number) where {Tsol, Tmsh, Tdim}

  Tres = promote_type(Tsol, Tmsh)
  gamma = params.gamma
  Un = zero(Tres)
  dA = zero(Tmsh)
  rhoLinv = 1/qL[1]

  pL = calcPressure(params, qL)
  aL = sqrt(gamma*pL*rhoLinv)  # speed of sound

  for i=1:Tdim
    Un += dir[i]*qL[i+1]*rhoLinv
    dA += dir[i]*dir[i]
  end

  dA2 = sqrt(dA)

  lambda_max = absvalue(Un) + dA2*aL
  rhoLinv_bar = zero(Tsol)
 

  # reverse sweep
  fac = Un > 0 ? 1 : -1
  Un_bar = fac*lambda_bar
#  dA2_bar = aL*lambda_bar
  aL_bar = dA2*lambda_bar

  
  # dA2 = sqrt(dA)
#  dA_bar = dA2_bar/(2*dA2)

  rhoLinv_bar = zero(Tsol)
  for i=1:Tdim
    qL_bar[i+1] += dir[i]*rhoLinv*Un_bar
    rhoLinv_bar += dir[i]*qL[i+1]*Un_bar
    # dir is not being differentated here
  end

  # aL
  pL_bar = gamma*rhoLinv*aL_bar/(2*aL)
  rhoLinv_bar += gamma*pL*aL_bar/(2*aL)

  calcPressure_revq(params, qL, qL_bar, pL_bar)

  qL_bar[1] += -rhoLinv*rhoLinv*rhoLinv_bar
  
  return lambda_max
end



