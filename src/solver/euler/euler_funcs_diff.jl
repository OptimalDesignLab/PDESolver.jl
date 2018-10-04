# differentiated version of functions in euler_funcs.jl
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
  fill!(flux_jac, 0.0); fill!(res_jac, 0.0)

  for i=1:mesh.numEl
    fill!(flux_jac, 0)
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
      # multiply by Minv if needed
      for q=1:mesh.numNodesPerElement
        for p=1:mesh.numNodesPerElement
          val = mesh.jac[p, i]/sbp.w[p]  # entry in Minv
          @simd for m=1:mesh.numDofPerNode
            @simd for n=1:mesh.numDofPerNode
              res_jac[n, m, p, q] *= val
            end
          end
        end
      end
    end

    # assemble element level jacobian into the residual
    assembleElement(assembler, mesh, i, res_jac)
    fill!(res_jac, 0.0)
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


  @assert eqn.params.use_Minv != 1  # use_Minv not supported

  data = eqn.params.calc_volume_integrals_data
  @unpack data flux_jac res_jac nrm
  fill!(flux_jac, 0.0); fill!(res_jac, 0.0)

  for i=1:mesh.numEl
    fill!(flux_jac, 0)
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
      # multiply by Minv if needed
      for q=1:mesh.numNodesPerElement
        for p=1:mesh.numNodesPerElement
          val = mesh.jac[p, i]/sbp.w[p]  # entry in Minv
          @simd for m=1:mesh.numDofPerNode
            @simd for n=1:mesh.numDofPerNode
              res_jac[n, m, p, q] *= val
            end
          end
        end
      end
    end

    # assemble element level jacobian into the residual
    assembleElement(assembler, mesh, i, res_jac)
    fill!(res_jac, 0.0)
    # flux_jac gets overwritten, so no need to zero it 

  end  # end loop i

  return nothing
end  # end function


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


