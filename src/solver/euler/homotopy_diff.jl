# differentiated version of homotopy.jl

"""
  Differentiated version of [`getLambdaMax`](@ref)

  **Inputs**

   * params: ParamType
   * qL: vector of conservative variables at a node
   * dir: direction vector (can be scaled)

  **Inputs/Outputs**

   * qL_dot: derivative of lambda max wrt qL
"""
function getLambdaMax_diff{Tsol, Tres, Tmsh, Tdim}(params::ParamType{Tdim}, 
                      qL::AbstractVector{Tsol},
                      dir::AbstractVector{Tmsh},
                      lambda_dot::AbstractVector{Tres})

  Tres = promote_type(Tsol, Tmsh)
  gamma = params.gamma
  Un = zero(Tres)
  dA = zero(Tmsh)
  rhoLinv = 1/qL[1]
  rhoLinv_dotL1 = -1/(rhoLinv*rhoLinv)

  p_dot = params.p_dot
  pL = calcPressure_ditt(params, qL, p_dot)
  aL = sqrt(gamma*pL*rhoLinv)  # speed of sound


  Un_dotL1 = dir[1]*q
  Un += dir[1]*qL[2]*rhoLinv

  Un += dir[2]*qL[3]*rhoLinv

  for i=1:Tdim
    Un += dir[i]*qL[i+1]*rhoLinv
    dA += dir[i]*dir[i]
  end

  lambda_max = absvalue(Un) + dA*aL

  return lambda_max
end


