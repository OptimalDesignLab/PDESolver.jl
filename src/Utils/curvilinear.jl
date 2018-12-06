# functions used for curvilinear calculations

"""
  Calculates the S matrices in the cartesian directions for a curvilinear
  element.

  Inputs:
    sbp: an AbstractOperator
    dxidx: the scaled mapping jacobian, dxidx/|J| at every node, Tdim x Txim
           x numNodesPerElement

  Inputs/Outputs:
    S: array to be populated with S, numNodesPerElement x numNodesPerElement x
       x Tdim

"""
function calcSCurvilinear(sbp::AbstractOperator, dxidx::AbstractArray{Tmsh, 3}, S::AbstractArray{T, 3}) where {T, Tmsh}
# calculate S in x, y, z directions

  Tdim = size(S, 3)
  nnodes = size(S, 1)
  fill!(S, 0.0)

  #TODO: test performance
  #      maybe doing 2 passed over the data in order would be better?
  for d1=1:Tdim  # cartesian direction
    for d2=1:Tdim  # parametric direction
      for i=1:nnodes
        for j=1:(i - 1)
          S[j, i, d1] += 0.5*dxidx[d2, d1, j ]*sbp.Q[j, i, d2] - 0.5*dxidx[d2, d1, i]*sbp.Q[i, j, d2]
          S[i, j, d1] = -S[j, i, d1]

        end
      end
    end
  end


  return nothing
end

"""
  Reverse mode of [`calcSCurvilinar`](@ref)

  **Inputs**

   * sbp
   * dxidx
   * S
   * S_bar: seed values for reverse mode

  **Outputs**

   * dxidx_bar: updated with reverse propigation of S_bar
"""
function calcSCurvilinear_rev(sbp::AbstractOperator, dxidx::AbstractArray{Tmsh, 3},
                               dxidx_bar::AbstractArray{Tmsh, 3},
                               S::AbstractArray{T, 3}, S_bar::AbstractArray{T, 3}) where {T, Tmsh}

  Tdim = size(S, 3)
  nnodes = size(S, 1)
  fill!(S, 0.0)

  for d1=1:Tdim  # cartesian direction
    for d2=1:Tdim  # parametric direction
      for i=1:nnodes
        for j=1:(i - 1)
          S[j, i, d1] += 0.5*dxidx[d2, d1, j ]*sbp.Q[j, i, d2] - 0.5*dxidx[d2, d1, i]*sbp.Q[i, j, d2]
          S[i, j, d1] = -S[j, i, d1]

          dxidx_bar[d2, d1, j] += 0.5*sbp.Q[j, i, d2]*S_bar[j, i, d1]
          dxidx_bar[d2, d1, i] -= 0.5*sbp.Q[i, j, d2]*S_bar[j, i, d1]

          dxidx_bar[d2, d1, j] -= 0.5*sbp.Q[j, i, d2]*S_bar[i, j, d1]
          dxidx_bar[d2, d1, i] += 0.5*sbp.Q[i, j, d2]*S_bar[i, j, d1]


        end
      end
    end
  end


  return nothing
end





"""
  Calculates the E matrices in the cartesian directions using the metric
  terms.  This is *not* the E used in the final curvilinear discretization.
  Note that this E is not symmetric

  Inputs:
    sbp: an AbstractOperator
    dxidx: the scaled mapping jacobian, dxidx/|J| at every node, Tdim x Txim
           x numNodesPerElement

  Inputs/Outputs:
    E: array to be populated with E, numNodesPerElement x numNodesPerElement x
       x Tdim

"""
function calcECurvilinear(sbp::AbstractOperator, dxidx::AbstractArray{Tmsh, 3}, E::AbstractArray{T, 3}) where {T, Tmsh}
# calculate E in x, y, z directions
# this isn't the E actually used by the code, it is Ex =  E_xi*metrics_xi_x + 
# E_eta*metrics_eta_x

  Tdim = size(E, 3)
  nnodes = size(E, 1)
  fill!(E, 0.0)
  for d1=1:Tdim  # cartesian direction
    for d2=1:Tdim  # parametric direction
      for i=1:nnodes
        for j=1:nnodes
          E_d2 = sbp.Q[j, i, d2] + sbp.Q[i, j, d2]
          E[j, i, d1] += E_d2*dxidx[d2, d1, i]
#          E[i, j, d1] = E[j, i, d1]
        end
      end
    end
  end

  return nothing
end

"""
  Calculates the D matrices in the cartesian directions using the metric
  terms.  This is *not* the D used in the final curvilinear discretization.

  Inputs:
    sbp: an AbstractOperator
    dxidx: the scaled mapping jacobian, dxidx/|J| at every node, Tdim x Txim
           x numNodesPerElement

  Inputs/Outputs:
    D: array to be populated with D, numNodesPerElement x numNodesPerElement x
       x Tdim

"""

function calcDCurvilinear(sbp::AbstractOperator, dxidx::AbstractArray{Tmsh, 3}, D::AbstractArray{T, 3}) where {T, Tmsh}
# calculate D in x, y, z directions
# this isn't the D actually used by the code, it is Dx =  D_xi*metrics_xi_x + 
# D_eta*metrics_eta_x

  Tdim = size(D, 3)
  nnodes = size(D, 1)
  fill!(D, 0.0)
  for d1=1:Tdim  # cartesian direction
    for d2=1:Tdim  # parametric direction
      for i=1:nnodes
        for j=1:nnodes
          h_j = 1./sbp.w[j]
          D[j, i, d1] += h_j*sbp.Q[j, i, d2]
        end
      end
    end
  end

  return nothing
end


