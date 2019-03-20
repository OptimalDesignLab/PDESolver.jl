# differentiated version of SBP cartesian

import SummationByParts: UnaryFunctor, Add, Subtract

#------------------------------------------------------------------------------
# Computing operator matrices

"""
  Computes the SBP differentiation operator in the Cartesian directions.

  **Inputs**

   * sbp
   * dxidx: scaled mapping jacobian, `dim` x `dim` x `numNodesPerElement`
   * jac: mapping jacobian determinant, `numNodesPerElement`

  **Inputs/Outputs**

   * Dx: `numNodesPerElement` x `numNodesPerElement` x `dim` array to be
         overwritten with the D operator in each direction

"""
function calcDx(sbp::AbstractOperator, dxidx::Abstract3DArray,
                jac::AbstractVector, Dx::Abstract3DArray)

  numNodesPerElement = size(sbp.Q, 1)
  dim = size(sbp.Q, 3)

  # Dx = |J| * (dxi/dx * 1/|J|) * 1/w * Qxi + |J| * (deta/dx * 1/|J|) * 1/w * Qeta
  fill!(Dx, 0)
  @simd for d1=1:dim  # cartesian direction

    @simd for d2=1:dim  # parametric direction
      @simd for j=1:numNodesPerElement
        @simd for i=1:numNodesPerElement
          Dx[i, j, d1] += dxidx[d2, d1, i]*sbp.Q[i, j, d2]
        end
      end
    end

    
    # apply |J|/w
    @simd for j=1:numNodesPerElement
      @simd for i=1:numNodesPerElement
        Dx[i, j, d1] *= jac[i]/sbp.w[i]
      end
    end
    
  end

  return nothing
end

"""
  Computes the transpose of `calcDx`, same arguments
"""
function calcDxTransposed(sbp::AbstractOperator, dxidx::Abstract3DArray,
                          jac::AbstractVector, Dx::Abstract3DArray)

  calcDx(sbp, dxidx, jac, Dx)
  transpose3D(Dx)

  return nothing
end

"""
  Calculates the SBP weak differentiation operator Q in the Cartesian
  directions.

  **Inputs**

   * sbp
   * dxidx: scaled mapping jacobian, `dim` x `dim` x `numNodesPerElement`

  **Inputs/Outputs**

   * Qx: `numNodesPerElement` x `numNodesPerElement` x `dim` array to be
         overwritten with the Q operator in each direction

"""
function calcQx(sbp::AbstractOperator, dxidx::Abstract3DArray, Qx::Abstract3DArray)

  numNodesPerElement = size(sbp.Q, 1)
  dim = size(sbp.Q, 3)

  # Qx =  (dxi/dx * 1/|J|) * Qxi +  (deta/dx * 1/|J|) * Qeta
  fill!(Qx, 0)
  @simd for d1=1:dim  # cartesian direction
    @simd for d2=1:dim  # parametric direction
      @simd for j=1:numNodesPerElement
        @simd for i=1:numNodesPerElement
          Qx[i, j, d1] += dxidx[d2, d1, i]*sbp.Q[i, j, d2]
        end
      end
    end
  end

  return nothing
end


"""
  Computes the transpose of `calcQx`, same arguments.
"""
function calcQxTransposed(sbp::AbstractOperator, dxidx::Abstract3DArray,
                          Qx::Abstract3DArray)

  calcQx(sbp, dxidx, Qx)
  transpose3D(Qx)

  return nothing
end


"""
  Internal function for computing the in-place transpose of a 3D array
  (transposed the first 2 dimensions).

  **Inputs/Outputs**

   * Dx: a 3D array
"""
function transpose3D(Dx::Abstract3DArray)

  numNodesPerElement = size(Dx, 1)
  dim = size(Dx, 3)

  @simd for d1=1:dim
    @simd for j=1:numNodesPerElement
      @simd for i=(j+1):numNodesPerElement
        tmp = Dx[j, i, d1]
        Dx[j, i, d1] = Dx[i, j, d1]
        Dx[i, j, d1] = tmp
      end
    end
  end

  return nothing
end


#------------------------------------------------------------------------------
# Applying operator matrices (inner function)


"""
  Internal method for 3D -> 5D

  **Inputs**

   * Dx: a 3D operator array
   * nodes: the list of nodes to compute the output for
   * t1_dot: the input Jacobian
   * zero_output: if true, t2_dot will be overwritten.  If false, it will be
                  updated to
   * op: a `UnaryFunctor`, determines if the result will be added or subtracted
         from `t2_dot`.  Note that this takes effect whether `zero_output`
         is either true of false.

  **Inputs/Outputs**

   * t2_dot: the output Jacobian
"""
function applyOperatorJac(Dx::Abstract3DArray, nodes::AbstractVector,
                          t1::AbstractArray{T, 3},
                          t2::AbstractArray{T2, 5}, zero_output::Bool=true,
                          op::UnaryFunctor=Add()) where {T, T2}
# t1 is numDofPerNode x numDofPerNode x numNodesPerElement (nodewise jacobian
# at every element).  This is similar to the 4D case, except t1[i, j, k, q] 
# = 0 when k != q

# compute Dx * t1, Dy * t1, Dz*t1

  numNodesPerElement = size(Dx, 1)
  dim = size(Dx, 3)
  s1 = size(t1, 1); s2 = size(t1, 2)
#=
  @assert size(t2, 1) == numDofPerNode
  @assert size(t2, 2) == numDofPerNode
  @assert size(t2, 3) == dim
  @assert size(t2, 4) == numNodesPerElement
  @assert size(t2, 5) == numNodesPerElement
=#

  @simd for q=1:numNodesPerElement
    @simd for p in nodes
      # zero out next tile
      # For arrays larger than the L1 cache, this is probably better than
      # zeroing out the full matrix and then doing a second pass to write the
      # final values
      if zero_output
        @simd for d=1:dim
          @simd for j=1:s2
            @simd for i=1:s1
              t2[i, j, d, p, q] = 0
            end
          end
        end
      end
      @simd for d=1:dim
        @simd for j=1:s2
          @simd for i=1:s1
            t2[i, j, d, p, q] += op(Dx[p, q, d]*t1[i, j, q])
          end
        end
      end
    end
  end

  return nothing
end


"""
  Internal method for 4D -> 5D
"""
function applyOperatorJac(Dx::Abstract3DArray, nodes::AbstractVector,
                          t1::AbstractArray{T, 4},
                          t2::AbstractArray{T2, 5}, zero_output::Bool=true,
                          op::UnaryFunctor=Add()) where {T, T2}
# 4D -> 5D
# compute Dx * t1, Dy * t1, Dz*t1

  numNodesPerElement = size(Dx, 1)
  dim = size(Dx, 3)
  s1 = size(t1, 1); s2 = size(t1, 2)
#=
  @assert size(t2, 1) == numDofPerNode
  @assert size(t2, 2) == numDofPerNode
  @assert size(t2, 3) == dim
  @assert size(t2, 4) == numNodesPerElement
  @assert size(t2, 5) == numNodesPerElement
=#

  @simd for q=1:numNodesPerElement
    @simd for p in nodes
      # zero out next tile
      # For arrays larger than the L1 cache, this is probably better than
      # zeroing out the full matrix and then doing a second pass to write the
      # final values
      if zero_output
        @simd for d=1:dim
          @simd for j=1:s2
            @simd for i=1:s1
              t2[i, j, d, p, q] = 0
            end
          end
        end
      end
      @simd for k=1:numNodesPerElement  # summed index
        @simd for d=1:dim
          @simd for j=1:s2
            @simd for i=1:s1
              t2[i, j, d, p, q] += op(Dx[p, k, d]*t1[i, j, k, q])
            end
          end
        end
      end
    end
  end

  return nothing
end

"""
  Internal method for 5D -> 4D
"""
function applyOperatorJac(Dx::Abstract3DArray, nodes::AbstractVector,
                          t1::AbstractArray{T, 5},
                          t2::AbstractArray{T2, 4}, zero_output::Bool=true,
                          op::UnaryFunctor=Add()) where {T, T2}
# 5D -> 4D
# compute Dx * t1, Dy * t1, Dz*t1

  numNodesPerElement = size(Dx, 1)
  dim = size(Dx, 3)
  s1 = size(t1, 1); s2 = size(t1, 2)
#=
  @assert size(t2, 1) == numDofPerNode
  @assert size(t2, 2) == numDofPerNode
  @assert size(t2, 3) == numNodesPerElement
  @assert size(t2, 4) == numNodesPerElement
=#
  @simd for q=1:numNodesPerElement
    @simd for p in nodes

      # zero out next tile
      if zero_output
        @simd for j=1:s2
          @simd for i=1:s1
            t2[i, j, p, q] = 0
          end
        end
      end

      @simd for k=1:numNodesPerElement  # summed index
        @simd for d=1:dim
          @simd for j=1:s2
            @simd for i=1:s1
              t2[i, j, p, q] += op(Dx[p, k, d]*t1[i, j, d, k, q])
            end
          end
        end
      end
    end
  end

  return nothing
end


#------------------------------------------------------------------------------
# User interface functions
# These functions don't care about array dimensionality, to make things simpler

"""
  Alias to make type signatures simpler.  This type can be any Jacobian
  (3D, 4D, 5D)
"""
const AnyJacobian{T} = Union{AbstractArray{T, 3}, AbstractArray{T, 4}, AbstractArray{T, 5}}

"""
  Computes the Jacobian of applying an operator matrix to an array. This method
  computes the result for all nodes of the element.  Methods are available for:

   * 4D -> 5D: computes Jacobian of t2 = [Dx * t1, Dy * t1, Dz*t1].
               The input Jacobian `t1_dot` is `numDofPerNode` x `numDofPerNode`
               x `numNodesPerElement` x `numNodesPerElement`.  The layout
               is dR[i, p]/dq[j, q] = t2_dot[i, j, p, q].
               The output
               Jacobian `t2_dot` is `numDofPerNode` x `numDofPerNode` x `dim`
               x `numNodesPerElement` x `numNodesPerElement`.  The result of
               `Dx * t1_dot` is stored in `t2_dot[:, :, 1, :, :]`, and similarly
               `t2_dot[:, :, 2, :, :]` contains `Dy * t1_dot`.
   * 3D -> 5D: similar to 4D -> 5D, but the input Jacobian is `numDofPerNode`
                x `numDofPerNode` x `numNodesPerElement`.  This array represents
                the `numDofPerNode` x `numDofPerNode` Jacobian at each node
                of the element (ie. the derivative of a quantity at node `i`
                wrt node `j` is zero, so the 4D array becomes 3D).
   * 5D -> 4D: computes the Jacobian of t2 = Dx * tx + Dy * ty + Dz * tz.
               The Jacobians have the same dimensions as the 4D -> 5D method.
               In this case `t2_dot[:, :, d, :, :] contains `tx_dot` when
               `d = 1`, and `ty_dot` when `d = 2`.

  **Inputs**

   * Dx: a 3D operator array
   * t1_dot: the input Jacobian
   * zero_output: if true, t2_dot will be overwritten.  If false, it will be
                  updated to
   * op: a `UnaryFunctor`, determines if the result will be added or subtracted
         from `t2_dot`.  Note that this takes effect whether `zero_output`
         is either true of false.

  **Inputs/Outputs**

   * t2_dot: the output Jacobian

"""
function applyOperatorJac(Dx::Abstract3DArray,
                          t1::AnyJacobian, t2::AnyJacobian,
                          zero_output::Bool=true,
                          op::UnaryFunctor=Add())

  nodes = 1:size(Dx, 1)
  applyOperatorJac(Dx, nodes, t1, t2, zero_output, op)

  return nothing
end

"""
  This method computes the result for only the nodes in the stencil of `R`
  (3rd dimension of a 4D Jacobian, 4th dimension of a 5D Jacobian).
  See the other method for details.

  **Inputs**

   * Dx
   * sbpface: an `AbstractFace` that describes the stencil of `R`
   * face: the local face number
   * t1_dot
   * zero_output
   * op

  **Inputs/Outputs**

   * t2_dot
"""
function applyOperatorJac(Dx::Abstract3DArray, sbpface::AbstractFace,
                          face::Integer,
                          t1::AnyJacobian, t2::AnyJacobian,
                          zero_output::Bool=true,
                          op::UnaryFunctor=Add())

  nodes = sview(sbpface.perm, :, face)
  applyOperatorJac(Dx, nodes, t1, t2, zero_output, op)

  return nothing
end


#------------------------------------------------------------------------------
# Reverse mode wrt q (the input)

const TwoOr3DArray{T} = Union{AbstractMatrix{T}, Abstract3DArray{T}}
"""
  Reverse mode of `applyTxTransposed` with respect to the input array
  `w`.  These methods work `w` as a 2D or 3D array, and similarly for `wx`

  **Inputs**

   * sbp
   * dxidx
   * wxi
   * wx_bar: the seed array for back-propigation
   * op: UnaryFunctor determining whether to add or subtract from the output

  **Inputs/Outputs**

   * w_bar: the array to be updated with the result of back-propigation
"""
function applyQxTransposed_revq(sbp, w::TwoOr3DArray, dxidx::Abstract3DArray,
                                wxi::Abstract3DArray, wx::TwoOr3DArray,
                                op::SummationByParts.UnaryFunctor=SummationByParts.Add())

  applyQx(sbp, wx, dxidx, wxi, w, op)

end


"""
  Reverse mode of `applyQx`.  See other function for details
"""
function applyQx_revq(sbp, w::TwoOr3DArray, dxidx::Abstract3DArray,
                 wxi::Abstract3DArray, wx::TwoOr3DArray,
                 op::SummationByParts.UnaryFunctor=SummationByParts.Add())

  applyQxTransposed(sbp, wx, dxidx, wxi, w, op)
end


"""
  Reverse mode of of `applyDxTransposed`.  See other function for details
"""
function applyDxTransposed_revq(sbp, w::TwoOr3DArray, dxidx::Abstract3DArray,
                                jac::AbstractVector,
                                wxi::Abstract3DArray, wx::TwoOr3DArray,
                                op::SummationByParts.UnaryFunctor=SummationByParts.Add())

  applyDx(sbp, wx, dxidx, jac, wxi, w, op)

end


"""
  Reverse mode of `applyDx`.  See other function for details.
"""
function applyDx_revq(sbp, w::TwoOr3DArray, dxidx::Abstract3DArray,
                 jac::AbstractVector,
                 wxi::Abstract3DArray, wx::TwoOr3DArray,
                 op::SummationByParts.UnaryFunctor=SummationByParts.Add())

  applyDxTransposed(sbp, wx, dxidx, jac, wxi, w, op)
end

#------------------------------------------------------------------------------
# revm, first method

"""
  Computes meverse mode of Qx^T * w with respect to the metrics in all
  Cartesian directions.

  **Inputs**

   * sbp: the SBP operator
   * w: `numDofPerNode` x `numNodesPerElement` array of values
   * dxidx: `dim` x `dim` x `numNodesPerElement` array containing the metrics
   * wx_bar: seed vector for reverse mode, same size as `wxi`
   * op: a UnaryFunctor that determines if values are added or subtracted to
         `wx`

  **Inputs/Outputs**

   * wxi: temporary array, `numDofPerNode` x `numNodesPerElement` x `dim`
   * dxidx_bar: array updated with result, same size as `dxidx`.
"""
function applyQxTransposed_revm(sbp, w::AbstractMatrix, dxidx::Abstract3DArray,
                             dxidx_bar::Abstract3DArray,
                             wxi::Abstract3DArray, wx_bar::Abstract3DArray,
                             op::SummationByParts.UnaryFunctor=SummationByParts.Add())

  # The idea is to compute dw/dxi first, and then use dxi/dx to rotate those
  # arrays to be d/dx

  # compute dw/dxi
  numDofPerNode, numNodesPerElement, dim = size(wxi)
  for d=1:dim
    smallmatmat!(w, sview(sbp.Q, :, :, d), sview(wxi, :, :, d))
  end

  # dw/dx = dw/dxi * dxi/dx + dw/dy * dy/dxi
  @simd for d1=1:dim
    @simd for i=1:numNodesPerElement
      @simd for d2=1:dim
        @simd for j=1:numDofPerNode
          #wx[j, i, d1] += op(wxi[j, i, d2]*dxidx[d2, d1, i])
          dxidx_bar[d2, d1, i] += op(wxi[j, i, d2]*wx_bar[j, i, d1])
        end
      end
    end
  end

  return nothing
end

"""
  Reverse mode of `applyDx`with respect to the metrics.  See
  [`applyQxTransposed_revm`](@ref) for details

  **Inputs**

   * sbp
   * w
   * dxidx
   * jac
   * wx_bar

  **Inputs/Outputs**

   * jac_bar: updated with result
   * dxidx_bar: updated with result
   * wxi: temporary array
"""
function applyDx_revm(sbp, w::AbstractMatrix, dxidx::Abstract3DArray,
                 dxidx_bar::Abstract3DArray,
                 jac::AbstractVector, jac_bar::AbstractVector,
                 wxi::Abstract3DArray, wx_bar::Abstract3DArray,
                 op::SummationByParts.UnaryFunctor=SummationByParts.Add())

  numDofPerNode, numNodesPerElement, dim = size(wxi)
  fill!(wxi, 0)
  for d=1:dim
    differentiateElement!(sbp, d, w, sview(wxi, :, :, d))
  end

  @simd for d1=1:dim
    @simd for j=1:numNodesPerElement
      @simd for d2=1:dim
#        fac = jac[j]*dxidx[d2, d1, j]
        @simd for k=1:numDofPerNode
          #wx[k, j, d1] += op(wxi[k, j, d2] * fac)
          jac_bar[j] += op(wxi[k, j, d2]*dxidx[d2, d1, j]*wx_bar[k, j, d1])
          dxidx_bar[d2, d1, j] += op(wxi[k, j, d2]*jac[j]*wx_bar[k, j, d1])
        end
      end
    end
  end

  return nothing
end

"""
  Reverse mode of `applyDxTransposed` with respect to the metrics.  See
  [`applyDx_revm`](@ref) for details.
"""
function applyDxTransposed_revm(sbp, w::AbstractMatrix, dxidx::Abstract3DArray,
                 dxidx_bar::Abstract3DArray,
                 jac::AbstractVector, jac_bar::AbstractVector,
                 wxi_bar::Abstract3DArray, wx_bar::Abstract3DArray,
                 op::SummationByParts.UnaryFunctor=SummationByParts.Add())


  numDofPerNode, numNodesPerElement, dim = size(wxi_bar)

  # Dx^T * w = (Dxi^T * |dxi/dx| * u + Deta^T * |deta/dx|* w)
  # so rotate w using mapping jacobian first, then apply D
  
  for d1=1:dim
    fill!(wxi_bar, 0)
    for d2=1:dim
#=
      for j=1:numNodesPerElement
        for k=1:numDofPerNode
          # we only need temporary storage in the shape of w, but to be
          # consistent with the other functions we use the wxi array
          wxi[k, j, d2] = w[k, j]*(dxidx[d2, d1, j]*jac[j])
        end
      end

      differentiateElement!(sbp, d2, ro_sview(wxi, :, :, d2), sview(wx, :, :, d1), op, true)
=#
      #--------------------
      # reverse sweep

      differentiateElement!(sbp, d2, ro_sview(wx_bar, :, :, d1),
                            sview(wxi_bar, :, :, d2), op, false)

      for j=1:numNodesPerElement
        for k=1:numDofPerNode
          dxidx_bar[d2, d1, j] += w[k, j]*jac[j]*wxi_bar[k, j, d2]
          jac_bar[j] += w[k, j]*dxidx[d2, d1, j]*wxi_bar[k, j, d2]
        end
      end
    end  # end d2
  end  # end d1

  return nothing
end

"""
  Reverse mode of `applyQx` with respect to the metrics.  See
  [`applyQxTransposed_revm`](@ref) for details.
"""
function applyQx_revm(sbp, w::AbstractMatrix, dxidx::Abstract3DArray,
                 dxidx_bar::Abstract3DArray,
                 wxi::Abstract3DArray, wx_bar::Abstract3DArray,
                 op::SummationByParts.UnaryFunctor=SummationByParts.Add())

  numDofPerNode, numNodesPerElement, dim = size(wxi)
  fill!(wxi, 0)
  for d=1:dim
    weakDifferentiateElement!(sbp, d, w, sview(wxi, :, :, d))
  end

  # reverse sweep
  @simd for d1=1:dim
    @simd for j=1:numNodesPerElement
      @simd for d2=1:dim
        @simd for k=1:numDofPerNode
          #wx[k, j, d1] += op(wxi[k, j, d2] * dxidx[d2, d1, j])
          dxidx_bar[d2, d1, j] += op(wxi[k, j, d2]*wx_bar[k, j, d1])
        end
      end
    end
  end

  return nothing
end

#------------------------------------------------------------------------------
# revm, second method


function applyQxTransposed_revm(sbp, w::Abstract3DArray, dxidx::Abstract3DArray,
                             dxidx_bar::Abstract3DArray,
                             wxi::Abstract3DArray, wx_bar::AbstractMatrix,
                             op::SummationByParts.UnaryFunctor=SummationByParts.Add())


  numDofPerNode, numNodesPerElement, dim = size(wxi)
  for d1=1:dim  # compute Q_d * w_d
    for d2=1:dim
      smallmatmat!(ro_sview(w, :, :, d1), ro_sview(sbp.Q, :, :, d2), sview(wxi, :, :, d2))
    end

    @simd for d2=1:dim
      @simd for j=1:numNodesPerElement
        @simd for k=1:numDofPerNode
          #wx[k, j] += op(dxidx[d2, d1, j]*wxi[k, j, d2])
          dxidx_bar[d2, d1, j] += op(wxi[k, j, d2]*wx_bar[k, j])
        end
      end
    end

  end  # end d1

  return nothing
end

function applyDx_revm(sbp, w::Abstract3DArray, dxidx::Abstract3DArray,
                 dxidx_bar::Abstract3DArray,
                 jac::AbstractVector, jac_bar::AbstractVector,
                 wxi::Abstract3DArray, wx_bar::AbstractMatrix,
                 op::SummationByParts.UnaryFunctor=SummationByParts.Add())

  numDofPerNode, numNodesPerElement, dim = size(w)

  @simd for d1=1:dim
    fill!(wxi, 0)
    for d2=1:dim
      differentiateElement!(sbp, d2, sview(w, :, :, d1), sview(wxi, :, :, d2))
    end

    @simd for j=1:numNodesPerElement
      @simd for d2=1:dim
        #fac = jac[j]*dxidx[d2, d1, j]
        @simd for k=1:numDofPerNode
          #wx[k, j] += op(wxi[k, j, d2] * fac)
          dxidx_bar[d2, d1, j] += op(wxi[k, j, d2]*jac[j]*wx_bar[k, j])
          jac_bar[j] += op(wxi[k, j, d2]*dxidx[d2, d1, j]*wx_bar[k, j])
        end
      end
    end

  end

  return nothing
end

function applyDxTransposed_revm(sbp, w::Abstract3DArray, dxidx::Abstract3DArray,
                 dxidx_bar::Abstract3DArray,
                 jac::AbstractVector, jac_bar::AbstractVector,
                 wxi_bar::Abstract3DArray, wx_bar::AbstractMatrix,
                 op::SummationByParts.UnaryFunctor=SummationByParts.Add())


  numDofPerNode, numNodesPerElement, dim = size(w)

  for d1=1:dim

#=
    for d2=1:dim
      for j=1:numNodesPerElement
        for k=1:numDofPerNode
          # we only need temporary storage in the shape of w, but to be
          # consistent with the other functions we use the wxi array
          wxi[k, j, d2] = w[k, j, d1]*(dxidx[d2, d1, j]*jac[j])
        end
      end

      differentiateElement!(sbp, d2, ro_sview(wxi, :, :, d2), wx, op, true)
=#
      # reverse sweep
    for d2=1:dim
      fill!(wxi_bar, 0)
      differentiateElement!(sbp, d2, wx_bar, sview(wxi_bar, :, :, d2), op, false)
      for j=1:numNodesPerElement
        for k=1:numDofPerNode
          # we only need temporary storage in the shape of w, but to be
          # consistent with the other functions we use the wxi array
          #wxi[k, j, d2] = w[k, j, d1]*(dxidx[d2, d1, j]*jac[j])
          dxidx_bar[d2, d1, j] += w[k, j, d1]*jac[j]*wxi_bar[k, j, d2]
          jac_bar[j] += w[k, j, d1]*dxidx[d2, d1, j]*wxi_bar[k, j, d2]
        end
      end
    end  # end d2

  end  # end d1

  return nothing
end

function applyQx_revm(sbp, w::Abstract3DArray, dxidx::Abstract3DArray,
                 dxidx_bar::Abstract3DArray,
                 wxi::Abstract3DArray, wx_bar::AbstractMatrix,
                 op::SummationByParts.UnaryFunctor=SummationByParts.Add())

  numDofPerNode, numNodesPerElement, dim = size(w)

  @simd for d1=1:dim
    fill!(wxi, 0)
    for d2=1:dim
      weakDifferentiateElement!(sbp, d2, sview(w, :, :, d1), sview(wxi, :, :, d2))
    end

    @simd for d2=1:dim
      @simd for j=1:numNodesPerElement
        @simd for k=1:numDofPerNode
          #wx[k, j] += op(wxi[k, j, d2] * dxidx[d2, d1, j])
          dxidx_bar[d2, d1, j] += op(wxi[k, j, d2]*wx_bar[k, j])
        end
      end
    end
  end

  return nothing
end


