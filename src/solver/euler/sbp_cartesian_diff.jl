# differentiated version of SBP cartesian

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
  Internal function @simd for computing the in-place transpose of a 3D array
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
# Applying operator matrices


function applyOperatorJac(Dx::Abstract3DArray, t1::AbstractArray{T, 4},
                          t2::AbstractArray{T2, 5}) where {T, T2}
# compute Dx * t1, Dy * t1, Dz*t1

  numNodesPerElement = size(Dx, 1)
  dim = size(Dx, 3)
  numDofPerNode = size(t1, 1)
#=
  @assert size(t2, 1) == numDofPerNode
  @assert size(t2, 2) == numDofPerNode
  @assert size(t2, 3) == dim
  @assert size(t2, 4) == numNodesPerElement
  @assert size(t2, 5) == numNodesPerElement
=#

  @simd for q=1:numNodesPerElement
    @simd for p=1:numNodesPerElement
      # zero out next tile
      # For arrays larger than the L1 cache, this is probably better than
      # zeroing out the full matrix and then doing a second pass to write the
      # final values
      @simd for d=1:dim
        @simd for i=1:numDofPerNode
          @simd for j=1:numDofPerNode
            t2[i, j, d, p, q] = 0
          end
        end
      end
      @simd for k=1:numNodesPerElement  # summed index
        @simd for d=1:dim
          @simd for i=1:numDofPerNode
            @simd for j=1:numDofPerNode
              t2[i, j, d, p, q] += Dx[p, k, d]*t1[i, j, k, q]
            end
          end
        end
      end
    end
  end

  return nothing
end


function applyOperatorJac(Dx::Abstract3DArray, t1::AbstractArray{T, 5},
                          t2::AbstractArray{T2, 4}) where {T, T2}
# compute Dx * t1, Dy * t1, Dz*t1

  numNodesPerElement = size(Dx, 1)
  dim = size(Dx, 3)
  numDofPerNode = size(t1, 1)
#=
  @assert size(t2, 1) == numDofPerNode
  @assert size(t2, 2) == numDofPerNode
  @assert size(t2, 3) == numNodesPerElement
  @assert size(t2, 4) == numNodesPerElement
=#
  @simd for q=1:numNodesPerElement
    @simd for p=1:numNodesPerElement

      # zero out next tile
      @simd for i=1:numDofPerNode
        @simd for j=1:numDofPerNode
          t2[i, j, p, q] = 0
        end
      end

      @simd for k=1:numNodesPerElement  # summed index
        @simd for d=1:dim
          @simd for i=1:numDofPerNode
            @simd for j=1:numDofPerNode
              t2[i, j, p, q] += Dx[p, k, d]*t1[i, j, d, k, q]
            end
          end
        end
      end
    end
  end

  return nothing
end

