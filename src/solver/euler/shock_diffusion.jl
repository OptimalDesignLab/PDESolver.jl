# the elementwise-constant diffusion model used for shock capturing

"""
  Applies the diffusion tensor to a given array, for an element 

  More specifically,

  ```
  for i=1:numNodesPerElement
    for d=1:mesh.dim
      for d2=1:mesh.dim
        flux[:, i, d] = C[:, :, d1, d2]*dx[:, i, d2]
      end
    end
  end
  ```
  where C[:, :, :, :] is the diffusion tensor computed at state w[:, i]

  **Inputs**

   * obj: the [`AbstractDiffusion`](@ref) object
   * w: the `numDofPerNode` x `numNodesPerElement` array of entropy variables
        for the element
   * i: the element number.  This is useful for diffusions where the coeffients
        are precomputed and stored in an array
   * dx: the values to multiply against, `numDofPerNode` x `numNodesPerElement`
         x `dim`

  **Inputs/Outputs**

   * flux: array to overwrite with the result, same size as `dx`

"""
function applyDiffusionTensor(obj::ShockDiffusion, w::AbstractMatrix,
                    i::Integer, dx::Abstract3DArray, flux::Abstract3DArray)

  # For shock capturing it is much easier because the viscoscity is diagonal
  # for C[:, :, j, j] and zero for the cross terms

  numDofPerNode, numNodesPerElement, dim = size(flux)
  ee = obj.ee[i]
  @simd for d=1:dim
    @simd for j=1:numNodesPerElement
      @simd for k=1:numDofPerNode
        flux[k, j, d] = ee*dx[k, j, d]
      end
    end
  end

  return nothing
end

#TODO: diff
