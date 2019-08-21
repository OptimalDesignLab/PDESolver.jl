# functions for interpolating between different SBP operators

#TODO: interpolate q or q_vec
"""
  This function interpolates a solution field from one SBP operator to
  another.  This requires the meshes be identical except for the difference
  in SBP operator.  In particular, the element numbering scheme must be the
  same for this to work.

  The field in question can be any field defined on the nodes of an SBP
  operator.  The interpolation is generally accurate to the degree of the
  SBP operator.

  **Inputs**

   * sbp_old: the SBP operator for the mesh with the solution to be interpolated
   * q_old: the 3D array containing the the solution to be interpolated
   * sbp_new: the SBP operator for the new mesh

  **Inputs/Outputs**

   * q_new: the array to put the new solution into (overwritten)
"""
function interpField(sbp_old::AbstractOperator, q_old::Abstract3DArray,
                        sbp_new::AbstractOperator, q_new::Abstract3DArray)

  @assert size(q_old, 3) == size(q_new, 3)
  @assert size(q_old, 1) == size(q_new, 1)

  # construct interpolation operator
  node_coords = calcnodes(sbp_new)
  interp_op = SummationByParts.buildinterpolation(sbp_old, node_coords)
  

  # apply interpolation operator
  applyInterpolation(interp_op, q_old, q_new)

  return nothing
end

"""
  This function interpolates a vector from one mesh to another.  The two
  meshes must have the elements numbered in the same order.

  Currently, the vector must have `mesh.numDofPerNode` degrees of freedom
  per node.

  **Inputs**

   * mesh_old: the old mesh
   * sbp_old: the old SBP operator
   * q_old: the vector containing the data, length `mesh_old.numDof`
   * mesh_new: the new mesh
   * sbp_new: the new SBP operator

  **Inputs/Outputs**
  
   * q_new: vector to overwrite with the interpolated values
"""
function interpField(mesh_old::AbstractMesh, sbp_old::AbstractOperator,
                     q_old::AbstractVector, mesh_new::AbstractMesh,
                     sbp_new::AbstractOperator, q_new::AbstractVector)

  @assert length(q_old) == mesh_old.numDof
  @assert length(q_new) == mesh_new.numDof
  # construct interpolation operator
  node_coords = calcnodes(sbp_new)
  interp_op = SummationByParts.buildinterpolation(sbp_old, node_coords)

  applyInterpolation(interp_op, mesh_old, q_old, mesh_new, q_new)

  return nothing
end


function applyInterpolation(interp_op::AbstractMatrix, q_old::Abstract3DArray,
                            q_new::Abstract3DArray)

  @assert size(q_old, 3) == size(q_new, 3)
  @assert size(q_old, 1) == size(q_new, 1)

  numEl = size(q_old, 3)

  for i=1:numEl
    qold_i = sview(q_old, :, :, i)
    qnew_i = sview(q_new, :, :, i)

    # to interpolate all dofs at the same time (in the layout of the q arrays)
    # do q_old[:, :, i]*interp_op.' = q_new[:, :, i]
    smallmatmatT!(qold_i, interp_op, qnew_i)
  end

  return nothing
end

function applyInterpolation(interp_op::AbstractMatrix,
                        mesh_old::AbstractMesh, q_old::AbstractVector{Told},
                        mesh_new::AbstractMesh, q_new::AbstractVector{Tnew}
                        ) where {Told, Tnew}

  q_el_old = zeros(Told, mesh_old.numDofPerNode, mesh_old.numNodesPerElement)
  q_el_new = zeros(Tnew, mesh_new.numDofPerNode, mesh_new.numNodesPerElement)
  for i=1:mesh_old.numEl

    @simd for j=1:mesh_old.numNodesPerElement
      @simd for k=1:mesh_old.numDofPerNode
        q_el_old[k, j] = q_old[mesh_old.dofs[k, j, i]]
      end
    end

    smallmatmatT!(q_el_old, interp_op, q_el_new)

    @simd for j=1:mesh_new.numNodesPerElement
      @simd for k=1:mesh_new.numDofPerNode
        q_new[mesh_new.dofs[k, j, i]] = q_el_new[k, j]
      end
    end

  end  # end i

  return nothing
end








