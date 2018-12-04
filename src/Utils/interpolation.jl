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
