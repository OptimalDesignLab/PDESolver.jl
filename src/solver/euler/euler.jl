# this file contains the functions to evaluate the right hand side of the weak form in the pdf

function evalVolumeIntegrals(mesh::AbstractMesh, operator::SBPOperator, eqn::EulerEquation, u::AbstractVector, u0::AbstractVector)
# evaluate all the integrals over the element (but not the boundary)
# does not do boundary integrals
# mesh : a mesh type, used to access information about the mesh
# operator : an SBP operator, used to get stiffness matricies and stuff
# eqn : a type that holds some constants needed to evaluate the equation
#	: also used for multiple dispatch
# u : solution vector to be populated (mesh.numDof entries)
# u0 : solution vector at previous timesteps (mesh.numDof entries)


# do calculations here


return nothing 

end

function evalBoundaryIntegrals(mesh::AbstractMesh, operator::SBPOperator, eqn::EulerEquation, u::AbstractVector, u0::AbstractVector)
# evaluate all the integrals over the boundary
# does not do boundary integrals
# mesh : a mesh type, used to access information about the mesh
# operator : an SBP operator, used to get stiffness matricies and stuff
# eqn : a type that holds some constants needed to evaluate the equation
#	: also used for multiple dispatch
# u : solution vector to be populated (mesh.numDof entries), partially populated by evalVolumeIntegrals
# u0 : solution vector at previous timesteps (mesh.numDof entries)


# do calculations here

return nothing

end

function applyMassMatrixInverse(mesh::AbstractMesh, operator::SBPOperator, eqn::EulerEquation, u::AbstractVector, u0::AbstractVector)
# apply the inverse of the mass matrix to the entire solution vector
# this is a not too memory inefficient way

numEl = getNumEl(mesh)
nnodes = operator.numnodes
dofpernode = getNumDofPerNode(mesh)
for i=1:numEl
  dofnums_i = getGlobalNodeNumbers(mesh, i)  # get dof nums for this element
  for j=1:nnodes
    for k=1:dofpernode
      dofnum_k = dofnums_i[k,j]
      u[dofnum_k] /= sbp.w[j]
    end
  end
end
  


return nothing 

end


