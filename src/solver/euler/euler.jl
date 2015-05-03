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
ICZero(mesh, sbp, eqn, u0)
x = zeros(Float64,(2,sbp.numnodes,getNumEl(mesh))); # nodal Coordinates of the marix
for i = 1:getNumEl(mesh)
  vtxcoord = getElementVertCoords(mesh, [i]);
  vtxcoord = squeeze(vtxcoord,3);
  vtxcoord = vtxcoord[1:2,:]
  vtxcoord = vtxcoord'
  x[:,:,i] = calcnodes(sbp, vtxcoord);
end

dxidx = zeros(Float64, (2,2,sbp.numnodes,getNumEl(mesh))); # Jacobian Matrix
jac = zeros(Float64, (sbp.numnodes,getNumEl(mesh))); # Determinant of the Jacobian Matrix
mappingjacobian!(sbp, x, dxidx, jac) # Get the Jocabian for transformation between actual and iso-parametric space

# Get The boundary faces

F1 = zeros(Float64, 2,sbp.numnodes,getNumEl(mesh))
F2 = zeros(Float64, 2,sbp.numnodes,getNumEl(mesh))
F3 = zeros(Float64, 2,sbp.numnodes,getNumEl(mesh))
F4 = zeros(Float64, 2,sbp.numnodes,getNumEl(mesh))

return nothing

end

function applyMassMatrixInverse(mesh::AbstractMesh, operator::SBPOperator, eqn::EulerEquation, u::AbstractVector, u0::AbstractVector)
# apply the inverse of the mass matrix to the entire solution vector
# this is a good, memory efficient implimentation

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

# some helper functions

function getF1(mesh::AbstractMesh, operator::SBPOperator, eqn::EulerEquation, u0::AbstractVector, element::Integer, f1::AbstractVector)
# gets the vector F1 (see weak form derivation) for a particular element
# for linear triangles, the size of F1 is 3 nodes * 4 dof per node = 12 entries
# element : number of element to fetch F1 for
# f1 : vector (of length 12) to populate with F1.  This vector is overwritten

println("entered getF1")
println("element number = ", element)

dofnums = getGlobalNodeNumbers(mesh, element)
println("dofnums = \n", dofnums)
u_vals = zeros(4)  # hold u0 values for a single node

for i=1:3  # loop over nodes
  println("at node ", i)
  u_vals = u0[dofnums[:,i]]  # get the u0 values
  println("u_vals = \n", u_vals)

  # calculate pressure
  internal_energy = u_vals[4]/u_vals[1] - 0.5*(u_vals[2]^2 + u_vals[3]^2)/(u_vals[1]^2)
  pressure = u_vals[1]*eqn.R*internal_energy/eqn.cv
  println("internal_energy = ", internal_energy, " , pressure = ", pressure)

  # calculate F1 for this node
  start_index = 4*(i-1) + 1
  f1[start_index] = u_vals[2]  # f1_1 (density equation)
  f1[start_index + 1] = (u_vals[2]^2)/u_vals[1] + pressure  # f1_2 (u1 equation)
  f1[start_index + 2] = (u_vals[2]*u_vals[3])/u_vals[1]  # f1_3 (u2 equation)
  f1[start_index + 3] = (u_vals[4] + pressure)*u_vals[2]/u_vals[1] # f1_4 (energy equation)

  print("\n")
end

println("F1 = \n", f1)

return nothing

end

function getF2(mesh::AbstractMesh, operator::SBPOperator, eqn::EulerEquation, u0::AbstractVector, element::Integer, f2::AbstractVector)
# gets the vector F2 (see weak form derivation) for a particular element
# for linear triangles, the size of F2 is 3 nodes * 4 dof per node = 12 entries
# element : number of element to fetch F2 for
# f2 : vector (of length 12) to populate with F2.  This vector is overwritten

println("entered getF2")
println("element number = ", element)

dofnums = getGlobalNodeNumbers(mesh, element)
println("dofnums = \n", dofnums)
u_vals = zeros(4)  # hold u0 values for a single node

for i=1:3  # loop over nodes
  println("at node ", i)
  u_vals = u0[dofnums[:,i]]  # get the u0 values
  println("u_vals = \n", u_vals)

  # calculate pressure
  internal_energy = u_vals[4]/u_vals[1] - 0.5*(u_vals[2]^2 + u_vals[3]^2)/(u_vals[1]^2)
  pressure = u_vals[1]*eqn.R*internal_energy/eqn.cv
  println("internal_energy = ", internal_energy, " , pressure = ", pressure)

  # calculate F1 for this node
  start_index = 4*(i-1) + 1
  f2[start_index] = u_vals[3]  # f2_1 (density equation)
  f2[start_index + 1] = (u_vals[2]*u_vals[3])/u_vals[1] # f2_2 (u1 equation)
  f2[start_index + 2] = (u_vals[3]^2)/u_vals[1] + pressure  # f2_3 (u2 equation)
  f2[start_index + 3] = (u_vals[4] + pressure)*u_vals[3]/u_vals[1] # f2_4 (energy equation)

  print("\n")
end

println("F2 = \n", f2)

return nothing

end

function assembleU(vec::AbstractVector, element::Integer, u::AbstractVector)
# assembles a vector vec (of size 12, coresponding to solution values for an element), into the global solution vector u
# element specifies which element number the number in vec belong to
println("entered assembleU")
println("element = ", element)
println("vec = \n", vec)
dofnums = getGlobalNodeNumbers(mesh, element)

for i=1:3  # loop over nodes
  dofnums_i = dofnums[:,i]
  start_index = 4*(i-1) + 1
  
  u[dofnums_i] += vec[start_index:(start_index+3)]
end

println("u = \n", u)

return nothing

end

function assembleU(vec::AbstractVector, element::Integer, component::Integer, u::AbstractVector)
# assembles a vector vec (of size 3, corresponding to the solution for one degree of freedom at each node) into the global solution vector u
# element specifies which element number the numbers in vec belong to
#  component specifies which dof of each node (1,2,3, or 4)

println("entered assembleU")
println("element = ", element)
println("component = ", component)
println("vec = \n", vec)
dofnums = getGlobalNodeNumbers(mesh, element)
println("dofnums = ", dofnums)

dofnums_comp = dofnums[component,:]
u[dofnums_comp.'] += vec

println("u = \n", u)

return nothing

end
