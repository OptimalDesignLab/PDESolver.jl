# this file contains the functions to evaluate the right hand side of the weak form in the pdf

function evalVolumeIntegrals(mesh::AbstractMesh, operator::SBPOperator, eqn::EulerEquation, u::AbstractVector, u0::AbstractVector)
# evaluate all the integrals over the elements (but not the boundary)
# does not do boundary integrals
# mesh : a mesh type, used to access information about the mesh
# operator : an SBP operator, used to get stiffness matricies and stuff
# eqn : a type that holds some constants needed to evaluate the equation
#	: also used for multiple dispatch
# u : solution vector to be populated (mesh.numDof entries)
# u0 : solution vector at previous timesteps (mesh.numDof entries)


println("Evaluating volume integrals")

numEl = getNumEl(mesh)

p = 1           # linear elements

u = zeros(Float64,mesh.numDof)
# u_el1 = zeros(Float64,2,

for element = 1:numEl

#   sbp = TriSBP{Float64}(degree=p)
  
  mappingjacobian!(sbp, x, dxidx, jac)

  dxi_dx = dxidx[1,1]
  dxi_dy = dxidx[1,2]
  deta_dx = dxidx[2,1]
  deta_dy = dxidx[2,2]

  getF1(mesh, sbp, eqn, u0, element, F1)
  getF2(mesh, sbp, eqn, u0, element, F2)

  # du/dt = inv(M)*(volumeintegrate!(f1/jac) + volumeintegrate!(f2/jac)

#   for node_ix = 1:numNodes

  src = 0
  f1 = zeros(3,1)
  f2 = zeros(3,1)
  f3 = zeros(3,1)
  f4 = zeros(3,1)

  # volumeintegrate! only does += to the result, so no need for res1, res2, etc
  # it is assumed that `u` is a rank-2 array, with the first dimension for the local-node index, and the second dimension for the element index.
  volumeintegrate!(sbp, f1/jac, res);
  volumeintegrate!(sbp, f2/jac, res);
  volumeintegrate!(sbp, f3/jac, res);
  volumeintegrate!(sbp, f4/jac, res);

#   u_el1 = volumeintegrate!(f1/jac) + volumeintegrate!(f2/jac) + volumeintegrate!(f3/jac) + volumeintegrate!(f4/jac)
  u_el1 = res

#   u_el2 = transpose(sbp.Q)*(F1*dxi_dx + F2*dxi_dy) + transpose(sbp.Q)*(F1*deta_dx + F2*deta_dy)
  u_el2 = eqn.bigQT_xi*(F1*dxi_dx + F2*dxi_dy) + eqn.bigQT_eta*(F1*deta_dx + F2*deta_dy)

  u_el = u_el1 + u_el2

  assembleU(u_el, element, u)

end

# function assembleU(vec::AbstractVector, element::Integer, u::AbstractVector)
# assembles a vector vec (of size 12, coresponding to solution values for an element), into the global solution vector u
# element specifies which element number the number in vec belong to



  

    



#------------------------------------------

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
#  internal_energy = u_vals[4]/u_vals[1] - 0.5*(u_vals[2]^2 + u_vals[3]^2)/(u_vals[1]^2)
#  pressure = u_vals[1]*eqn.R*internal_energy/eqn.cv
#  println("internal_energy = ", internal_energy, " , pressure = ", pressure)
  pressure = calcPressure(u_vals, eqn)

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
#  internal_energy = u_vals[4]/u_vals[1] - 0.5*(u_vals[2]^2 + u_vals[3]^2)/(u_vals[1]^2)
#  pressure = u_vals[1]*eqn.R*internal_energy/eqn.cv
#  println("internal_energy = ", internal_energy, " , pressure = ", pressure)
  pressure = calcPressure(u_vals, eqn)

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

function calcPressure(u_vals::AbstractVector, eqn::EulerEquation)
  # calculate pressure for a node
  # u is a vector of length 4

  internal_energy = u_vals[4]/u_vals[1] - 0.5*(u_vals[2]^2 + u_vals[3]^2)/(u_vals[1]^2)
  pressure = u_vals[1]*eqn.R*internal_energy/eqn.cv
  println("internal_energy = ", internal_energy, " , pressure = ", pressure)

  return pressure
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
