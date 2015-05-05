# this file contains the functions to evaluate the right hand side of the weak form in the pdf
# SL stands for semi-linear form, contents inside the inv(M)(...) parenthesis on pg. 7 of Jared's derivation

function evalVolumeIntegrals(mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerEquation, SL::AbstractVector, SL0::AbstractVector)
# evaluate all the integrals over the elements (but not the boundary)
# does not do boundary integrals
# mesh : a mesh type, used to access information about the mesh
# sbp : an SBP operator, used to get stiffness matricies and stuff
# eqn : a type that holds some constants needed to evaluate the equation
#	: also used for multiple dispatch
# SL : solution vector to be populated (mesh.numDof entries)
# SL0 : solution vector at previous timesteps (mesh.numDof entries)


  println("Evaluating volume integrals")
  
  numEl = getNumEl(mesh)
  
  p = 1           # linear elements
  
  # this is commented out because it gets zeroed out earlier, and we don't want to throw out other work
#   SL = zeros(Float64,mesh.numDof)
  
  for element = 1:numEl
  
    dofnums = getGlobalNodeNumbers(mesh, element)

  #   sbp = TriSBP{Float64}(degree=p)

    # From PUMI, PdePumiInterface.jl, line ~142
    x_3D = getElementVertCoords(mesh,[element])

    # need x_3D -> x_2D
    # required by mappingjacobian, dimensions:
    #   1: coord
    #   2: node
    #   3: elem
    x = sub(x_3D, 1:2, 1:3, :)

    # required by mappingjacobian, dimensions:
    #   1: ref coord
    #   2: phys coord
    #   3: node
    #   4: elem
    dxidx = zeros(Float64, 2, 2, 3, 1)

    # required by mappingjacobian, dimensions:
    #   1: node
    #   2: elem
    jac = zeros(Float64, 3, 1)
    
    mappingjacobian!(sbp, x, dxidx, jac)
  
    # cheating a little, constant across all nodes since linear, so use node #1 for all
    dxi_dx = dxidx[1,1,1,1]
    dxi_dy = dxidx[1,2,1,1]
    deta_dx = dxidx[2,1,1,1]
    deta_dy = dxidx[2,2,1,1]
  
    F1 = zeros(Float64, 12)
    F2 = zeros(Float64, 12)

    getF1(mesh, sbp, eqn, SL0, element, F1)
    getF2(mesh, sbp, eqn, SL0, element, F2)
  
    # du/dt = inv(M)*(volumeintegrate!(f1/jac) + volumeintegrate!(f2/jac)
  
  
    # ------------- setting source terms
    src_rho = 0.0
    src_rhou = 0.0
    src_rhov = 0.0
    src_E = 0.0

    # f's need 2 indices for 2D, and 3 entries corresponding to the vertices of a tri element
    # little f's: source terms; f1, for ex, source @ all 3 nodes for 1st conserved variable
    f1 = zeros(Float64,3,1)
    f2 = zeros(Float64,3,1)
    f3 = zeros(Float64,3,1)
    f4 = zeros(Float64,3,1)

    # ------------- populating source terms
    for node = 1:3
      f1[node,1] = src_rho
      f2[node,1] = src_rhou
      f3[node,1] = src_rhov
      f4[node,1] = src_E
    end

    result1 = zeros(size(f1))
    result2 = zeros(size(f2))
    result3 = zeros(size(f3))
    result4 = zeros(size(f4))

    # volumeintegrate! only does += to the result, so no need for result1, result2, etc
    # it is assumed that `u` is a rank-2 array, with the first dimension for the local-node index, and the second dimension for the element index.
    volumeintegrate!(sbp, f1/jac[1,1], result1)
    volumeintegrate!(sbp, f2/jac[1,1], result2)
    volumeintegrate!(sbp, f3/jac[1,1], result3)
    volumeintegrate!(sbp, f4/jac[1,1], result4)
  
    source_result = zeros(Float64,12)
    i = 1
    for node = 1:3
      source_result[i] = result1[node,1]
      source_result[i+1] = result2[node,1]
      source_result[i+2] = result3[node,1]
      source_result[i+3] = result4[node,1]
      i = i+4
    end
     
    flux_result = eqn.bigQT_xi*(F1*dxi_dx + F2*dxi_dy) + eqn.bigQT_eta*(F1*deta_dx + F2*deta_dy)

    println("F1hat: \n",F1)
    println("F2hat: \n",F2)
    println("bigQT_xi: \n",eqn.bigQT_xi)
    println("bigQT_eta: \n",eqn.bigQT_eta)
    println("dofnums: \n",dofnums)
    println("dxi_dx: \n",dxi_dx)
    println("dxi_dy: \n",dxi_dy)
    println("deta_dx: \n",deta_dx)
    println("deta_dy: \n",deta_dy)
  
    F1 = zeros(Float64, 12)
    F2 = zeros(Float64, 12)
  
    SL_el = source_result + flux_result
  
    # function assembleU(vec::AbstractVector, element::Integer, u::AbstractVector)
    # assembles a vector vec (of size 12, coresponding to solution values for an element), into the global solution vector u
    # element specifies which element number the number in vec belong to
    assembleSL(SL_el, element, SL)
  
  end
  
  return nothing 

end
#------------- end of evalVolumeIntegrals


function evalBoundaryIntegrals(mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerEquation, u::AbstractVector, u0::AbstractVector)
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
#------------- end of evalBoundaryIntegrals

function addEdgeStabilize(mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerEquation, SL::AbstractVector, SL0::AbstractVector)

  # alpha calculated like in edgestabilize! documentation
  # stabscale (U+a)*gamma*h^2 where U=u*n, where u is the velocity 
  #   (remember to scale by rho) and n is the unit normal vector, from nrm->dxidx, then scaled by length
  # ifaces needs to be calculated
  # x needs to be passed
  # need to clarify u vs res. maybe change the u variable name to semilinear 

  # u argument here is SL in a different format
#   edgestabilize!(sbp, ifaces, u, x, dxidx, jac, alpha, stabscale, res)



  return nothing

end



function applyMassMatrixInverse(mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerEquation, SL::AbstractVector, SL0::AbstractVector)
# apply the inverse of the mass matrix to the entire solution vector
# this is a good, memory efficient implimentation

numEl = getNumEl(mesh)
nnodes = sbp.numnodes
dofpernode = getNumDofPerNode(mesh)
for i=1:numEl
  dofnums_i = getGlobalNodeNumbers(mesh, i)  # get dof nums for this element
  for j=1:nnodes
    for k=1:dofpernode
      dofnum_k = dofnums_i[k,j]
      SL[dofnum_k] /= sbp.w[j]
    end
  end
end
  
return nothing 

end

# some helper functions

function getF1(mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerEquation, SL0::AbstractVector, element::Integer, f1::AbstractVector)
# gets the vector F1 (see weak form derivation) for a particular element
# for linear triangles, the size of F1 is 3 nodes * 4 dof per node = 12 entries
# element : number of element to fetch F1 for
# f1 : vector (of length 12) to populate with F1.  This vector is overwritten

println("entered getF1")
println("element number = ", element)

dofnums = getGlobalNodeNumbers(mesh, element)
println("dofnums = \n", dofnums)
SL_vals = zeros(4)  # hold SL0 values for a single node

for i=1:3  # loop over nodes
  println("at node ", i)
  SL_vals = SL0[dofnums[:,i]]  # get the SL0 values
  println("SL_vals = \n", SL_vals)

  # calculate pressure
#  internal_energy = SL_vals[4]/SL_vals[1] - 0.5*(SL_vals[2]^2 + SL_vals[3]^2)/(SL_vals[1]^2)
#  pressure = SL_vals[1]*eqn.R*internal_energy/eqn.cv
#  println("internal_energy = ", internal_energy, " , pressure = ", pressure)
  pressure = calcPressure(SL_vals, eqn)

  # calculate F1 for this node
  start_index = 4*(i-1) + 1
  f1[start_index] = SL_vals[2]  # f1_1 (density equation)
  f1[start_index + 1] = (SL_vals[2]^2)/SL_vals[1] + pressure  # f1_2 (momentum-x equation)
  f1[start_index + 2] = (SL_vals[2]*SL_vals[3])/SL_vals[1]  # f1_3 (momentum-y equation)
  f1[start_index + 3] = (SL_vals[4] + pressure)*SL_vals[2]/SL_vals[1] # f1_4 (energy equation)

  print("\n")
end

println("F1 = \n", f1)

return nothing

end

function getF2(mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerEquation, SL0::AbstractVector, element::Integer, f2::AbstractVector)
# gets the vector F2 (see weak form derivation) for a particular element
# for linear triangles, the size of F2 is 3 nodes * 4 dof per node = 12 entries
# element : number of element to fetch F2 for
# f2 : vector (of length 12) to populate with F2.  This vector is overwritten

println("entered getF2")
println("element number = ", element)

dofnums = getGlobalNodeNumbers(mesh, element)
println("dofnums = \n", dofnums)
SL_vals = zeros(4)  # hold SL0 values for a single node

for i=1:3  # loop over nodes
  println("at node ", i)
  SL_vals = SL0[dofnums[:,i]]  # get the SL0 values
  println("SL_vals = \n", SL_vals)

  # calculate pressure
#  internal_energy = SL_vals[4]/SL_vals[1] - 0.5*(SL_vals[2]^2 + SL_vals[3]^2)/(SL_vals[1]^2)
#  pressure = SL_vals[1]*eqn.R*internal_energy/eqn.cv
#  println("internal_energy = ", internal_energy, " , pressure = ", pressure)
  pressure = calcPressure(SL_vals, eqn)

  # calculate F1 for this node
  start_index = 4*(i-1) + 1
  f2[start_index] = SL_vals[3]  # f2_1 (density equation)
  f2[start_index + 1] = (SL_vals[2]*SL_vals[3])/SL_vals[1] # f2_2 (momentum-x equation)
  f2[start_index + 2] = (SL_vals[3]^2)/SL_vals[1] + pressure  # f2_3 (momentum-y equation)
  f2[start_index + 3] = (SL_vals[4] + pressure)*SL_vals[3]/SL_vals[1] # f2_4 (energy equation)

  print("\n")
end

println("F2 = \n", f2)

return nothing

end

function calcPressure(SL_vals::AbstractVector, eqn::EulerEquation)
  # calculate pressure for a node
  # SL is a vector of length 4

  internal_energy = SL_vals[4]/SL_vals[1] - 0.5*(SL_vals[2]^2 + SL_vals[3]^2)/(SL_vals[1]^2)
  pressure = SL_vals[1]*eqn.R*internal_energy/eqn.cv
  println("internal_energy = ", internal_energy, " , pressure = ", pressure)

  return pressure
end



function assembleSL(vec::AbstractVector, element::Integer, SL::AbstractVector)
# assembles a vector vec (of size 12, coresponding to solution values for an element), into the global solution vector SL
# element specifies which element number the number in vec belong to
println("entered assembleU")
println("element = ", element)
println("vec = \n", vec)
dofnums = getGlobalNodeNumbers(mesh, element)

for i=1:3  # loop over nodes
  dofnums_i = dofnums[:,i]
  start_index = 4*(i-1) + 1
  
  SL[dofnums_i] += vec[start_index:(start_index+3)]
end

println("SL = \n", SL)

return nothing

end

function assembleSL(vec::AbstractVector, element::Integer, component::Integer, SL::AbstractVector)
# assembles a vector vec (of size 3, corresponding to the solution for one degree of freedom at each node) into the global solution vector SL
# element specifies which element number the numbers in vec belong to
#  component specifies which dof of each node (1,2,3, or 4)

println("entered assembleU")
println("element = ", element)
println("component = ", component)
println("vec = \n", vec)
dofnums = getGlobalNodeNumbers(mesh, element)
println("dofnums = ", dofnums)

dofnums_comp = dofnums[component,:]
SL[dofnums_comp.'] += vec

println("SL = \n", SL)

return nothing

end
