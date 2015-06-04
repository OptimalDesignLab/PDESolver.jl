# this file contains the functions to evaluate the right hand side of the weak form in the pdf
# SL stands for semi-linear form, contents inside the inv(M)(...) parenthesis on pg. 7 of Jared's derivation


using SummationByParts
using PdePumiInterface
using Equation

function dataPrep(mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerEquation, SL::AbstractVector, SL0::AbstractVector)
# gather up all the data needed to do vectorized operatinos on the mesh
# linear elements only

# need: u (previous timestep solution), x (coordinates), dxidx, jac, res, array of interfaces

#println("Entered dataPrep()")

  u = eqn.q
  F_xi = eqn.F_xi

  F_eta = eqn.F_eta
  fill!(eqn.res, 0.0)


  for i=1:mesh.numEl  # loop over elements
    for j = 1:sbp.numnodes
      for k=1:4
	dofnum_k = mesh.dofs[k, j, i]
	u[k, j, i] = SL0[dofnum_k]
      end
    end
  end


  # calculate fluxes
  getEulerFlux(eqn, eqn.q, mesh.dxidx, F_xi, F_eta)
#  println("getEulerFlux @time printed above")




#  println("finished dataPrep()")
  return nothing
end # end function dataPrep

function evalVolumeIntegrals(mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerEquation, SL::AbstractVector, SL0::AbstractVector)
# evaluate all the integrals over the elements (but not the boundary)
# does not do boundary integrals
# mesh : a mesh type, used to access information about the mesh
# sbp : an SBP operator, used to get stiffness matricies and stuff
# eqn : a type that holds some constants needed to evaluate the equation
#	: also used for multiple dispatch
# SL : solution vector to be populated (mesh.numDof entries)
# SL0 : solution vector at previous timesteps (mesh.numDof entries)

#  println("size eqn.F_xi = ", size(eqn.F_xi), " size eqn.res = ", size(eqn.res), " sbp.numnodes = ", sbp.numnodes)
  weakdifferentiate!(sbp, 1, eqn.F_xi, eqn.res, trans=true)
  weakdifferentiate!(sbp, 2, eqn.F_eta, eqn.res, trans=true)


  # need source term here



end
#------------- end of evalVolumeIntegrals


function evalBoundaryIntegrals(mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerEquation, SL::AbstractVector, SL0::AbstractVector)
# evaluate all the integrals over the boundary
# does not do boundary integrals
# mesh : a mesh type, used to access information about the mesh
# operator : an SBP operator, used to get stiffness matricies and stuff
# eqn : a type that holds some constants needed to evaluate the equation
#	: also used for multiple dispatch
# u : solution vector to be populated (mesh.numDof entries), partially populated by evalVolumeIntegrals
# u0 : solution vector at previous timesteps (mesh.numDof entries)


#  println("====== start of evalBoundaryIntegrals ======")
#boundaryintegrate!(sbp, bndryfaces, u, x, dxidx, isentropicVortexBC, result)

boundaryintegrate2!(sbp, mesh.bndryfaces, eqn.q, mesh.coords, mesh.dxidx, isentropicVortexBC, eqn.res)
#boundaryintegrate!(sbp, bndryfaces, u, x, dxidx, rho1Energy2BC, result)


#  println("==== end of evalBoundaryIntegrals ====")


  return nothing

end
#------------- end of evalBoundaryIntegrals

# This function adds edge stabilization to a residual using Prof. Hicken's edgestabilize! in SBP
function addEdgeStabilize(mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerEquation, SL::AbstractVector, SL0::AbstractVector)

#  println("==== start of addEdgeStabilize ====")
  # alpha calculated like in edgestabilize! documentation
  # stabscale (U+a)*gamma*h^2 where U=u*n, where u is the velocity 
  #   (remember to scale by rho) and n is the unit normal vector, from nrm->dxidx, then scaled by length
  # ifaces needs to be calculated
  # x needs to be passed
  # need to clarify u vs res. maybe change the u variable name to semilinear 

#   println("====== Entering edge stabilize ======")
  

    # dxidx dimensions:
    #   1: ref coord
    #   2: phys coord
    #   3: node
    #   4: elem

#  u, x, dxidx, jac, result, interfaces = dataPrep(mesh, sbp, eqn, SL, SL0)
  numEl = getNumEl(mesh)

  alpha = zeros(Float64,2,2,sbp.numnodes,numEl)
  dxidx = mesh.dxidx # referency only, for code compatability
  jac = mesh.jac  # reference only

  # calculating alpha, required by edgestabilize!
  # this canbe compuated apriori
  for k = 1:numEl
    for i = 1:sbp.numnodes                                                                              
      for di1 = 1:2                                                                                     
        for di2 = 1:2                                                                                   
          alpha[di1,di2,i,k] = (dxidx[di1,1,i,k].*dxidx[di2,1,i,k] + dxidx[di1,2,i,k].*dxidx[di2,2,i,k])*jac[i,k]
        end                                                                                             
      end                                                                                               
    end                                                                                                 
  end


  function stabscale(u, dxidx, nrm)

#     println("==== entering stabscale ====")

    # grabbing conserved variables
    rho = u[1]
    vel_x = u[2]/rho
    vel_y = u[3]/rho
    Energy = u[4]

    # from JC's code below, eqn should still be in scope
    pressure = calcPressure(u, eqn)

    # solved eqn for e: E = rho*e + (1/2)*rho*u^2
    vel_squared = vel_x^2 + vel_y^2
    energy = Energy/rho - (1/2)*vel_squared

    # gamma stored in EulerEquation type
    gamma = eqn.gamma

#     println("pressure: ",pressure)
#     println("gamma: ",gamma)
#     println("rho: ",rho)
    # ideal gas law
    speed_sound = sqrt((gamma*pressure)/rho)

    # choice for edge stabilization constant: 
    #   refer to email from JH, 20150504:
    #   Anthony: there is little guidance in the literature regarding 
    #     gamma for the Euler equations.  I suggest starting with 
    #     gamma = 0.01.  If that fails (with a cfl of 1.0), then decrease 
    #     it by an order of magnitude at at time until RK is stable.  
    #     Once you find a value that works, try increasing it slowly.
    edge_stab_gamma = -0.01  # default
#     edge_stab_gamma = 0.0 
#     edge_stab_gamma = 0.00001

    # edge lengths component wise
    h_x = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
    h_y = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]

    # edge length
    h = sqrt(h_x^2 + h_y^2)

    # scaled velocity scalar
#     U = vel_x*(nrm[1]/h) + vel_y*(nrm[2]/h)
    U = vel_x*(h_x/h) + vel_y*(h_y/h)

#     return (U + speed_sound)*edge_stab_gamma*h^2
    return (abs(U) + speed_sound)*edge_stab_gamma*h^2

  end

  # u argument here is SL in a different format
  edgestabilize!(sbp, mesh.interfaces, eqn.q, mesh.coords, mesh.dxidx, mesh.jac, alpha, stabscale, eqn.res)

#  println("==== end of addEdgeStabilize ====")



  return nothing

end



# this function is deprecated
function applyDissipation(mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerEquation, SL::AbstractVector, SL0::AbstractVector)
# apply sketchy dissipation scheme

  # get some data
  u, x, dxidx, jac, res_xi, interfaces = dataPrep(mesh, sbp, eqn, SL, SL0)

  res2_xi = zeros(res_xi)


  res_eta= zeros(res_xi)
  res2_eta = zeros(res_xi)

  weakdifferentiate!(sbp, 1, u, res_xi)

  mass_matrix = sbp.w

  # apply inverse mass matrix twice
  res_xi[:, 1, :] /= mass_matrix[1]*mass_matrix[1]
  res_xi[:, 2, :] /= mass_matrix[2]*mass_matrix[2]
  res_xi[:, 3, :] /= mass_matrix[3]*mass_matrix[3]
 

   weakdifferentiate!(sbp, 1, res_xi, res2_xi; trans=true)

   # assemble
for element = 1:mesh.numEl
  for node = 1:sbp.numnodes
    vec = -res2_xi[:,node,element]
    assembleSLNode(vec, element, node, SL)
  end
#   println("- element #: ",element,"   result[:,:,element]:",result[:,:,element])
end


  # now do it in the eta direction

 weakdifferentiate!(sbp, 2, u, res_eta)

  mass_matrix = sbp.w

  # apply inverse mass matrix twice
  res_eta[:, 1, :] /= mass_matrix[1]*mass_matrix[1]
  res_eta[:, 2, :] /= mass_matrix[2]*mass_matrix[2]
  res_eta[:, 3, :] /= mass_matrix[3]*mass_matrix[3]
 
   weakdifferentiate!(sbp, 2, res_eta, res2_eta; trans=true)

   # assemble
for element = 1:mesh.numEl
  for node = 1:sbp.numnodes
    vec = -res2_eta[:,node,element]
    assembleSLNode(vec, element, node, SL)
  end
#   println("- element #: ",element,"   result[:,:,element]:",result[:,:,element])
end




  return nothing

end





function applyMassMatrixInverse(eqn::EulerEquation, SL::AbstractVector)
# apply the inverse mass matrix stored eqn to SL

#  SL .*= eqn.Minv  # this gives wrong answer


  ndof = length(SL)
  for i=1:ndof
    SL[i] *= eqn.Minv[i]
  end

  return nothing
end

#=
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
=#
# some helper functions

function getEulerFlux{T}(eqn::EulerEquation, q::AbstractArray{T,1}, dir::AbstractArray{T,1},  F::AbstractArray{T,1})
# calculates the Euler flux in a particular direction at a point
# eqn is the equation type
# q is the vector (of length 4), of the conservative variables at the point
# dir is a vector of length 2 that specifies the direction
# F is populated with the flux (is a vector of length 4)

# once the Julia developers fix slice notation and speed up subarrays, we can make a faster
# vectorized version of this

  press = calcPressure(q, eqn)
  U = (q[2]*dir[1] + q[3]*dir[2])/q[1]
  F[1] = q[1]*U
  F[2] = q[2]*U + dir[1]*press
  F[3] = q[3]*U + dir[2]*press
  F[4] = (q[4] + press)*U
 
  return nothing

end


function getEulerFlux{T}(eqn::EulerEquation, q::AbstractArray{T,3}, dxidx::AbstractArray{T,4},  F_xi::AbstractArray{T,3}, F_eta::AbstractArray{T,3})
# calculates the Euler flux for every node in the xi and eta directions
# eqn is the equation type
# q is the 3D array (4 by nnodes per element by nel), of the conservative variables
# dxidx is the 4D array (2 by 2 x nnodes per element by nel) that specifies the direction of xi and eta in each element (output from mappingjacobian!)
# F_xi is populated with the flux in xi direction (same shape as q)
# F_eta is populated with flux in eta direction

# once the Julia developers fix slice notation and speed up subarrays, we won't have to 
# vectorize like this (can calculate flux one node at a time inside a dedicated function

(ncomp, nnodes, nel) = size(q)  # get sizes of things

  for i=1:nel  # loop over elements
    for j=1:nnodes  # loop over nodes within element
      # get direction vector components (xi direction)
      nx = dxidx[1, 1, j, i]
      ny = dxidx[1, 2, j, i]
      # calculate pressure 
      press = (eqn.gamma-1)*(q[4, j, i] - 0.5*(q[2, j, i]^2 + q[3, j, i]^2)/q[1, j, i])

      # calculate flux in xi direction
      # hopefully elements of q get stored in a register for reuse in eta direction
      U = (q[2, j, i]*nx + q[3, j, i]*ny)/q[1, j, i]
      F_xi[1, j, i] = q[1, j, i]*U
      F_xi[2, j, i] = q[2, j, i]*U + nx*press
      F_xi[3, j, i] = q[3, j, i]*U + ny*press
      F_xi[4, j, i] = (q[4, j, i] + press)*U

      # get direction vector components (eta direction)
      nx = dxidx[2, 1, j, i]
      ny = dxidx[2, 2, j, i]

      # calculate xi flux
      U = (q[2, j, i]*nx + q[3, j, i]*ny)/q[1, j, i]
      F_eta[1, j, i] = q[1, j, i]*U
      F_eta[2, j, i] = q[2, j, i]*U + nx*press
      F_eta[3, j, i] = q[3, j, i]*U + ny*press
      F_eta[4, j, i] = (q[4, j, i] + press)*U
    end
  end


 
  return nothing

end



#=
function getF1(mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerEquation, SL0::AbstractVector, element::Integer, f1::AbstractVector)
# gets the vector F1 (see weak form derivation) for a particular element
# for linear triangles, the size of F1 is 3 nodes * 4 dof per node = 12 entries
# element : number of element to fetch F1 for
# f1 : vector (of length 12) to populate with F1.  This vector is overwritten

#   println("entered getF1")
#   println("element number = ", element)

  dofnums = getGlobalNodeNumbers(mesh, element)
#   println("dofnums = \n", dofnums)
  SL_vals = zeros(4)  # hold SL0 values for a single node

  for i=1:3  # loop over nodes
#     println("at node ", i)
    SL_vals = SL0[dofnums[:,i]]  # get the SL0 values
#     println("SL_vals = \n", SL_vals)

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

#     print("\n")
  end

#   println("F1 = \n", f1)

  return nothing

end

function getF2(mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerEquation, SL0::AbstractVector, element::Integer, f2::AbstractVector)
# gets the vector F2 (see weak form derivation) for a particular element
# for linear triangles, the size of F2 is 3 nodes * 4 dof per node = 12 entries
# element : number of element to fetch F2 for
# f2 : vector (of length 12) to populate with F2.  This vector is overwritten

#   println("entered getF2")
#   println("element number = ", element)

  dofnums = getGlobalNodeNumbers(mesh, element)
#   println("dofnums = \n", dofnums)
  SL_vals = zeros(4)  # hold SL0 values for a single node

  for i=1:3  # loop over nodes
#     println("at node ", i)
    SL_vals = SL0[dofnums[:,i]]  # get the SL0 values
#     println("SL_vals = \n", SL_vals)

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

#     print("\n")
  end

#   println("F2 = \n", f2)

  return nothing

end
=#
function calcPressure(q::AbstractVector, eqn::EulerEquation)
  # calculate pressure for a node
  # q is a vector of length 4 of the conservative variables

#  internal_energy = SL_vals[4]/SL_vals[1] - 0.5*(SL_vals[2]^2 + SL_vals[3]^2)/(SL_vals[1]^2)
#  pressure = SL_vals[1]*eqn.R*internal_energy/eqn.cv

  return  (eqn.gamma-1)*(q[4] - 0.5*(q[2]^2 + q[3]^2)/q[1])
  
#   println("internal_energy = ", internal_energy, " , pressure = ", pressure)


end


function assembleSolution(mesh::AbstractMesh, eqn::EulerEquation, SL::AbstractVector)


for i=1:mesh.numEl  # loop over elements
#  dofnums = getGlobalNodeNumbers(mesh, i)

  for j=1:mesh.numNodesPerElement
    for k=1:4  # loop over dofs on the node
      dofnum_k = mesh.dofs[k, j, i]
      SL[dofnum_k] += eqn.res[k,j,i]
#      SL[dofnums[k, j]] += eqn.res[k, j, i]
    end
  end
end
 


  return nothing
end

#=
function assembleSL(vec::AbstractVector, element::Integer, SL::AbstractVector)
# assembles a vector vec (of size 12, coresponding to solution values for an element), into the global solution vector SL
# element specifies which element number the number in vec belong to
#   println("entered assembleU")
#   println("element = ", element)
#   println("vec = \n", vec)
  dofnums = getGlobalNodeNumbers(mesh, element)

  start_index = 1
  for i=1:3  # loop over nodes
#    dofnums_i = dofnums[:,i]
#    start_index = 4*(i-1) + 1
 
    for j=1:4
      SL[j,i] += vec[start_index]
    end
  end

#   println("SL = \n", SL)

  return nothing

end

function assembleSL(vec::AbstractVector, element::Integer, component::Integer, SL::AbstractVector)
# assembles a vector vec (of size 3, corresponding to the solution for one degree of freedom at each node) into the global solution vector SL
# element specifies which element number the numbers in vec belong to
#  component specifies which dof of each node (1,2,3, or 4)

#   println("entered assembleU")
#   println("element = ", element)
#   println("component = ", component)
#   println("vec = \n", vec)
  dofnums = getGlobalNodeNumbers(mesh, element)
#   println("dofnums = ", dofnums)

#  dofnums_comp = dofnums[component,:]
#  SL[dofnums_comp.'] += vec

  for i=1:3
    SL[dofnums[component, i]] = vec[i]
  end

#   println("SL = \n", SL)

  return nothing

end

function assembleSLNode(vec::AbstractVector, element::Integer, node::Integer, SL::Vector)
# vec is vector of length 4 (conservative variables for a single node
# element is the element number that the node belongs to
# node is the node on the element (1, 2 or 3 for linear triangles)
# SL is global solution vector

  dofnums = getGlobalNodeNumbers(mesh, element)
#  dofnums_n = dofnums[:, node]

#  SL[dofnums_n] += vec

  for i=1:4
    SL[dofnums[i, node]] += vec[i]
  end
  
  return nothing

end  # end function
=#

