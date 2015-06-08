# this file contains the functions to evaluate the right hand side of the weak form in the pdf
# SL stands for semi-linear form, contents inside the inv(M)(...) parenthesis on pg. 7 of Jared's derivation


#using SummationByParts
#using PdePumiInterface
#include("sbp_interface.jl")
#include("stabilization.jl")

export evalEuler


# this function is what the timestepper calls
function evalEuler(t, SL0, SL, extra_args)
# SL is popualted with du/dt
# SL0 is q at previous timestep
# t is current timestep
# extra_args is unpacked into object needed to evaluation equation


mesh = extra_args[1]
sbp = extra_args[2]
eqn = extra_args[3]


@time dataPrep(mesh, sbp, eqn, SL, SL0)
println("dataPrep @time printed above")
@time evalVolumeIntegrals(mesh, sbp, eqn, SL, SL0)
println("volume integral @time printed above")
@time evalBoundaryIntegrals(mesh, sbp, eqn, SL, SL0)
println("boundary integral @time printed above")



@time addEdgeStabilize(mesh, sbp, eqn, SL, SL0)
println("edge stabilizing @time printed above")


@time assembleSolution(mesh, eqn, SL)
println("assembly @time printed above")

@time applyMassMatrixInverse(eqn, SL)
println("Minv @time printed above")

#applyDissipation(mesh, sbp, eqn, SL, SL0)



#print current error
err_norm = norm(SL)/mesh.numDof
print(" ", err_norm)


return nothing
#return SL

end  # end evalEuler




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

boundaryintegrate2!(sbp, mesh.bndryfaces, eqn.q, mesh.coords, mesh.dxidx, isentropicVortexBC, eqn.res, mesh, eqn)
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

#   println("====== Entering edge stabilize ======")
  


  # u argument here is SL in a different format
  edgestabilize!(sbp, mesh.interfaces, eqn.q, mesh.coords, mesh.dxidx, mesh.jac, eqn.edgestab_alpha, stabscale, eqn.res, mesh, eqn)

#  println("==== end of addEdgeStabilize ====")



  return nothing

end



# this function is deprecated
function applyDissipation(mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerEquation, SL::AbstractVector, SL0::AbstractVector)
# apply sketchy dissipation scheme
# this doesn't work now that the code has been reorganized

  # get some data
#  u, x, dxidx, jac, res_xi, interfaces = dataPrep(mesh, sbp, eqn, SL, SL0)

  u = mesh.q
  x = mesh.coords
  dxidx = mesh.dxidx
  jac = mesh.jac
  res_xi = zeros(eqn.res)

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
    for j=1:mesh.numNodesPerElement
      for k=1:4  # loop over dofs on the node
	dofnum_k = mesh.dofs[k, j, i]
	SL[dofnum_k] += eqn.res[k,j,i]
      end
    end
  end
   
  return nothing
end



