# this file contains the functions to evaluate the right hand side of the weak form in the pdf
# SL stands for semi-linear form, contents inside the inv(M)(...) parenthesis on pg. 7 of Jared's derivation


# Rules: 
# 1. function that take a composite type with abstract fields *cannot* use those
# fields directly (they can pass those fields to other functions however
# This leads directly to a two level structure for the code: high level function
# that take in composite types and low level function that take in arrays and
# perform calculations on them

# at this time, no composite types have abstract fields (all AbstractArrays
# were turned into arrays)

# the reason for this is that the compiler does not compile new version of the 
# function based
# on the types of the fields of a composite type. Passing the fields of the typ
# e to other functions fixes this problem because the fields are now arguments,
# so the compiler can specialize the code

# 2.  Arrays should not be returned from functions.  The caller should allocate
# and array and pass it into the function

# this prevents repeatedly allocating new arrays


#using SummationByParts
#using PdePumiInterface
#include("sbp_interface.jl")
#include("stabilization.jl")

export evalEuler

# Rules for paramaterization:
# Tmsh = mesh data type
# Tsbp = SBP operator data type
# Tsol = equation, SL, SL0 data type


# this function is what the timestepper calls
function evalEuler(t, mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerEquation,  SL0, SL)
# SL is popualted with du/dt
# SL0 is q at previous timestep
# t is current timestep
# extra_args is unpacked into object needed to evaluation equation


@time dataPrep(mesh, sbp, eqn, SL0)
println("dataPrep @time printed above")
evalVolumeIntegrals(mesh, sbp, eqn)
#println("after evalVolumeIntegrals, isnan: ", isnan(eqn.res))
#println("volume integral @time printed above")
evalBoundaryIntegrals(mesh, sbp, eqn)
#println("boundary integral @time printed above")



addEdgeStabilize(mesh, sbp, eqn)
#println("edge stabilizing @time printed above")


assembleSolution(mesh, eqn, SL)
#println("assembly @time printed above")

applyMassMatrixInverse(eqn, SL)
#println("Minv @time printed above")

#applyDissipation(mesh, sbp, eqn, SL, SL0)



#print current error
#err_norm = norm(SL)/mesh.numDof
#print(" ", err_norm)


return nothing
#return SL

end  # end evalEuler



function dataPrep{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator{Tsbp}, eqn::EulerEquation{Tsol}, SL0::AbstractVector{Tsol})
# gather up all the data needed to do vectorized operatinos on the mesh
# linear elements only
# disassembles SL0 into eqn.q
# calculates all mesh wide quantities in eqn

# need: u (previous timestep solution), x (coordinates), dxidx, jac, res, array of interfaces

#println("Entered dataPrep()")

#  println("eqn.q = ", eqn.q)

  u = eqn.q
  F_xi = eqn.F_xi

  F_eta = eqn.F_eta
  fill!(eqn.res, 0.0)

  # disassemble SL0 into eqn.
  for i=1:mesh.numEl  # loop over elements
    for j = 1:sbp.numnodes
      for k=1:4
	dofnum_k = mesh.dofs[k, j, i]
	u[k, j, i] = SL0[dofnum_k]
      end
    end
  end

  # disassmble SL0 into eqn.q
#  disassembleSolution(mesh, mesh.dofs, eqn.q, SL0)



  # calculate fluxes
  getEulerFlux(eqn, eqn.q, mesh.dxidx, F_xi, F_eta)
  isentropicVortexBC(mesh, sbp, eqn)
  stabscale(mesh, sbp, eqn)
#  println("getEulerFlux @time printed above")




#  println("finished dataPrep()")
  return nothing
end # end function dataPrep

function evalVolumeIntegrals{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator{Tsbp}, eqn::EulerEquation{Tsol})
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


function evalBoundaryIntegrals{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator{Tsbp}, eqn::EulerEquation{Tsol})
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
function addEdgeStabilize{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator{Tsbp}, eqn::EulerEquation{Tsol})

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
function applyDissipation(mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerEquation)
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





function applyMassMatrixInverse{Tsol}(eqn::EulerEquation{Tsol}, SL::AbstractVector{Tsol})
# apply the inverse mass matrix stored eqn to SL

#  SL .*= eqn.Minv  # this gives wrong answer


  ndof = length(SL)
  for i=1:ndof
    SL[i] *= eqn.Minv[i]
  end

  return nothing
end

# some helper functions


function getEulerJac_wrapper{T}(q::AbstractArray{T,1}, F::AbstractArray{T,1})
  dir = [1.0, 0.0]
#  F = zeros(T, 4)
  sbp = TriSBP{Float64}()
  mesh = PumiMesh2{Float64}(".null", "../../mesh_files/quarter_vortex3l.smb", 1, sbp; dofpernode=4)
  eqn = EulerEquation1{T, T}(mesh, sbp)

  @time getEulerFlux(eqn, q, dir, F)

  return F

end


function getEulerFlux{Tmsh, Tsol}(eqn::EulerEquation{Tsol}, q::AbstractArray{Tsol,1}, dir::AbstractArray{Tmsh,1},  F::AbstractArray{Tsol,1})
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




function getEulerFlux{Tmsh, Tsol}(eqn::EulerEquation{Tsol}, q::AbstractArray{Tsol,3}, dxidx::AbstractArray{Tmsh,4},  F_xi::AbstractArray{Tsol,3}, F_eta::AbstractArray{Tsol,3})
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

fluxJac = forwarddiff_jacobian!(getEulerJac_wrapper, Float64, fadtype=:dual; n=4, m=4)



function calcPressure{Tsol}(q::AbstractArray{Tsol,1}, eqn::EulerEquation{Tsol})
  # calculate pressure for a node
  # q is a vector of length 4 of the conservative variables

#  internal_energy = SL_vals[4]/SL_vals[1] - 0.5*(SL_vals[2]^2 + SL_vals[3]^2)/(SL_vals[1]^2)
#  pressure = SL_vals[1]*eqn.R*internal_energy/eqn.cv

  return  (eqn.gamma-1)*(q[4] - 0.5*(q[2]^2 + q[3]^2)/q[1])
  
#   println("internal_energy = ", internal_energy, " , pressure = ", pressure)


end


function disassembleSolution{Tmsh, Tsol}(mesh::AbstractMesh{Tmsh}, dofs::AbstractArray{Int, 3}, u::AbstractArray{Tsol,3}, SL0::AbstractArray{Tsol, 1})

  # disassemble SL0 into eqn.
  for i=1:mesh.numEl  # loop over elements
    for j = 1:mesh.numNodesPerElement
      for k=1:4
	dofnum_k = dofs[k, j, i]
	u[k, j, i] = SL0[dofnum_k]
      end
    end
  end

  return nothing
end


function assembleSolution{Tmsh, Tsol}(mesh::AbstractMesh{Tmsh}, eqn::EulerEquation{Tmsh}, SL::AbstractArray{Tsol,1})

  println("in assembleSolution")

#  println("size(mesh.dofs) = ", size(mesh.dofs))
#  println("size(eqn.res = ", size(eqn.res))

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



