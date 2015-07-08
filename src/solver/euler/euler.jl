# this file contains the functions to evaluate the right hand side of the weak form in the pdf
# SL stands for semi-linear form, contents inside the inv(M)(...) parenthesis on pg. 7 of Jared's derivation

# Organization:
# The code is organized into 3 levels of functions.  High level functions
# only call mid level functions, passing the Mesh, Equation, and SBP objects
# to them.  At this point, there shouldn't be any need for more high level
# functions.
#
# The mid level functions either call SBP functions or have loops over
# the arrays stored in the Equation, calling low level functions
# to calculate quantities at each node.  Use the ArrayViews package to pass
# portions of a larger array to the low level functions.
#
# Low level functions calculate quantities for a single node.  Although
# these functions will likely take in an ArrayView, their argument
# types should be AbstractArrays{T, N}, specifying T according to 
# what type the original array is (from either the Equation or Mesh object)

# The proper order for arguments is (mesh, sbp, eqn).


# Rules: 
# 1. function that take a composite type with abstract fields *cannot* use those
# fields directly (they can pass those fields to other functions however
# This leads directly to a two level structure for the code: high level function
# that take in composite types and low level function that take in arrays and
# perform calculations on them

# at this time, no composite types have abstract fields (all AbstractArrays
# were turned into arrays), so this is a non-issue

# the reason for this is that the compiler does not compile new version of the 
# function based
# on the types of the fields of a composite type. Passing the fields of the typ
# e to other functions fixes this problem because the fields are now arguments,
# so the compiler can specialize the code

# 2.  Arrays should not be returned from functions.  The caller should allocate
# and array and pass it into the function

# this allows reusing the same array during a loop (rather than 
# allocating a new array)



export evalEuler

# Rules for paramaterization:
# Tmsh = mesh data type
# Tsbp = SBP operator data type
# Tsol = equation, SL, SL0 data type


# this function is what the timestepper calls
# high level function
function evalEuler(mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerEquation,  SL0, SL, t=0.0)
# SL is popualted with du/dt
# SL0 is q at previous timestep
# t is current timestep
# extra_args is unpacked into object needed to evaluation equation

dataPrep(mesh, sbp, eqn, SL0)
#println("dataPrep @time printed above")
evalVolumeIntegrals(mesh, sbp, eqn)
#println("volume integral @time printed above")

evalBoundaryIntegrals(mesh, sbp, eqn)
#println("boundary integral @time printed above")



addStabilization(mesh, sbp, eqn)
#println("edge stabilizing @time printed above")


assembleSolution(mesh, eqn, SL)
#println("assembly @time printed above")

applyMassMatrixInverse(eqn, SL)
#println("Minv @time printed above")

#applyDissipation(mesh, sbp, eqn, SL, SL0)



#print current error
#err_norm = norm(SL)/mesh.numDof
#print(" ", err_norm)


#print("\n")

return nothing
#return SL

end  # end evalEuler


# high level functions
function dataPrep{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator{Tsbp}, eqn::AbstractEulerEquation{Tsol}, SL0::AbstractVector{Tsol})
# gather up all the data needed to do vectorized operatinos on the mesh
# disassembles SL0 into eqn.q
# calculates all mesh wide quantities in eqn

# need: u (previous timestep solution), x (coordinates), dxidx, jac, res, array of interfaces

#println("Entered dataPrep()")

#  println("eqn.q = ", eqn.q)

  u = eqn.q
  F_xi = eqn.F_xi

  # zero out res
  fill!(eqn.res, 0.0)
  
  # disassemble SL0 into eqn.q
  disassembleSolution(mesh, eqn, SL0)
  # disassmble SL0 into eqn.q


  checkDensity(eqn)
#  println("checkDensity @time printed above")

  checkPressure(eqn)
#  println("checkPressure @time printed above")

  # calculate fluxes
#  getEulerFlux(eqn, eqn.q, mesh.dxidx, view(F_xi, :, :, :, 1), view(F_xi, :, :, :, 2))
  getEulerFlux(mesh, eqn)
  getIsentropicVortexBoundaryFlux(mesh, sbp, eqn)
#  isentropicVortexBC(mesh, sbp, eqn)
  stabscale(mesh, sbp, eqn)
#  println("getEulerFlux @time printed above")




#  println("finished dataPrep()")
  return nothing
end # end function dataPrep


function checkDensity(eqn::EulerEquation)
# check that density is positive

(ndof, nnodes, numel) = size(eqn.q)

for i=1:numel
  for j=1:nnodes
    @assert( real(eqn.q[1, j, i]) > 0.0)
  end
end

return nothing

end

function checkPressure(eqn::EulerEquation)
# check that density is positive

(ndof, nnodes, numel) = size(eqn.q)

for i=1:numel
  for j=1:nnodes
    q_vals = view(eqn.q, :, j, i)
    press = calcPressure(q_vals, eqn)
    @assert( press > 0.0)
  end
end

return nothing

end


# mid level function
function evalVolumeIntegrals{Tmsh, Tsbp, Tsol, Tdim}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator{Tsbp}, eqn::EulerEquation{Tsol, Tdim})
# evaluate all the integrals over the elements (but not the boundary)
# does not do boundary integrals
# mesh : a mesh type, used to access information about the mesh
# sbp : an SBP operator, used to get stiffness matricies and stuff
# eqn : a type that holds some constants needed to evaluate the equation
#	: also used for multiple dispatch
# SL : solution vector to be populated (mesh.numDof entries)
# SL0 : solution vector at previous timesteps (mesh.numDof entries)

#  println("size eqn.F_xi = ", size(eqn.F_xi), " size eqn.res = ", size(eqn.res), " sbp.numnodes = ", sbp.numnodes)

  for i=1:Tdim
    weakdifferentiate!(sbp, i, view(eqn.F_xi, :, :, :, i), eqn.res, trans=true)
  end


  # need source term here



end
#------------- end of evalVolumeIntegrals


function evalBoundaryIntegrals{Tmsh, Tsbp, Tsol, Tdim}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator{Tsbp}, eqn::EulerEquation{Tsol, Tdim})
# evaluate all the integrals over the boundary
# mesh : a mesh type, used to access information about the mesh
# sbp : an SBP operator, used to get stiffness matricies and stuff
# eqn : a type that holds some constants needed to evaluate the equation
#	: also used for multiple dispatch


#  println("====== start of evalBoundaryIntegrals ======")
#boundaryintegrate!(sbp, bndryfaces, u, x, dxidx, isentropicVortexBC, result)

#boundaryintegrate2!(sbp, mesh.bndryfaces, eqn.q, mesh.coords, mesh.dxidx, isentropicVortexBC, eqn.res, mesh, eqn)

boundaryintegrate!(sbp, mesh.bndryfaces, eqn.bndryflux, eqn.res)
#boundaryintegrate!(sbp, bndryfaces, u, x, dxidx, rho1Energy2BC, result)


#  println("==== end of evalBoundaryIntegrals ====")


  return nothing

end
#------------- end of evalBoundaryIntegrals

# This function adds edge stabilization to a residual using Prof. Hicken's edgestabilize! in SBP
# mid level function
function addStabilization{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator{Tsbp}, eqn::EulerEquation{Tsol})

#  println("==== start of addStabilization ====")
  # alpha calculated like in edgestabilize! documentation
  # stabscale (U+a)*gamma*h^2 where U=u*n, where u is the velocity 
  #   (remember to scale by rho) and n is the unit normal vector, from nrm->dxidx, then scaled by length
  # ifaces needs to be calculated
  # x needs to be passed

  


  # u argument here is SL in a different format
#  edgestabilize!(sbp, mesh.interfaces, eqn.q, mesh.coords, mesh.dxidx, mesh.jac, eqn.edgestab_alpha, stabscale, eqn.res, mesh, eqn)

  edgestabilize!(sbp, mesh.interfaces, eqn.q, mesh.coords, mesh.dxidx, mesh.jac, eqn.edgestab_alpha, eqn.stabscale, eqn.res)

#  println("==== end of addStabilization ====")



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



# mid level function
function getEulerFlux{Tmsh, Tsol, Tdim}(mesh::AbstractMesh{Tmsh}, eqn::EulerEquation{Tsol, Tdim})
# calculate Euler flux in parametric coordinate directions, stores it in eqn.F_xi

  nrm = zeros(Tmsh, 2)
  for i=1:mesh.numEl  # loop over elements
    for j=1:mesh.numNodesPerElement  # loop over nodes on current element
      q_vals = unsafe_view(eqn.q, :, j, i)
      # put this loop here (rather than outside) to we don't have to fetch
      # q_vals twice, even though writing to the flux vector is slower
      # it might be worth copying the normal vector rather than
      # doing an unsafe_view
      for k=1:Tdim  # loop over dimensions  
	# this will dispatch to the proper calcEulerFlux
	nrm[1] = mesh.dxidx[k, 1, j, i]
	nrm[2] = mesh.dxidx[k, 2, j, i]
#        nrm = unsafe_view(mesh.dxidx, k, :, j, i) # this causes a type stability problem
        calcEulerFlux(eqn, q_vals, nrm, unsafe_view(eqn.F_xi, :, j, i, k))
      end
    end
  end

  return nothing
end

# this function is deprecated in factor of getEulerFlux()
# useful for benchmarking purposes
function getEulerFlux2{Tmsh, Tsol}( mesh::AbstractMesh{Tmsh}, eqn::EulerEquation{Tsol})
# calculates the Euler flux for every node in the xi and eta directions
# eqn is the equation type
# q is the 3D array (4 by nnodes per element by nel), of the conservative variables
# dxidx is the 4D array (2 by 2 x nnodes per element by nel) that specifies the direction of xi and eta in each element (output from mappingjacobian!)
# F_xi is populated with the flux in xi direction (same shape as q)
# F_eta is populated with flux in eta direction

# once the Julia developers fix slice notation and speed up subarrays, we won't have to 
# vectorize like this (can calculate flux one node at a time inside a dedicated function

q = eqn.q
dxidx = mesh.dxidx
F_xi = view(eqn.F_xi, :, :, :, 1)
F_eta = view(eqn.F_xi, :, :, :, 2)

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



# mid level function (although it doesn't really need to Tdim)
function applyMassMatrixInverse{Tsol, Tdim}(eqn::EulerEquation{Tsol, Tdim}, SL::AbstractVector{Tsol})
# apply the inverse mass matrix stored eqn to SL

#  SL .*= eqn.Minv  # this gives wrong answer


  ndof = length(SL)
  for i=1:ndof
    SL[i] *= eqn.Minv[i]
  end

  return nothing
end





# mid level function (although it doesn't need Tdim)
function disassembleSolution{Tmsh, Tsol, Tdim}(mesh::AbstractMesh{Tmsh}, eqn::EulerEquation{Tsol, Tdim}, SL0::AbstractArray{Tsol, 1})
  # disassemble SL0 into eqn.
  for i=1:mesh.numEl  # loop over elements
    for j = 1:mesh.numNodesPerElement
      for k=1:(Tdim+2)
	dofnum_k = mesh.dofs[k, j, i]
	eqn.q[k, j, i] = SL0[dofnum_k]
      end
    end
  end

  return nothing
end

# mid level function (although it doesn't need Tdim)
function assembleSolution{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh}, eqn::EulerEquation{Tsol}, SL::AbstractArray{Tres,1})

#  println("in assembleSolution")

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


# low level function
function calcEulerFlux{Tmsh, Tsol}(eqn::EulerEquation{Tsol, 2}, q::AbstractArray{Tsol,1}, dir::AbstractArray{Tmsh},  F::AbstractArray{Tsol,1})
# calculates the Euler flux in a particular direction at a point
# eqn is the equation type
# q is the vector (of length 4), of the conservative variables at the point
# dir is a vector of length 2 that specifies the direction
# F is populated with the flux (is a vector of length 4)
# 2D  only


  press = calcPressure(q, eqn)
  U = (q[2]*dir[1] + q[3]*dir[2])/q[1]
  F[1] = q[1]*U
  F[2] = q[2]*U + dir[1]*press
  F[3] = q[3]*U + dir[2]*press
  F[4] = (q[4] + press)*U
 
  return nothing

end

# low level function
function calcEulerFlux{Tmsh, Tsol}(eqn::EulerEquation{Tsol, 3}, q::AbstractArray{Tsol,1}, dir::AbstractArray{Tmsh},  F::AbstractArray{Tsol,1})
# calculates the Euler flux in a particular direction at a point
# eqn is the equation type
# q is the vector (of length 5), of the conservative variables at the point
# dir is a vector of length 3 that specifies the direction
# F is populated with the flux (is a vector of length 5)
# 3D  only

# once the Julia developers fix slice notation and speed up subarrays, we can make a faster
# vectorized version of this

  press = calcPressure(q, eqn)
  U = (q[2]*dir[1] + q[3]*dir[2] + q[4]*dir[3])/q[1]
  F[1] = q[1]*U
  F[2] = q[2]*U + dir[1]*press
  F[3] = q[3]*U + dir[2]*press
  F[4] = q[4]*U + dir[3]*press
  F[5] = (q[5] + press)*U
 
  return nothing

end




# low level function
function calcPressure{Tsol}(q::AbstractArray{Tsol,1}, eqn::EulerEquation{Tsol, 2})
  # calculate pressure for a node
  # q is a vector of length 4 of the conservative variables

#  internal_energy = SL_vals[4]/SL_vals[1] - 0.5*(SL_vals[2]^2 + SL_vals[3]^2)/(SL_vals[1]^2)
#  pressure = SL_vals[1]*eqn.R*internal_energy/eqn.cv

  return  (eqn.gamma_1)*(q[4] - 0.5*(q[2]*q[2] + q[3]*q[3])/q[1])
  
#   println("internal_energy = ", internal_energy, " , pressure = ", pressure)


end

# low level function
function calcPressure{Tsol}(q::AbstractArray{Tsol,1}, eqn::EulerEquation{Tsol, 3})
  # calculate pressure for a node
  # q is a vector of length 5 of the conservative variables

#  internal_energy = SL_vals[4]/SL_vals[1] - 0.5*(SL_vals[2]^2 + SL_vals[3]^2)/(SL_vals[1]^2)
#  pressure = SL_vals[1]*eqn.R*internal_energy/eqn.cv

  return  (eqn.gamma_1)*(q[5] - 0.5*(q[2]*q[2] + q[3]*q[3] + q[4]*q[4])/q[1])
  
#   println("internal_energy = ", internal_energy, " , pressure = ", pressure)


end


