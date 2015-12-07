# this file contains the functions to evaluate the right hand side of the weak form in the pdf
# res_vec stands for semi-linear form, contents inside the inv(M)(...) parenthesis on pg. 7 of Jared's derivation

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

# The proper order for arguments is (mesh, sbp, eqn, opts).


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

@doc """
This module is organized into 3 levels of functions: high, middle, and low.
The high level functions take the mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerEquation, and opts (options dictionary), They do not know the types their 
arguments are paramaterized on. There is only one method for each high level function.  All they do is call mid level functions.

Mid level functions take the same arguments as high level functions but know
the types they are paramaterized on and do the right thing for all paramters.
There is only one method for each mid level function.  Mid level functions
typically pass pieces of the large arrays stored in the eqn and mesh objects
to low level functions that actually do computation.

Low level functions usually take in vectors of whatever values they need, do
a computation, and store the result in an argument vector.  These function 
are usually passed an ArrayView, so writing to the argument vector writes 
directly to the large array stored in the eqn object.  Thus there is no 
loss of efficiency by using low level functions.
"""
export evalEuler, init

# Rules for paramaterization:
# Tmsh = mesh data type
# Tsbp = SBP operator data type
# Tsol = equation, res_vec, q_vec data type


@doc """
### EulerEquationMod.evalEuler

  This function drives the evaluation of the EulerEquations.
  It is agnostic to the dimension of the equation. and the types the arguments
  are paramaterized on.

  The function calls only high level functions, all of which take the same
  four arguments.  Mid level function also take the same arguments.

  The input/output variables are eqn.q and eqn.res, respectively.
  eqn.q_vec and eqn.res_vec exist as reusable storage

  The function disassembleSolution takes q_vec and puts it into eqn.q
  The function assembleSolution takes eqn.res and puts it into res_vec

  Arguments:
    * mesh  : a mesh object
    * sbp   : SBP operator object
    * eqn   : an EulerData object
    * opts  : options dictionary

  The optional time argument is used for unsteady equations

"""->
# this function is what the timestepper calls
# high level function
function evalEuler(mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerData, opts, t=0.0)
  # res_vec is popualted with du/dt
  # q_vec is q at previous timestep
  # t is current timestep
  # extra_args is unpacked into object needed to evaluation equation
  
  dataPrep(mesh, sbp, eqn, opts)
  # println("\n eqn.bndryflux = \n", eqn.bndryflux)

  #println("dataPrep @time printed above")
  evalVolumeIntegrals(mesh, sbp, eqn, opts)
  #println("volume integral @time printed above")
  
  #----------------------------------------------------------------------------
  #=
  bndryfluxPhysical = zeros(eqn.bndryflux)
  getPhysBCFluxes(mesh, sbp, eqn, opts, bndryfluxPhysical)
  #println("bndryfluxPhysical = \n", bndryfluxPhysical)
  #println("eqn.bndryflux = \n", eqn.bndryflux)
  bndryfluxPhysical = -1*bndryfluxPhysical
  boundaryintegrate!(sbp, mesh.bndryfaces, bndryfluxPhysical, eqn.res)
  =#

  if opts["use_GLS"]
    GLS(mesh,sbp,eqn)
  end
  
  #=
  bndryfluxPhysical = -1*bndryfluxPhysical
  boundaryintegrate!(sbp, mesh.bndryfaces, bndryfluxPhysical, eqn.res)
  =#
  #----------------------------------------------------------------------------

  evalBoundaryIntegrals(mesh, sbp, eqn)
  # println("\n eqn.res = \n", eqn.res)
  #println("boundary integral @time printed above")


  addStabilization(mesh, sbp, eqn, opts)
  #println("edge stabilizing @time printed above")


  # no more assembling solution
  #assembleSolution(mesh, sbp, eqn, res_vec)
  #println("assembly @time printed above")

  #applyMassMatrixInverse(eqn, res_vec)
  #println("Minv @time printed above")

  #println("after Minv, sum(res) = ", sum(eqn.res))
  #println("    sum(res_vec) = ", sum(res_vec))
  #applyDissipation(mesh, sbp, eqn, res_vec, q_vec)
 
  #print current error
  #err_norm = norm(res_vec)/mesh.numDof
  #print(" ", err_norm)


  #print("\n")

  return nothing
end  # end evalEuler



@doc """
### EulerEquationMod.init

  This function performs any operations that need to be performed before
  the first residual evaluation.
  Any operations that need to be performed at the beginning of *each* 
  residual evaluation should go in dataPrep
"""
# high level functions
function init{Tmsh, Tsol}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, eqn::AbstractEulerData{Tsol}, opts)

  println("\nInitializing Euler module")
  # get BC functors
  getBCFunctors(mesh, sbp, eqn, opts)


  return nothing
end


@doc """
### EulerEquationMod.dataPrep

  This function calculates all the quantities that are stored over the entire 
  mesh, and performs a few sanity check such as ensuring density and pressure
  are positive.

  This is a high level function
"""
# high level function
function dataPrep{Tmsh,  Tsol}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, eqn::AbstractEulerData{Tsol}, opts)
# gather up all the data needed to do vectorized operatinos on the mesh
# disassembles q_vec into eqn.q
# calculates all mesh wide quantities in eqn

# need: u (previous timestep solution), x (coordinates), dxidx, jac, res, array of interfaces

#println("Entered dataPrep()")

#  println("eqn.q = ", eqn.q)

  # apply filtering to input
  if eqn.params.use_filter
    applyFilter(mesh, sbp, eqn, eqn.q, opts)
  end


  u = eqn.q
  flux_parametric = eqn.flux_parametric

  # zero out res
  fill!(eqn.res, 0.0)
  
  # disassemble q_vec into eqn.q
#  disassembleSolution(mesh, sbp, eqn, opts, q_vec)
  # disassmble q_vec into eqn.q

  getAuxVars(mesh, eqn)
#  println("getAuxVars @time printed above")


  checkDensity(eqn)
#  println("checkDensity @time printed above")

  checkPressure(eqn)
#  println("checkPressure @time printed above")

  # calculate fluxes
#  getEulerFlux(eqn, eqn.q, mesh.dxidx, view(flux_parametric, :, :, :, 1), view(flux_parametric, :, :, :, 2))
  getEulerFlux(mesh, sbp,  eqn, opts)
#  println("getEulerFlux @time printed above")



#  getIsentropicVortexBoundaryFlux(mesh, sbp, eqn)
   getBCFluxes(mesh, sbp, eqn, opts)
   # println("\n eqn.bndryflux = \n", eqn.bndryflux)
#   println("getBCFluxes @time printed above")
#  isentropicVortexBC(mesh, sbp, eqn)
  stabscale(mesh, sbp, eqn)
#  println("stabscale @time printed above")




#  println("finished dataPrep()")
  return nothing
end # end function dataPrep


@doc """
### EulerEquationMod.getBCFluxes

  This function calls other functions to calculate the boundary fluxes, passing
  them pieces of the array needed.  This populates eqn.bndryflux.  It also
  calles writeBoundary()

  This is a mid level function
"""->
# this is a mid level function
function getBCFluxes(mesh, sbp, eqn, opts)
  #get all the fluxes for all the boundary conditions and save them in eqn.bndryflux

  #println("mesh.bndry_funcs = ", mesh.bndry_funcs)

  for i=1:mesh.numBC
  #  println("computing flux for boundary condition ", i)
    functor_i = mesh.bndry_funcs[i]
    start_index = mesh.bndry_offsets[i]
    end_index = mesh.bndry_offsets[i+1]
    bndry_facenums_i = view(mesh.bndryfaces, start_index:(end_index - 1))
    bndryflux_i = view(eqn.bndryflux, :, :, start_index:(end_index - 1))
  
    calcBoundaryFlux(mesh, sbp, eqn, functor_i, bndry_facenums_i, bndryflux_i)
  end

  writeBoundary(mesh, sbp, eqn, opts)

  return nothing
end

@doc """
### EulerEquationMod.writeBoundary 

  This function writes information about the boundary faces and fluxes to files.
  It is controlled by the input argument writeboundary, of type Bool.

  It generates the files:
    * boundaryfaces.dat : writes mesh.bndryfaces, an array with eltype Boundary
                          to a file, one element per line
    * boundaryflux.dat  : writes the element, local node number and boundary 
                          flux to a line in a human readable format
    * boundaryflux2.dat : writes the real part ofmesh.bndryflux to space 
                          delimited file

   This is a high level function.
"""->
function writeBoundary(mesh, sbp, eqn, opts)

  if !eqn.params.writeboundary
    return nothing
  end

    face_name = "boundaryfaces.dat"
    flux_name = "boundaryflux.dat"
    flux_dlm = "boundaryflux2.dat"

    rmfile(face_name)
    rmfile(flux_name)
    rmfile(flux_dlm)


  # write boundaryfaces.dat
  f = open(face_name, "a+")
  for i=1:length(mesh.bndryfaces)
    println(f, mesh.bndryfaces[i])
  end
  close(f)

  # write boundaryflux.dat
  f = open(flux_name, "a+")
  for i=1:mesh.numBoundaryEdges
    el = mesh.bndryfaces[i].element
    face = mesh.bndryfaces[i].face
    for j=1:sbp.numfacenodes
      jb = sbp.facenodes[j, face]
      println(f, "el ", el, ", node_index ", jb, ", flux = ", real(eqn.bndryflux[:, j, i]))
    end
  end
  close(f)

  # write boundaryflux2.dat
  writedlm(flux_dlm, real(eqn.bndryflux))
  
  return nothing
end



@doc """
### EulerEquationMod.checkDensity

  This function checks that the density is positive. If not, an error is thrown.
  Because density is stored in the eqn.q array, performing this check takes 
  very little time.

  Arguments:
    * EulerData

  This is a low level function.
"""->
# low level function
function checkDensity(eqn::EulerData)
# check that density is positive

(ndof, nnodes, numel) = size(eqn.q)

for i=1:numel
  for j=1:nnodes
    @assert( real(eqn.q[1, j, i]) > 0.0, "element $i, node $j. Density < 0")
  end
end

return nothing

end


@doc """
### EulerEquationMod.checkPressure

  This function checks that pressure is positive.  If not, an error is thrown.
  Because pressure is stored in the auxiliary variables array, this check 
  takes very little time

  Arguments:
    * EulerData

  This is a low level function
"""->
function checkPressure(eqn::EulerData)
# check that density is positive

(ndof, nnodes, numel) = size(eqn.q)

for i=1:numel
  for j=1:nnodes
#    q_vals = view(eqn.q, :, j, i)
#    press = calcPressure(q_vals, eqn)
    aux_vars = view(eqn.aux_vars,:, j, i)
    press = @getPressure(aux_vars)
#    press = getPressure(aux_vars)
#    println("press = ", press)
    @assert( real(press) > 0.0, "element $i, node $j")
  end
end

return nothing

end

@doc """
### EulerEquationMod.evalVolumeIntegrals

  This function evaluates the volume integrals of the Euler equations by
  calling the appropriate SBP function on the Euler flux stores in eqn.
  This function knows the dimension of the equation but does the right
  thing for all dimensions.  eqn.res is updated with the result

  This is a mid level function.
"""
# mid level function
function evalVolumeIntegrals{Tmsh,  Tsol, Tdim}(mesh::AbstractMesh{Tmsh}, 
                             sbp::SBPOperator, eqn::EulerData{Tsol, Tdim}, opts)
# evaluate all the integrals over the elements (but not the boundary)
# does not do boundary integrals
# mesh : a mesh type, used to access information about the mesh
# sbp : an SBP operator, used to get stiffness matricies and stuff
# eqn : a type that holds some constants needed to evaluate the equation
#	: also used for multiple dispatch
# res_vec : solution vector to be populated (mesh.numDof entries)
# q_vec : solution vector at previous timesteps (mesh.numDof entries)

#  println("size eqn.flux_parametric = ", size(eqn.flux_parametric), " size eqn.res = ", size(eqn.res), " sbp.numnodes = ", sbp.numnodes)
  
  if opts["Q_transpose"] == true
    for i=1:Tdim
      weakdifferentiate!(sbp, i, view(eqn.flux_parametric, :, :, :, i), eqn.res, trans=true)
    end
  else
    for i=1:Tdim
      weakdifferentiate!(sbp, i, -1*view(eqn.flux_parametric, :, :, :, i), eqn.res, trans=false)
    end
  end  # end if

  # artificialViscosity(mesh, sbp, eqn) 

  # hAverage = AvgMeshSize(mesh, eqn)
  # println("Average Mesh Size = ", hAverage)
  # need source term here



end
#------------- end of evalVolumeIntegrals


#=
@doc """
  This function evaluates the advective terms of the strong form.
  eqn.res is updates with the result

"""->
function evalAdvectiveStrong{Tmsh, Tsol, Tdim}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, eqn::EulerData{Tsol, Tdim}, opts)

  for i=1:Tdim
    differentiate!(sbp, i, 

=#
@doc """
  This function evaluates the boundary integrals in the Euler equations by 
  calling the appropriate SBP function on eqn.bndryflux, which must be populated
  before calling this function.  eqn.res is updated with the result

  This is a mid level function

"""->
# mid level function
function evalBoundaryIntegrals{Tmsh,  Tsol, Tdim}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, eqn::EulerData{Tsol, Tdim})
# evaluate all the integrals over the boundary
# mesh : a mesh type, used to access information about the mesh
# sbp : an SBP operator, used to get stiffness matricies and stuff
# eqn : a type that holds some constants needed to evaluate the equation
#	: also used for multiple dispatch


#  println("====== start of evalBoundaryIntegrals ======")
#boundaryintegrate!(sbp, bndryfaces, u, x, dxidx, isentropicVortexBC, result)

#boundaryintegrate2!(sbp, mesh.bndryfaces, eqn.q, mesh.coords, mesh.dxidx, isentropicVortexBC, eqn.res, mesh, eqn)

#println("eqn.bndryflux = ")
#println(eqn.bndryflux)


#idx = mesh.bndry_offsets[2]
#el = mesh.bndryfaces[idx].element
#println("bndry element = ", el)
#println("boundary flux = ", eqn.bndryflux[:, :, el])


boundaryintegrate!(sbp, mesh.bndryfaces, eqn.bndryflux, eqn.res)

#boundaryintegrate!(sbp, bndryfaces, u, x, dxidx, rho1Energy2BC, result)


#  println("==== end of evalBoundaryIntegrals ====")


  return nothing

end
#------------- end of evalBoundaryIntegrals


@doc """
### EulerEquationMod.addStabilization

  This function add whatever form of stabilization opts specifies by calling
  the appropriate function.  Some arrays may need to be populated before 
  calling this function, depending on the type of stabilization used.

  This is a mid level function
"""->
# This function adds edge stabilization to a residual using Prof. Hicken's edgestabilize! in SBP
# mid level function
function addStabilization{Tmsh,  Tsol}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, eqn::EulerData{Tsol}, opts)

#  println("==== start of addStabilization ====")
  # alpha calculated like in edgestabilize! documentation
  # stabscale (U+a)*gamma*h^2 where U=u*n, where u is the velocity 
  #   (remember to scale by rho) and n is the unit normal vector, from nrm->dxidx, then scaled by length
  # ifaces needs to be calculated
  # x needs to be passed

  


  # u argument here is res_vec in a different format
#  edgestabilize!(sbp, mesh.interfaces, eqn.q, mesh.coords, mesh.dxidx, mesh.jac, eqn.edgestab_alpha, stabscale, eqn.res, mesh, eqn)

  if eqn.params.use_edgestab
#    println("applying edge stabilization")
    edgestabilize!(sbp, mesh.interfaces, eqn.q, mesh.coords, mesh.dxidx, mesh.jac, eqn.edgestab_alpha, eqn.stabscale, eqn.res)
  end

  if eqn.params.use_res_filter
#    println("applying residual filter")
    applyFilter(mesh, sbp, eqn, eqn.res, opts, trans=true)
  end

  if eqn.params.use_dissipation
#    println("applying artificial dissipation")
    applyDissipation(mesh, sbp, eqn, eqn.q, opts)
  end

#  println("==== end of addStabilization ====")



  return nothing

end


# some helper functions
function getEulerJac_wrapper{T}(q::AbstractArray{T,1}, F::AbstractArray{T,1})
  dir = [1.0, 0.0]
#  F = zeros(T, 4)
  sbp = TriSBP{Float64}()
  mesh = PumiMesh2{Float64}(".null", "../../mesh_files/quarter_vortex3l.smb", 1, sbp; dofpernode=4)
  eqn = EulerData1{T, T}(mesh, sbp)

  @time getEulerFlux(eqn, q, dir, F)

  return F

end

@doc """
### EulerEquationMod.getAuxVars

  This function calculates any extra variables that are stored across the mesh
  using the conservative variables eqn.q.

  Thi is a mid level function
"""->
# mid level function
function getAuxVars{Tmsh, Tsol, Tdim}(mesh::AbstractMesh{Tmsh}, eqn::EulerData{Tsol, Tdim})
# calculate all auxiliary variables

  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      q_vals = view(eqn.q, :, j, i)

      # calculate pressure
      press = calcPressure(q_vals, eqn.params)
      @setPressure(eqn.aux_vars, j, i, press)
    end
  end

  return nothing
end

@doc """ 
### EulerEquationMod.getEulerFlux

  This function calculates the Euler flux across the entire mesh by passing
  pieces of the eqn.q, eqn.aux_vars, eqn.f_xi and eqn.params to a low level
  function.  The flux is calculated in the xi and eta directions, scaled (mulitiplied)
  by the Jacobian (so that when performing the integral we don't have to explictly
  divide by the jacobian, it just cancels out with the jacobian factor introduced
  here.

  This is a mid level function
"""->
# mid level function
function getEulerFlux{Tmsh, Tsol, Tdim}(mesh::AbstractMesh{Tmsh}, 
                                        sbp::SBPOperator,  
                                        eqn::EulerData{Tsol, Tdim}, opts)
# calculate Euler flux in parametric coordinate directions, stores it in eqn.flux_parametric

  nrm = zeros(Tmsh, 2)
  for i=1:mesh.numEl  # loop over elements
    for j=1:mesh.numNodesPerElement  # loop over nodes on current element
      q_vals = view(eqn.q, :, j, i)
      aux_vars = view(eqn.aux_vars, :, j, i)
      # put this loop here (rather than outside) to we don't have to fetch
      # q_vals twice, even though writing to the flux vector is slower
      # it might be worth copying the normal vector rather than
      # doing an view
      for k=1:Tdim  # loop over dimensions  
	# this will dispatch to the proper calcEulerFlux
	      nrm[1] = mesh.dxidx[k, 1, j, i]
	      nrm[2] = mesh.dxidx[k, 2, j, i]
#	nrm_mag = sqrt(nrm[1]*nrm[1] + nrm[2]*nrm[2])
#	nrm[1] /= nrm_mag
#	nrm[2] /= nrm_mag
#        nrm = view(mesh.dxidx, k, :, j, i) # this causes a type stability problem
        flux = view(eqn.flux_parametric, :, j, i, k)
        calcEulerFlux(eqn.params, q_vals, aux_vars, nrm, flux)
      end
    end
  end


  writeFlux(mesh, sbp, eqn, opts)

  return nothing
end

@doc """
### EulerEquationMod.writeFlux

  This function writes the real part of Euler flux to a file named Fxi.dat, space delimited, controlled by the input options writeflux, of type Bool.

  This is a high level function.
"""->
function writeFlux(mesh, sbp, eqn, opts)

   if !eqn.params.writeflux
     return nothing
   end
 
   fname = "Fxi.dat"
   rmfile(fname)
   writedlm(fname, real(eqn.flux_parametric))

   return nothing
end


@doc """
### EulerEquationMod.getEulerFlux2

  This function calcules the euler flux over the entire mesh directly (ie.
  does not call a low level function.  This function is deprecated, although
  useful for benchmarking purposes.  2D only.

  This is a mid level function
"""->
# this function is deprecated in factor of getEulerFlux()
# useful for benchmarking purposes
function getEulerFlux2{Tmsh, Tsol}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator,
                                   eqn::EulerData{Tsol}, opts)
# calculates the Euler flux for every node in the xi and eta directions
# eqn is the equation type
# q is the 3D array (4 by nnodes per element by nel), of the conservative variables
# dxidx is the 4D array (2 by 2 x nnodes per element by nel) that specifies the direction of xi and eta in each element (output from mappingjacobian!)
# flux_parametric is populated with the flux in xi direction (same shape as q)
# F_eta is populated with flux in eta direction

# once the Julia developers fix slice notation and speed up subarrays, we won't have to 
# vectorize like this (can calculate flux one node at a time inside a dedicated function

q = eqn.q
dxidx = mesh.dxidx
flux_parametric = view(eqn.flux_parametric, :, :, :, 1)
F_eta = view(eqn.flux_parametric, :, :, :, 2)

(ncomp, nnodes, nel) = size(q)  # get sizes of things

  for i=1:nel  # loop over elements
    for j=1:nnodes  # loop over nodes within element
      # get direction vector components (xi direction)
      nx = dxidx[1, 1, j, i]
      ny = dxidx[1, 2, j, i]
      # calculate pressure 
      press = (eqn.params.gamma-1)*(q[4, j, i] - 0.5*(q[2, j, i]^2 + q[3, j, i]^2)/q[1, j, i])

      # calculate flux in xi direction
      # hopefully elements of q get stored in a register for reuse in eta direction
      U = (q[2, j, i]*nx + q[3, j, i]*ny)/q[1, j, i]
      flux_parametric[1, j, i] = q[1, j, i]*U
      flux_parametric[2, j, i] = q[2, j, i]*U + nx*press
      flux_parametric[3, j, i] = q[3, j, i]*U + ny*press
      flux_parametric[4, j, i] = (q[4, j, i] + press)*U

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


@doc """
### EulerEquationMod.applyMassMatrixInverse

  This function multiplies eqn.res_vec (the residual in vector form  by eqn.Minv,
  the diagonal mass matrix.  This is a very fast function because all values
  are precomputed and stores linearly in memory.

  This is a mid level function, and does the correct thing regardless of the
  dimension of the equation.
"""->
# mid level function (although it doesn't really need to Tdim)
function applyMassMatrixInverse{Tsol, Tdim}(eqn::EulerData{Tsol, Tdim}, 
                                            res_vec::AbstractVector{Tsol})
  # apply the inverse mass matrix stored eqn to res_vec

  ndof = length(res_vec)
  for i=1:ndof
    res_vec[i] *= eqn.Minv[i]
  end

  return nothing
end




#=
@doc """
### EulerEquationMod.disassembleSolution

  This takes eqn.q_vec (the initial state), and disassembles it into eqn.q, the
  3 dimensional array of conservative variables.  This function uses mesh.dofs
  to speed the process.

  This is a mid level function, and does the right thing regardless of equation
  dimension.
"""->
=#
# mid level function (although it doesn't need Tdim)
function disassembleSolution{Tmsh, Tsol, Tdim}(mesh::AbstractMesh{Tmsh},
                             sbp, eqn::EulerData{Tsol, Tdim}, opts, 
                             q_vec::AbstractArray{Tsol, 1})
  # disassemble q_vec into eqn.
  for i=1:mesh.numEl  # loop over elements
    for j = 1:mesh.numNodesPerElement
      for k=1:(Tdim+2)
	      dofnum_k = mesh.dofs[k, j, i]
	      eqn.q[k, j, i] = q_vec[dofnum_k]
      end
    end
  end

  writeQ(mesh, sbp, eqn, opts)

  return nothing
end

@doc """
### EulerEquationMod.writeQ

  This function writes the real part of the solution variables eqn.q to a space 
  delimited file called q.dat, controlled by the input options writeq, of type bool

  This is a high level function.
"""->
function writeQ(mesh, sbp, eqn, opts)

  if !eqn.params.writeq
    return nothing
  end

  fname = "q.dat"
  rmfile(fname)
  writedlm(fname, real(eqn.q))

  return nothing
end



@doc """
### EulerEquationMod.assembleSolution

  This function takes the 3D array of variables in arr and 
  reassmbles is into the vector res_vec.  Note that
  This is a reduction operation and requires eqn.res_vec to be zerod before
  calling this function.

  This is a mid level function, and does the right thing regardless of
  equation dimension
"""->
# mid level function (although it doesn't need Tdim)
function assembleSolution{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, eqn::EulerData{Tsol}, opts, arr, res_vec::AbstractArray{Tres,1})
# arr is the array to be assembled into res_vec

#  println("in assembleSolution")

#  println("size(mesh.dofs) = ", size(mesh.dofs))
#  println("size(eqn.res = ", size(eqn.res))


#  println("before assembly, sum(res) = ", sum(eqn.res))
#  println("    sum(res_vec) = ", sum(res_vec))

  for i=1:mesh.numEl  # loop over elements
    for j=1:mesh.numNodesPerElement
      for k=1:4  # loop over dofs on the node
	      dofnum_k = mesh.dofs[k, j, i]
	      res_vec[dofnum_k] += arr[k,j,i]
      end
    end
  end
  

#  println("after assembly, sum(res) = ", sum(eqn.res))
#  println("    sum(res_vec) = ", sum(res_vec))

  return nothing
end

@doc """
### EulerEquationMod.calcEulerFlux

   This function calculates the Euler flux from the conservative variables at
   a single node in a particular direction.  2D only.

   Inputs:
   params  : ParamaterType{2}
   q  : vector of conservative variables
   aux_vars : vector of auxiliary variables
   dir :  unit vector in direction to calculate the flux

   Inputs/Outputs:
   F  : vector to populate with the flux

   The Tdim paramater of params determine whether this method or the 3D 
   version is called.

   This is a low level function
"""->
# low level function
function calcEulerFlux{Tmsh, Tsol, Tres}(params::ParamType{2}, 
                      q::AbstractArray{Tsol,1}, 
                      aux_vars::AbstractArray{Tres, 1}, 
                      dir::AbstractArray{Tmsh},  F::AbstractArray{Tsol,1})
# calculates the Euler flux in a particular direction at a point
# eqn is the equation type
# q is the vector (of length 4), of the conservative variables at the point
# aux_vars is the vector of auxiliary variables at the point
# dir is a vector of length 2 that specifies the direction
# F is populated with the flux (is a vector of length 4)
# 2D  only


#  press = calcPressure(q, params)
  press = @getPressure(aux_vars)
  U = (q[2]*dir[1] + q[3]*dir[2])/q[1]
  F[1] = q[1]*U
  F[2] = q[2]*U + dir[1]*press
  F[3] = q[3]*U + dir[2]*press
  F[4] = (q[4] + press)*U
 
  return nothing

end

@doc """
### EulerEquationMod.calcEulerFlux
  This is the 3D method.  All arguments are same as the 2D version.
"""->
# low level function
function calcEulerFlux{Tmsh, Tsol}(params::ParamType{3}, q::AbstractArray{Tsol,1}, dir::AbstractArray{Tmsh},  F::AbstractArray{Tsol,1})
# calculates the Euler flux in a particular direction at a point
# eqn is the equation type
# q is the vector (of length 5), of the conservative variables at the point
# dir is a vector of length 3 that specifies the direction
# F is populated with the flux (is a vector of length 5)
# 3D  only

# once the Julia developers fix slice notation and speed up subarrays, we can make a faster
# vectorized version of this

  press = calcPressure(q, params)
  U = (q[2]*dir[1] + q[3]*dir[2] + q[4]*dir[3])/q[1]
  F[1] = q[1]*U
  F[2] = q[2]*U + dir[1]*press
  F[3] = q[3]*U + dir[2]*press
  F[4] = q[4]*U + dir[3]*press
  F[5] = (q[5] + press)*U
 
  return nothing

end



@doc """
### EulerEquationMod.calcPressure

  This function calculates the pressure from the conservative variables at a
  node in 2D.  It returns a single value.

  Inputs:
    q  : vector of conservative variables
    params : ParamType{2}

  The paramater of params determines whether the 2D or 3D method is dispatched.

  This is a low level function.
"""->
# low level function
function calcPressure{Tsol}(q::AbstractArray{Tsol,1}, params::ParamType{2})
  # calculate pressure for a node
  # q is a vector of length 4 of the conservative variables

#  internal_energy = res_vec_vals[4]/res_vec_vals[1] - 0.5*(res_vec_vals[2]^2 + res_vec_vals[3]^2)/(res_vec_vals[1]^2)
#  pressure = res_vec_vals[1]*eqn.R*internal_energy/eqn.cv

  return  (params.gamma_1)*(q[4] - 0.5*(q[2]*q[2] + q[3]*q[3])/q[1])
  
#   println("internal_energy = ", internal_energy, " , pressure = ", pressure)


end

@doc """
### EulerEquationMod.calcPressure

  3D method.  See 2D method documentation
"""->
# low level function
function calcPressure{Tsol}(q::AbstractArray{Tsol,1}, params::ParamType{3})
  # calculate pressure for a node
  # q is a vector of length 5 of the conservative variables

#  internal_energy = res_vec_vals[4]/res_vec_vals[1] - 0.5*(res_vec_vals[2]^2 + res_vec_vals[3]^2)/(res_vec_vals[1]^2)
#  pressure = res_vec_vals[1]*eqn.R*internal_energy/eqn.cv

  return  (params.gamma_1)*(q[5] - 0.5*(q[2]*q[2] + q[3]*q[3] + q[4]*q[4])/q[1])
  
#   println("internal_energy = ", internal_energy, " , pressure = ", pressure)


end


