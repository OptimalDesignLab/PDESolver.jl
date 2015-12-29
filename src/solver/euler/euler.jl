# this file contains the functions to evaluate the right hand side of the weak form in the pdf
# res_vec stands for semi-linear form, contents inside the inv(M)(...) parenthesis on pg. 7 of Jared's derivation

# Organization:
# The code is organized into 3 levels of functions.  High level functions
# only call mid level functions, passing the Mesh, Equation, and SBP objects
# to them.  At this point, there shouldn't be any need for more high level
# functions.  High level functions should use AbstractMesh and 
# AbstractEulerData as the types of the mesh and equation objects, to signify 
# that they are implementation independent.
#
# The mid level functions either call SBP functions or have loops over
# the arrays stored in the Equation, calling low level functions
# to calculate quantities at each node.  Use the ArrayViews package to pass
# portions of a larger array to the low level functions.  These functions 
# should use AbstractMesh and EulerData as the types of the mesh and equation
# objects, to signify they are independent of the mesh implementation but 
# do interact with the fields of the equation object.
#
# Low level functions calculate quantities for a single node.  Although
# these functions will likely take in an ArrayView, their argument
# types should be AbstractArrays{T, N}, specifying T according to 
# what type the original array is (from either the Equation or Mesh object)
# They should also take the eqn.params object as an argument, which has
# all of the static parameters as the the equation object, so it can be
# used to dispatch to the correct method based on the static paremters.


# For high and mid level functions, the proper order for arguments 
# is (mesh, sbp, eqn, opts).
# For low level functions, the order of arguments should be:  eqn.params first, 
# then all input (read-only) arguments, then any output arguments.

# Rules: 
# 1. function that take a composite type with abstract fields *cannot* use those
# fields directly (they can pass those fields to other functions however
# This leads directly to a two level structure for the code: high level function
# that take in composite types and low level function that take in arrays and
# perform calculations on them

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
### EulerEquationMod General Description
This module is organized into 3 levels of functions: high, middle, and low.

The high level functions take the mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerEquation, and opts (options dictionary), They do not know the types their 
arguments are paramaterized on. There is only one method for each high level function.  All they do is call mid level functions.

Mid level functions take the same arguments as high level functions but know
the types they are paramaterized on and do the right thing for all parameters.
There is only one method for each mid level function.  Mid level functions
typically pass pieces of the large arrays stored in the eqn and mesh objects
to low level functions that actually do computation.

Low level functions usually take in vectors of whatever values they need, do
a computation, and store the result in an argument vector.  These function 
are usually passed an ArrayView, so writing to the argument vector writes 
directly to the large array stored in the eqn object.  Thus there is no 
loss of efficiency by using low level functions.
The Params type is also passed to low level function, which has all the 
type parameters as the EulerEquation object, so it can be used for dispatch.
"""
export evalEuler, init, convertToEntropy_, convertToConservative_

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
  eqn.q_vec and eqn.res_vec exist for reusable storage *outside* the residual
  evaluation.  They should never be used inside the residual evaluation.

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

  eqn.params.t = t  # record t to params
  
  dataPrep(mesh, sbp, eqn, opts)
#println("dataPrep @time printed above")


  evalVolumeIntegrals(mesh, sbp, eqn, opts)
  #println("volume integral @time printed above")

  # delete this if unneeded or put it in a function.  It doesn't belong here,
  # in a high level function.
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
  #println("boundary integral @time printed above")


  addStabilization(mesh, sbp, eqn, opts)
  #println("edge stabilizing @time printed above")


  
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
function init{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, eqn::AbstractEulerData{Tsol, Tres}, opts, pmesh=mesh)

  println("\nInitializing Euler module")
  # get BC functors
  getBCFunctors(mesh, sbp, eqn, opts)
  getBCFunctors(pmesh, sbp, eqn, opts)


  return nothing
end


function majorIterationCallback(itr::Integer, mesh::AbstractMesh, sbp::SBPOperator, eqn::AbstractEulerData, opts)

#  println("Performing major Iteration Callback")
#=
  if itr == 0
    #----------------------------------------------------------------------------
    # Storing the initial density value at all the nodes
    global const vRho_act = zeros(mesh.numNodes)
    k = 1
    for counter1 = 1:4:length(eqn.q_vec)
      vRho_act[k] = eqn.q_vec[counter1]
      k += 1
    end
    println("Actual Density value succesfully extracted")
    #--------------------------------------------------------------------------
  else
    #--------------------------------------------------------------------------
    # Calculate the error in density
    vRho_calc = zeros(vRho_act)
    k = 1
    for counter1 = 1:4:length(eqn.q_vec)
      vRho_calc[k] = eqn.q_vec[counter1]
      k += 1
    end
    ErrDensityVec = vRho_calc - vRho_act
    # ErrDensity1 = norm(ErrDensityVec, 1)/mesh.numNodes
    # ErrDensity2_discrete = norm(ErrDensityVec, 2)/mesh.numNodes
    # println("DensityErrorNormL1 = ", ErrDensity1) 
    ErrDensity2 = 0.0
    k = 1
    for counter1 = 1:length(ErrDensityVec)
      ErrDensity2 += real(ErrDensityVec[counter1])*eqn.M[k]*real(ErrDensityVec[counter1])
      k += 4
    end
    ErrDensity2 = sqrt(ErrDensity2)
    println("DensityErrorNormL2 = ", ErrDensity2)
    # println("Discrete density error norm L2 = ", ErrDensity2_discrete)



    #--------------------------------------------------------------------------
  end
=#
  if opts["write_entropy"]
    # calculate the entropy norm
    val = zero(Float64)
    for i=1:mesh.numDofPerNode:mesh.numDof
      q_vals = view(eqn.q_vec, i:(i+3))
      s = calcEntropy(eqn.params, q_vals)
      val += real(s)*eqn.M[i]*real(s)
    end
    f = open(opts["write_entropy_fname"], "a+")
    println(f, itr, " ", eqn.params.t, " ",  val)
    close(f)
  end


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
function dataPrep{Tmsh,  Tsol, Tres}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, eqn::AbstractEulerData{Tsol, Tres}, opts)
# gather up all the data needed to do vectorized operatinos on the mesh
# calculates all mesh wide quantities in eqn

#println("Entered dataPrep()")


  # apply filtering to input
  if eqn.params.use_filter
    applyFilter(mesh, sbp, eqn, eqn.q, opts)
  end

  # zero out res
  fill!(eqn.res, 0.0)
  fill!(eqn.res_edge, 0.0)
 
  getAuxVars(mesh, eqn)
#  println("getAuxVars @time printed above")

  if opts["check_density"]
    checkDensity(eqn)
    # println("checkDensity @time printed above")
  end

  if opts["check_pressure"]
    checkPressure(eqn)
    # println("checkPressure @time printed above")
  end

  # calculate fluxes
#  getEulerFlux(eqn, eqn.q, mesh.dxidx, view(flux_parametric, :, :, :, 1), view(flux_parametric, :, :, :, 2))
  getEulerFlux(mesh, sbp,  eqn, opts)
#  println("getEulerFlux @time printed above")



   getBCFluxes(mesh, sbp, eqn, opts)
#   println("getBCFluxes @time printed above")
  
   stabscale(mesh, sbp, eqn)
#  println("stabscale @time printed above")

#  println("finished dataPrep()")
  return nothing
end # end function dataPrep


@doc """
### EulerEquationMod.getBCFluxes

  This function calls other functions to calculate the boundary fluxes, passing
  them pieces of the array needed.  This populates eqn.bndryflux.  It also
  calls writeBoundary() to do any requested output.

  This is a mid level function
"""->
# this is a mid level function
function getBCFluxes(mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerData, opts)
  #get all the fluxes for all the boundary conditions and save them in eqn.bndryflux

  #println("mesh.bndry_funcs = ", mesh.bndry_funcs)

  for i=1:mesh.numBC
  #  println("computing flux for boundary condition ", i)
    functor_i = mesh.bndry_funcs[i]
    start_index = mesh.bndry_offsets[i]
    end_index = mesh.bndry_offsets[i+1]
    bndry_facenums_i = view(mesh.bndryfaces, start_index:(end_index - 1))
    bndryflux_i = view(eqn.bndryflux, :, :, start_index:(end_index - 1))
 
    # call the function that calculates the flux for this boundary condition
    # passing the functor into another function avoid type instability
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



#------------------------------------------------------------------------------
# sanity check functions
#------------------------------------------------------------------------------
@doc """
### EulerEquationMod.checkDensity

  This function checks that the density is positive. If not, an error is thrown.
  Because density is stored in the eqn.q array, performing this check takes 
  very little time.

  Arguments:
    * EulerData

  This is a mid level function.
"""->
# mid level function
function checkDensity{Tsol}(eqn::EulerData{Tsol})
# check that density is positive

(ndof, nnodes, numel) = size(eqn.q)
q_cons = zeros(Tsol, ndof)  # conservative variables
for i=1:numel
  for j=1:nnodes
    convertToConservative(eqn.params, view(eqn.q, :, j, i), q_cons)
    @assert( real(q_cons[1]) > 0.0, "element $i, node $j. Density < 0")
  end
end

return nothing

end


@doc """
### EulerEquationMod.checkPressure

  This function checks that pressure is positive.  If not, an error is thrown.
  Because pressure is stored in the auxiliary variables array, this check 
  takes very little time.

  Arguments:
    * EulerData

  This is a mid level function
"""->
function checkPressure(eqn::EulerData)
# check that density is positive

(ndof, nnodes, numel) = size(eqn.q)

for i=1:numel
  for j=1:nnodes
    aux_vars = view(eqn.aux_vars,:, j, i)
    press = @getPressure(aux_vars)
    @assert( real(press) > 0.0, "element $i, node $j")
  end
end

return nothing

end


#------------------------------------------------------------------------------
# functions that evaluate terms in the weak form
#------------------------------------------------------------------------------
@doc """
### EulerEquationMod.evalVolumeIntegrals

  This function evaluates the volume integrals of the Euler equations by
  calling the appropriate SBP function on the Euler flux stored in
  eqn.flux_parametric.
  This function knows the dimension of the equation and does the right
  thing for all dimensions.  eqn.res is updated with the result

  This is a mid level function.
"""
# mid level function
function evalVolumeIntegrals{Tmsh,  Tsol, Tdim}(mesh::AbstractMesh{Tmsh}, 
                             sbp::SBPOperator, eqn::EulerData{Tsol, Tdim}, opts)
  
  if opts["Q_transpose"] == true
    for i=1:Tdim
      weakdifferentiate!(sbp, i, view(eqn.flux_parametric, :, :, :, i), eqn.res, trans=true)
    end
  else
    for i=1:Tdim
      #TODO: do this more efficiently
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

  boundaryintegrate!(sbp, mesh.bndryfaces, eqn.bndryflux, eqn.res)


  return nothing

end  # end evalBoundaryIntegrals


@doc """
### EulerEquationMod.addStabilization

  This function add whatever form of stabilization opts specifies by calling
  the appropriate function.  Some arrays may need to be populated before 
  calling this function, depending on the type of stabilization used.

  This is a mid level function
"""->
# mid level function
function addStabilization{Tmsh,  Tsol}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, eqn::EulerData{Tsol}, opts)

#  println("==== start of addStabilization ====")

  # ----- Edge Stabilization -----#
  if eqn.params.use_edgestab
#    println("applying edge stabilization")
    if opts["use_edge_res"]
      edgestabilize!(mesh, sbp, eqn, mesh.interfaces, eqn.q, mesh.coords, 
                     mesh.dxidx, mesh.jac, eqn.edgestab_alpha, eqn.stabscale, 
                     eqn.res, eqn.res_edge)
    else
      edgestabilize!(sbp, mesh.interfaces, eqn.q, mesh.coords, mesh.dxidx, 
                     mesh.jac, eqn.edgestab_alpha, eqn.stabscale, eqn.res)
    end
  end

  # ----- Filtering -----
  if eqn.params.use_res_filter
#    println("applying residual filter")
    applyFilter(mesh, sbp, eqn, eqn.res, opts, trans=true)
  end

  # ----- Artificial Dissipation -----
  if eqn.params.use_dissipation
#    println("applying artificial dissipation")
    applyDissipation(mesh, sbp, eqn, eqn.q, opts)
  end

#  println("==== end of addStabilization ====")



  return nothing

end


# some helper functions
# this is a test for an algorithmic differentiation package
function getEulerJac_wrapper{T}(q::AbstractArray{T,1}, F::AbstractArray{T,1})
  dir = [1.0, 0.0]
#  F = zeros(T, 4)
  sbp = TriSBP{Float64}()
  mesh = PumiMesh2{Float64}(".null", "../../mesh_files/quarter_vortex3l.smb", 1, sbp; dofpernode=4)
  eqn = EulerData1{T, T}(mesh, sbp)

  @time getEulerFlux(eqn, q, dir, F)

  return F

end

fluxJac = forwarddiff_jacobian!(getEulerJac_wrapper, Float64, fadtype=:dual; n=4, m=4)


@doc """
### EulerEquationMod.getAuxVars

  This function calculates any extra variables that are stored across the mesh
  using the conservative variables eqn.q.  Currently only calculates pressure.

  Thi is a mid level function
"""->
# mid level function
function getAuxVars{Tmsh, Tsol, Tdim}(mesh::AbstractMesh{Tmsh}, eqn::EulerData{Tsol, Tdim})
# calculate all auxiliary variables

  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      q_vals = view(eqn.q, :, j, i)

      # calculate pressure
      press = calcPressure(eqn.params, q_vals)
      @setPressure(eqn.aux_vars, j, i, press)
    end
  end

  return nothing
end


@doc """ 
### EulerEquationMod.getEulerFlux

  This function calculates the Euler flux across the entire mesh by passing
  pieces of the eqn.q, eqn.aux_vars, eqn.f_xi and eqn.params to a low level
  function.  The flux is calculated in the xi and eta directions, 
  scaled (mulitiplied) by the mapping jacobian (so that when performing the 
  integral we don't have to explictly divide by the jacobian, it just cancels 
  out with the jacobian factor introduced here.

  Calls writeFlux to do any requested output.

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
        # don't do an array view because strided views are type-unstable
        nrm[1] = mesh.dxidx[k, 1, j, i]
        nrm[2] = mesh.dxidx[k, 2, j, i]
        flux = view(eqn.flux_parametric, :, j, i, k)

	# this will dispatch to the proper calcEulerFlux
        calcEulerFlux(eqn.params, q_vals, aux_vars, nrm, flux)
      end
    end
  end


  writeFlux(mesh, sbp, eqn, opts)

  return nothing
end


@doc """
### EulerEquationMod.writeFlux

  This function writes the real part of Euler flux to a file named Fxi.dat, 
  space delimited, controlled by the input options 'writeflux', of type Bool.

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


@doc """
### EulerEquationMod.applyMassMatrixInverse

  This function multiplies eqn.res_vec (the residual in vector form  by eqn.Minv,
  the diagonal mass matrix.  This is a very fast function because all values
  are precomputed and stored linearly in memory.

  This is a mid level function, and does the correct thing regardless of the
  dimension of the equation.

  Aliasing restrictions: none
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




@doc """
### EulerEquationMod.disassembleSolution

  This takes eqn.q_vec (the initial state), and disassembles it into eqn.q, the
  3 dimensional array.  This function uses mesh.dofs
  to speed the process.

  This function also calls writeQ to do any requested output.

  Inputs:
    mesh
    sbp
    eqn
    opts

  This is a mid level function, and does the right thing regardless of equation
  dimension.

  Aliasing restrictions: none
"""->
# mid level function (although it doesn't need Tdim)
function disassembleSolution{Tmsh, Tsol, Tdim, T}(mesh::AbstractMesh{Tmsh}, sbp,
                             eqn::EulerData{Tsol, Tdim}, opts, 
                             q_arr::AbstractArray{T, 3}, 
                             q_vec::AbstractArray{T, 1})
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
  delimited file called q.dat, controlled by the input options 'writeq', of type bool

  This is a high level function.
"""->
function writeQ(mesh, sbp, eqn, opts)

  if !eqn.params.writeq
    return nothing
  end

  fname = "q.dat"
  rmfile(fname)
  writedlm(fname, eqn.q)

  return nothing
end



@doc """
### EulerEquationMod.assembleSolution

  This function takes the 3D array of variables in arr and 
  reassmbles is into the vector res_vec.  Note that
  This is a reduction operation and zeros res_vec before performing the 
  operation, unless zero_res is set to false

  This is a mid level function, and does the right thing regardless of
  equation dimension
"""->
# mid level function (although it doesn't need Tdim)
function assembleSolution{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh}, 
                         sbp::SBPOperator, eqn::EulerData{Tsol}, opts, 
                         arr::Abstract3DArray, res_vec::AbstractArray{Tres,1}, 
                         zero_resvec=true)
# arr is the array to be assembled into res_vec

#  println("in assembleSolution")

  if zero_resvec
    fill!(res_vec, 0.0)
  end


  for i=1:mesh.numEl  # loop over elements
    for j=1:mesh.numNodesPerElement
      for k=1:4  # loop over dofs on the node
        dofnum_k = mesh.dofs[k, j, i]
        res_vec[dofnum_k] += arr[k,j,i]
      end
    end
  end
  
  return nothing
end


#------------------------------------------------------------------------------
# convert variables without checking what type they currently are (unsafe)
#------------------------------------------------------------------------------
# convert to entropy variables without checking if we are already in them
@doc """
### EulerEquationMod.convertToEntropy_

  Converts the conservative variables at a node to the entropy variables
  regardless of the static parameter var_type.

  Inputs:
  params  : ParamType{2} used to dispatch to the proper method
  qe  : vector (of length 4) of conservative variables
  
  Inputs/outputs
  qc : vector (of length 4) of entropy variables.  Contents of vector are
       overwritten

  Aliasing: none (qc and qe *can* be the same vector)

  Ths is a low level function.
"""->

function convertToEntropy_{Tsol}(params::ParamType{2}, qc::AbstractArray{Tsol,1},
                             qe::AbstractArray{Tsol, 1})

  gamma = params.gamma
  gamma_1 = params.gamma_1

  k1 = 0.5*(qc[2]^2 + qc[3]^2)/qc[1]
  rho_int = qc[4] -k1
  s = log(gamma_1*rho_int/(qc[1]^gamma))
  fac = 1.0/rho_int
  tmp1 = -qc[1]*fac
  qe[1] = (rho_int*(gamma + 1 -s) - qc[4])*fac
  qe[2] = qc[2]*fac
  qe[3] = qc[3]*fac
  qe[4] = tmp1

  return nothing
end

#=
# define vector version used in startup.jl
function convertToEntropy_{Tmsh, Tsol, Tdim, Tres}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, eqn::EulerData{Tsol, Tdim, Tres}, opts, q_vec::AbstractArray{Tsol, 1})

  for i=1:mesh.numDofPerNode:mesh.numDof
    q_view = view(q_vec, i:(i+mesh.numDofPerNode-1))
    convertToEntropy_(eqn.params, q_view, q_view)
  end

  return nothing
end
=#

@doc """
### EulerEquationMod.convertToConservative_

  Converts  the entropy variables at a node to the conservative variables
  regardless of the static parameter var_type.

  Inputs:
  params  : ParamType{2} used to dispatch to the proper method
  qe  : vector (of length 4) of entropy variables
  
  Inputs/outputs
  qc : vector (of length 4) of conservative variables.  Contents of vector are
       overwritten

  Aliasing: none (qc and qe *can* be the same vector)

  Ths is a low level function.
"""->
function convertToConservative_{Tsol}(params::ParamType{2}, 
                  qe::AbstractArray{Tsol,1}, qc::AbstractArray{Tsol, 1})

  gamma = params.gamma
  gamma_1 = params.gamma_1
  k1 = 0.5*(qe[2]^2 + qe[3]^2)/qe[4]
  s = gamma - qe[1] + k1
  rho_int = exp(-s/gamma_1)*(gamma_1/(-qe[4])^gamma)^(1/gamma_1)
  qc[1] = -qe[4]*rho_int
  qc[2] = qe[2]*rho_int
  qc[3] = qe[3]*rho_int
  qc[4] = (1.0 - k1)*rho_int
end




#------------------------------------------------------------------------------
# converting to conservative variables safely
#------------------------------------------------------------------------------
@doc """
### EulerEquationMod.convertToConservative

  Converts  the entropy variables at a node to the conservative variables.

  Inputs:
  params  : ParamType{2, :entropy} used to dispatch to the proper method
  qe  : vector (of length 4) of entropy variables
  
  Inputs/outputs
  qc : vector (of length 4) of conservative variables.  Contents of vector are
       overwritten

  Aliasing: none (qc and qe *can* be the same vector)

  Ths is a low level function.
"""->
function convertToConservative{Tsol}(params::ParamType{2, :entropy}, 
                  qe::AbstractArray{Tsol,1}, qc::AbstractArray{Tsol, 1})

  convertToConservative_(params, qe, qc)

end

@doc """
# low level function

  Converts conservative variables to conservative variables (ie. it
  copies the input to the output).  This method exists to values can be 
  converted without knowing whether they are conservative or entropy.

  Aliasing restrictions: none (qc and qe *can* be the same vector)

"""->
function convertToConservative{Tsol}(params::ParamType{2, :conservative}, 
                  qe::AbstractArray{Tsol,1}, qc::AbstractArray{Tsol, 1})
# maybe the compiler will be smart enough to make this a no-op?
  for i=1:length(qe)
    qc[i] = qe[i]
  end

  return nothing
end


@doc """
# mid level function

  Converts the array (3D form) of entropy variables to conservative variables 
  in place.  If the array is already in conservative variables this is a no-op.
  Other methods exist if q_arr is a vector.

  Inputs:
    mesh
    sbp
    eqn
    opts

  Inputs/Outputs:
    q_arr: array of values (3D) to be converted

  Aliasing: no restrictions
"""->
#3D array entropy -> conservative
function convertToConservative{Tmsh, Tsol, Tdim, Tres}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, eqn::EulerData{Tsol, Tdim, Tres, :entropy}, opts, q_arr::AbstractArray{Tsol, 3})

  for i=1:mesh.numEl  # loop over elements
    for j=1:mesh.numNodesPerElement
	q_view = view(q_arr, :, j, i)  # reuse memory
	convertToConservative(eqn.params, q_view, q_view)
    end
  end

  return nothing
end

# 3D array conservative -> conservative
function convertToConservative{Tmsh, Tsol, Tdim, Tres}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, eqn::EulerData{Tsol, Tdim, Tres, :conservative}, opts, q_arr::AbstractArray{Tsol, 3})

  return nothing
end


# q_vec conversion entropy -> conservative
function convertToConservative{Tmsh, Tsol, Tdim, Tres}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, eqn::EulerData{Tsol, Tdim, Tres, :entropy}, opts, q_vec::AbstractArray{Tsol, 1})

  for i=1:mesh.numDofPerNode:mesh.numDof
    q_view = view(q_vec, i:(i+mesh.numDofPerNode-1))
    convertToConservative(eqn.params, q_view, q_view)
  end

  return nothing
end

# q_vec conversion conservative -> conservative
function convertToConservative{Tmsh, Tsol, Tdim, Tres}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, eqn::EulerData{Tsol, Tdim, Tres, :conservative}, opts, q_arr::AbstractArray{Tsol, 1})

  return nothing
end



#------------------------------------------------------------------------------
# converting to entropy variables
#------------------------------------------------------------------------------

# converting to entropy variables with checking
@doc """
# low level function
  Converts the conservative variables at a node to entropy variables

  Input:
  params : ParamType{s, :conservative}
  qc  : vector (of length 4) of conservative variables

  Outputs:
  qe : vector (of length 4) of conservative variables.  Contents of vector are
       overwritten

  Aliasing restrictions: none (qc and qe *can* be the same vector)
"""->
function convertToEntropy{Tsol}(params::ParamType{2, :conservative}, 
               qc::AbstractArray{Tsol,1}, qe::AbstractArray{Tsol,1})

  convertToEntropy_(params, qc, qe)

end


@doc """
# low level function
  Converts the entropy variables to entropy variables (ie. it copies the 
  input to the output).  This method exists so variables can be converted 
  to entropy variables without knowing what type they are.

  Aliasing restrictions: none (qc and qe *can* be the same vector)

"""->
function convertToEntropy{Tsol}(params::ParamType{2, :entroyp}, 
               qc::AbstractArray{Tsol,1}, qe::AbstractArray{Tsol,1})

  for i=1:length(qc)
    qe[i] = qc[i]
  end

  return nothing
end



@doc """
# mid level function

  Converts the array (3D form) of conservative variables to entropy variables 
  in place.  If the array is already in entropy variables this is a no-op

  Methods also exist for the 1D form.

  Inputs:
    mesh
    sbp
    eqn
    opts

  Inputs/Outputs:
    q_arr: array of values (3D) to be converted

  Aliasing: no restrictions
"""->
# 3D array conservative -> entropy
function convertToEntropy{Tmsh, Tsol, Tdim, Tres}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, eqn::EulerData{Tsol, Tdim, Tres, :conservative}, opts, q_arr::AbstractArray{Tsol, 3})
i


  for i=1:mesh.numEl  # loop over elements
    for j=1:mesh.numNodesPerElement
	q_view = view(q_arr, :, j, i)  # reuse memory
	convertToEntropy(eqn.params, q_view, q_view)
    end
  end

  return nothing
end


# 3D array entropy -> entropy
function convertToEntropy{Tmsh, Tsol, Tdim, Tres}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, eqn::EulerData{Tsol, Tdim, Tres, :entropy}, opts, q_arr::AbstractArray{Tsol, 3})

  return nothing
end



# q_vec conversion conservative -> entropy
function convertToEntropy{Tmsh, Tsol, Tdim, Tres}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, eqn::EulerData{Tsol, Tdim, Tres, :conservative}, opts, q_vec::AbstractArray{Tsol, 1})

  for i=1:mesh.numDofPerNode:mesh.numDof
    q_view = view(q_vec, i:(i+mesh.numDofPerNode-1))
    convertToEntropy_(eqn.params, q_view, q_view)
  end

  return nothing
end

# q_vec conversion entropy -> entropy
function convertToEntropy{Tmsh, Tsol, Tdim, Tres}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, eqn::EulerData{Tsol, Tdim, Tres, :entropy}, opts, q_arr::AbstractArray{Tsol, 1})

  return nothing
end



# calculating the Euler flux
#------------------------------------------------------------------------------
@doc """
### EulerEquationMod.calcEulerFlux

   This function calculates the Euler flux from the conservative variables at
   a single node in a particular direction.  2D only.

   Inputs:
   params  : ParamaterType{2, :conservative}
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
function calcEulerFlux{Tmsh, Tsol, Tres}(params::ParamType{2, :conservative}, 
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
# low level function
    Calculates the Euler flux from entropy variables

    Inputs:
    params : ParameterType{2, :entropy}
    q : vector of entropy variables
    aux_vars : vector of auxiliary variables
    dir : vector specifying the direction to caculate the flux

    Inputs/Outputs:
    F  : vector to populate with the flux
 
    This is a low level function.  The second static parameter of 
    the ParameterType is used to dispatch to the right method for
    entropy or conservative variables
"""->
function calcEulerFlux{Tmsh, Tsol, Tres}(params::ParamType{2, :entropy}, 
                       q::AbstractArray{Tsol,1}, 
                       aux_vars::AbstractArray{Tres, 1}, 
                       dir::AbstractArray{Tmsh},  F::AbstractArray{Tsol,1})

  gamma = params.gamma
  gamma_1 = params.gamma_1

  # calculate some intermediate quantities
  k1 = 0.5*(q[2]^2 + q[3]^2)/q[4]  # a constant from Hughes' paper
  s = gamma - q[1] + k1    # entropy
    # internal energy (rho*i in Hughes) - not specific internal energy e
  rho_int = exp(-s/gamma_1)*(gamma_1/((-q[4])^gamma))^(1/gamma_1)
  U = q[2]*dir[1] + q[3]*dir[2]
  fac = rho_int/q[4]

  # now we can actually calculate the flux
  F[1] = q[4]*U*fac
  F[2] = (dir[1]*gamma_1*q[4] - q[2]*U)*fac
  F[3] = (dir[2]*gamma_1*q[4] - q[3]*U)*fac
  F[4] = U*(k1 - gamma)*fac

  return nothing
end

@doc """
### EulerEquationMod.calcEulerFlux
  This is the 3D method.  All arguments are same as the 2D version.
"""->
# low level function
function calcEulerFlux{Tmsh, Tsol}(params::ParamType{3, :conservative}, 
                       q::AbstractArray{Tsol,1}, dir::AbstractArray{Tmsh},  
                       F::AbstractArray{Tsol,1})
# calculates the Euler flux in a particular direction at a point
# eqn is the equation type
# q is the vector (of length 5), of the conservative variables at the point
# dir is a vector of length 3 that specifies the direction
# F is populated with the flux (is a vector of length 5)
# 3D  only

# once the Julia developers fix slice notation and speed up subarrays, we can make a faster
# vectorized version of this

  press = calcPressure(params, q)
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
    params : ParamType{2, :conservative}
    q  : vector of conservative variables

  The parameter of params determines whether the 2D or 3D method is dispatched.

  This is a low level function.

  Aliasing restrictions: none
"""->
# low level function
function calcPressure{Tsol}(params::ParamType{2, :conservative}, 
                            q::AbstractArray{Tsol,1} )
  # calculate pressure for a node
  # q is a vector of length 4 of the conservative variables

  return  (params.gamma_1)*(q[4] - 0.5*(q[2]*q[2] + q[3]*q[3])/q[1])
  
end


@doc """
### EulerEquationMod.calcPressure
  This function calculates pressure using the entropy variables.

  Inputs:
    params : ParamType{2, :entropy}
    q  : vector of entropy varaibles

  returns pressure

  Aliasing restrictions: none
"""->
function calcPressure{Tsol}(params::ParamType{2, :entropy}, 
                            q::AbstractArray{Tsol,1})

  gamma = params.gamma
  gamma_1 = params.gamma_1
  
  k1 = 0.5*(q[2]*q[2] + q[3]*q[3])/q[4]  # a constant from Hughes' paper
  s = gamma - q[1] + k1    # entropy
  rho_int = exp(-s/gamma_1)*(gamma_1/((-q[4])^gamma))^(1/gamma_1)
  return gamma_1*rho_int
end


@doc """
### EulerEquationMod.calcPressure

  3D method.  See 2D method documentation
"""->
# low level function
function calcPressure{Tsol}(params::ParamType{3, :conservative}, 
                            q::AbstractArray{Tsol,1} )
  # calculate pressure for a node
  # q is a vector of length 5 of the conservative variables

  return  (params.gamma_1)*(q[5] - 0.5*(q[2]*q[2] + q[3]*q[3] + q[4]*q[4])/q[1])
  

end



@doc """
### EulerEquationMod.calcSpeedofSound

  This function calculates the speed of sound at a node and returns it.
  Methods are available for both conservative and entropy variables.

  Inputs:
    params:  ParamType{2, var_type}
    q  vector of solution variables at a node

  Returns: speed of sound

  This is a low level function

  Aliasing restrictions: none
"""->

function calcSpeedofSound{Tsol}(params::ParamType{2, :conservative}, q::AbstractArray{Tsol, 1})
# calculates teh speed of sond at a node
  pressure = calcPressure(params, q)
  return sqrt((params.gamma*pressure)/q[1])

end



function calcSpeedofSound{Tsol}(params::ParamType{2, :entropy}, 
                                q::AbstractArray{Tsol, 1})
# calculate speed of sound using the same formula as conservative variables,
# just rewriting all variables in entropy variables

#  printbacktrace()
#  println("q = ", q)

  gamma = params.gamma
  gamma_1 = params.gamma_1
  k1 = 0.5*(q[2]^2 + q[3]^2)/q[4]  # a constant from Hughes' paper
  pressure = calcPressure(params, q)
  s = gamma - q[1] + k1    # entropy

  rho_int = exp(-s/gamma_1)*(gamma_1/((-q[4])^gamma))^(1/gamma_1)
  rho = -q[4]*rho_int

#  println("rho = ", rho)
  return sqrt((params.gamma*pressure)/rho)
end


@doc """
### EulerEquationMod.calcEntropy

  This function calculates the entropy at a node and returns it.  Method are 
  available for conservative and entropy variables

  Inputs:
    params: ParamType{2, var_type}, used to dispatch to the right method.
    q: vector of solution variables at a node.

  Returns: entropy

  This is a low level function

  Aliasing Restrictions: none

"""->
function calcEntropy{Tsol}(params::ParamType{2, :conservative}, 
                           q::AbstractArray{Tsol,1} )

  gamma = params.gamma
  gamma_1 = params.gamma_1

  rho_int = q[4] - 0.5*(q[2]*q[2] + q[3]*q[3])/q[1]
  return log(gamma_1*rho_int/(q[1]^gamma))
end

function calcEntropy{Tsol}(q::AbstractArray{Tsol,1}, params::ParamType{2, :entropy})

  gamma = params.gamma
  gamma_1 = params.gamma_1

  return gamma - q[1] + 0.5*(q[2]*q[2] + q[3]*q[3])/q[4]
end


@doc """
### EulerEquationMod.calcA0

  This function calculates the A0 (ie. dq/dv, where q are the conservative 
  and v are the entropy variables) for a node, and stores it A0

  The formation of A0 is given in Hughes

  Inputs:
    params : ParamType{2, :entropy}
    q  : vector of entropy variables, length 4


  Inputs/Outputs:
  A0 : 4x4 matrix populated with A0.  Overwritten

  Aliasing restrictions: none
"""->
function calcA0{Tsol}(params::ParamType{2, :entropy}, q::AbstractArray{Tsol,1},
                      A0::AbstractArray{Tsol, 2})


  gamma = params.gamma
  gamma_1 = params.gamma_1

  k1 = 0.5*(q[2]^2 + q[3]^2)/q[4]  # a constant from Hughes' paper
  k2 = k1 - gamma
  k3 = k1*k1 - 2*gamma*k1 + gamma
#    k4 = k2 - gamma_1
  s = gamma - q[1] + k1    # entropy

  rho_int = exp(-s/gamma_1)*(gamma_1/((-q[4])^gamma))^(1/gamma_1)

  fac = rho_int/(gamma_1*q[4])

  # calculate the variables used in Hughes A.1
  c1 = gamma_1*q[4] - q[2]*q[2]
  c2 = gamma_1*q[4] - q[3]*q[3]

  d1 = -q[2]*q[3]

  e1 = q[2]*q[4]
  e2 = q[3]*q[4]

  # populate the matrix
  # the matrix is symmetric, but we don't use it because I think populating
  # the matrix will be faster if the matrix is write-only
  A0[1,1] = -q[4]*q[4]*fac
  A0[2,1] = e1*fac
  A0[3,1] = e2*fac
  A0[4,1] = q[4]*(1-k1)*fac
  A0[1,2] = e1*fac  # symmetric
  A0[2,2] = c1*fac
  A0[3,2] = d1*fac
  A0[4,2] = q[2]*k2*fac
  A0[1,3] = e2*fac  # symmetric
  A0[2,3] = d1*fac  # symmetric
  A0[3,3] = c2*fac
  A0[4,3] = q[3]*k2*fac
  A0[1,4] = q[4]*(1-k1)*fac  # symmetric
  A0[2,4] = q[2]*k2*fac   # symmetric
  A0[3,4] = q[3]*k2*fac  # symmetric
  A0[4,4] = -k3*fac

    return nothing
end

@doc """
# EulerEquationMod.calcA0Inv

  Calculates inv(A0), where A0 = dq/dv, where q are the conservative variables
  at a node and v are the entropy varaibles at a node, using the entropy 
  variables.  

  Inputs:
    params : ParamType{2, :entropy}
    q  : vector (length 4) of entropy variables at a node

  Inputs/Outputs:
    A0inv : matrix to be populated with inv(A0).  Overwritten.

  Aliasing restrictions: none
"""->
function calcA0Inv{Tsol}(params::ParamType{2, :entropy},
                   q::AbstractArray{Tsol,1},  
                   A0inv::AbstractArray{Tsol, 2})
  gamma = params.gamma
  gamma_1 = params.gamma_1

  k1 = 0.5*(q[2]^2 + q[3]^2)/q[4]  # a constant from Hughes' paper
   s = gamma - q[1] + k1    # entropy

  rho_int = exp(-s/gamma_1)*(gamma_1/((-q[4])^gamma))^(1/gamma_1)

  fac = -1/(rho_int*q[4])

  d1 = -q[2]*q[3]
  e1 = q[2]*q[4]
  e2 = q[3]*q[4]


  A0inv[1,1] = (k1*k1 + gamma)*fac
  A0inv[2,1] = k1*q[2]*fac
  A0inv[3,1] = k1*q[3]*fac
  A0inv[4,1] = (k1 + 1)*q[4]*fac
  A0inv[1,2] = k1*q[2]*fac  # symmetry
  A0inv[2,2] = (q[2]*q[2] - q[4])*fac
  A0inv[3,2] = -d1*fac
  A0inv[4,2] = e1*fac
  A0inv[1,3] = k1*q[3]*fac  # symmetry
  A0inv[2,3] = -d1*fac  # symmetry
  A0inv[3,3] = (q[3]*q[3] - q[4])*fac
  A0inv[4,3] = e2*fac
  A0inv[1,4] = (k1 + 1)*q[4]*fac  # symmetry
  A0inv[2,4] = e1*fac  # symmetry
  A0inv[3,4] = e2*fac  # symmetry
  A0inv[4,4] = q[4]*q[4]*fac

    return nothing
end

@doc """
### EulerEquationMod.matVecA0inv

  This function takes a 3D array and multiplies it in place by the inv(A0) 
  matrix (calculated at each node), inplace, (where A0 =dq/dv, where q are the 
  conservative variables and v are some other variables), inplace.
  Methods are available for conservative and entropy variables.

  For conservative variables, A0 is the identity matrix and this is a no-op

  Inputs:
    mesh
    sbp
    eqn
    opts

  Inputs/Outputs:
    res_arr: the array to multiply 

  Aliasing restrictions: none
"""->
function matVecA0inv{Tmsh, Tsol, Tdim, Tres}(mesh::AbstractMesh{Tmsh}, 
                     sbp::SBPOperator, 
                     eqn::EulerData{Tsol, Tdim, Tres, :entropy}, opts, 
                     res_arr::AbstractArray{Tsol, 3})
# multiply a 3D array by inv(A0) in-place, useful for explicit time stepping
# res_arr *can* alias eqn.q safely
  A0inv = Array(Tsol, mesh.numDofPerNode, mesh.numDofPerNode)
  res_vals = Array(Tsol, mesh.numDofPerNode)
  q_vals = Array(Tsol, mesh.numDofPerNode)
  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      # copy values into workvec
      for k=1:mesh.numDofPerNode
	q_vals[k] = eqn.q[k, j, i]
	res_vals[k] = res_arr[k, j, i]
      end

      res_view = view(res_arr, :, j, i)
      # get A0Inv for this node
      calcA0Inv(eqn.params, q_vals, A0inv)

      smallmatvec!(A0inv, res_vals, res_view)
    end
  end

  return nothing
end

# no-op, because for conservative variables this is A0inv is the identity matrix
function matVecA0inv{Tmsh, Tsol, Tdim, Tres}(mesh::AbstractMesh{Tmsh}, 
                     sbp::SBPOperator, 
                     eqn::EulerData{Tsol, Tdim, Tres, :conservative}, 
                     opts, res_arr::AbstractArray{Tsol, 3})

  return nothing
end

@doc """
### EulerEquationMod.matVecA0

  This function is the same as matVecA0inv, except it multiplies by A0 not
  A0inv.  See its documention.

"""->
function matVecA0{Tmsh, Tsol, Tdim, Tres}(mesh::AbstractMesh{Tmsh}, 
                  sbp::SBPOperator, eqn::EulerData{Tsol, Tdim, Tres, :entropy},
                  opts, res_arr::AbstractArray{Tsol, 3})
# multiply a 3D array by inv(A0) in-place, useful for explicit time stepping
# res_arr *can* alias eqn.q safely
# a non-alias tolerant implimention wold avoid copying q_vals
  A0 = Array(Tsol, mesh.numDofPerNode, mesh.numDofPerNode)
  res_vals = Array(Tsol, mesh.numDofPerNode)
  q_vals = Array(Tsol, mesh.numDofPerNode)
  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      # copy values into workvec
      for k=1:mesh.numDofPerNode
	q_vals[k] = eqn.q[k, j, i]
	res_vals[k] = res_arr[k, j, i]
      end

      res_view = view(res_arr, :, j, i)
      # get A0Inv for this node
      calcA0(eqn.params, q_vals, A0)

      smallmatvec!(A0, res_vals, res_view)
    end
  end

  return nothing
end

#no-op
function matVecA0{Tmsh, Tsol, Tdim, Tres}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, eqn::EulerData{Tsol, Tdim, Tres, :conservative}, opts, res_arr::AbstractArray{Tsol, 3})

  return nothing
end




@doc """
### EulerEquationMod.calcA1

  This function calculates the A1 (ie. dF1/dv, where F1 is the first column of 
  theEuler flux and v are the entropy variables) for a node

  The formation of A1 is given in Hughes paper

  Inputs:
    params : ParamType{2, :entropy}
    q  : vector of entropy variables, length 4

  Inputs/Outputs:
  A1 : 4x4 matrix to be populated.  Overwritten


"""->
function calcA1{Tsol}(params::ParamType{2, :entropy}, q::AbstractArray{Tsol,1}, 
                      A1::AbstractArray{Tsol, 2})


  gamma = params.gamma
  gamma_1 = params.gamma_1

  k1 = 0.5*(q[2]^2 + q[3]^2)/q[4]  # a constant from Hughes' paper
  k2 = k1 - gamma
  k3 = k1*k1 - 2*gamma*k1 + gamma
  k4 = k2 - gamma_1
  k5 = k2*k2 - gamma_1*(k1 + k2)

  s = gamma - q[1] + k1    # entropy

  rho_int = exp(-s/gamma_1)*(gamma_1/((-q[4])^gamma))^(1/gamma_1)

  fac = rho_int/(gamma_1*q[4]*q[4])

  # calculate the variables used in Hughes A.1
  c1 = gamma_1*q[4] - q[2]*q[2]
  c2 = gamma_1*q[4] - q[3]*q[3]

  d1 = -q[2]*q[3]

  e1 = q[2]*q[4]
#  e2 = q[3]*q[4]

  println("fac = ", fac)
  println("e1 = ", e1)
  println("q[4] = ", q[4])
  # populate the matrix
  # the matrix is symmetric, but we don't use it because I think populating
  # the matrix will be faster if the matrix is write-only
  A1[1,1] = e1*q[4]*fac
  A1[2,1] = c1*q[4]*fac
  A1[3,1] = d1*q[4]*fac
  A1[4,1] = k2*e1*fac
  A1[1,2] = c1*q[4]*fac  # symmetric
  A1[2,2] = -(c1 + 2*gamma_1*q[4])*q[2]*fac
  A1[3,2] = -c1*q[3]*fac
  A1[4,2] = (c1*k2 + gamma_1*q[2]*q[2])*fac
  A1[1,3] = d1*q[4]*fac  # symmetric
  A1[2,3] = -c1*q[3]*fac  # symmetric
  A1[3,3] = -c2*q[2]*fac
  A1[4,3] = k4*d1*fac
  A1[1,4] = k2*e1*fac  # symmetric
  A1[2,4] = (c1*k2 + gamma_1*q[2]*q[2])*fac   # symmetric
  A1[3,4] = k4*d1*fac  # symmetric
  A1[4,4] = k5*q[2]*fac

  println("A1[1,1] = ", A1[1,1])
    return nothing
end


@doc """
### EulerEquationMod.calcA2

  This function calculates the A2 (ie. dF2/dv, where F2 is the second column of the
  Euler flux.
  and v are the entropy variables) for a node

  The formation of A2 is given in Hughes

  Inputs:
    params : ParamType{2, :entropy},
    q  : vector of entropy variables, length 4
  Inputs/Outputs:
  A2 : 4x4 matrix to be populated.  Overwritten


"""->
function calcA2{Tsol}(params::ParamType{2, :entropy}, q::AbstractArray{Tsol,1}, 
                      A2::AbstractArray{Tsol, 2})


  gamma = params.gamma
  gamma_1 = params.gamma_1

  k1 = 0.5*(q[2]^2 + q[3]^2)/q[4]  # a constant from Hughes' paper
  k2 = k1 - gamma
  k3 = k1*k1 - 2*gamma*k1 + gamma
  k4 = k2 - gamma_1
  k5 = k2*k2 - gamma_1*(k1 + k2)

  s = gamma - q[1] + k1    # entropy

  rho_int = exp(-s/gamma_1)*(gamma_1/((-q[4])^gamma))^(1/gamma_1)

  fac = rho_int/(gamma_1*q[4]*q[4])

  # calculate the variables used in Hughes A.1
  c1 = gamma_1*q[4] - q[2]*q[2]
  c2 = gamma_1*q[4] - q[3]*q[3]

  d1 = -q[2]*q[3]

#  e1 = q[2]*q[4]
  e2 = q[3]*q[4]

  # populate the matrix
  # the matrix is symmetric, but we don't use it because I think populating
  # the matrix will be faster if the matrix is write-only
  A2[1,1] = e2*q[4]*fac
  A2[2,1] = d1*q[4]*fac
  A2[3,1] = c2*q[4]*fac
  A2[4,1] = k2*e2*fac
  A2[1,2] = d1*q[4]*fac  # symmetric
  A2[2,2] = -c1*q[3]*fac
  A2[3,2] = -c2*q[2]*fac
  A2[4,2] = k4*d1*fac
  A2[1,3] = c2*q[4]*fac  # symmetric
  A2[2,3] = -c2*q[2]*fac  # symmetric
  A2[3,3] = -(c2 + 2*gamma_1*q[4])*q[3]*fac
  A2[4,3] = (c2*k2 + gamma_1*q[3]*q[3])*fac
  A2[1,4] = k2*e2*fac  # symmetric
  A2[2,4] = k4*d1*fac   # symmetric
  A2[3,4] = (c2*k2 + gamma_1*q[3]*q[3])*fac  # symmetric
  A2[4,4] = k5*q[3]*fac

    return nothing
end

@doc """
### EulerEquationMod.calcMaxWaveSpeed

  This function calculates the maximum wave speed (ie. acoustic wave speed)
  present in the domain and returns it.  Methods are available for conservative and entropy
  variables.

  This is a mid level function
"""->
function calcMaxWaveSpeed{Tsol, Tdim, Tres}(mesh, sbp, 
                          eqn::EulerData{Tsol, Tdim, Tres, :conservative}, opts)
# calculate the maximum wave speed (ie. characteristic speed) on the mesh
# uses solution vector q, not array
  q = eqn.q
  max_speed = zero(eltype(q))
  for i=1:4:length(q)
    q_i = view(q, i:(i+3))
    a = calcSpeedofSound(eqn.params, q_i)
    ux = q_i[2]/q_i[1]
    uy = q_i[3]/q_i[1]
    u_norm = sqrt(ux*ux + uy*uy)
    speed = a + u_norm

    if speed > max_speed
      max_speed = speed
    end  # end if statement
  end  # end loop over vector q

  return max_speed
end  # end function


function calcMaxWaveSpeed{Tsol, Tdim, Tres}(mesh, sbp, 
                          eqn::EulerData{Tsol, Tdim, Tres, :entropy}, opts)
# calculate the maximum wave speed (ie. characteristic speed) on the mesh
# uses solution vector q, not array
  q = eqn.q
  gamma = eqn.params.gamma
  gamma_1 = eqn.params.gamma_1
 
  max_speed = zero(eltype(q))
  for i=1:4:length(q)
    q_i = view(q, i:(i+3))
    a = calcSpeedofSound(eqn.params, q_i)

   
    k1 = 0.5*(q_i[2]*q_i[2] + q_i[3]*q_i[3])/q_i[4]  # a constant from Hughes' paper
    s = gamma - q_i[1] + k1    # entropy
    rho_int = exp(-s/gamma_1)*(gamma_1/((-q_i[4])^gamma))^(1/gamma_1)
    rho = -rho_int*q_i[4]
    ux = rho_int*q[2]/rho
    uy = rho_int*q[3]/rho
    u_norm = sqrt(ux*ux + uy*uy)
    speed = a + u_norm

    if speed > max_speed
      max_speed = speed
    end  # end if statement
  end  # end loop over vector q

  return max_speed
end  # end function
