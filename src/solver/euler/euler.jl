# this file contains the functions to evaluate the right hand side of the weak form in the pdf

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
function evalEuler(mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerData, opts, 
                   t=0.0)

  eqn.params.t = t  # record t to params
  
  dataPrep(mesh, sbp, eqn, opts)
#  println("dataPrep @time printed above")


  evalVolumeIntegrals(mesh, sbp, eqn, opts)
#  println("volume integral @time printed above")

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
#  println("boundary integral @time printed above")


  addStabilization(mesh, sbp, eqn, opts)
#  println("stabilizing @time printed above")


  
#  print("\n")

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
function init{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, 
              eqn::AbstractEulerData{Tsol, Tres}, opts, pmesh=mesh)

  println("\nInitializing Euler module")
  # get BC functors
  getBCFunctors(mesh, sbp, eqn, opts)
  getBCFunctors(pmesh, sbp, eqn, opts)


  return nothing
end


function majorIterationCallback(itr::Integer, mesh::AbstractMesh, 
                                sbp::SBPOperator, eqn::AbstractEulerData, opts)

#  println("Performing major Iteration Callback")


    if opts["write_vis"] && ((itr % opts["output_freq"])) == 0 || itr == 1
      vals = abs(real(eqn.q_vec))  # remove unneded imaginary part
      saveSolutionToMesh(mesh, vals)
      fname = string("solution_", itr)
      writeVisFiles(mesh, fname)
    end
 

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
function dataPrep{Tmsh,  Tsol, Tres}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator,
                                     eqn::AbstractEulerData{Tsol, Tres}, opts)
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
function evalVolumeIntegrals{Tmsh,  Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, 
                             sbp::SBPOperator, eqn::EulerData{Tsol, Tres, Tdim}, opts)
  
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



end  # end evalVolumeIntegrals


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
function evalBoundaryIntegrals{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, 
                               sbp::SBPOperator, eqn::EulerData{Tsol, Tres, Tdim})

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
function addStabilization{Tmsh,  Tsol}(mesh::AbstractMesh{Tmsh}, 
                          sbp::SBPOperator, eqn::EulerData{Tsol}, opts)

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
    applyDissipation(mesh, sbp, eqn, opts, eqn.q)
  end

  if opts["use_GLS2"]
     applyGLS2(mesh, sbp, eqn, opts)
  end

#  println("==== end of addStabilization ====")



  return nothing

end


#------------------------------------------------------------------------------
# other functions
#------------------------------------------------------------------------------

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
function disassembleSolution{Tmsh, Tsol, Tres, Tdim, T}(mesh::AbstractMesh{Tmsh}, sbp,
                             eqn::EulerData{Tsol, Tres, Tdim}, opts, 
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

# mid level function (although it doesn't need Tdim)
function assembleArray{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh}, 
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
      for k=size(arr, 1)  # loop over dofs on the node

        dofnum_k = mesh.dofs[k, j, i]
        dofnum_k1 = div(dofnum_k, mesh.numDofPerNode) + 1 # get node number

        res_vec[dofnum_k1] = arr[k,j,i]
      end
    end
  end
  
  return nothing
end




