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

The high level functions take the mesh::AbstractMesh, sbp::AbstractSBP, eqn::EulerEquation, and opts (options dictionary), They do not know the types their 
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

# Rules for paramaterization:
# Tmsh = mesh data type
# Tsbp = SBP operator data type
# Tsol = equation, res_vec, q_vec data type

import PDESolver.evalResidual

@doc """
### EulerEquationMod.evalResidual

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
function evalResidual(mesh::AbstractMesh, sbp::AbstractSBP, eqn::EulerData, 
                     opts::Dict, t=0.0)

#  println("----- entered evalResidual -----")

  time = eqn.params.time
  eqn.params.t = t  # record t to params
  myrank = mesh.myrank

#  println("entered evalResidual")
#  println("q1319-3 = ", eqn.q[:, 3, 1319])
  time.t_send += @elapsed if opts["parallel_type"] == 1
    startDataExchange(mesh, opts, eqn.q,  eqn.q_face_send, eqn.q_face_recv, eqn.params.f)
    #  println(" startDataExchange @time printed above")
  end
 

  time.t_dataprep += @elapsed dataPrep(mesh, sbp, eqn, opts)
#  println("dataPrep @time printed above")


  time.t_volume += @elapsed if opts["addVolumeIntegrals"]
    evalVolumeIntegrals(mesh, sbp, eqn, opts)
#    println("volume integral @time printed above")
  end

  # delete this if unneeded or put it in a function.  It doesn't belong here,
  # in a high level function.
  #----------------------------------------------------------------------------
  #=
  bndryfluxPhysical = zeros(eqn.bndryflux)
  getPhysBCFluxes(mesh, sbp, eqn, opts, bndryfluxPhysical)
  #println("bndryfluxPhysical = \n", bndryfluxPhysical)
  #println("eqn.bndryflux = \n", eqn.bndryflux)
  bndryfluxPhysical = -1*bndryfluxPhysical
  boundaryintegrate!(mesh.sbpface, mesh.bndryfaces, bndryfluxPhysical, eqn.res, SummationByParts.Subtract())
  =#

  if opts["use_GLS"]
    GLS(mesh,sbp,eqn)
  end
  
  #=
  bndryfluxPhysical = -1*bndryfluxPhysical
  boundaryintegrate!(mesh.sbpface, mesh.bndryfaces, bndryfluxPhysical, eqn.res, SummationByParts.Subtract())
  =#
  #----------------------------------------------------------------------------

  time.t_bndry += @elapsed if opts["addBoundaryIntegrals"]
   evalBoundaryIntegrals(mesh, sbp, eqn)
   #println("boundary integral @time printed above")
  end

  time.t_stab += @elapsed if opts["addStabilization"]
    addStabilization(mesh, sbp, eqn, opts)
#    println("stabilizing @time printed above")
  end

  time.t_face += @elapsed if mesh.isDG && opts["addFaceIntegrals"]
    evalFaceIntegrals(mesh, sbp, eqn, opts)
#    println("face integral @time printed above")
  end

  time.t_sharedface += @elapsed if mesh.commsize > 1
    evalSharedFaceIntegrals(mesh, sbp, eqn, opts)
    #println("evalSharedFaceIntegrals @time printed above")
  end

  time.t_source += @elapsed evalSourceTerm(mesh, sbp, eqn, opts)
  #println("source integral @time printed above")

  # apply inverse mass matrix to eqn.res, necessary for CN
  if opts["use_Minv"]
    applyMassMatrixInverse3D(mesh, sbp, eqn, opts, eqn.res)
  end

  return nothing
end  # end evalResidual



@doc """
### EulerEquationMod.init

  This function performs any operations that need to be performed before
  the first residual evaluation.
  Any operations that need to be performed at the beginning of *each* 
  residual evaluation should go in dataPrep
"""
# high level functions
function init{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP, 
              eqn::AbstractEulerData{Tsol, Tres}, opts, pmesh=mesh)

#  println("\nInitializing Euler module")

  initMPIStructures(mesh, opts)
  # get BC functors
  getBCFunctors(mesh, sbp, eqn, opts)
  getBCFunctors(pmesh, sbp, eqn, opts)

  getSRCFunctors(mesh, sbp, eqn, opts)
  if mesh.isDG
    getFluxFunctors(mesh, sbp, eqn, opts)
    getFaceElementFunctors(mesh, sbp, eqn, opts)
  end



  return nothing
end


function majorIterationCallback{Tmsh, Tsol, Tres, Tdim}(itr::Integer, 
                               mesh::AbstractMesh{Tmsh}, 
                               sbp::AbstractSBP, 
                               eqn::EulerData{Tsol, Tres, Tdim}, opts, f::IO)

#  println("Performing major Iteration Callback")

  myrank = mesh.myrank
#  println("eqn.q = \n", eqn.q)
  # undo multiplication by inverse mass matrix
  res_vec_orig = eqn.M.*copy(eqn.res_vec)
  res_orig = reshape(res_vec_orig, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)

    if opts["write_vis"] && (((itr % opts["output_freq"])) == 0 || itr == 1)
      vals = real(eqn.q_vec)  # remove unneded imaginary part

#      println("writing vtk, q1319-3 = ", eqn.q[:, 3, 1319])
#      dofs = mesh.dofs[:, 3, 1319]
#      println("writing vtk, q_vec1319-3 = ", vals[dofs])
      saveSolutionToMesh(mesh, vals)
      fname = string("solution_", itr)
      writeVisFiles(mesh, fname)
#=
      # DEBUGGING: write error to file
      q_exact = zeros(eqn.q_vec)
      ex_func = ICDict[opts["IC_name"]]
      ex_func(mesh, sbp, eqn, opts, q_exact)
      q_err = real(eqn.q_vec) - q_exact
      saveSolutionToMesh(mesh, q_err)
      fname = string("error_", itr)
      writeVisFiles(mesh, fname)
=#
    end
 
    # add an option on control this or something.  Large blocks of commented
    # out code are bad
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
  if opts["write_entropy"] && (itr % opts["write_entropy_freq"] == 0)
    # calculate the entropy norm
    val = calcEntropyIntegral(mesh, sbp, eqn, opts, eqn.q_vec)

    # compute w^T * res_vec
    val2 = contractResEntropyVars(mesh, sbp, eqn, opts, eqn.q_vec, res_vec_orig)

    # DEBUGGING: compute the potential flux from q
    #            directly, to verify the boundary terms are the problem
#    val3 = calcInterfacePotentialFlux(mesh, sbp, eqn, opts, eqn.q)
#    val3 += calcVolumePotentialFlux(mesh, sbp, eqn, opts, eqn.q)

    @mpi_master begin
      f = open(opts["write_entropy_fname"], "a+")
      println(f, itr, " ", eqn.params.t, " ",  val, " ", val2)
      close(f)
    end
  end

  if opts["write_integralq"]
    integralq_vals = integrateQ(mesh, sbp, eqn, opts, eqn.q_vec)
    f = open(opts["write_integralq_fname"], "a+")
    print(f, itr, " ", eqn.params.t)
    for i=1:length(integralq_vals)
      print(f, " ", integralq_vals[i])
    end
    print(f, "\n")
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
function dataPrep{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP,
                                     eqn::AbstractEulerData{Tsol, Tres}, opts)
# gather up all the data needed to do vectorized operatinos on the mesh
# calculates all mesh wide quantities in eqn
# this is almost the exact list of everything we *shouldn't* be storing, but
# rather recalculating on the fly

#println("Entered dataPrep()")

#  println("typeof(eqn) = ", typeof(eqn))
#  println("typeof(eqn.params) = ", typeof(eqn.params))

  # apply filtering to input
  if eqn.params.use_filter
    applyFilter(mesh, sbp, eqn, eqn.q, opts)
  end

  # zero out res
  fill!(eqn.res, 0.0)
  fill!(eqn.res_edge, 0.0)
 
  getAuxVars(mesh, eqn)
#  println("  getAuxVars @time printed above")
  
  if opts["check_density"]
    checkDensity(eqn)
#    println("  checkDensity @time printed above")
  end

  if opts["check_pressure"]
#    throw(ErrorException("I'm done"))
    checkPressure(eqn)
#    println("  checkPressure @time printed above")
  end

  # calculate fluxes
#  getEulerFlux(eqn, eqn.q, mesh.dxidx, sview(flux_parametric, :, :, :, 1), sview(flux_parametric, :, :, :, 2))
  getEulerFlux(mesh, sbp,  eqn, opts)
#  println("  getEulerFlux @time printed above")


  if mesh.isDG
    fill!(eqn.q_bndry, 0.0)
    fill!(eqn.q_face, 0.0)
    fill!(eqn.flux_face, 0.0)
    interpolateBoundary(mesh, sbp, eqn, opts, eqn.q, eqn.q_bndry)

    if opts["face_integral_type"] == 1
      interpolateFace(mesh, sbp, eqn, opts, eqn.q, eqn.q_face)
      calcFaceFlux(mesh, sbp, eqn, eqn.flux_func, mesh.interfaces, eqn.flux_face)
    end
  end
#  println("  DG dataPrep @time printed above")
  fill!(eqn.bndryflux, 0.0)
  getBCFluxes(mesh, sbp, eqn, opts)
#   println("  getBCFluxes @time printed above")
 
  # is this needed for anything besides edge stabilization?
  if eqn.params.use_edgestab 
    stabscale(mesh, sbp, eqn)
  end
#  println("  stabscale @time printed above")

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
    convertToConservative(eqn.params, sview(eqn.q, :, j, i), q_cons)
    if real(q_cons[1]) < 0.0
      println("q_conservative = ", q_cons)
    end
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
    q = sview(eqn.q, :, j, i)
    aux_vars = sview(eqn.aux_vars,:, j, i)
    press = @getPressure(aux_vars)
    @assert( real(press) > 0.0, "element $i, node $j, q = $q, press = $press")
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
                             sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim}, opts)

  integral_type = opts["volume_integral_type"]
  if integral_type == 1
    if opts["Q_transpose"] == true
      for i=1:Tdim
        weakdifferentiate!(sbp, i, sview(eqn.flux_parametric, :, :, :, i), eqn.res, trans=true)
      end
    else
      for i=1:Tdim
        weakdifferentiate!(sbp, i, sview(eqn.flux_parametric, :, :, :, i), eqn.res, SummationByParts.Subtract(), trans=false)
      end
    end  # end if
  elseif integral_type == 2
    calcVolumeIntegralsSplitForm(mesh, sbp, eqn, opts, eqn.volume_flux_func)
  else
    throw(ErrorException("Unsupported volume integral type = $integral_type"))
  end

  # artificialViscosity(mesh, sbp, eqn) 

end  # end evalVolumeIntegrals


@doc """
  This function evaluates the boundary integrals in the Euler equations by 
  calling the appropriate SBP function on eqn.bndryflux, which must be populated
  before calling this function.  eqn.res is updated with the result

  This is a mid level function

"""->
# mid level function
function evalBoundaryIntegrals{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, 
                               sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim})

  #TODO: remove conditional
  if mesh.isDG
    boundaryintegrate!(mesh.sbpface, mesh.bndryfaces, eqn.bndryflux, eqn.res, SummationByParts.Subtract())
  else
    boundaryintegrate!(mesh.sbpface, mesh.bndryfaces, eqn.bndryflux, eqn.res, SummationByParts.Subtract())
  end


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
                          sbp::AbstractSBP, eqn::EulerData{Tsol}, opts)

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
     applyGLS3(mesh, sbp, eqn, opts)
#    test_GLS(mesh, sbp, eqn, opts)
  end

#  println("==== end of addStabilization ====")



  return nothing

end

@doc """
### EulerEquationMod.evalFaceIntegrals

  This function evaluates the face integrals in a DG formulation and
  updates the residual.  The array eqn.flux_face must already be populated
  with the face flux.

  Inputs:
    mesh: an AbstractDGMesh
    sbp: an SBP operator
    eqn: an EulerData object
    opts: the options dictonary

  Outputs:
    none

"""->
function evalFaceIntegrals{Tmsh, Tsol}(mesh::AbstractDGMesh{Tmsh}, 
                           sbp::AbstractSBP, eqn::EulerData{Tsol}, opts)

  face_integral_type = opts["face_integral_type"]
  if face_integral_type == 1
#    println("calculating regular face integrals")
    interiorfaceintegrate!(mesh.sbpface, mesh.interfaces, eqn.flux_face, eqn.res, SummationByParts.Subtract())

  elseif face_integral_type == 2
#    println("calculating ESS face integrals")
   
    getFaceElementIntegral(mesh, sbp, eqn, eqn.face_element_integral_func,  
                           eqn.flux_func, mesh.sbpface, mesh.interfaces)

  else
    throw(ErrorException("Unsupported face integral type = $face_integral_type"))
  end

  # do some output here?
  return nothing

end
                            
#=
@doc """
### EulerEquationMod.sendParallelData

  This function interpolates the data into the send buffer and post
  the Isends and Irecvs.  It does not wait for them to finish

  Inputs:
    mesh
    sbp
    eqn
    opts
"""->
function sendParallelData(mesh::AbstractDGMesh, sbp, eqn, opts)

  for i=1:mesh.npeers
    # interpolate
    mesh.send_waited[i] = getSendData(mesh, opts, eqn.q, mesh.bndries_local[i], eqn.q_face_send[i], mesh.send_reqs[i], mesh.send_waited[i])
  end

  exchangeFaceData(mesh, opts, eqn.q_face_send, eqn.q_face_recv)

  return nothing
end
=#
@doc """
### EulerEquationMod.evalSharedFaceIntegrals

  This function does the computation that needs the parallel
  communication to have finished already, namely the face integrals
  for the shared faces

  Inputs:
    mesh
    sbp
    eqn
    opts
"""->
function evalSharedFaceIntegrals(mesh::AbstractDGMesh, sbp, eqn, opts)

#  println(eqn.params.f, "evaluating shared face integrals")
  face_integral_type = opts["face_integral_type"]
  if face_integral_type == 1

    if opts["parallel_data"] == "face"
      calcSharedFaceIntegrals(mesh, sbp, eqn, opts, eqn.flux_func)
    elseif opts["parallel_data"] == "element"
      calcSharedFaceIntegrals_element(mesh, sbp, eqn, opts, eqn.flux_func)
    else
      throw(ErrorException("unsupported parallel data type"))
    end

  elseif face_integral_type == 2

    getSharedFaceElementIntegrals_element(mesh, sbp, eqn, opts, eqn.face_element_integral_func,  eqn.flux_func)
  else
    throw(ErrorException("unsupported face integral type = $face_integral_type"))
  end

  return nothing
end

@doc """
### EulerEquationMod.evalSourceTerm

  This function performs all the actions necessary to update eqn.res
  with the source term.  The source term is stored in eqn.src_func.  It is
  an abstract field, so it cannot be accessed (performantly) directly, so
  it is passed to an inner function.

  Inputs:
    mesh : Abstract mesh type
    sbp  : Summation-by-parts operator
    eqn  : Euler equation object
    opts : options dictonary

  Outputs: none

  Aliasing restrictions: none

"""->
function evalSourceTerm{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh},
                     sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim}, 
                     opts)


  # placeholder for multiple source term functionality (similar to how
  # boundary conditions are done)
  if opts["use_src_term"]
    applySourceTerm(mesh, sbp, eqn, opts, eqn.src_func)
  end

  return nothing
end  # end function

@doc """
### EulerEquationMod.applyMassMatrixInverse3D

  This function applies the 3D inverse mass matrix to an array. 
    The array passed in should always be eqn.res

  Inputs:
    mesh: mesh object, needed for numEl and numDofPerNode fields
    sbp: sbp object, needed for numnodes field
    eqn: equation object, needed for Minv3D field
    opts
    arr: the 3D array to have the 3D mass matrix inverse applied to it

"""->
function applyMassMatrixInverse3D(mesh, sbp, eqn, opts, arr)

  for i = 1:mesh.numEl
    for j = 1:sbp.numnodes
      for k = 1:mesh.numDofPerNode
        arr[k, j, i] = eqn.Minv3D[k, j, i] * arr[k, j, i]
      end
    end
  end

  return arr
end

