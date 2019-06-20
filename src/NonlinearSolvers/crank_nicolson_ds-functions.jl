#------------------------------------------------------------------------------------------------------------
# important functions for crank_nicolson_ds
# 
# Sections:
#   1. calcLinearOperator & modifyJacCN for CNDS LO's
#   2. modifyCNJacForMatFreeCheck & modifyCNJacForMatFreeCheck_reverse, used for mat-free jac-vec prod verification
#   3. checkpointing setup for CNDS


"""
  All CN LOs that have matrices
"""
const CNDSMatLO = Union{CNDSDenseLO, CNDSSparseDirectLO, CNDSPetscMatLO}
const CNDSHasMat = Union{CNDSDenseLO, CNDSSparseDirectLO, CNDSPetscMatLO}


"""
  calcLinearOperator function

  Mostly a copy of calcLinearOperator(lo::CNMatLO...)
  defined in crank_nicolson.jl.

  Only difference is debugging output and the inclusion of the CNDS stabilization.
"""
function calcLinearOperator(lo::CNDSMatLO, mesh::AbstractMesh,
                            sbp::AbstractSBP, eqn::AbstractSolutionData,
                            opts::Dict, ctx_residual, t)

  # println(BSTDOUT, "    entered cLO(lo::CNDSMatLO...) in crank_nicolson_ds-functions.jl")

  stabilize_this_LOupdate = ctx_residual[9]
  # println(BSTDOUT, " stabilize_this_LOupdate: ", stabilize_this_LOupdate)

  # println(BSTDOUT, "     typeof(lo): ", typeof(lo))
  # println(BSTDOUT, "     typeof(lo.lo_inner): ", typeof(lo.lo_inner))
  # println(BSTDOUT, "     calling inner cLO in cLO(lo::CNDSMatLO...) in crank_nicolson_ds-functions.jl")
  calcLinearOperator(lo.lo_inner, mesh, sbp, eqn, opts, ctx_residual, t)

  ########################################################################################
  # Here is where we stabilize
  ########################################################################################
  if opts["stabilize_v"] && stabilize_this_LOupdate == true
    stabilizeCNDSLO(lo, mesh, sbp, eqn, opts, ctx_residual, t)
  end

  # println(BSTDOUT, "     calling modifyJacCN from cLO()")

  lo_innermost = getBaseLO(lo)
  # writedlm("lo_innermost_A-before_modifyJacCN.dat", lo_innermost.A)   # can't do this with Petsc matrices, will hang
  modifyJacCN(lo, mesh, sbp, eqn, opts, ctx_residual, t)
  # writedlm("lo_innermost_A-after_modifyJacCN.dat", lo_innermost.A)    # can't do this with Petsc matrices, will hang

  # println(BSTDOUT, "    leaving cLO(lo::CNDSMatLO...) in crank_nicolson_ds-functions.jl")

  return nothing
end

"""
  Takes the Jacobian of the physics and modifies it to be the Crank-Nicolson
  Jacobian.

  CN_Jac = I - dt/2 * physicsJac


  **Inputs**

   * lo: a CNHasMat
   * mesh
   * sbp
   * eqn
   * opts
   * ctx_residual: delta_t must be the 3rd element
   * t
"""
function modifyJacCN(lo::CNDSHasMat, mesh, sbp, eqn, opts, ctx_residual, t)


  # println(BSTDOUT, "      modifyJacCN(lo::CNDSHasMat) called")

  lo2 = getBaseLO(lo)
  h = ctx_residual[3]

  assembly_begin(lo2.A, MAT_FINAL_ASSEMBLY)
  assembly_end(lo2.A, MAT_FINAL_ASSEMBLY)

  # scale jac by -delta_t/2
#  scale_factor = h*-0.5
  petsc_scale_factor = PetscScalar(-h*0.5)
  scale!(lo2.A, petsc_scale_factor)

  # add the identity
  diagonal_shift!(lo2.A, 1)

  return nothing
end


#------------------------------------------------------------------------------
# other functions: modifyCNJac functions for mat-free check

"""
  Takes a Crank-Nicolson Jacobian and modifies it to be of the form needed to 
  verify a mat-free Jacobian-vector product.

  Input:
    CN_Jac = I - dt/2 * physicsJac
  Output:
    CN_Jac_mod = - I - dt/2 * physicsJac


  **Inputs**

   * lo: a CNHasMat
   * mesh
   * sbp
   * eqn
   * opts
   * ctx_residual: delta_t must be the 3rd element
   * t
"""
function modifyCNJacForMatFreeCheck(lo::CNDSHasMat, mesh, sbp, eqn, opts, ctx_residual, t)

  println(BSTDOUT, " modifyCNJacForMatFreeCheck(lo::CNHasMat...) called")

  lo2 = getBaseLO(lo)
  h = ctx_residual[3]

  assembly_begin(lo2.A, MAT_FINAL_ASSEMBLY)
  assembly_end(lo2.A, MAT_FINAL_ASSEMBLY)

  # add a -I twice
  diagonal_shift!(lo2.A, -1)
  diagonal_shift!(lo2.A, -1)

  return nothing

end   # end function modifyCNJacForMatFreeCheck

"""
  Takes a CN Jacobian of type needed to verify a mat-free Jacobian-vector product,
    and reverses it so that it is a standard CN Jacobian afterwards.

  Input:
    CN_Jac_mod = - I - dt/2 * physicsJac
  Output:
    CN_Jac = I - dt/2 * physicsJac

  **Inputs**

   * lo: a CNHasMat
   * mesh
   * sbp
   * eqn
   * opts
   * ctx_residual: delta_t must be the 3rd element
   * t
"""
function modifyCNJacForMatFreeCheck_reverse(lo::CNDSHasMat, mesh, sbp, eqn, opts, ctx_residual, t)

  println(BSTDOUT, " modifyCNJacForMatFreeCheck_reverse(lo::CNHasMat...) called")

  lo2 = getBaseLO(lo)
  h = ctx_residual[3]

  assembly_begin(lo2.A, MAT_FINAL_ASSEMBLY)
  assembly_end(lo2.A, MAT_FINAL_ASSEMBLY)

  # add a I twice
  diagonal_shift!(lo2.A, 1)
  diagonal_shift!(lo2.A, 1)

  return nothing

end   # end function modifyCNJacForMatFreeCheck_reverse

#------------------------------------------------------------------------------
# Auxiliary functions
#------------------------------------------------------------------------------

"""
  Converts an array of size numEl to one of size numDofs.
  It does this by simply assigning the value corresponding
  to each element to all the dofs in the element.

  Inputs:
    mesh: usual
    array: array of length mesh.numEl to populate a mesh.numDof length array

  Input/Output:
    vec: array of length mesh.numDof that will be overwritten with the spread-out
         data from 'array'
"""
function convertArrayOfElsToDofs(mesh, array::Array, vec::Array)

  @assert length(array) == mesh.numEl
  @assert length(vec) == mesh.numDof

  dofsPerEl = div(mesh.numDof, mesh.numEl)    # == numDofPerNode*numNodesPerEl

  ctr = 1

  for el_ix = 1:length(array)
    for dof_ix = 1:dofsPerEl
      vec[ctr] = array[el_ix]
      ctr += 1
    end
  end

  return nothing
end

#------------------------------------------------------------------------------
# Checkpointing functions
#------------------------------------------------------------------------------

"""
  This type stores all the data needed for:
    1) CN to restart (time index), and
    2) the direct sensitivity evolution to continue

"""
mutable struct CNDSCheckpointData <: AbstractCheckpointData
  i::Int    # current time step
  i_test::Int
  numDof::Int   # number of DOFs in the mesh, needed for v_vec sizing
  v_vec::Array{Float64,1}   # storing the direct sensitivity
  drag_array::Array{Float64, 1}
  term23::Float64
end

""" 
"""
function CNDSCheckpointData(chkpointer::Checkpointer, comm_rank::Integer)

  chkpoint_data = readLastCheckpointData(chkpointer, comm_rank)

  return chkpoint_data::CNDSCheckpointData
end

function CNDS_checkpoint_setup(mesh, opts, myrank, finaliter)
  is_restart = opts["is_restart"]
  ncheckpoints = opts["ncheckpoints"]

  if ! is_restart
    # this is a new simulation, create all the stuff needed to checkpoint
    # note that having zero checkpoints is valid
    istart = 2
    i_test = istart*10
    numDof = mesh.numDof
    v_vec = zeros(Float64, numDof)
    drag_array = zeros(Float64, finaliter)
    term23 = 0.0
    chkpointdata = CNDSCheckpointData(istart, i_test, numDof, v_vec, drag_array, term23)
    chkpointer = Checkpointer(myrank, ncheckpoints)
    skip_checkpoint = false
  else
    # this is a restart, load existing data
    # using default value of 0 checkpoints is ok
    chkpointer = Checkpointer(opts, myrank)
    chkpointdata = CNDSCheckpointData(chkpointer, myrank)
    skip_checkpoint = true  # when restarting, don't immediately write a checkpoint
  end

  return chkpointer, chkpointdata, skip_checkpoint
end

