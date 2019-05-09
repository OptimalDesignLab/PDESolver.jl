function getCNDSPCandLO(mesh, sbp, eqn, opts)
  
  jac_type = opts["jac_type"]

  if jac_type <= 2
    pc = PCNone(mesh, sbp, eqn, opts)
  elseif opts["use_volume_preconditioner"]
    pc = CNVolumePC(mesh, sbp, eqn, opts)
  else
    pc = CNMatPC(mesh, sbp, eqn, opts)
  end

  # if jac_type == 1    # Dense jac
    # lo = CNDSDenseLO(pc, mesh, sbp, eqn, otps)
  # elseif jac_type == 2    # Sparse jac, julia matrix
    # lo = CNDSSparseDirectLO(pc, mesh, sbp, eqn, opts)
  # elseif jac_type == 3
    # lo = CNDSPetscMatLO(pc, mesh, sbp, eqn, opts)
  # elseif jac_type == 4
    # lo = CNDSPetscMatFreeLO(pc, mesh, sbp, eqn, opts)
  # end
  if jac_type == 1    # Dense jac
    lo = CNDenseLO(pc, mesh, sbp, eqn, otps)
  elseif jac_type == 2    # Sparse jac, julia matrix
    lo = CNSparseDirectLO(pc, mesh, sbp, eqn, opts)
  elseif jac_type == 3
    lo = CNPetscMatLO(pc, mesh, sbp, eqn, opts)
  elseif jac_type == 4
    lo = CNPetscMatFreeLO(pc, mesh, sbp, eqn, opts)
  end



  return pc, lo
end

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
function modifyCNJacForMatFreeCheck(lo::CNHasMat, mesh, sbp, eqn, opts, ctx_residual, t)

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
function modifyCNJacForMatFreeCheck_reverse(lo::CNHasMat, mesh, sbp, eqn, opts, ctx_residual, t)

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
  new_res_vec_Maimag::Array{Complex{Float64},1}  # storing the last time step's new_res_vec_Maimag (needed in dRdM calcs)
end

""" 
"""
function CNDSCheckpointData(chkpointer::Checkpointer, comm_rank::Integer)

  chkpoint_data = readLastCheckpointData(chkpointer, comm_rank)

  return chkpoint_data::CNDSCheckpointData
end

function CNDS_checkpoint_setup(mesh, opts, myrank)
  is_restart = opts["is_restart"]
  ncheckpoints = opts["ncheckpoints"]

  if ! is_restart
    # this is a new simulation, create all the stuff needed to checkpoint
    # note that having zero checkpoints is valid
    istart = 2
    i_test = istart*10
    numDof = mesh.numDof
    v_vec = zeros(Float64, numDof)
    new_res_vec_Maimag = zeros(Complex{Float64}, numDof)
    chkpointdata = CNDSCheckpointData(istart, i_test, numDof, v_vec, new_res_vec_Maimag)
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



# #------------------------------------------------------------------------------
# # Linear operators

# #TODO: make a macro to generate these definitions

# mutable struct CNDSDenseLO <: AbstractDenseLO
  # lo_inner::NewtonDenseLO
# end

# function CNDSDenseLO(pc::PCNone, mesh::AbstractMesh,
                    # sbp::AbstractSBP, eqn::AbstractSolutionData, opts::Dict)

  # lo_inner = NewtonDenseLO(pc, mesh, sbp, eqn, opts)

  # return CNDSDenseLO(lo_inner)
# end

# mutable struct CNDSSparseDirectLO <: AbstractSparseDirectLO
  # lo_inner::NewtonSparseDirectLO
# end

# function CNDSSparseDirectLO(pc::PCNone, mesh::AbstractMesh,
                    # sbp::AbstractSBP, eqn::AbstractSolutionData, opts::Dict)

  # lo_inner = NewtonSparseDirectLO(pc, mesh, sbp, eqn, opts)

  # return CNDSSparseDirectLO(lo_inner)
# end


# mutable struct CNDSPetscMatLO <: AbstractPetscMatLO
  # lo_inner::NewtonPetscMatLO
# end

# function CNDSPetscMatLO(pc::AbstractPetscPC, mesh::AbstractMesh,
                    # sbp::AbstractSBP, eqn::AbstractSolutionData, opts::Dict)

  # lo_inner = NewtonPetscMatLO(pc, mesh, sbp, eqn, opts)

  # return CNDSPetscMatLO(lo_inner)
# end

# """
  # All CN LOs that have matrices
# """
# const CNDSMatLO = Union{CNDSDenseLO, CNDSSparseDirectLO, CNDSPetscMatLO}

# function calcLinearOperator(lo::CNDSMatLO, mesh::AbstractMesh,
                            # sbp::AbstractSBP, eqn::AbstractSolutionData,
                            # opts::Dict, ctx_residual, t)

  # calcLinearOperator(lo.lo_inner, mesh, sbp, eqn, opts, ctx_residual, t)

  # modifyJacCN(lo, mesh, sbp, eqn, opts, ctx_residual, t)

  # return nothing
# end

# # applyLinearOperator & applyLinearOperatorTranspose --- needed?

