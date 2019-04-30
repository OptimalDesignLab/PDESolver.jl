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

