
# Note: passed in for the eqn argument is 'eqn_nextstep'
function stabilizeCNDSLO(lo, mesh, sbp, eqn, opts, ctx_residual, t)

  # ctx_residual = (f, eqn, h, newton_data, stab_A, stab_assembler, clipJacData)
  f =               ctx_residual[1]
  eqn_old =         ctx_residual[2]
  h =               ctx_residual[3]
  newton_data =     ctx_residual[4]
  stab_A =          ctx_residual[5]
  stab_assembler =  ctx_residual[6]
  clipJacData =     ctx_residual[7]

  println(BSTDOUT, "        stabilizeCNDSLO called")

  evalJacobianStrong(mesh, sbp, eqn, opts, stab_assembler, t)

  # filterDiagJac
  #   location: jacobian_diag.jl
  #
  #   The third argument is q_vec.
  #   It is used as part of computing the quadprog stabilization (findStablePerturbation!).
  #   It is not used when opts["stabilization_method"] is "clipJac" or "clipJacFast".
  #   For the explicit stabilization what was passed in was 'real(tmp_imag)'.
  #   TODO: figure out from derivation what should be passed in for this arg.
  filterDiagJac(mesh, opts, eqn.q_vec, clipJacData, stab_A, eigs_to_remove="neg")

  # TODO: here is where we then add the stabilizer to the LO

  return numberOfEigchgs

end

#=
        if opts["stabilize_v"]

          # NOTE: we are now stabilizing inside calcLinearOperator
          error("no evalJacobianStrong call. remove this error only when derivation vs implementation is complete.")

          # Recalculate dRdq
          # filterDiagJac(mesh, opts, real(tmp_imag), clipJacData, stab_A, eigs_to_remove="neg")
          # filterDiagJac(mesh, opts, real(tmp_imag), clipJacData, stab_A, eigs_to_remove="pos")

          # loop over blocks
          # blocksize is set above (during DiagJac init) as mesh.numDofPerNode*mesh.numNodesPerElement
          nblocks = size(stab_A.A, 3)       # third dimension of our block diag Jac is the block index
          ix_petsc_row = zeros(PetscInt, blocksize)
          ix_petsc_col = zeros(PetscInt, blocksize)
          block_to_add = zeros(PetscScalar, blocksize, blocksize)
          for block_ix = 1:nblocks

            # TODO: no offsets present: this may not be correct in parallel
            for row_ix = 1:length(ix_petsc_row)
              # set the row indicies that we will insert into
              ix_petsc_row[row_ix] = blocksize*(block_ix-1)+row_ix
            end
            for col_ix = 1:length(ix_petsc_col)
              ix_petsc_col[col_ix] = blocksize*(block_ix-1)+col_ix
            end

            # println(BSTDOUT, "\n ix_petsc_row: ", ix_petsc_row)
            # println(BSTDOUT, " ix_petsc_col: ", ix_petsc_col)

            for row_ix = 1:length(ix_petsc_row)
              for col_ix = 1:length(ix_petsc_col)
                block_to_add[row_ix, col_ix] = stab_A.A[row_ix, col_ix, block_ix]
              end
            end

            # We should be subtracting, so we should scale block_to_add by -1.0
            scale!(block_to_add, -1.0)

            # now subtract the filtered DiagJac to the actual Jacobian, which will remove the positive eigenvalues of
            #   the strong Jacobian from the actual Jacobian
        
            # Add the negated block to the existing Jac inside the ls_ds LO object
            set_values1!(lo_ds_innermost.A, ix_petsc_row, ix_petsc_col, block_to_add, ADD_VALUES)

          end

          MatAssemblyBegin(lo_ds_innermost.A, MAT_FINAL_ASSEMBLY)
          MatAssemblyEnd(lo_ds_innermost.A, MAT_FINAL_ASSEMBLY)

        end   # end if opts["stabilize_v"]
=#
