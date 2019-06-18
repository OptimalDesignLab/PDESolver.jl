
# Note: passed in for the eqn argument is 'eqn_nextstep'
function stabilizeCNDSLO(lo_ds, mesh, sbp, eqn, opts, ctx_residual, t)

  # ctx_residual = (f, eqn, h, newton_data, stab_A, stab_assembler, clipJacData)
  f =               ctx_residual[1]
  eqn_old =         ctx_residual[2]
  h =               ctx_residual[3]
  newton_data =     ctx_residual[4]
  stab_A =          ctx_residual[5]
  stab_assembler =  ctx_residual[6]
  clipJacData =     ctx_residual[7]
  v_vec =           ctx_residual[8]

  myrank = eqn.myrank

  # get i from t
  i = round(Int, t/h + 2)
  lo_ds_innermost = getBaseLO(lo_ds)

  # Note: attempting to print eigenvalues of lo_ds_innermost.A here 
  #   with eigvals(full(lo_ds_innermost.A)) takes TOO MUCH MEMORY & TIME

  # We have to zero out the DiagJac, as assembleElement inside evalJacobianStrong accumulates.
  MatZeroEntries(stab_assembler.A)

  if opts["stabilize_on_which_dFdq"] == "Minv"
  elseif opts["stabilize_on_which_dFdq"] == "noMinv"
    eqn.params.use_Minv = 0
  end

  # stores the strong Jacobian (volume Jacobian) into stab_assembler.A
  evalJacobianStrong(mesh, sbp, eqn, opts, stab_assembler, t)

  @mpi_master if (i % opts["output_freq"]) == 0
    if opts["write_strongJac_eigvals"]
      # eigenvalue plotting, strong Jac, before any filtering
      eigs_strongJac_before_stab = eigvals(stab_assembler.A)
      filename = string("i", i,"-eigs_strongJac_before_stab.dat")
      writedlm(filename, eigs_strongJac_before_stab)

      # copy it into a separate DiagJac so the stabilization can be analyzed
      # type notes:
      #   typeof(stab_assembler): NonlinearSolvers.AssembleDiagJacData{Complex{Float64}}
      #   typeof(stab_assembler.A): NonlinearSolvers.DiagJac{Complex{Float64}}
      #   typeof(stab_assembler.A.A): Array{Complex{Float64},3}
      diagJac_for_comparison = deepcopy(stab_assembler.A.A)
    end
  end

  # filterDiagJac
  #   location: jacobian_diag.jl
  #
  #   The third argument is q_vec in the fn signature.
  #   It is used as part of computing the quadprog stabilization (findStablePerturbation!).
  #   It is not used when opts["stabilization_method"] is "clipJac" or "clipJacFast".
  #   For the explicit stabilization what was passed in was 'real(tmp_imag)'.
  #   For implicit: It should be v_vec. Of the previous time step.
  #     See the AFOSR report's section on quadprog; it uses the adjoint.
  #     I suppose it could be the next time step, but that would require some implicit solving 
  #     and KSP iterations? Maybe room for investigation later. (future work)
  #     Because right after this stabilizeCNDSLO is called, 
  #     linearSolve is called to find v_vec^(n+1)
  eigs_to_remove = opts["eigs_to_remove"]
  # TODO: make sure the fourth argument is v_vec!!    TODO
  # eqn.q_vec allows it to proceed past a SingularException, but it is incorrect.

  # v_vec is zeros at i == 2; causes a SingularException
  if i == 2 && opts["stabilization_method"] == "quadprog"
    println(BSTDOUT, " breaking from stabilizeCNDSLO; i == 2 and quadprog selected but v_vec is zeros")
    # need to reset eqn.params.use_Minv! or else it'll be wrong for the next primal evalResidual
    if opts["stabilize_on_which_dFdq"] == "noMinv"
      eqn.params.use_Minv = 1
    end
    return
  end

  numEigChgsAllEls, numEigChgs_arrayEls = filterDiagJac(mesh, eqn, opts, v_vec, clipJacData, 
                                                        stab_A, eigs_to_remove=eigs_to_remove)
  @mpi_master println(BSTDOUT, " numEigChgsAllEls: ", numEigChgsAllEls)

  @mpi_master if (i % opts["output_freq"]) == 0

    if opts["write_numEigChgs_arrayEls"]
      numEigChgs_vec = zeros(mesh.numDof)
      convertArrayOfElsToDofs(mesh, numEigChgs_arrayEls, numEigChgs_vec)
      saveSolutionToMesh(mesh, numEigChgs_vec)
      fname = string("numEigChgs_arrayEls_", i)
      writeVisFiles(mesh, fname)
    end

    if opts["write_strongJac_eigvals"]
      # eigenvalue plotting, strong Jac, after filtering
      eigs_whatisremovedfromJac = eigvals(stab_assembler.A)
      filename = string("i", i,"-eigs_whatisremovedfromJac.dat")
      writedlm(filename, eigs_whatisremovedfromJac)

      # subtract the stabilization term from the diagJac_for_comparison
      # This is done for eigcluster analysis only.
      # The stabilization term is intended to be subtracted from the full Jac below.
      diagJac_for_comparison = diagJac_for_comparison - stab_assembler.A.A
      eigs_strongJac_after_stab = eigvals(diagJac_for_comparison)
      filename = string("i", i,"-eigs_strongJac_after_stab.dat")
      writedlm(filename, eigs_strongJac_after_stab)
    end

  end

  # Now add each block of the stabilized strong jacobian to the full Jacobian
  # We are converting between the 2D element Jacobian in each block of the DiagJac
  #   to the 4D form required by assembleElement.
  # DiagJac dims: (blocksize, blocksize, numEl)
  # res_jac dims: (numDofPerNode, numDofPerNode, numNodesPerElement, numNodesPerElement)
  assembler = _AssembleElementData(lo_ds_innermost.A, mesh, sbp, eqn, opts)
  blocksize = mesh.numDofPerNode*mesh.numNodesPerElement
  this_res_jac = zeros(Complex{Float64}, mesh.numDofPerNode, mesh.numDofPerNode, 
                       mesh.numNodesPerElement, mesh.numNodesPerElement)

  if opts["stabilize_on_which_dFdq"] == "noMinv"
    eqn.params.use_Minv = 1
  end

  for el_ix = 1:mesh.numEl

    # This does not make a difference.
    if opts["zeroout_this_res_jac"] == true
      fill!(this_res_jac, 0.0)
    end

    for q = 1:mesh.numNodesPerElement
      for p = 1:mesh.numNodesPerElement

        if opts["stabilize_on_which_dFdq"] == "noMinv"
          if opts["stab_Minv_val"] == "one"
            Minv_val = 1.0        # TODO TODO why does this work, and not the proper Minv?
          elseif opts["stab_Minv_val"] == "actual"
            Minv_val = mesh.jac[p, el_ix]/sbp.w[p]  # entry in Minv
          else
            error("Minv_val option set improperly.")
          end
          # 20190612
          # mistakenly using i (timestep index) instead of el_ix caused some stabilization
        end
        @simd for j = 1:mesh.numDofPerNode
          @simd for i = 1:mesh.numDofPerNode

            # Within each DiagJac block, it is indexed along one dimension as:
            #   all the dofs on node 1, all the dofs on node 2, etc. 
            # This permits the following conversion:
            i1 = i + (p-1)*mesh.numDofPerNode
            j1 = j + (q-1)*mesh.numDofPerNode

            this_res_jac[i, j, p, q] = stab_A.A[i1, j1, el_ix]

            if opts["stabilize_on_which_dFdq"] == "noMinv"
              this_res_jac[i, j, p, q] *= Minv_val
            end

          end
        end

      end   # end loop over p
    end   # end loop over q

    # this_res_jac should contain all the positive eigs, so if we subtract, 
    #   we are left with only negative and zero eigenvalues.
    # MUST scale this_res_jac by -1.0 for clipJac
    # if opts["stabilization_method"] != "quadprog"     # TODO: figure out if this is necessary
      scale!(this_res_jac, -1.0)
    # end

    assembleElement(assembler, mesh, el_ix, this_res_jac)
    # This is calling function assembleElement(helper::_AssembleElementData{PetscMat}, 
    #                                          mesh::AbstractMesh,
    #                                          elnum::Integer, jac::AbstractArray{T, 4}) where T
    # in jacobian.jl. Line 888 or so

  end   # end loop over elements
  
  # Note: attempting to print eigenvalues of lo_ds_innermost.A here 
  #   with eigvals(full(lo_ds_innermost.A)) takes TOO MUCH MEMORY & TIME

  return nothing

end

