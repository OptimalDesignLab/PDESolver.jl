
export cnAdjLoadChkpt, cnAdjCalcdRdu, cnAdjDirect

"""
function cnAdjLoadChkpt:
  loads a checkpoint from disk based on desired forward sweep index

  returns eqn object, which has only the necessary correct fields for the unsteady adjoint
"""
# function cnAdjLoadChkpt(mesh, sbp, opts, adj, physics_func, i_fwd, t)
function cnAdjLoadChkpt(mesh, sbp, opts, adj, physics_func, i_fwd)

  # initialize dummy eqn object for jacobian calculation use
  # Note: can't just call this:
  #   eqn_dummy = AdvectionData_{Tsol, Tres, Tdim, Tmsh}(mesh, sbp, opts)
  # because it needs AdvectionEquationMod loaded, which needs PDESolver loaded, which can't be loaded more than once.
  eqn_dummy = eqn_deepcopy(mesh, sbp, adj, opts)                     # allocate a dummy eqn object

  # println("in cnAdjLoadChkpt: pointer(adj.q_vec): ", pointer(adj.q_vec))
  # println("in cnAdjLoadChkpt: pointer(adj.res_vec): ", pointer(adj.res_vec))
  # println("in cnAdjLoadChkpt: pointer(eqn_dummy.q_vec):   ", pointer(eqn_dummy.q_vec))
  # println("in cnAdjLoadChkpt: pointer(eqn_dummy.res_vec): ", pointer(eqn_dummy.res_vec))

  qvec_filename = string("qvec_for_adj-", i_fwd, ".dat")
  println("  Calculating Jac using forward sweep data from: ", qvec_filename)
  q_vec_with_complex = readdlm(qvec_filename)
  eqn_dummy.q_vec = q_vec_with_complex[:,1]     # because readdlm gives a field to the zero-valued complex part
  array1DTo3D(mesh, sbp, eqn_dummy, opts, eqn_dummy.q, eqn_dummy.q_vec)

  check_q_qvec_consistency(mesh, sbp, eqn_dummy, opts)
  println("  eqn_dummy.q_vec loaded. i_fwd: ", i_fwd)
  # print_qvec_coords(mesh, sbp, eqn_dummy, opts)


  vis_filename = string("solution_loadedfromdisk_ifwd-", i_fwd)
  saveSolutionToMesh(mesh, real(eqn_dummy.q_vec))
  writeVisFiles(mesh, vis_filename)

  # sync up eqn_dummy.q and eqn_dummy.q_vec

  # TODO: is this needed here? YES. 
  #     explain why here
  # TODO new: I don't think so now that eqn_deepcopy is properly implemented
  # eqn_dummy.q = reshape(eqn_dummy.q_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
  # eqn_dummy.res = reshape(eqn_dummy.res_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)

  #=
  # make eqn object consistent
  evalResidual(mesh, sbp, eqn_dummy, opts)
  array3DTo1D(mesh, sbp, eqn_dummy, opts, eqn_dummy.res, eqn_dummy.res_vec)

  # TODO: fix ordering & repetition of these
  eqn_dummy.q = reshape(eqn_dummy.q_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
  eqn_dummy.res = reshape(eqn_dummy.res_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
  =#
  # looks like the above commented section doesn't make a difference

  # use startDataExchange to sync up q/q_vec and q_face_send/recv
  if opts["parallel_type"] == 2 && mesh.npeers > 0
    startDataExchange(mesh, opts, eqn_dummy.q, eqn_dummy.q_face_send, eqn_dummy.q_face_recv, eqn_dummy.params.f)
  end

  return eqn_dummy

end   # end function cnAdjLoadChkpt

"""
function cnAdjCalcdRdu:
  calculates the Jacobian dRdu at an arbitrary time t given an associated eqn_dummy object

  returns a jac
"""
function cnAdjCalcdRdu(mesh, sbp, opts, eqn_dummy, physics_func, i_fwd, t)
  # Note: probably don't need to allocate another jac, but this should cause no problems aside from allocation cost
  newton_data_discard, jac, rhs_vec_discard = setupNewton(mesh, mesh, sbp, eqn_dummy, opts, physics_func, () )

  # make sure we're doing complex step! since res_copy is zeros, it would mess up the FD calc
  assert(opts["jac_method"] == 2)
  epsilon = opts["epsilon"]
  pert = complex(0, epsilon)

  # checked: eqn_dummy is loaded properly
  # checked: cjc args are correct
  # checked: norm of eqn_dummy.q_vec in cJC matches it here

  println("  calculating jac for loaded eqn_dummy. i_fwd = ", i_fwd,", t = ", t)
  calcJacobianComplex(newton_data_discard, mesh, sbp, eqn_dummy, opts, physics_func, pert, jac, t)

  println(" ::::::::::::: in cnAdjCalcdRdu, after calcJacobianComplex, norm(jac): ", norm(jac, 1))

  # filename = string("jac_from_cnAdjCalcdRdu-i_fwd-",i_fwd,".dat")
  # writedlm(filename, jac)

  return jac
end   # end function cnAdjCalcdRdu

""" 
function cnAdjDirect:
  performs a direct solve of the linear adjoint PDE

  (*)_(i+1) corresponds to t
  (*)_i corresponds to t_nextstep

  given psi_(i+1)

  returns psi_i
"""
function cnAdjDirect(mesh, sbp, opts, adj, physics_func, jac, i_fwd, h, t_steps, t, dRdu_global_rev)
  # Note: t is passed in (instead of t_nextstep), but t_nextstep is calc'd & used below
  #
  # Note: adj is not altered in this function. It is only used for passing into cnAdjLoadChkpt for an eqn_deepcopy

# adj_nextstep.q_vec = cnAdjDirect(mesh, sbp, opts, adj_nextstep, jac, i_fwd, t_nextstep)

  t_nextstep = t - h
  t_nextstep = negativeZeroCheck(t_nextstep)   # ensure negative zero is changed to zero

  # Load checkpoint at time t, aka i+1
  eqn_dummy = cnAdjLoadChkpt(mesh, sbp, opts, adj, physics_func, i_fwd)   # adj is only used for a eqn_deepcopy within
                                                                          #   cnAdjLoadCheckpt
  # checked: eqn_dummy is loaded properly
  # Use time t's eqn.q_vec and t_nextstep, aka i, to compute the jac (dRdu) at t_nextstep aka i
  jac = cnAdjCalcdRdu(mesh, sbp, opts, eqn_dummy, physics_func, i_fwd, t_nextstep)
  dRdu_i = transpose(jac)       # see derivation for transpose justification

  # TODO: double check that there is not an off by one error on dRdu:
  #       even though calculated at t_nextstep, am I loading q_vec at i+1 instead of i?

  I = eye(length(adj.q_vec))

  if opts["uadj_global"]
    ### GLOBAL dRdu section
    B1_pos_for_dRdu_global =  1.0*I - 0.5*h*dRdu_i
    B1_neg_for_dRdu_global = -1.0*I - 0.5*h*dRdu_i

    blksz = mesh.numDof
    # blksz = 3   # testing 44
    # B1_pos_for_dRdu_global = ones(blksz, blksz)*44
    # B1_neg_for_dRdu_global = ones(blksz, blksz)*55

    println(" GLOBAL: forming global dRdu (rev), i_fwd = $i_fwd")
    println(" GLOBAL: size of B1_pos_for_dRdu_global: ", size(B1_pos_for_dRdu_global))
    row_ix_start = (i_fwd-2)*blksz + 1
    row_ix_end = (i_fwd-2)*blksz + blksz
    col_ix_start = (i_fwd-2)*blksz + 1
    col_ix_end = (i_fwd-2)*blksz + blksz
    println(" GLOBAL: diagonal block ix's: [",row_ix_start,":",row_ix_end,", ",col_ix_start,":", col_ix_end,"]")
    dRdu_global_rev[row_ix_start:row_ix_end, col_ix_start:col_ix_end] = B1_pos_for_dRdu_global

    if i_fwd != t_steps + 1
      col_ix_start = (i_fwd-2)*blksz + blksz + 1
      col_ix_end = (i_fwd-2)*blksz + 2*blksz
      row_ix_start = (i_fwd-2)*blksz + 1
      row_ix_end = (i_fwd-2)*blksz + blksz
      println(" GLOBAL: off-diagonal block ix's: [",row_ix_start,":",row_ix_end,", ",col_ix_start,":", col_ix_end,"]")

      dRdu_global_rev[row_ix_start:row_ix_end, col_ix_start:col_ix_end] = B1_neg_for_dRdu_global
    end   # end of "if i_fwd != t_steps + 1"
    ### end GLOBAL
  end # end of "if opts["uadj_global"]

  dJdu_i = calcdJdu_CS(mesh, sbp, eqn_dummy, opts, h, t)  # obtain dJdu at time step i

  # 20170831. Looks like there's a transpose missing here on B1. 
  # See notes, LHS, 08/30/2017, for rigorous explanation.
  # B1 = I - 0.5*h*dRdu_i
  # B2 = adj.q_vec + 0.5*h*dRdu_i*adj.q_vec - dJdu_i
  B1 = transpose(I - 0.5*h*dRdu_i)
  B2 = adj.q_vec + 0.5*h*dRdu_i*adj.q_vec - dJdu_i

  nextstep_q_vec = zeros(adj.q_vec)
  nextstep_q_vec = B1\B2

  # filename = string("eqn_dummyqvec_from_cnAdjDirect-i_fwd-",i_fwd,".dat")
  # writedlm(filename, eqn_dummy.q_vec)

  return nextstep_q_vec, jac

end     # end function cnAdjDirect
