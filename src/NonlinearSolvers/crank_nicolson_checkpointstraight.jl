
export cnAdjLoadChkpt, cnAdjCalcdRdu

function cnAdjLoadChkpt(mesh, sbp, opts, adj, physics_func, i_actual, t)

  # initialize dummy eqn object for jacobian calculation use
  # Note: can't just call this:
  #   eqn_dummy = AdvectionData_{Tsol, Tres, Tdim, Tmsh}(mesh, sbp, opts)
  # because it needs AdvectionEquationMod loaded, which needs PDESolver loaded, which can't be loaded more than once.
  eqn_dummy = deepcopy(adj)                     # allocate a dummy eqn object

  filename = string("qvec_for_adj-", i_actual, ".dat")
  println("Setting IC for reverse sweep with file: ", filename)
  q_vec_with_complex = readdlm(filename)
  eqn_dummy.q_vec = q_vec_with_complex[:,1]     # because readdlm gives a field to the zero-valued complex part

  # sync up eqn_dummy.q and eqn_dummy.q_vec
  disassembleSolution(mesh, sbp, eqn_dummy, opts, eqn_dummy.q, eqn_dummy.q_vec)

  # use startDataExchange to sync up q/q_vec and q_face_send/recv
  if opts["parallel_type"] == 2 && mesh.npeers > 0
    startDataExchange(mesh, opts, eqn_dummy.q, eqn_dummy.q_face_send, eqn_dummy.q_face_recv, eqn_dummy.params.f)
  end

  return eqn_dummy

end   # end function cnAdjLoadChkpt

function cnAdjCalcdRdu(mesh, sbp, opts, eqn_dummy, physics_func, t)
  # Note: probably don't need to allocate another jac, but this should cause no problems aside from allocation cost
  newton_data_discard, jac, rhs_vec_discard = setupNewton(mesh, mesh, sbp, eqn_dummy, opts, physics_func)

  # make sure we're doing complex step! since res_copy is zeros, it would mess up the FD calc
  assert(opts["jac_method"] == 2)
  epsilon = opts["epsilon"]
  pert = complex(0, epsilon)

  calcJacobianComplex(newton_data_discard, mesh, sbp, eqn_dummy, opts, physics_func, pert, jac, t)

  return jac
end   # end function cnAdjCalcdRdu
