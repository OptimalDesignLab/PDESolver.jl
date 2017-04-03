
export cnAdjLoadChkpt, cnAdjCalcdRdu, cnAdjDirect

"""
function cnAdjLoadChkpt:
  loads a checkpoint from disk based on desired forward sweep index

  returns eqn object, which has only the necessary correct fields for the unsteady adjoint
"""
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

"""
function cnAdjCalcdRdu:
  calculates the Jacobian dRdu at an arbitrary time t given an associated eqn_dummy object

  returns a jac
"""
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

""" 
function cnAdjDirect:
  performs a direct solve of the linear adjoint PDE

  (*)_(i+1) corresponds to t
  (*)_i corresponds to t_nextstep

  given psi_(i+1)

  returns psi_i
"""
function cnAdjDirect(mesh, sbp, opts, adj, physics_func, jac, i_actual, h, t)
# adj_nextstep.q_vec = cnAdjDirect(mesh, sbp, opts, adj_nextstep, jac, i_actual, t_nextstep)

  # TODO: verify that passing t (instead of t_nextstep) into this is ok wrt cnAdjRhs & cnAdjJac

  t_nextstep = t - h

  eqn_dummy = cnAdjLoadChkpt(mesh, sbp, opts, adj, physics_func, i_actual, t_nextstep)
  jac = cnAdjCalcdRdu(mesh, sbp, opts, eqn_dummy, physics_func, t_nextstep)
  dRdu_i = transpose(jac)     

  # TODO: off by one error on dRdu? even though calculated at t_nextstep, am I loading q_vec at i+1 instead of i?

  dJdu_i = calcdJdu_CS(mesh, sbp, eqn_dummy, opts)  # obtain dJdu at time step i

  I = eye(length(adj.q_vec))

  B1 = I - 0.5*h*dRdu_i

  B2 = adj.q_vec + 0.5*h*dRdu_i*adj.q_vec - dJdu_i

  nextstep_q_vec = zeros(adj.q_vec)
  nextstep_q_vec = B1\B2

  return nextstep_q_vec


end     # end function cnAdjDirect
