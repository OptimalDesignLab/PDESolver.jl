# This entire code block is called with an include().
# It is wrapped in:  if opts["perturb_Ma_CN"]
# Such that it should only be included/evaluated when opts["perturb_Ma_CN"] is set

if i != 1     # can't do F(q^(n+1)) + F(q^(n)) until we have both

  # Store both q_vec's & both res_vec's. Will be put back after all DS calcs.
  for ix_dof = 1:mesh.numDof
    beforeDS_eqn_q_vec[ix_dof] = eqn.q_vec[ix_dof]
    beforeDS_eqn_nextstep_q_vec[ix_dof] = eqn_nextstep.q_vec[ix_dof]

    beforeDS_eqn_res_vec[ix_dof] = eqn.res_vec[ix_dof]
    beforeDS_eqn_nextstep_res_vec[ix_dof] = eqn_nextstep.res_vec[ix_dof]
  end


  #------------------------------------------------------------------------------
  # This section calculates dR/dM

  Ma_pert = opts["perturb_Ma_magnitude"]
  pert = complex(0, Ma_pert)
  Ma += pert

  # now evalResidual to store into F(q^(n+1))
  f(mesh, sbp, eqn_nextstep, opts)      # F(q^(n+1)) -- with Ma perturbation, now in eqn_nextstep.res_vec

  for ix_dof = 1:mesh.numDof
    old_res_vec_Maimag[ix_dof] = new_res_vec_Maimag[ix_dof]        # F(q^(n)) -- with Ma perturbation
    new_res_vec_Maimag[ix_dof] = eqn_nextstep.res_vec[ix_dof]      # F(q^(n+1)) -- with Ma perturbation

    # form unsteady residual (res_hat) with F(q^(n)) and F(q^(n+1))
    res_hat_vec[ix_dof] = -0.5*eqn.Minv[ix_dof]*dt * (new_res_vec_Maimag[ix_dof] + old_res_vec_Maimag[ix_dof])
  end

  # obtain dR/dM using the complex step method
  for ix_dof = 1:mesh.numDof
    dRdM_vec[ix_dof] = imag(dRdM_vec[ix_dof])/Ma_pert_mag         # TODO: alloc dRdM_vec
  end

  Ma -= pert
  #------------------------------------------------------------------------------

  #------------------------------------------------------------------------------
  # This section calculates dR/dq * v^(n)

  # v_vec currently holds v at timestep n: v^(n)

  # need to add epsilon*v_vec*im (evi) to q_vec, then re-evaluate res_hat at the updated q
  # so: need to add evi to old_q_vec too. Because to form res_hat, need F(q^(n)) and F(q^(n+1))
  for ix_dof = 1:mesh.numDof
    eqn.q_vec[ix_dof]          += Ma_pert_mag*im*v_vec[ix_dof]
    eqn_nextstep.q_vec[ix_dof] += Ma_pert_mag*im*v_vec[ix_dof]
  end


  f(mesh, sbp, eqn, opts)             # F(q^(n) + evi) now in eqn.res_vec
  f(mesh, sbp, eqn_nextstep, opts)    # F(q^(n+1) + evi) now in eqn_nextstep.res_vec

  for ix_dof = 1:mesh.numDof
    
    # form unsteady residual (res_hat) with F(q^(n) + evi) and F(q^(n+1) + evi)
    res_hat_vec[ix_dof] = -0.5*eqn.Minv[ix_dof]*dt * (eqn_nextstep.res_vec[ix_dof] + eqn.res_vec[ix_dof])

    # calc dRdq * v^(n) by doing matrix-free complex step
    dRdq_vn_prod[ix_dof] = imag(res_hat_vec_tmpforprod[ix_dof])/Ma_pert_mag

    # remove the imaginary component from q used for matrix_dof-free product    # TODO: necessary?
    # eqn.q_vec[ix_dof] = complex(real(eqn.q_vec[ix_dof]))
    # eqn_nextstep.q_vec[ix_dof] = complex(real(eqn_nextstep.q_vec[ix_dof]))

    # combine (-dRdM - dRdq * v^(n)) into b     TODO: alloc b_vec
    b_vec[ix_dof] = - dRdM_vec[ix_dof] - dRdq_vn_prod[ix_dof]

  end

  #------------------------------------------------------------------------------
  # Now the calculation of v_ix at n+1
  #--------

  #--------

  #------------------------------------------------------------------------------
  # Restore q_vec & res_vec for both eqn & eqn_nextstep, as they were
  #   before the DS calcs
  for ix_dof = 1:mesh.numDof
    eqn.q_vec[ix_dof] = beforeDS_eqn_q_vec[ix_dof]
    eqn_nextstep.q_vec[ix_dof] = beforeDS_eqn_nextstep_q_vec[ix_dof]

    eqn.res_vec[ix_dof] = beforeDS_eqn_res_vec[ix_dof]
    eqn_nextstep.res_vec[ix_dof] = beforeDS_eqn_nextstep_res_vec[ix_dof]
  end


  #------------------------------------------------------------------------------
  # getting from v_vec to term23
  fill!(dDdu, 0.0)     # TODO: alloc dDdu
  evalFunctionalDeriv(mesh, sbp, eqn, opts, objective, dDdu)    # dDdu is func_deriv_arr
  fill!(dDdu_vec, 0.0)     # not sure this is necessary
  array3DTo1D(mesh, sbp, eqn, opts, dDdu, dDdu_vec)      # dDdu -> dDdu_vec

  for v_ix = 1:length(v_vec)    # TODO: alloc v_vec
    # this accumulation occurs across all dofs and all time steps.
    term23 += quad_weight * dDdu_vec[v_ix] * v_vec[v_ix]
  end


end     # end 'if i != 1'

