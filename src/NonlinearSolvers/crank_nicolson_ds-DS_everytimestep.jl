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

  Ma_pert_mag = opts["perturb_Ma_magnitude"]
  pert = complex(0, Ma_pert_mag)
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
    dRdM_vec[ix_dof] = imag(dRdM_vec[ix_dof])/Ma_pert_mag
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

    # combine (-dRdM - dRdq * v^(n)) into b
    b_vec[ix_dof] = - dRdM_vec[ix_dof] - dRdq_vn_prod[ix_dof]

  end

  #------------------------------------------------------------------------------
  # Now the calculation of v_ix at n+1
  # Solve:
  #   dR/dq^(n+1) v^(n+1) = b
  #--------

  fill!(v_vec, 0.0)

  # Recalculate dRdq
  calcLinearOperator(ls_ds, mesh, sbp, eqn, opts, ctx_residual, t)

  # linearSolve: solves Ax=b for x
  #   ls::StandardLinearSolver
  #   b::AbstractVector     -> RHS
  #   x::AbstractVector     -> what is solved for
  #   verbose=5
  linearSolve(ls_ds, b_vec, v_vec, verbose)      
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

  #### The above is all new CSR code

  #------------------------------------------------------------------------------
  # getting from v_vec to term23
  fill!(dDdu, 0.0)
  evalFunctionalDeriv(mesh, sbp, eqn, opts, objective, dDdu)    # dDdu is func_deriv_arr
  fill!(dDdu_vec, 0.0)     # not sure this is necessary
  array3DTo1D(mesh, sbp, eqn, opts, dDdu, dDdu_vec)      # dDdu -> dDdu_vec

  for v_ix = 1:length(v_vec)
    # this accumulation occurs across all dofs and all time steps.
    term23 += quad_weight * dDdu_vec[v_ix] * v_vec[v_ix]
  end

  #------------------------------------------------------------------------------
  # here is where we should be calculating the 'energy' to show that it is increasing over time
  #   'energy' = L2 norm of the solution
  # JEH: So, at each time step, evaluate: sum_{i,j,k} q[i,j,k]*q[i,j,k]*sbp.w[j]*jac[j,k]
  #      (here I assume jac is proportional to the element volume)
  if opts["write_L2vnorm"]
    L2_v_norm = calcNorm(eqn, v_vec)
    @mpi_master println(f_L2vnorm, i, "  ", L2_v_norm)

    if (i % opts["output_freq"]) == 0
      @mpi_master flush(f_L2vnorm)
    end
  end



end     # end 'if i != 1'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# start of old DS code (the above is mostly new CSR stuff (except for the dDdu & term23 stuff))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# for visualization of element level DS energy: end of all stages
if opts["write_L2vnorm"]
  # special case for calculating v_energy:
  #   Need to calculate v_vec after the first stage for use in the first-stage-only v_energy calculation.
  # TODO: why only first-stage?
  # TODO: stage parameter divide
  # TODO: store old_q_vec
  #   Previously, this would only be done at the end of all the stages
  # for v_ix = 1:length(v_vec)
    # v_vec[v_ix] = imag(q_vec[v_ix])/Ma_pert_mag         # v_vec alloc'd outside timestep loop
  # end
  # RK4: the v_vec assignment is commented out because it should be set inside perturb_Ma just above
  fill!(R_stab, 0.0)
  for j = 1:length(eqn.q_vec)
    # R_stab[j] = (q_vec[j] - old_q_vec[j])/(fac*delta_t)     # This was LSERK's
    R_stab[j] = (eqn.q_vec[j] - old_q_vec[j])/(dt)
  end
  fill!(v_energy, 0.0)
  for j = 1:length(eqn.q_vec)
    v_energy[j] = v_vec[j]*eqn.M[j]*imag(R_stab[j])/Ma_pert_mag
  end

  if (i % output_freq) == 0
    saveSolutionToMesh(mesh, v_energy)
    fname = string("v_energy_stageall_", i)
    writeVisFiles(mesh, fname)
  end

  # i  calcNorm    avg   min   max
  v_energy_mean_local = mean(v_energy)
  # v_energy_min_local = minimum(v_energy)      # not doing abs on purpose
  # v_energy_max_local = maximum(v_energy)

  v_energy_mean = MPI.Allreduce(v_energy_mean_local, MPI.SUM, mesh.comm)
  # v_energy_min = MPI.Allreduce(v_energy_min_local, MPI.SUM, mesh.comm)
  # v_energy_max = MPI.Allreduce(v_energy_max_local, MPI.SUM, mesh.comm)

  v_energy_norm = calcNorm(eqn, v_energy)

  # @mpi_master println(f_v_energy, i, "  ", real(v_energy_norm), "  ", real(v_energy_mean), "  ", real(v_energy_min), "  ", real(v_energy_max))
  @mpi_master println(f_v_energy_stageall, i, "  ", real(v_energy_norm))
  if (i % 500) == 0
    @mpi_master flush(f_v_energy_stageall)
  end

end
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
