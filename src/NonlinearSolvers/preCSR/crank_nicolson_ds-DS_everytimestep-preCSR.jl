
    #------------------------------------------------------------------------------
    # direct sensitivity of Cd wrt M : calculation each time step
    if opts["perturb_Ma"]

      # v is the direct sensitivity, du/dM
      # Ma has been perturbed during setup, in types.jl when eqn.params is initialized
      for v_ix = 1:length(v_vec)
        v_vec[v_ix] = imag(eqn.q_vec[v_ix])/Ma_pert_mag         # v_vec alloc'd outside timestep loop
      end

      # dDdu is the partial deriv of the functional wrt the state: dCd/du
      fill!(dDdu, 0.0)     # initialized before timestepping loop
      # evalFunctional calls disassembleSolution, which puts q_vec into q
      # should be calling evalFunctional, not calcFunctional.
      #     disassemble isn't getting called. but it shouldn't matter b/c DG
      # EulerEquationMod.evalFunctionalDeriv(mesh, sbp, eqn, opts, objective, dDdu)    # dDdu is func_deriv_arr
      evalFunctionalDeriv(mesh, sbp, eqn, opts, objective, dDdu)    # dDdu is func_deriv_arr
      println(" >>>> i: ", i, "  quad_weight: ", quad_weight, "  dDdu: ", vecnorm(dDdu), "  v_vec: ", vecnorm(v_vec))

      # do the dot product of the two terms, and save
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

    end   # end if opts["perturb_Ma"]

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
