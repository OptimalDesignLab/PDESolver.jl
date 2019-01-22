
  if opts["perturb_Ma"]

    @mpi_master if opts["stabilize_v"]
      close(f_stabilize_v)
    end

    @mpi_master close(f_drag)

    @mpi_master println(" eqn.params.Ma: ", eqn.params.Ma)
    @mpi_master println(" Ma_pert: ", Ma_pert)
    eqn.params.Ma -= Ma_pert      # need to remove perturbation now
    @mpi_master println(" pert removed from Ma")
    @mpi_master println(" eqn.params.Ma: ", eqn.params.Ma)

    finaliter = calcFinalIter(t_steps, itermax)
    Cd, dCddM = calcDragTimeAverage(mesh, sbp, eqn, opts, dt, finaliter)   # will use eqn.params.Ma
    term23 = term23 * 1.0/t     # final step of time average: divide by total time
    global_term23 = MPI.Allreduce(term23, MPI.SUM, mesh.comm)
    total_dCddM = dCddM + global_term23

    # Cd calculations
    @mpi_master begin
      f_total_dCddM = open("total_dCddM.dat", "w")
      println(f_total_dCddM, " Cd: ", Cd)
      println(f_total_dCddM, " dCd/dM: ", dCddM)
      println(f_total_dCddM, " global_term23: ", global_term23)
      println(f_total_dCddM, " total dCd/dM: ", total_dCddM)
      flush(f_total_dCddM)
      close(f_total_dCddM)
      println(" Cd: ", Cd)
      println(" dCd/dM: ", dCddM)
      println(" global_term23: ", global_term23)
      println(" total dCd/dM: ", total_dCddM)
    end

  end   # end if opts["perturb_Ma"]

  @mpi_master begin
    f_Ma = open("Ma.dat", "w")
    println(f_Ma, eqn.params.Ma)
    close(f_Ma)
    f_dt = open("delta_t.dat", "w")
    println(f_dt, dt)
    close(f_dt)

    println(" ")
    println(" run parameters that were used:")
    if opts["perturb_Ma"]
      println("    Ma: ", eqn.params.Ma + Ma_pert)
    else
      println("    Ma: ", eqn.params.Ma)
    end
    println("    aoa: ", eqn.params.aoa)
    println("    dt: ", dt)
    println("    a_inf: ", eqn.params.a_free)
    println("    rho_inf: ", eqn.params.rho_free)
    println("    c: ", 1.0)
    println("    mesh.coord_order: ", mesh.coord_order)
    println(" ")
    println("    opts[stabilization_method]: ", opts["stabilization_method"])
    println("    opts[output_freq]: ", opts["output_freq"])
    println("    opts[use_itermax]: ", opts["use_itermax"])
    println("    opts[itermax]: ", opts["itermax"])
    println("    opts[use_checkpointing]: ", opts["use_checkpointing"])
    println("    opts[checkpoint_freq]: ", opts["checkpoint_freq"])
    println("    opts[ncheckpoints]: ", opts["ncheckpoints"])
    println(" ")
  end

  if opts["write_L2vnorm"]
    @mpi_master close(f_L2vnorm)
    # @mpi_master close(f_v_energy)
    @mpi_master close(f_v_energy_stageall)
  end
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
