
  # Note: mesh, sbp, eqn, opts already in scope here due to function signature

  #------------------------------------------------------------------------------
  # direct sensitivity of Cd wrt M : setup
  if opts["perturb_Ma"]
    term23 = zero(eqn.params.Ma)      # type stable version of 'term23 = 0.0'
    Ma_pert_mag = opts["perturb_Ma_magnitude"]
    Ma_pert = complex(0, Ma_pert_mag)

    dDdu = zeros(eqn.q)
    dDdu_vec = zeros(Complex128, mesh.numDofPerNode * mesh.numNodesPerElement * mesh.numEl,)
    @mpi_master println("Direct sensitivity setup done.")
  end   # end if opts["perturb_Ma"]

  #------------------------------------------------------------------------------
  # capture direct sensitivity at the IC
  # v is the direct sensitivity, du/dM
  # Ma has been perturbed during setup, in types.jl when eqn.params is initialized
  if opts["write_drag"]
    objective = createFunctional(mesh, sbp, eqn, opts, 1)    # 1 is the functional num
    drag = real(evalFunctional(mesh, sbp, eqn, opts, objective))
    # note about drag writing: file_dict populated and file opened in src/solver/euler/types.jl
    @mpi_master f_drag = eqn.file_dict[opts["write_drag_fname"]]
    @mpi_master println(f_drag, 1, " ", drag)
    println("i: ", 1, "  myrank: ", myrank,"  drag: ", drag)
    @mpi_master flush(f_drag)
  end
  if opts["write_L2vnorm"]
    @mpi_master f_L2vnorm = eqn.file_dict[opts["write_L2vnorm_fname"]]
  end

  if opts["write_L2vnorm"]
    # for visualization of element level DS energy
    old_q_vec = zeros(eqn.q_vec)
    R_stab = zeros(eqn.q_vec)
    v_energy = zeros(eqn.q_vec)

    @mpi_master f_v_energy_stageall = open("v_energy_data_stageall.dat", "w")
  end
  if opts["perturb_Ma"]

    # this is the IC, so it gets the first time step's quad_weight
    i = 1       # note that timestep loop below starts at i = 2
    finaliter = calcFinalIter(t_steps, itermax)
    quad_weight = calcQuadWeight(i, dt, finaliter)

    #------------------------------------------------------------------------------
    # allocation of objects for stabilization routine
    if opts["stabilize_v"]

      stab_A = DiagJac(Complex128, mesh.numDofPerNode*mesh.numNodesPerElement, mesh.numEl)
      stab_assembler = AssembleDiagJacData(mesh, sbp, eqn, opts, stab_A)
      clipJacData = ClipJacData(mesh.numDofPerNode*mesh.numNodesPerElement)

      # Bv = zeros(Float64, length(q_vec), )
      Bv = zeros(Complex128, length(eqn.q_vec))
      tmp_imag = zeros(Float64, length(eqn.q_vec))
      dqimag_vec = zeros(Bv)

      @mpi_master f_stabilize_v = open("stabilize_v_updates.dat", "w")        # TODO: buffered IO

    end

    # Note: no stabilization of q_vec at the IC
    #   no evalResidual yet, and hasn't entered the time-stepper yet
    #   so no appropriate scaling factors like delta_t or fac or anything

    v_vec = zeros(eqn.q_vec)      # direct sensitivity vector       This is at the IC
    # for v_ix = 1:length(v_vec)
      # v_vec[v_ix] = imag(eqn.q_vec[v_ix])/Ma_pert_mag
    # end
    fill!(v_vec, 0.0)       # v at IC is all 0's  ++++ new CS'ing R method (CSR) <- CSR comment on all lines new to CS'ing R method

    q_vec_save_imag = Array{Float64}(length(eqn.q_vec))     # set up array for saving imag component of q (CSR)

    dDdu = zeros(eqn.q)      # First allocation of dDdu. fill! used below, during timestep loop
    # evalFunctional calls disassembleSolution, which puts q_vec into q
    # should be calling evalFunctional, not calcFunctional.
    evalFunctionalDeriv(mesh, sbp, eqn, opts, objective, dDdu)    # dDdu is func_deriv_arr
    println(" >>>> i: ", i, "  quad_weight: ", quad_weight, "  dDdu: ", vecnorm(dDdu), "  v_vec: ", vecnorm(v_vec))    # matches rk4_ds to 1e-5 or so, consistent with drag variation

    # do the dot product of the two terms, and save
    # this dot product is: dJdu*dudM
    dDdu_vec = zeros(Complex128, mesh.numDofPerNode * mesh.numNodesPerElement * mesh.numEl,)
    array3DTo1D(mesh, sbp, eqn, opts, dDdu, dDdu_vec)      # dDdu -> dDdu_vec

    for v_ix = 1:length(v_vec)
      # this accumulation occurs across all dofs and all time steps.
      term23 += quad_weight * dDdu_vec[v_ix] * v_vec[v_ix]
    end

    if opts["write_L2vnorm"]
      L2_v_norm = calcNorm(eqn, v_vec)
      @mpi_master println(f_L2vnorm, i, "  ", L2_v_norm)
    end

  end   # end if opts["perturb_Ma"]

