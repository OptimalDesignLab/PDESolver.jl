# implementation of 5 stage, 4th order Low Storage Explicit Runge Kutta
# from Carpenter and Kennedy's Fourth-Order 2N-Storage Runge Kutta Schemes
# NASA Technical Memorandum 109112

# this code borrows some infrastructure (pde_pre_func, pde_post_func) from
# the classical rk4 file
# Also borrows the RK$CheckpointData


"""
  This function implements the 5 stage, 4th order Low Storage Explicit Runge Kutta scheme of
  Carpenter and Kennedy

  Arguments:
    f: a function that evalutes du/dt = f(q, t)
    delta_t: the time step
    t_max: the maximum time value
    q_vec: vector (of length numDof) containing initial solution.  Will contain final solution 
           at exit
    res_vec: vector to store the residual in during evaluation of f.  The contents of this vector
             at exit is undefined
    pre_func: function to call after new values are written into q_vec but before f is called
    post_func: function to call after f is called but before res_vec is accessed
    ctx: tuple arguments of f (ie. f = f(ctx...))
    opts: options dictionary
    timing: a Timings object
    
  Keyword Arguments:
    majorIterationCallback: function to call after first function evaluation of each time step, ie. 
                            when q_vec and res_vec have been updated.  Useful for logging.  Defaults
                            to no-op
    res_tol: stopping tolerance for residual (useful for pseudo-timestepping), default -1.0
    real_time: whether or not to advance time (ie. pseudo timestepping or not) default faulse

  See the documentation for rk4.
"""

function lserk54(f::Function, delta_t::AbstractFloat, t_max::AbstractFloat, 
             q_vec::AbstractVector, res_vec::AbstractVector, pre_func, 
             post_func, ctx, opts, timing::Timings=Timings(); 
             majorIterationCallback=((a...) -> (a...)), 
             res_tol = -1.0, real_time=false)

  # LSERK coefficients
  const a_coeffs = [0; 
                    -567301805773.0/1357537059087.0; 
                    -2404267990393.0/2016746695238.0; 
                    -3550918686646.0/2091501179385.0; 
                    -1275806237668.0/842570457699.0]

  const b_coeffs = [1432997174477.0/9575080441755.0; 
                    5161836677717.0/13612068292357.0; 
                    1720146321549.0/2090206949498.0; 
                    3134564353537.0/4481467310338.0; 
                    2277821191437.0/14882151754819.0]

  const c_coeffs = [0; 
                    1432997174477/9575080441755; 
                    2526269341429/6820363962896; 
                    2006345519317/3224310063776; 
                    2802321613138/2924317926251]

  myrank = MPI.Comm_rank(MPI.COMM_WORLD)
#  MPI.Barrier(MPI.COMM_WORLD)
  if myrank == 0
    println(BSTDOUT, "\nEntered lserk54")
    println(BSTDOUT, "res_tol = ", res_tol)
  end
#  flush(BSTDOUT)
#  MPI.Barrier(MPI.COMM_WORLD)
# res_tol is alternative stopping criteria


  # unpack options
  output_freq = opts["output_freq"]::Int
  use_itermax = opts["use_itermax"]::Bool
  if use_itermax
    itermax = opts["itermax"]
  end

  use_checkpointing = opts["use_checkpointing"]::Bool
  chkpoint_freq = opts["checkpoint_freq"]::Int
  ncheckpoints = opts["ncheckpoints"]::Int


  t = 0.0  # timestepper time
  treal = 0.0  # real time (as opposed to pseudo-time)
  t_steps = round(Int, t_max/delta_t)
  println(BSTDOUT, "t_steps: ",t_steps)
  println(BSTDOUT, "delta_t = ", delta_t)

  (m,) = size(q_vec)

  # allocate storage
  # this is actually a 3N scheme because the function f cannot overwrite its
  # input
  dq_vec = zeros(q_vec)

  if myrank == 0
    _f1 = open("convergence.dat", "a")
    f1 = BufferedIO(_f1)
  end

  # setup all the checkpointing related data
  chkpointer, chkpointdata, skip_checkpoint = explicit_checkpoint_setup(opts, myrank)
  istart = chkpointdata.i

  #------------------------------------------------------------------------------
  # direct sensitivity of Cd wrt M : setup
  if opts["perturb_Ma"]
    term23 = 0.0
    Ma_pert = opts["perturb_Ma_magnitude"]
    quad_weight = delta_t
  end   # end if opts["perturb_Ma"]

  # Initialize quadrature weight for trapezoidal rule
  #   This will be adjusted within the loop for the first & final time steps

  #------------------------------------------------------------------------------

  flush(BSTDOUT)
  #------------------------------------------------------------------------------
  # Main timestepping loop
  timing.t_timemarch += @elapsed for i=istart:(t_steps + 1)

    if opts["perturb_Ma"]
      if i == istart || i == (t_steps + 1)
        quad_weight = delta_t/2.0             # first & last time step, trapezoid rule quadrature weight
      end
    end   # end if opts["perturb_Ma"]

    t = (i-2)*delta_t

    @mpi_master if i % output_freq == 0
       println(BSTDOUT, "\ntimestep ",i)
       if i % output_freq == 0
         flush(BSTDOUT)
       end
    end

    if use_checkpointing && i % chkpoint_freq == 0
      if skip_checkpoint    # skip only the first checkpoint
        skip_checkpoint = false
      else
        @mpi_master println(BSTDOUT, "Saving checkpoint at timestep ", i)
        skip_checkpoint = false
        # save all needed variables to the chkpointdata
        chkpointdata.i = i

        if countFreeCheckpoints(chkpointer) == 0
          freeOldestCheckpoint(chkpointer)  # make room for a new checkpoint
        end

        # save the checkpoint
        saveNextFreeCheckpoint(chkpointer, ctx..., opts, chkpointdata)
      end   # end of if skip_checkpoint check
    end   # end of if use_checkpointing check


    #--------------------------------------------------------------------------
    # stage 1
#    f(params, u, F_vals, t_i)

    pre_func(ctx..., opts)
    if real_time treal = t end
    timing.t_func += @elapsed f(ctx..., opts, treal)
    sol_norm = post_func(ctx..., opts)
 
    #--------------------------------------------------------------------------
    # callback and logging
    timing.t_callback += @elapsed majorIterationCallback(i, ctx..., opts, BSTDOUT) # dirsens note: here is where drag is written

    # logging
    @mpi_master if i % 1 == 0
      println(f1, i, " ", sol_norm)
    end
    
    @mpi_master if i % output_freq == 0
      println(BSTDOUT, "flushing convergence.dat to disk")
      flush(f1)
    end

    # check stopping conditions
    if (sol_norm < res_tol) && !real_time
      if myrank == 0
        println(BSTDOUT, "breaking due to res_tol, res norm = $sol_norm")
        close(f1)
        flush(BSTDOUT)
      end
      break
    end

    if use_itermax && i > itermax
      if myrank == 0
        println(BSTDOUT, "breaking due to itermax")
        close(f1)
        flush(BSTDOUT)
      end
      break
    end


    #--------------------------------------------------------------------------
    # remaining stages

    # stage 1 update

    fac = b_coeffs[1]
    for j=1:length(q_vec)
      dq_vec[j] = delta_t*res_vec[j]
      q_vec[j] += fac*dq_vec[j]
    end

    # loop over remaining stages
    for stage=2:5
      pre_func(ctx..., opts) 
      if real_time
        treal = t + c_coeffs[stage]*delta_t 
      end
      timing.t_func += @elapsed f( ctx..., opts, treal)
      post_func(ctx..., opts, calc_norm=false)

      # update
      fac = a_coeffs[stage]
      fac2 = b_coeffs[stage]
      for j=1:length(q_vec)
        dq_vec[j] = fac*dq_vec[j] + delta_t*res_vec[j]
        q_vec[j] += fac2*dq_vec[j]
      end
    end  # end loop over stages

    #------------------------------------------------------------------------------
    # direct sensitivity of Cd wrt M : calculation each time step
    if opts["perturb_Ma"]

      # v is the direct sensitivity, du/dM
      # Ma has been perturbed during setup, in types.jl when eqn.params is initialized
      v_vec = zeros(q_vec)      # direct sensitivity vector
      for v_ix = 1:length(v_vec)
        v_vec[v_ix] = imag(q_vec[v_ix])/norm(Ma_pert)
      end

      # term2 is the partial deriv of the functional wrt the state: dCd/du
      #=
      # this is what needs to be adapted for fwd mode
      term2 = zeros(q)
      calcFunctionalDeriv(mesh, sbp, eqn, opts, functionalData, term2)    # term2 is func_deriv_arr
      =#

      term2 = zeros(q)
      calcFunctionalDeriv(mesh, sbp, eqn, opts, functionalData, term2)    # term2 is func_deriv_arr

      # TODO: need mesh, functionalData
      #   mesh: easy pack into ctx
      #   objective: figure out this; it's created every evalResidual call otherwise. Maybe call it again? So
      #               a new objective would be created every time step

      # new method: get dCd/du analytically
      # eqn: dCd/du = 4*D/(u_inf^3*rho_inf*c)
      # rho_free: eqn.params.rho_free
      # chord: 1.0      DOUBLE CHECK THIS
      # figure out how to get drag. it's being calculated within majorIterationCallback
      # how to get u: for each dof: q[2]/q[1] - x component, q[3]/q[1] - y component
      # wait---- this is u_inf?
      #   if so, Ma = u/a, so u = Ma*a, so u_inf = a_free*Ma

      # (mesh, sbp, eqn) = ctx...
      mesh = ctx[1]   # fastest way to grab mesh from ctx?
      sbp = ctx[2]   # fastest way to grab mesh from ctx?
      eqn = ctx[3]   # fastest way to grab mesh from ctx?

      #=
      # this is wrong. u_inf is not correct.
      Ma_unpert = real(eqn.params.Ma)
      u_inf = eqn.params.a_free*Ma_unpert
      chord = 1.0
      rho_inf = eqn.params.rho_free

      objective = EulerEquationMod.createObjectiveFunctionalData(mesh, sbp, eqn, opts)
      drag = real(evalFunctional(mesh, sbp, eqn, opts, objective))

      term2 = (4*drag) / (u_inf^3 * rho_inf * chord)
      =#

      # do the dot product of the two terms, and save
      for v_ix = 1:length(v_vec)
        term23 += quad_weight * term2[v_ix] * v_vec[v_ix]     # this accumulation occurs across all dofs and all time steps.
      end
    end   # end if opts["perturb_Ma"]

  end  # end loop over timesteps




  #------------------------------------------------------------------------------
  # LSERK end of time step stuff

  # final update
  t += delta_t

  @mpi_master println("---------------------------------------------")
  @mpi_master println("   LSERK: final time step reached. t = $t")
  @mpi_master println("---------------------------------------------")

  if opts["perturb_Ma"]
    @mpi_master f_term23 = open("term23.dat", "w")
    @mpi_master println(f_term23, i, " ", term23)
    @mpi_master flush(f_term23)
    @mpi_master close(f_term23)
    @mpi_master f_Ma = open("Ma.dat", "w")
    @mpi_master println(f_Ma, i, " ", real(eqn.params.Ma))
    @mpi_master close(f_Na)
  end   # end if opts["perturb_Ma"]


  if myrank == 0
    close(f1)
  end

  flush(BSTDOUT)

  return t
end  # end lserk54

"""
  See rk4 method with same signature
"""
function lserk54(f::Function, h::AbstractFloat, t_max::AbstractFloat, 
             q_vec::AbstractVector, res_vec::AbstractVector, ctx, opts,
             timing::Timings=Timings(); 
             majorIterationCallback=((a...) -> (a...)), res_tol=-1.0, 
             real_time=false)

  t = lserk54(f::Function, h::AbstractFloat, t_max::AbstractFloat,
              q_vec::AbstractVector, res_vec::AbstractVector,
              pde_pre_func, pde_post_func, ctx, opts, timing; 
              majorIterationCallback=majorIterationCallback, res_tol=res_tol,
              real_time=real_time)

        return t
end
