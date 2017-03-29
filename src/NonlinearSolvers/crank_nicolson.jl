# crank_nicolson.jl
# Crank-Nicolson implicit solver for PDEs

export crank_nicolson, cnResidual

push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/euler"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Debugging"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))

"""
crank_nicolson

  This function performs Crank-Nicolson implicit time solution, using a function 
  of the form du/dt = f(u, t)
  
  Arguments:
    * f  : function evaluation, must have signature (ctx..., opts, t)
    * h  : time step size
    * t_max: time value to stop time stepping (time starts at 0)
    * q_vec: vector of the u values
    * res_vec: vector of du/dt values (the output of the function f)
    * pre_func: function to to be called after the new u values are put into
                q_vec but before the function f is evaluated.  Must have
                signature: pre_func(ctx..., opts)
    * post_func: function called immediately after f is called.  The function
                 must have the signature res_norm = post_func(ctx..., opts, 
                 calc_norm=true),
                 where res_norm is a norm of res_vec, and calc_norm determines
                 whether or not to calculate the norm.
    * ctx: a tuple (or any iterable container) of the objects needed by
           f, pre_func, and post func.  The tuple is splatted before being
           passed to the functions.
    * opts : options dictionary

    Keyword Arguments:
    * majorIterationCallback: a callback function called after the first
                              stage, useful to do output or logging
    * res_tol : keyword arg, residual topping tolerance
    * neg_time : step through time in the negative direction,
                 starting at t_max, stepping with h, and ending at 0.0.


   The eqn.q_vec should hold the whichever variables (conservative or
   entropy) that the simulation should use.

   The idea behind pre_func, post func, and ctx is that they abstract what kind
   of system CN is timestepping. CN only needs to know about q_vec and
   res_vec.

   For physics modules, ctx should be (mesh, sbp, eqn) and q_vec and res_vec 
   should be eqn.q_vec and eqn.res_vec.

   TODO: fully document eqn/eqn_nextstep
"""
function crank_nicolson{Tmsh, Tsol}(f::Function, h::AbstractFloat, t_max::AbstractFloat,
                        mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP, eqn::AbstractSolutionData{Tsol},
                        opts, res_tol=-1.0; neg_time=false, obj_fn=obj_zero, store_u_to_disk=false)
                        # NEWNEW: neg_time, obj_fn
  #----------------------------------------------------------------------
#   throw(ErrorException("Crank-Nicolson is in development. Exiting."))

  myrank = MPI.Comm_rank(MPI.COMM_WORLD)
  fstdout = BufferedIO(STDOUT)
  if myrank == 0
    println(fstdout, "\nEntered Crank-Nicolson")
    println(fstdout, "res_tol = ", res_tol)
  end

  output_freq = opts["output_freq"]::Int
  write_vis = opts["write_vis"]::Bool
  use_itermax = opts["use_itermax"]::Bool
  jac_type = opts["jac_type"]
  if use_itermax
    itermax = opts["itermax"]
  end

  if jac_type == 4
    throw(ErrorException("CN not implemented for matrix-free ops. (jac_type cannot be 4)"))
  end

  if myrank == 0
    _f1 = open("convergence.dat", "a+")
    f1 = BufferedIO(_f1)
  end
 
  if neg_time == false
    # start time at 0.0
    t = 0.0
  else    # negative time for unsteady adjoint
    # start time at t_max
    t = t_max
  end

  # calculate t_steps, the number of time steps that CN will take
  t_steps = round(Int, t_max/h)

  # make a copy of the eqn object for storage of t_(n+1) information
  eqn_nextstep = deepcopy(eqn)
  # TODO: copyForMultistage does not give correct values.
  #     deepcopy works for now, but uses more memory than copyForMultistage, if it worked
#   eqn_nextstep = copyForMultistage(eqn)
  eqn_nextstep.q = reshape(eqn_nextstep.q_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
  eqn_nextstep.res = reshape(eqn_nextstep.res_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)

  # Initialize adjoint pde eqn objects
  if neg_time == true
    adj = deepcopy(eqn)
    adj_nextstep = deepcopy(eqn)
  end

  @debug1 println("============ In CN ============")

  #-------------------------------------------------------------------------------
  # allocate Jac outside of time-stepping loop
  println("===== neg_time: ", neg_time, " =====")
  if neg_time == false
    newton_data, jac, rhs_vec = setupNewton(mesh, mesh, sbp, eqn, opts, f)
  else
    newton_data, jac, rhs_vec = setupNewton(mesh, mesh, sbp, adj, opts, f)
  end

  # Setting IC for reverse sweep
  if neg_time == true

    #----------------
    # this section obtains q_vec at n
    i_actual = t_steps + 1  # index during forward sweep of the n'th q_vec. +1 instead of +3-i because the loop adds 2
    filename = string("qvec_for_adj-", i_actual, ".dat")
    q_vec_with_complex = readdlm(filename)
    eqn_dummy = deepcopy(adj)
    eqn_dummy.q_vec = q_vec_with_complex[:,1]

    # sync up eqn_dummy.q and eqn_dummy.q_vec
    disassembleSolution(mesh, sbp, eqn_dummy, opts, eqn_dummy.q, eqn_dummy.q_vec)

    # use startDataExchange to sync up q/q_vec and q_face_send/recv
    if opts["parallel_type"] == 2 && mesh.npeers > 0
      startDataExchange(mesh, opts, eqn_dummy.q, eqn_dummy.q_face_send, eqn_dummy.q_face_recv, eqn_dummy.params.f)
    end

    newton_data_discard, jac, rhs_vec_discard = setupNewton(mesh, mesh, sbp, eqn_dummy, opts, physics_func)

    # make sure we're doing complex step! since res_copy is zeros, it would mess up the FD calc
    assert(opts["jac_method"] == 2)
    epsilon = opts["epsilon"]
    pert = complex(0, epsilon)

    calcJacobianComplex(newton_data_discard, mesh, sbp, eqn_dummy, opts, physics_func, pert, jac, t)

    dRdu_n = jac    # dRdu_i: we don't need dRdu_(i-1), see derivation
    #----------------

    # need dJdu
    dJdu = calcdJdu_CS(mesh, sbp, eqn_dummy, opts)
    # need dRdu_n
    I = eye(length(eqn_dummy.q_vec))
    B = (I - (h/2) * (dRdu_n))^T
    psi = transpose(B)\(-dJdu)
    adj.q_vec = psi
  end

  for i = 2:(t_steps + 1)

    @debug1 println(eqn.params.f, "====== CN: at the top of time-stepping loop, t = $t, i = $i")
    @debug1 flush(eqn.params.f)

    #----------------------------
    # zero out Jac
    #   this works for both PETSc and Julia matrices.
    #   when jac is a Julia matrix, this is effectively wrapping: fill!(jac, 0.0)
    PetscMatZeroEntries(jac)

    # NOTE:
    # majorIterationCallback: called before every step of Newton's method. signature: 
    #   majorIterationCallback(itr, mesh::AbstractMesh, sbp::AbstractSBP, eqn::AbstractEulerData, opts)

    # NOTE: Must include a comma in the ctx tuple to indicate tuple
    # f is the physics function, like evalEuler

    #-------------
    # objective function section
    # 1. read option to indicate which obj fun
    # 2. call it, complex step it, and store it in dJdu
    if neg_time == true
      dJdu = zeros(Tsol, length(eqn.q_vec))
      #dJdu = calcdJdu(mesh, sbp, eqn, opts)
      # dJdu = calcObjectiveFn(mesh, sbp, adj, opts; isDeriv=true)
      J = calcObjectiveFn(mesh, sbp, adj, opts; isDeriv=false)
      dJdu = calcdJdu_CS(mesh, sbp, adj, opts)
    end

      


    if neg_time == false
      ctx_residual = (f, eqn, h, newton_data)
    else
      # i is the time step index in the reverse sweep.
      #   It moves forward from 2 to (t_steps+1) even though the adjoint is going backwards in time.
      #   The adjoint calculation requires data from the forward sweep to be loaded from disk.
      #   At a given time step in the adjoint solve, although the index i corresponds to this
      #     loop's time step, the index i does not correspond to the same i of the forward sweep.
      #   The adjustment is not just (t_steps - i) because the loop starts at 2 and ends at t_steps + 1.
      i_actual = t_steps + 3 - i
      ctx_residual = (f, adj, h, newton_data, i_actual, dJdu)
    end

    @debug1 println(fstdout, "in CN: before call to newtonInner")

    # time step update: h is passed in as argument to crank_nicolson
    if neg_time == false
      # need to add h in the forward time usage
      t_nextstep = t + h
    else
      # need to add h in the forward time usage
      t_nextstep = t - h
    end

    if neg_time == false
      # allow for user to select CN's internal Newton's method. Only supports dense FD Jacs, so only for debugging
      if opts["cleansheet_CN_newton"]
        # cnNewton: in cnNewton.jl
        cnNewton(mesh, sbp, opts, h, f, eqn, eqn_nextstep, t)
      else
        @time newtonInner(newton_data, mesh, sbp, eqn_nextstep, opts, cnRhs, cnJac, jac, rhs_vec, ctx_residual, t)
      end
    else      # call newtonInner using cnAdjJac and cnAdjRhs
      @time newtonInner(newton_data, mesh, sbp, adj_nextstep, opts, cnAdjRhs, cnAdjJac, jac, rhs_vec, ctx_residual, t)
      # TODO: contents of ctx_residual? needs to have adj instead of eqn
    end

    # This allows the solution to be updated from _nextstep without a deepcopy.
    #   There are two memory locations used by eqn & eqn_nextstep, 
    #   and this flips the location of eqn & eqn_nextstep every time step
    if neg_time == false
      eqn_temp = eqn
      eqn = eqn_nextstep
      eqn_nextstep = eqn_temp
    else
      adj_temp = adj
      adj = adj_nextstep
      adj_nextstep = adj_temp
    end

    # for adjoint_straight option: stores every time step's q to disk
    # TODO: decide between storing eqn or just eqn.q
    if store_u_to_disk == true
      filename = string("qvec_for_adj-", i, ".dat")
      writedlm(filename, eqn.q_vec)
    end

    # Note: we now need to copy the updated q over for the initial newton guess
    if neg_time == false
      for i = 1:mesh.numDof
        eqn_nextstep.q_vec[i] = eqn.q_vec[i]
      end
      disassembleSolution(mesh, sbp, eqn_nextstep, opts, eqn_nextstep.q, eqn_nextstep.q_vec)
    else
      for i = 1:mesh.numDof
        adj_nextstep.q_vec[i] = adj.q_vec[i]
      end
      disassembleSolution(mesh, sbp, adj_nextstep, opts, adj_nextstep.q, adj_nextstep.q_vec)
    end


    t = t_nextstep        # update time step

  end   # end of t step loop

  # depending on how many timesteps we do, this may or may not be necessary
  #   usage: copy!(dest, src)
  copy!(eqn, eqn_nextstep)      # copying eqn_nextstep to eqn
  writedlm("solution_final_inCN.dat", real(eqn.q_vec))


  if jac_type == 3      # if jac is a Petsc matrix, it needs to be freed when we're done using it
    # contents of ctx_newton: (jacp, x, b, ksp)
    NonlinearSolvers.destroyPetsc(jac, newton_data.ctx_newton...)
  end

  @debug1 println("============= end of CN: t = $t ===============")
  return t

end   # end of crank_nicolson function

function calcdJdu_CS{Tmsh, Tsol}(mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP, eqn::AbstractSolutionData{Tsol}, opts)

  # complex step it
  pert = complex(0, 1e-20)

  integrand_deriv = zeros(Tsol, length(eqn.q_vec))

  for i = 1:length(eqn.q_vec)
    eqn.q_vec[i] += pert

    J_arr = calcObjectiveFn(mesh, sbp, eqn, opts)
    J = J_arr[1]
    # println("=== in dJdu_CS: typeof(J_arr): ", typeof(J_arr))
    # println("=== in dJdu_CS: typeof(J): ", typeof(J))
    # println("=== in dJdu_CS: typeof(integrand_deriv): ", typeof(integrand_deriv))
    # println("=== in dJdu_CS: typeof(pert): ", typeof(pert))
    integrand_deriv[i] = imag(J)/norm(pert)
    eqn.q_vec[i] -= pert
  end

  return integrand_deriv

end

# function calcObjectiveFn{Tmsh, Tsol}(mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP, eqn::AbstractSolutionData{Tsol}, opts)
function calcObjectiveFn{Tmsh, Tsol}(mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP, eqn::AbstractSolutionData{Tsol}, opts; isDeriv=false)

  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
  if mesh.isDG
    boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
  end

  # TODO: get functional edges in a non BS way
  functional_edges = 3
  nDof = 1

  local_functional_val = zeros(Tsol, nDof)
  println("===dJdu=== size(local_functional_val): ", size(local_functional_val))
  println("===dJdu=== size(eqn.q_vec): ", size(eqn.q_vec))
  # eqn.q_vec: (96,)
  # local_functional_val: (96,)


  for itr = 1:length(functional_edges)
    g_edge_number = functional_edges[itr]
    itr2 = 0
    for itr2 = 1:mesh.numBC
      if findfirst(mesh.bndry_geo_nums[itr2], g_edge_number) > 0
        break
      end
    end

    start_index = mesh.bndry_offsets[itr2]
    end_index = mesh.bndry_offsets[itr2+1]
    idx_range = start_index:(end_index-1)
    bndry_facenums = sview(mesh.bndryfaces, idx_range)

    nfaces = length(bndry_facenums)

    integrand = zeros(Tsol, 1, mesh.sbpface.numnodes, nfaces)
    println("===dJdu=== size(integrand): ", size(integrand))
    println("===dJdu=== nfaces: ", nfaces)
    println("===dJdu=== mesh.sbpface.numnodes: ", mesh.sbpface.numnodes)

    for i = 1:nfaces
      bndry_i = bndry_facenums[i]
      global_facenum = idx_range[i]

      for j = 1:mesh.sbpface.numnodes
        #q = sview(eqn.q_bndry, :, j, global_facenum)
        q = eqn.q_bndry[:, j, global_facenum]
        # println("====== type of q: ", typeof(q))
        # println("====== size of q: ", size(q))
        # println("====== q: ", q)
        # convertToConservative(eqn.params, q, q2)

        # replaces calcBoundaryFunctionalIntegrand
        # integrand = zeros(Tsol, ndof, mesh.sbpface.numnodes, nfaces)    # dims?
        # integrand[1, j, i] = q.^2
        # TODO: figure out why [1]'s are required
        if isDeriv == false                 # calculates J = int(u^2)
          integrand[1, j, i] = q[1]*q[1]
        else                                # calculates dJdu = deriv(int(u^2)) = 2*u
          integrand[1, j, i] = 2*q[1]
        end

        # TODO: how to get the analytical derivative outside of integral


      end   # end of loop: j = 1:mesh.sbpfacenumnodes


      val_per_geom_edge = zeros(Tsol, 1)
      println("===dJdu=== size(val_per_geom_edge): ", size(val_per_geom_edge))

      # use integratefunctional, not boundaryintegrate: why?
      integratefunctional!(mesh.sbpface, mesh.bndryfaces[idx_range], integrand, val_per_geom_edge)
#       boundaryintegrate!(mesh.sbpface, mesh.bndryfaces[idx_range], integrand, val_per_geom_edge)

#       println(" size of val_per_geom_edge: ", size(val_per_geom_edge))

      local_functional_val[:] += val_per_geom_edge[:]


      # TODO:
      # serial: print out local_functional_val, compare with analytical
      # parallel: mpi all reduce then do the same


    end   # end of loop: i = 1:nfaces

  end   # end of loop: itr = 1:length(functional_edges)

  return local_functional_val

end

"""
  obj_zero

  Inputs: none
  Outputs: 0

  Zero-valued objective function
"""
function obj_zero()
  return 0.0
end
