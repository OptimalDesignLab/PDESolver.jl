
# the goal is to replace newton.jl.
# this will go into CN in the time-stepping loop
function cnNewton(mesh, sbp, opts, h, physics_func, eqn, eqn_nextstep, t)
  println("++++++++++++++++ clean sheet Newton being run ++++++++++")

  println("---- physics_func: ",physics_func)

  # Jac on eqn or eqn_nextstep?

  epsilon = 1e-8
  t_nextstep = t + h

  jac = zeros(mesh.numDof, mesh.numDof)

  # emulates physicsJac
  # so we need to step through the jacobian column wise.
  #   d res[1]/d q[1]   d res[1]/d q[2]   d res[1]/d q[3] ...
  #   d res[2]/d q[1]   d res[2]/d q[2]   d res[2]/d q[3] ...
  #   d res[3]/d q[1]   d res[3]/d q[2]   d res[3]/d q[3] ...
  #   ...               ...               ...

  newton_itermax = 2
  delta_q_vec = zeros(eqn_nextstep.q_vec)

  # newton_loop starting here?
  for newton_i = 1:newton_itermax

    #--------------------------
    # emulates physicsJac
    unperturbed_q_vec = copy(eqn_nextstep.q_vec)

    physics_func(mesh, sbp, eqn_nextstep, opts, t_nextstep)
    # needed b/c physics_func only updates eqn.res
    assembleSolution(mesh, sbp, eqn_nextstep, opts, eqn_nextstep.res, eqn_nextstep.res_vec)
    # Comment here about mass matrix inv multiplication TODO
    applyMassMatrixInv(mesh, eqn_nextstep, eqn_nextstep.res_vec)
    unperturbed_res_vec = copy(eqn_nextstep.res_vec)

    for i = 1:mesh.numDof
      eqn_nextstep.q_vec[i] = eqn_nextstep.q_vec[i] + epsilon

      physics_func(mesh, sbp, eqn_nextstep, opts, t_nextstep)
      assembleSolution(mesh, sbp, eqn_nextstep, opts, eqn_nextstep.res, eqn_nextstep.res_vec)
      applyMassMatrixInv(mesh, eqn, eqn_nextstep.res_vec)

      jac[:,i] = (eqn_nextstep.res_vec - unperturbed_res_vec)/epsilon

      eqn_nextstep.q_vec[i] = unperturbed_q_vec[i]

    end

    #--------------------------
    # emulates cnJac
    scale!(jac, -0.5*h)
    for i = 1:mesh.numDof
      jac[i,i] += 1
    end

    #--------------------------
    # emulates cnRhs
    #   what this is doing:
    #   u_(n+1) - 0.5*dt* (del dot G_(n+1)) - u_n - 0.5*dt* (del dot G_n)
    #=
    physics_func(mesh, sbp, eqn, opts, t)
    assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
    physics_func(mesh, sbp, eqn_nextstep, opts, t_nextstep)
    assembleSolution(mesh, sbp, eqn_nextstep, opts, eqn_nextstep.res, eqn_nextstep.res_vec)

    rhs_vec = zeros(eqn.q_vec)

    for i = 1:mesh.numDof
      rhs_vec[i] = eqn_nextstep.q_vec[i] - h*0.5*eqn_nextstep.res_vec[i] - eqn.q_vec[i] - h*0.5*eqn.res_vec[i]
    end
    =#
    physics_func(mesh, sbp, eqn, opts, t)
    assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
    applyMassMatrixInv(mesh, eqn, eqn.res_vec)
    current_t_step_contribution = zeros(eqn.q_vec)
    for i = 1:mesh.numDof
      current_t_step_contribution[i] = - eqn.q_vec[i] - h*0.5*eqn.res_vec[i]
    end

    # Test for 3D Minv results
    # this works!
#     res_vec_control = deepcopy(eqn.res_vec)
#     res_vec_test = deepcopy(eqn.res_vec)
#     res_control = deepcopy(eqn.res)
#     res_test = deepcopy(eqn.res)
#     applyMassMatrixInv3D(mesh, sbp, eqn, res_test)
#     assembleSolution(mesh, sbp, eqn, opts, res_test, res_vec_test)
#     println("=+=+=+ norm of diff btwn res_vec_test & res_vec_control: ", norm(res_vec_test - res_vec_control))

    physics_func(mesh, sbp, eqn_nextstep, opts, t_nextstep)
    assembleSolution(mesh, sbp, eqn_nextstep, opts, eqn_nextstep.res, eqn_nextstep.res_vec)
    applyMassMatrixInv(mesh, eqn_nextstep, eqn_nextstep.res_vec)
    next_t_step_contribution = zeros(eqn.q_vec)
    for i = 1:mesh.numDof
      next_t_step_contribution[i] = eqn_nextstep.q_vec[i] - h*0.5*eqn_nextstep.res_vec[i] 
    end

    rhs_vec = zeros(eqn.q_vec)

    for i = 1:mesh.numDof
      rhs_vec[i] = current_t_step_contribution[i] + next_t_step_contribution[i]
    end

    # TODO: check these args
    rhs_norm = calcNorm(eqn, rhs_vec, strongres=true)

    #--------------------------
    # start of actual Newton
    neg_rhs = scale(rhs_vec, -1.0)

    fill!(delta_q_vec, 0.0)
    delta_q_vec = jac\neg_rhs
    fill!(jac, 0.0)

    for i = 1:mesh.numDof
      eqn_nextstep.q_vec[i] += delta_q_vec[i]
    end

    rhs_norm_tol = 1e-6
    if rhs_norm < rhs_norm_tol
      println("=== cnNewton converged with rhs_norm under $rhs_norm_tol -- newton iters: $newton_i ===")
      return nothing
    end

  end   # end of newton iterations

  println("=== cnNewton did not converge ===")
  return nothing


end

# TODO: comment here
function applyMassMatrixInv(mesh, eqn, vec)

  for k = 1:mesh.numDof
    vec[k] = eqn.Minv[k] * vec[k]
  end

  return vec
end

# TODO: comment here
function applyMassMatrixInv3D(mesh, sbp, eqn, arr)

  for i = 1:mesh.numEl
    for j = 1:sbp.numnodes
      for k = 1:mesh.numDofPerNode
        arr[k, j, i] = eqn.Minv3D[k, j, i] * arr[k, j, i]
      end
    end
  end

  return arr
end


