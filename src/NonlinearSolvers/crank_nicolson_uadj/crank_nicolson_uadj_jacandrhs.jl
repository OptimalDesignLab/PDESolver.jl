# crank_nicolson_jacandrhs.jl
#
# Contains Jacobian and RHS calculation functions,
#   for both forward sweep CN equation and reverse sweep adjoint CN equation
#----------------------------------------------------------------------------------
export cnJac_uadj, cnRhs_uadj
export cnAdjJac, cnAdjRhs

function cnAdjJac(newton_data, mesh, sbp, adj_nextstep, opts, jac, ctx, t)
  # adj_nextstep contains psi_i
  # ctx will pass in adj, which contains psi_(i+1)
  # eqn contains q at time i (check not i+1)

  myrank = MPI.Comm_rank(MPI.COMM_WORLD)

  physics_func = ctx[1]
  adj = ctx[2]
  h = ctx[3]
  newton_data = ctx[4]
  i_fwd = ctx[5]

  # Forming the CN Adjoint Jacobian:
  #   call physicsJac with eqn_nextstep & t_nextstep
  #   then form cnAdjJac = I - dt/2 * physics_Jac

  #println(" === typeof eqn: ", typeof(eqn))
  # TODO: why is eqn Array{Float64,2}
  @debug1 println(" === typeof adj: ", typeof(adj))
  # Note: setupNewton shouldn't be here, bc Newton can't properly iterate on it if the jac is cleared. 
  #       jac is passed in as arg to cnAdjJac

  # TODO: need to double check that t_nextstep is used here, not t. I believe it should be t_nextstep
  t_nextstep = t - h
  t_nextstep = negativeZeroCheck(t_nextstep)   # ensure negative zero is changed to zero

  # Fixed 20170330: physicsJac needs to be calculated at time i, which corresponds to adj_nextstep in the reverse sweep
  NonlinearSolvers.physicsJac(newton_data, mesh, sbp, adj_nextstep, opts, jac, ctx, t_nextstep)

  # need to flush assembly cache before performing the scale operation.
  #   These are defined for Julia matrices; they're just noops
  MatAssemblyBegin(jac, PETSc.MAT_FINAL_ASSEMBLY)
  MatAssemblyEnd(jac, PETSc.MAT_FINAL_ASSEMBLY)

  #--------------------------
  # applying dt/2 to jac
  # Jacobian is always 2D
  scale_factor = h*-0.5
  # make this into a petsc_scale_factor
  petsc_scale_factor = PetscScalar(scale_factor)

  # PetscMatScale is defined for all jac types, PETSc and Julia.
  #   When jac is julia array, this effectively does: scale!(jac, scale_factor)
  PetscMatScale(jac, petsc_scale_factor)

  # PetscMatAssembly___ not necessary here; PetscMatScale is provably local so it doesn't cache its stuff

  #--------------------------
  # adding identity
  ix_petsc_row = zeros(PetscInt, 1)
  ix_petsc_col = zeros(PetscInt, 1)
  value_to_add = zeros(PetscScalar, 1, 1)
  value_to_add[1,1] = 1.0
  flag = PETSc.ADD_VALUES

  for i = 1:mesh.numDof
    ix_petsc_row[1,1] = i + mesh.dof_offset       # jac diag index + local to global offset
    ix_petsc_col[1,1] = i + mesh.dof_offset       # jac diag index + local to global offset

    # PETSc function: set_values1!(Jac, [2], [2], Jac[2,2] + 1)
    # args: array, row (as an array of len 1), col (as an array of len 1), new values (as a 2D array of len 1)
    #   set_values1! has different methods for both PETSc matrices and Julia matrices
    #   for Julia dense & sparse arrays, set_values1! effectively does this: jac[i,i] += 1
    #   in serial, mesh.dof_offset is set to 0 automatically
    set_values1!(jac, ix_petsc_row, ix_petsc_col, value_to_add, flag)
  end

  # set_values1! only caches the results; need to be assembled. This happens in petscSolve in petsc_funcs.jl
  #   (The assemble funcs are defined for Julia matrices; they're just noops)

  # jac is now I - dt/2 * physics_jac

end     # end of function cnAdjJac


"""
NonlinearSolvers.cnJac_uadj

  Jac of the CN calculation.
  Effectively a wrapper for physicsJac, because the CN Jac is:
    CN_Jac = I - dt/2 * physicsJac

  ctx:    
    physics func must be the first element, i.e. evalEuler
    eqn_nextstep must be the second element
    h must be the third element
    newton_data must be the fourth element
"""
function cnJac_uadj(newton_data, mesh, sbp, eqn_nextstep, opts, jac, ctx, t)

  myrank = MPI.Comm_rank(MPI.COMM_WORLD)

  physics_func = ctx[1]
  eqn = ctx[2]
  h = ctx[3]
  newton_data = ctx[4]

  jac_type = opts["jac_type"]

  t_nextstep = t + h

  # Forming the CN Jacobian:
  #   call physicsJac with eqn_nextstep & t_nextstep
  #   then form CN_Jac = I - dt/2 * physics_Jac

  NonlinearSolvers.physicsJac(newton_data, mesh, sbp, eqn_nextstep, opts, jac, ctx, t_nextstep)

  # need to flush assembly cache before performing the scale operation.
  #   These are defined for Julia matrices; they're just noops
  PetscMatAssemblyBegin(jac, PETSc.MAT_FINAL_ASSEMBLY)
  PetscMatAssemblyEnd(jac, PETSc.MAT_FINAL_ASSEMBLY)

  #--------------------------
  # applying dt/2 to jac
  # Jacobian is always 2D
  scale_factor = h*-0.5
  # make this into a petsc_scale_factor
  petsc_scale_factor = PetscScalar(scale_factor)

  # PetscMatScale is defined for all jac types, PETSc and Julia
  # when jac is julia array, this effectively does: scale!(jac, scale_factor)
  PetscMatScale(jac, petsc_scale_factor)

  # PetscMatAssembly___ not necessary here; PetscMatScale is provably local so it doesn't cache its stuff

  #--------------------------
  # adding identity
  ix_petsc_row = zeros(PetscInt, 1)
  ix_petsc_col = zeros(PetscInt, 1)
  value_to_add = zeros(PetscScalar, 1, 1)
  value_to_add[1,1] = 1.0
  flag = PETSc.ADD_VALUES

  for i = 1:mesh.numDof
    ix_petsc_row[1,1] = i + mesh.dof_offset       # jac diag index + local to global offset
    ix_petsc_col[1,1] = i + mesh.dof_offset       # jac diag index + local to global offset

    # PETSc function: set_values1!(Jac, [2], [2], Jac[2,2] + 1)
    # args: array, row (as an array of len 1), col (as an array of len 1), new values (as a 2D array of len 1)
    #   set_values1! has different methods for both PETSc matrices and Julia matrices
    #   for Julia dense & sparse arrays, set_values1! effectively does this: jac[i,i] += 1
    #   in serial, mesh.dof_offset is set to 0 automatically
    set_values1!(jac, ix_petsc_row, ix_petsc_col, value_to_add, flag)
  end

  # set_values1! only caches the results; need to be assembled. This happens in petscSolve in petsc_funcs.jl
  #   (The assemble funcs are defined for Julia matrices; they're just noops)

  # jac is now I - dt/2 * physics_jac

  return nothing

end     # end of function cnJac_uadj

"""
NonlinearSolvers.cnAdjRhs

  RHS of the CN unsteady adjoint calculation

  cnAdjRhs counts on the fact that crank_nicolson above will call it stepping through time in the negative direction

"""
function cnAdjRhs(mesh::AbstractMesh, sbp::AbstractSBP, adj_nextstep::AbstractSolutionData, opts, rhs_vec, ctx, t)

  physics_func = ctx[1]
  adj = ctx[2]
  h = ctx[3]
  newton_data = ctx[4]
  i_fwd = ctx[5]
  dJdu = ctx[6]

  # TODO: need to double check that t_nextstep is used here, not t. I believe it should be t_nextstep
  t_nextstep = t - h    # adjoint going backwards in time
  t_nextstep = negativeZeroCheck(t_nextstep)   # ensure negative zero is changed to zero

  eqn_dummy = cnAdjLoadChkpt(mesh, sbp, opts, adj, physics_func, i_fwd, t_nextstep)
  jac = cnAdjCalcdRdu(mesh, sbp, opts, eqn_dummy, physics_func, t_nextstep)
  dRdu_i = transpose(jac)    # dRdu_i: we don't need dRdu_(i-1), see derivation

  # Fix 20170330: the dRdu_i used below in forming rhs_vec needs to be transposed
  dRdu_i = transpose(jac)    # dRdu_i: we don't need dRdu_(i-1), see derivation

  # println(" mesh.numDof: ", mesh.numDof)

  # println(" size dJdu: ", size(dJdu))
  # println(" size rhs_vec: ", size(rhs_vec))
  # println(" size adj.q_vec: ", size(adj.q_vec))
  # println(" size adj_nextstep.q_vec: ", size(adj_nextstep.q_vec))
  # println(" size dRdu_i: ", size(dRdu_i))

  # println(" typeof dJdu: ", typeof(dJdu))
  # println(" typeof rhs_vec: ", typeof(rhs_vec))
  # println(" typeof adj.q_vec: ", typeof(adj.q_vec))
  # println(" typeof adj_nextstep.q_vec: ", typeof(adj_nextstep.q_vec))
  # println(" typeof dRdu_i: ", typeof(dRdu_i))

  for i = 1:mesh.numDof

    # the Jacobian vector product needs to be jac * a vector, so use q_vec, not q in:
    #   dRdu_i*adj_nextstep.q_vec

    # adj_nextstep corresponds to psi_i
    # adj corresponds to psi_i+1

    rhs_vec[i] = dJdu[i] 
                  + adj_nextstep.q_vec[i]
                  - 0.5*h*(dRdu_i*adj_nextstep.q_vec)[i] 
                  - adj.q_vec[i]
                  - 0.5*h*(dRdu_i*adj.q_vec)[i]


    # note: something about how I'm multiplying dRdu_i & q_vec

    # Note: as written above: the q_vec field is being used for psi in the unsteady adj eqn

  end     # end of loop: i = 1:mesh.numDof


end     # end of function cnAdjRhs

"""
NonlinearSolvers.cnRhs_uadj

  RHS of the CN calculation

  ctx:    
    physics func must be the first element, i.e. evalEuler
    eqn_nextstep must be the second element
    h must be the third element

"""
# function cnRhs_uadj(mesh::AbstractMesh, sbp::AbstractSBP, eqn_nextstep::AbstractSolutionData, opts, rhs_vec, ctx, t)
function cnRhs_uadj(mesh, sbp, eqn_nextstep, opts, rhs_vec, ctx, t)

  # eqn comes in through ctx_residual, which is set up in CN before the newtonInner call

  physics_func = ctx[1]
  eqn = ctx[2]
  h = ctx[3]

  t_nextstep = t + h

  # assembleSolution needs to be called after physics_func TODO: explain why
  physics_func(mesh, sbp, eqn_nextstep, opts, t_nextstep)
  assembleSolution(mesh, sbp, eqn_nextstep, opts, eqn_nextstep.res, eqn_nextstep.res_vec)

  if opts["parallel_type"] == 2 && mesh.npeers > 0
    startDataExchange(mesh, opts, eqn.q, eqn.q_face_send, eqn.q_face_recv, eqn.params.f)
  end
  physics_func(mesh, sbp, eqn, opts, t)
  assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)

  #   what this is doing:
  #   u_(n+1) - 0.5*dt* (del dot G_(n+1)) - u_n - 0.5*dt* (del dot G_n)
  for i = 1:mesh.numDof

    temp1 = eqn_nextstep.q_vec[i] - 0.5*h*eqn_nextstep.res_vec[i]
    temp2 = eqn.q_vec[i] + 0.5*h*eqn.res_vec[i]

    rhs_vec[i] = temp1 - temp2 

    # NOTE: question: is there a sign problem here? should rhs_vec = -rhs_vec ?
    #     NO. this negative gets applied in newton.jl, where res_0[i] = -res_0[i]

  end

  return nothing

end     # end of function cnRhs_uadj
