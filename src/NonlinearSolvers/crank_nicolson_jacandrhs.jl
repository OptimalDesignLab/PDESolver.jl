# crank_nicolson_jacandrhs.jl
#
# Contains Jacobian and RHS calculation functions,
#   for both forward sweep CN equation and reverse sweep adjoint CN equation
#----------------------------------------------------------------------------------

function cnAdjJac(newton_data, mesh, sbp, adj_nextstep, opts, eqn, ctx, t)
  # adj_nextstep contains psi_i
  # ctx will pass in adj, which contains psi_(i+1)
  # eqn contains q at time i (check not i+1)

  myrank = MPI.Comm_rank(MPI.COMM_WORLD)

  physics_func = ctx[1]
  adj = ctx[2]
  h = ctx[3]
  newton_data = ctx[4]
  i_actual = ctx[5]

  # Forming the CN Adjoint Jacobian:
  #   call physicsJac with eqn_nextstep & t_nextstep
  #   then form cnAdjJac = I - dt/2 * physics_Jac

  # instead of allocating another cnJac, modify this jac
  newton_data, jac, rhs_vec = setupNewton(mesh, mesh, sbp, eqn, opts, f)
  # get jacobian from eqn here
  NonlinearSolvers.physicsJac(newton_data, mesh, sbp, eqn, opts, jac, ctx, t)

  # need to flush assembly cache before performing the scale operation.
  #   These are defined for Julia matrices; they're just noops
  PetscMatAssemblyBegin(jac, PETSc.PETSC_MAT_FINAL_ASSEMBLY)
  PetscMatAssemblyEnd(jac, PETSc.PETSC_MAT_FINAL_ASSEMBLY)

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
  flag = PETSc.PETSC_ADD_VALUES

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

end


"""
NonlinearSolvers.cnJac

  Jac of the CN calculation.
  Effectively a wrapper for physicsJac, because the CN Jac is:
    CN_Jac = I + dt/2 * physicsJac

  ctx:    
    physics func must be the first element, i.e. evalEuler
    eqn_nextstep must be the second element
    h must be the third element
    newton_data must be the fourth element
"""
function cnJac(newton_data, mesh::AbstractMesh, sbp::AbstractSBP,
               eqn_nextstep::AbstractSolutionData, opts, jac, ctx, t)

  myrank = MPI.Comm_rank(MPI.COMM_WORLD)

  physics_func = ctx[1]
  eqn = ctx[2]
  h = ctx[3]
  newton_data = ctx[4]

  jac_type = opts["jac_type"]

  t_nextstep = t + h

  # Forming the CN Jacobian:
  #   call physicsJac with eqn_nextstep & t_nextstep
  #   then form CN_Jac = I + dt/2 * physics_Jac

  NonlinearSolvers.physicsJac(newton_data, mesh, sbp, eqn_nextstep, opts, jac, ctx, t_nextstep)

  # need to flush assembly cache before performing the scale operation.
  #   These are defined for Julia matrices; they're just noops
  PetscMatAssemblyBegin(jac, PETSc.PETSC_MAT_FINAL_ASSEMBLY)
  PetscMatAssemblyEnd(jac, PETSc.PETSC_MAT_FINAL_ASSEMBLY)

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
  flag = PETSc.PETSC_ADD_VALUES

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

  # jac is now I + dt/2 * physics_jac

  return nothing

end

"""
NonlinearSolvers.cnAdjRhs

  RHS of the CN unsteady adjoint calculation

  cnAdjRhs counts on the fact that crank_nicolson above will call it stepping through time in the negative direction


"""
function cnAdjRhs(mesh, sbp, adj_nextstep, opts, rhs_vec, ctx, t)

  physics_func = ctx[1]
  adj = ctx[2]
  h = ctx[3]
  newton_data = ctx[4]
  i_actual = ctx[5]
  dJdu = ctx[6]

  # TODO need to put actual dJdu here.
#   dJdu = zeros(mesh.numDof)

  filename = string("eqn_for_adj-", i_actual, ".dat")
  eqn = readdlm(filename)

  newton_data_discard, jac, rhs_vec_discard = setupNewton(mesh, mesh, sbp, eqn, opts, physics_func)

  # TODO: check newton_data
  #       actually check all args
  calcJacobianComplex(newton_data, mesh, sbp, eqn, ops, func, res_copy, pert, jac, t)

  t_nextstep = t - h    # adjoint going backwards in time

  for i = 1:mesh.numDof

    # TODO adj vs adj_vec in dRdu_i * adj.adj____
    rhs_vec[i] = dJdu[i] + adj_nextstep.adj_vec[i] - 0.5*h*dRdu_i*adj_nextstep.adj[i] - adj.adj_vec[i] - 0.5*h*dRdu_i*adj.adj[i]

  end


end

"""
NonlinearSolvers.cnRhs

  RHS of the CN calculation

  ctx:    
    physics func must be the first element, i.e. evalEuler
    eqn_nextstep must be the second element
    h must be the third element

"""
# function cnRhs(mesh::AbstractMesh, sbp::AbstractSBP, eqn_nextstep::AbstractSolutionData, opts, rhs_vec, ctx, t)
function cnRhs(mesh, sbp, eqn_nextstep, opts, rhs_vec, ctx, t)

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

end     # end of function cnRhs

"""
NonlinearSolvers.pde_pre_func

  The pre-function for solving partial differential equations with a physics
  module.  The only operation it performs is disassembling eqn.q_vec into
  eqn.q

  Inputs:
    mesh
    sbp
    eqn
    opts
"""
function pde_pre_func(mesh, sbp, eqn, opts)
  
  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
end


"""
NonlinearSolvers.pde_post_func

  The post-function for solving partial differential equations with a physics
  module.  This function multiplies by A0inv, assembles eqn.res into
  eqn.res_vec, multiplies by the inverse mass matrix, and calculates
  the SBP approximation to the integral L2 norm

  Inputs:
    mesh
    sbp
    eqn
    opts

"""
function pde_post_func(mesh, sbp, eqn, opts; calc_norm=true)
  eqn.multiplyA0inv(mesh, sbp, eqn, opts, eqn.res)
  eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
  for j=1:length(eqn.res_vec) eqn.res_vec[j] = eqn.Minv[j]*eqn.res_vec[j] end
  if calc_norm
    local_norm = calcNorm(eqn, eqn.res_vec)
    eqn.params.time.t_allreduce += @elapsed global_norm = MPI.Allreduce(local_norm*local_norm, MPI.SUM, mesh.comm)
    return sqrt(global_norm)
  end

   return nothing
end
