# crank_nicolson_jacandrhs.jl
#
# Contains Jacobian and RHS calculation functions,
#   for both forward sweep CN equation and reverse sweep adjoint CN equation
#----------------------------------------------------------------------------------
export cnJac, cnRhs
export cnAdjJac, cnAdjRhs

# function cnJac(newton_data, mesh, sbp, eqn_nextstep, opts, jac, ctx, t)
# function cnAdjJac(newton_data, mesh, sbp, adj_nextstep, opts, eqn, ctx, t)
# TODO: why was eqn in the function signature instead of jac
function cnAdjJac(newton_data, mesh, sbp, adj_nextstep, opts, jac, ctx, t)
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
  #println(" === typeof eqn: ", typeof(eqn))
  # TODO: why is eqn Array{Float64,2}
  println(" === typeof adj: ", typeof(adj))
  newton_data, jac, rhs_vec = setupNewton(mesh, mesh, sbp, adj, opts, physics_func)
  # TODO: setupNewton shouldn't be here? because Newton can't properly iterate on this
  # get jacobian from eqn here
  NonlinearSolvers.physicsJac(newton_data, mesh, sbp, adj, opts, jac, ctx, t)

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
function cnJac(newton_data, mesh, sbp, eqn_nextstep, opts, jac, ctx, t)

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
function cnAdjRhs(mesh::AbstractMesh, sbp::AbstractSBP, adj_nextstep::AbstractSolutionData, opts, rhs_vec, ctx, t)

  # this doesn't work inside here
  # using PDESolver
  # using AdvectionEquationMod

  #push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/advection"))
  #importall AdvectionEquationMod

  physics_func = ctx[1]
  adj = ctx[2]
  h = ctx[3]
  newton_data = ctx[4]
  i_actual = ctx[5]
  dJdu = ctx[6]

  # TODO need to put actual dJdu here.
#   dJdu = zeros(mesh.numDof)

#   Tsol = opts["Tsol"]
#   Tres = opts["Tres"]
#   Tmsh = opts["Tmsh"]
#   Tdim = opts["Tdim"]

  # initialize dummy eqn object for jacobian calculation use
#   eqn_dummy = AdvectionData_{Tsol, Tres, Tdim, Tmsh}(mesh, sbp, opts)
#   eqn_dummy = AdvectionEquationMod.AdvectionData_{Tsol, Tres, Tdim, Tmsh}(mesh, sbp, opts)
  eqn_dummy = deepcopy(adj)

  # load needed q_vec checkpoint file into eqn_dummy
  filename = string("qvec_for_adj-", i_actual, ".dat")
#   eqn_dummy.q_vec = readdlm(filename)
  q_vec_with_complex = readdlm(filename)
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

  dRdu_i = jac

  t_nextstep = t - h    # adjoint going backwards in time

  println(" mesh.numDof: ", mesh.numDof)

  println(" size dJdu: ", size(dJdu))
  println(" size rhs_vec: ", size(rhs_vec))
  println(" size adj.q_vec: ", size(adj.q_vec))
  println(" size adj_nextstep.q_vec: ", size(adj_nextstep.q_vec))
  println(" size dRdu_i: ", size(dRdu_i))

  println(" typeof dJdu: ", typeof(dJdu))
  println(" typeof rhs_vec: ", typeof(rhs_vec))
  println(" typeof adj.q_vec: ", typeof(adj.q_vec))
  println(" typeof adj_nextstep.q_vec: ", typeof(adj_nextstep.q_vec))
  println(" typeof dRdu_i: ", typeof(dRdu_i))

  for i = 1:mesh.numDof

    # the Jacobian vector product needs to be jac * a vector, so use q_vec, not q in:
    #   dRdu_i*adj_nextstep.q_vec
#     rhs_vec[i] = dJdu[i] 
    rhs_vec[i] = dJdu[1]
                  + adj_nextstep.q_vec[i]
                  - 0.5*h*(dRdu_i*adj_nextstep.q_vec)[i] 
                  - adj.q_vec[i]
                  - 0.5*h*(dRdu_i*adj.q_vec)[i]

    # note: something about how I'm multiplying dRdu_i & q_vec

    # TODO: actually q_vec instead of adj_vec ?
    # NOTE: as written above: I think I'm using the q_vec field for psi in the unsteady adj eqn

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
