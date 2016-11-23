# crank_nicolson.jl
# Crank-Nicolson implicit solver for PDEs

export crank_nicolson, cnResidual

push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/euler"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Debugging"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))

# DEBUG = false
DEBUG = true


@doc """
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
    * real_time : do actual time marching, not pseudo-time marching

   The eqn.q_vec should hold the whichever variables (conservative or
   entropy) that the simulation should use.

   The idea behind pre_func, post func, and ctx is that they abstract what kind
   of system rk4 is timestepping.  rk4 only needs to know about q_vec and
   res_vec.

   For physics modules, ctx should be (mesh, sbp, eqn) and q_vec and res_vec 
   should be eqn.q_vec and eqn.res_vec.
"""->
function crank_nicolson(f::Function, h::AbstractFloat, t_max::AbstractFloat,
                        mesh::AbstractMesh, sbp::AbstractSBP, eqn::AbstractSolutionData,
                        opts, res_tol=-1.0, real_time=true)

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

  if myrank == 0
    _f1 = open("convergence.dat", "a+")
    f1 = BufferedIO(_f1)
  end

  t = 0.0
  t_steps = round(Int, t_max/h)

  eqn_nextstep = deepcopy(eqn)
  # TODO TODO TODO: copyForMultistage does not give correct values.
  #     deepcopy works for now, but uses more memory than copyForMultistage, if it worked
#   eqn_nextstep = copyForMultistage(eqn)
  eqn_nextstep.q = reshape(eqn_nextstep.q_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
  eqn_nextstep.res = reshape(eqn_nextstep.res_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)

  # for the number of times eqn data is flipped btwn one or the other memory locations
  nflips_eqn = 0

  q_vec_old_DEBUG = zeros(eqn.q_vec)

  println("============ In CN ============")

  #-------------------------------------------------------------------------------
  # allocate Jac outside of time-stepping loop
  # NOTE 20161103: supplying eqn_nextstep does not work for x^2 + t^2 case, need to use eqn
  newton_data, jac, rhs_vec = setupNewton(mesh, mesh, sbp, eqn, opts)

  for i = 2:(t_steps + 1)

    println("CN: at the top of time-stepping loop, t = $t, i = $i")
    q_vec_old_DEBUG = deepcopy(eqn.q_vec)

    # zero out Jac
    fill!(jac, 0.0)

    # TODO: Allow for some kind of stage loop: ES-Dirk

    # TODO: output freq

    DEBUG = false
    if DEBUG
      q_file = "q$i.dat"
      writedlm(q_file, eqn.q_vec)
    end
    if DEBUG
      res_file = "res$i.dat"
      writedlm(res_file, eqn.res_vec)
    end

    # NOTE:
    # majorIterationCallback: called before every step of Newton's method
    # majorIterationCallback(itr, mesh::AbstractMesh, sbp::AbstractSBP, eqn::AbstractEulerData, opts)

#     if use_itermax && i > itermax
#       if myrank == 0
#         println(fstdout, "breaking due to itermax")
#         close(f1)
#         flush(fstdout)
#       end
#       break
#     end
#     println("mark3")

    # NOTE: Must include a comma in the ctx tuple to indicate tuple
    # f is the physics function, like evalEuler

    # NOTE: eqn_nextstep changed to eqn 20161013
#     ctx_residual = (f, eqn_nextstep, h, newton_data)
    ctx_residual = (f, eqn, h, newton_data)

    println(fstdout, "in CN: before call to newtonInner")

    t_nextstep = t + h

#     @time newtonInner(newton_data, mesh, sbp, eqn_nextstep, opts, cnRhs, cnJac, jac, rhs_vec, ctx_residual, t)
    cnNewton(mesh, h, f, eqn, eqn_nextstep, t)

    # This allows the solution to be updated from _nextstep without a deepcopy.
    # There are two memory locations used by eqn & eqn_nextstep, 
    #   and this flips the location of eqn & eqn_nextstep every time step
    nflips_eqn += 1
    eqn_temp = eqn
    eqn = eqn_nextstep
    eqn_nextstep = eqn_temp

    # 20161116
    # Note: we now need to copy the updated q over for the initial newton guess
    for i = 1:mesh.numDof
      eqn_nextstep.q_vec[i] = eqn.q_vec[i]
    end
    disassembleSolution(mesh, sbp, eqn_nextstep, opts, eqn_nextstep.q, eqn_nextstep.q_vec)

    delta_q_vec_DEBUG = eqn.q_vec - q_vec_old_DEBUG
    println("============+++++++++ norm(q_vec): ", norm(eqn.q_vec))
    println("============+++++++++ norm(delta_q_vec_DEBUG): ", norm(delta_q_vec_DEBUG))

#     if (nflips_eqn % 2) == 0
#       eqn = eqn1
#       eqn_nextstep = eqn2
#     else
#       eqn = eqn2
#       eqn_nextstep = eqn1
#     end

    # Old way
#     eqn = deepcopy(eqn_nextstep)
#     eqn.q = reshape(eqn.q_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
#     eqn.res = reshape(eqn.res_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
  
    t = t_nextstep

  end   # end of t step loop

  # depending on how many timesteps we do, this may or may not be necessary
  #   usage: copy!(dest, src)
  copy!(eqn, eqn_nextstep)
  writedlm("solution_final_inCN.dat", real(eqn.q_vec))


  if jac_type == 3
    # contents of ctx_newton: (jacp, x, b, ksp)
    NonlinearSolvers.destroyPetsc(jac, newton_data.ctx_newton...)
  end

  println("============= end of CN: t = $t ===============")
  return t

end   # end of crank_nicolson function

@doc """
###NonlinearSolvers.cnJac

  Jac of the CN calculation.
  Effectively a wrapper for physicsJac, because the CN Jac is:
    CN_Jac = I + dt/2 * physicsJac

  ctx:    
    physics func must be the first element, i.e. evalEuler
    eqn_nextstep must be the second element
    h must be the third element
    newton_data must be the fourth element
"""->
function cnJac(newton_data, mesh::AbstractMesh, sbp::AbstractSBP,
               eqn_nextstep::AbstractSolutionData, opts, jac, ctx, t)

  myrank = MPI.Comm_rank(MPI.COMM_WORLD)

  physics_func = ctx[1]
  # NOTE: eqn instead of eqn_nextstep, 20161013
  eqn = ctx[2]
#   eqn_nextstep = ctx[2]
  h = ctx[3]
  newton_data = ctx[4]

  t_nextstep = t + h

  # Orig CN: 
  #   call physicsJac with eqn_nextstep & t_nextstep
  #   then form CN_Jac = I + dt/2 * physics_Jac

  NonlinearSolvers.physicsJac(newton_data, mesh, sbp, eqn_nextstep, opts, jac, ctx, t_nextstep)
# TODO: CN dismantle 20161116 NonlinearSolvers.physicsJac(newton_data, mesh, sbp, eqn_nextstep, opts, jac, ctx, t)


  # Jacobian is always 2D
  scale!(jac, h*-0.5)

  # adding identity
  for i = 1:mesh.numDof

    jac[i,i] += 1

  end

#   println("==== in cnJac: h: $h ====")

#   writedlm("jac_inside_cnJac.dat", full(jac))

  return nothing

end

@doc """
###NonlinearSolvers.cnRhs

  RHS of the CN calculation

  ctx:    
    physics func must be the first element, i.e. evalEuler
    eqn_nextstep must be the second element
    h must be the third element

"""->
function cnRhs(mesh::AbstractMesh, sbp::AbstractSBP, eqn_nextstep::AbstractSolutionData, opts, rhs_vec, ctx, t)

  # eqn comes in through ctx_residual, which is set up in CN before the newtonInner call

  physics_func = ctx[1]
  # NOTE: changed to eqn 20161013
#   eqn_nextstep = ctx[2]
  eqn = ctx[2]
  h = ctx[3]

  t_nextstep = t + h

  physics_func(mesh, sbp, eqn, opts, t)
  println("---- physics_func: ",physics_func)
  assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)

# TODO: CN dismantle 20161116 physics_func(mesh, sbp, eqn_nextstep, opts, t)
# Orig CN: physics_func called with eqn_nextstep & t_nextstep
# dismantled: physics_func called with eqn_nextstep & t
  physics_func(mesh, sbp, eqn_nextstep, opts, t_nextstep)
  assembleSolution(mesh, sbp, eqn_nextstep, opts, eqn_nextstep.res, eqn_nextstep.res_vec)

  #   what this is doing:
  #   u_(n+1) - 0.5*dt* (del dot G_(n+1)) - u_n - 0.5*dt* (del dot G_n)
  for i = 1:mesh.numDof
# #     rhs_vec[i] = eqn_nextstep.q_vec[i] - 0.5*h*eqn_nextstep.res_vec[i] - eqn.q_vec[i] - 0.5*h*eqn.res_vec[i]

# TODO: CN dismantle 20161116 rhs_vec[i] = eqn_nextstep.res_vec[i]
    temp1 = eqn_nextstep.q_vec[i] - 0.5*h*eqn_nextstep.res_vec[i]
    temp2 = eqn.q_vec[i] + 0.5*h*eqn.res_vec[i]
# 
    rhs_vec[i] = temp1 - temp2 

    # NOTE: question: is there a sign problem here? should rhs_vec = -rhs_vec ?
    #     NO. this negative gets applied in newton.jl, where res_0[i] = -res_0[i]

    # local dof index, node, el
    loc_dof, node, el = ind2sub(mesh.coords, findfirst(mesh.dofs, i))

    DEBUG = false
    if DEBUG
      # mesh.coords indices: [dimension index, node num, el num]
      # 20161118: this threw an error in the vortex case???
      x = mesh.coords[1, node, el]
      y = mesh.coords[2, node, el]

      println("x: $x    y: $y    t: $t    t_nextstep: $t_nextstep")
      println("dof_ix: $i    rhs_vec[i]: ",rhs_vec[i])
      println("dof_ix: $i    eqn.q_vec[i]: ",eqn.q_vec[i])
      println("dof_ix: $i    eqn.res_vec[i]: ",eqn.res_vec[i])
      println("dof_ix: $i    eqn_nextstep.q_vec[i]: ",eqn_nextstep.q_vec[i])
      println("dof_ix: $i    eqn_nextstep.res_vec[i]: ",eqn_nextstep.res_vec[i])
      println("-")
    end

  end

#   println("==== in cnRhs: h: $h ====")

#   writedlm("rhs_inside_cnRhs.dat", full(rhs_vec))

  return nothing

end     # end of function cnRhs

@doc """
### NonlinearSolvers.pde_pre_func

  The pre-function for solving partial differential equations with a physics
  module.  The only operation it performs is disassembling eqn.q_vec into
  eqn.q

  Inputs:
    mesh
    sbp
    eqn
    opts
"""->
function pde_pre_func(mesh, sbp, eqn, opts)
  
  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
end


@doc """
### NonlinearSolvers.pde_post_func

  The post-function for solving partial differential equations with a physics
  module.  This function multiplies by A0inv, assembles eqn.res into
  eqn.res_vec, multiplies by the inverse mass matrix, and calculates
  the SBP approximation to the integral L2 norm

  Inputs:
    mesh
    sbp
    eqn
    opts

"""->
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

# the goal is to replace newton.jl.
# this will go into CN in the time-stepping loop
function cnNewton(mesh, h, physics_func, eqn, eqn_nextstep, t)

  println("---- physics_func: ",physics_func)

  # Jac on eqn or eqn_nextstep?

  epsilon = 1e-8

  jac = zeros(mesh.numDof, mesh.numDof)

  # emulates physicsJac
  # so we need to step through the jacobian column wise.
  #   d res[1]/d q[1]   d res[1]/d q[2]   d res[1]/d q[3] ...
  #   d res[2]/d q[1]   d res[2]/d q[2]   d res[2]/d q[3] ...
  #   d res[3]/d q[1]   d res[3]/d q[2]   d res[3]/d q[3] ...
  #   ...               ...               ...

  newton_itermax = 14
  delta_q_vec = zeros(eqn_nextstep.q_vec)

  # newton_loop starting here?
  for newton_i = 1:newton_itermax

    #--------------------------
    # emulates physicsJac
    unperturbed_q_vec = copy(eqn_nextstep.q_vec)
  
    physics_func(mesh, sbp, eqn_nextstep, opts, t_nextstep)
    # needed b/c physics_func only updates eqn.res
    assembleSolution(mesh, sbp, eqn_nextstep, opts, eqn_nextstep.res, eqn_nextstep.res_vec)
    unperturbed_res_vec = copy(eqn_nextstep.res_vec)
  
    for i = 1:mesh.numDof
      eqn_nextstep.q_vec[i] = eqn_nextstep.q_vec[i] + epsilon
  
      physics_func(mesh, sbp, eqn_nextstep, opts, t_nextstep)
      assembleSolution(mesh, sbp, eqn_nextstep, opts, eqn_nextstep.res, eqn_nextstep.res_vec)
  
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
    physics_func(mesh, sbp, eqn, opts, t)
    assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
    physics_func(mesh, sbp, eqn_nextstep, opts, t_nextstep)
    assembleSolution(mesh, sbp, eqn_nextstep, opts, eqn_nextstep.res, eqn_nextstep.res_vec)
  
    rhs_vec = zeros(eqn.q_vec)
  
    for i = 1:mesh.numDof
      rhs_vec[i] = eqn_nextstep.q_vec[i] - h*0.5*eqn_nextstep.res_vec[i] - eqn.q_vec[i] - h*0.5*eqn.res_vec[i]
    end
  
    # TODO: check these args
    rhs_norm = calcNorm(eqn, rhs_vec, strongres=true)
    
    #--------------------------
    # start of actual Newton
    neg_rhs = scale(rhs, -1.0)

    fill!(delta_q_vec, 0.0)
    # TODO: do these need to be column wise? why in newton.jl is it [:]?
    delta_q_vec = jac\neg_rhs
    fill!(jac, 0.0)

#     eqn_nextstep.q_vec += delta_q_vec

    for i = 1:mesh.numDof
      eqn_nextstep.q_vec[i] += delta_q_vec[i]
    end

    rhs_norm_tol = 1e-6
    if rhs_norm < rhs_norm_tol
      println("=== cnNewton converged with rhs_norm under $rhs_norm_tol ===")
      return nothing
    end


  end   # end of newton iterations


  println("=== cnNewton did not converge ===")
  return nothing


end


