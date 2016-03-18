# rk4.jl
# Runge Kutta 4th order solver for ODEs
# Anthony Ashley


export rk4

# base RK4 method:
# dxdt = f(t,x)



# Inputs:
#   f:      function, that accepts input: (scalar t, vector x_old, vector x_new)
#   h:      delta t
#   x_ic:   initial condition for x
#   t_max:  length of time to step through
# Outputs:
#   x:      solved x at t_max

@doc """
rk4

  This function does 4th order Runge Kutta time stepping

  Arguments:
    * f  : function evaluation, must have signature (ctx..., opts, t), must have signature (ctx..., opts, t)
    * h  : time step size
    * t_max : time value to stop time stepping (time starts at 0)
    * mesh : AbstractMesh
    * sbp : AbstractSBP 
    * eqn : AbstractSolutionData
    * opts : options dictionary
    * res_tol : keyword arg, residual topping tolerance
    * real_time : do actual time marching, not pseudo-time marching

   The eqn.q_vec should hold the whichever variables (conservative or
   entropy) that the simulation should use.
"""->
function rk4(f::Function, h::AbstractFloat, t_max::AbstractFloat, q_vec::AbstractVector, res_vec::AbstractVector, pre_func, post_func, ctx, opts; majorIterationCallback=((a...) -> (a...)), res_tol = -1.0, real_time=false)
#function rk4(f::Function, h::AbstractFloat, t_max::AbstractFloat, mesh::AbstractMesh, sbp, eqn::AbstractSolutionData, opts; res_tol = -1.0, real_time=false) 
#function rk4(f, h, x_new, x_ic, t_max, extra_args)

  println("\nEntered rk4")
# res_tol is alternative stopping criteria

  # unpack options
  output_freq = opts["output_freq"]::Int
  write_vis = opts["write_vis"]::Bool

#  q_vec = eqn.q_vec
#  res_vec = eqn.res_vec
#  extra_args = (mesh, sbp, eqn)

  t = 0.0  # timestepper time
  treal = 0.0  # real time (as opposed to pseudo-time)
  t_steps = round(Int, t_max/h)
  println("t_steps: ",t_steps)

  (m,) = size(q_vec)
#   x = Array(Float64,3,t_steps+2)
#  x = Array(Float64,m,t_steps+1)

  f1 = open("convergence.dat", "a+")

#  x[:,1] = x_ic

#  x_old = q_vec
  x_old = zeros(q_vec)
  x_old[:] = q_vec
  k1 = zeros(x_old)
  k2 = zeros(x_old)
  k3 = zeros(x_old)
  k4 = zeros(x_old)

  x2 = zeros(x_old)
  x3 = zeros(x_old)
  x4 = zeros(x_old)


  for i=2:(t_steps + 1)


#    update_msg = string("RK4 i: ",i,"\n")
#    write(STDERR,update_msg)
#    print("\nRK4 i : ", i)
  if i % output_freq == 0
     println("\ntimestep ",i)
  end

#    println("in rk4, i = ", i)
#    println("in rk4, t = ", t)

#    x_old = x[:,i-1]

#    println("eqn.q_vec = ", eqn.q_vec)

#     println("\nrk4 stage 1")
 #   eqn.q_vec = x_old
#    eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    pre_func(ctx..., opts)
    if real_time treal = t end
    f( ctx..., opts, treal)
    sol_norm = post_func(ctx..., opts)
#    fill!(eqn.res_vec, 0.0)
#    eqn.multiplyA0inv(mesh, sbp, eqn, opts, eqn.res)
#    eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
#    for j=1:length(eqn.res_vec) k1[j] = eqn.Minv[j]*eqn.res_vec[j] end
    for j=1:m
      k1[j] = res_vec[j]
    end
    x2[:] = x_old + (h/2)*k1  # because the Euler equations are written 
                              # in conservative variables, the residual
			      # is r(q) not r(v)


#    sol_norm = norm(eqn.res_vec)/mesh.numDof
     majorIterationCallback(i, ctx..., opts)
   
#     sol_norm = calcNorm(eqn, k1)
    if i % 1 == 0
#=
      # Calculate the error in density
      vRho_calc = zeros(vRho_act)
      k = 1
      for i = 1:4:length(eqn.q_vec)
        vRho_calc[k] = eqn.q_vec[i]
        k += 1
      end

      ErrDensity = norm(vRho_calc - vRho_act)/mesh.numNodes 
      println("DensityErrorNorm = ", ErrDensity)
      println("Solution Residual Norm = ", sol_norm)
=#
#      println("writing to convergence.dat")
      println(f1, i, " ", sol_norm)
    end
    
    if i % output_freq == 0
      println("flushing convergence.dat to disk")
#      close(f1)
#      f1 = open("convergence.dat", "a+")
      flush(f1)
    end



    if (sol_norm < res_tol)
      println("breaking due to res_tol")
      flush(f1)
     # put solution into q_vec
#     fill!(eqn.q_vec, 0.0)
#     eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)

      # put residual into eqn.res_vec
#      eqn.res_vec[:] = res_0
 
      break
    end

    # stage 2
#    println("\nrk4 stage 2")
    q_vec[:] = x2
    pre_func(ctx..., opts) 
#    eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    # convert to entropy variables here
    if real_time  treal = t + h/2 end
    f( ctx..., opts, treal)
    post_func(ctx..., opts)
#    fill!(eqn.res_vec, 0.0)
#    eqn.multiplyA0inv(mesh, sbp, eqn, opts, eqn.res)
#    eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
#    k2[:] = eqn.res_vec
#    for j=1:length(eqn.res_vec) k2[j] = eqn.Minv[j]*eqn.res_vec[j] end
    for j=1:m
      k2[j] = res_vec[j]
    end

    x3[:] = x_old + (h/2)*k2

    # stage 3
#    println("\nrk4 stage 3")
    q_vec[:] = x3
    pre_func(ctx..., opts)
    if real_time treal= t + h/2 end
#    eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    f( ctx..., opts, treal)
    post_func(ctx..., opts)

#    fill!(eqn.res_vec, 0.0)
#    eqn.multiplyA0inv(mesh, sbp, eqn, opts, eqn.res)
#    eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
#    k3[:] = eqn.res_vec
#    for j=1:length(eqn.res_vec) k3[j] = eqn.Minv[j]*eqn.res_vec[j] end
    for j=1:m
      k3[j] = res_vec[j]
    end

    # stage 4
#    println("\nrk4 stage 4")
    x4[:] = x_old + h*k3
    q_vec[:] = x4
    pre_func(ctx..., opts)
#   eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    if real_time treal = t + h end
    f( ctx..., opts, treal)
    post_func(ctx..., opts)
#    fill!(eqn.res_vec, 0.0)
#    eqn.multiplyA0inv(mesh, sbp, eqn, opts, eqn.res)
#    eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
#    for j=1:length(eqn.res_vec) k4[j] = eqn.Minv[j]*eqn.res_vec[j] end
#    k4 = eqn.res_vec[:]
    for j=1:m
      k4[j] = res_vec[j]
    end
#=
    println("k1 = \n", k1)
    println("k2 = \n", k2)
    println("k3 = \n", k3)
    println("k4 = \n", k4)

    println("q_old = \n", x_old)
=#
    # update
    x_old[:] = x_old + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
    q_vec[:] = x_old


#    println("q_new = \n", q_vec)

    fill!(k1, 0.0)
    fill!(k2, 0.0)
    fill!(k3, 0.0)
    fill!(k4, 0.0)


#    x[:,i] = x_old + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
#    println("==== RK4 ==== i: ",i)
#     println("x[:,i]: ",x[:,i])
#    println("k1: ",k1)
    t = t + h

  end

  close(f1)

  # put solution into q_veca
#  fill!(eqn.q_vec, 0.0)
#  eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)


  # put residual into eqn.res_vec
#  eqn.res_vec[:] = k1
 
#=
  # final result needs to be returned in a different variable for AD
  println("coping x_old to res_vec")
  println("x_old = ", x_old)
  for i = 1:length(x_old)
    res_vec[i] = x_old[i]
  end
=#
#  println("eqn.res_vec = ", eqn.res_vec)
#  println("res_vec = ", res_vec)

#  writedlm("rk4_output.dat",x,",")
#   writecsv("rk4_output.dat",x," ")
#  return x[:, t_steps+1], x
  return t

end

# this is the version for solving PDEs
# it uses the pde_pre_func and pde_post_func below
function rk4(f::Function, h::AbstractFloat, t_max::AbstractFloat, q_vec::AbstractVector, res_vec::AbstractVector, ctx, opts; majorIterationCallback=((a...) -> (a...)), res_tol=-1.0, real_time=false)
    rk4(f::Function, h::AbstractFloat, t_max::AbstractFloat, q_vec::AbstractVector, res_vec::AbstractVector, pde_pre_func, pde_post_func, ctx, opts; majorIterationCallback=majorIterationCallback, res_tol =res_tol, real_time=real_time)

end

function rk4(f::Function, h::AbstractFloat, t_max::AbstractFloat, mesh, sbp, eqn, opts; res_tol=-1.0, real_time=false)

  rk4(f, h, t_max, eqn.q_vec, eqn.res_vec, pde_pre_func, pde_post_func, (mesh, sbp, eqn), opts; majorIterationCallback=eqn.majorIterationCallback, res_tol=res_tol, real_time=real_time)
end


function pde_pre_func(mesh, sbp, eqn, opts)
  
  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
end

function pde_post_func(mesh, sbp, eqn, opts)
  eqn.multiplyA0inv(mesh, sbp, eqn, opts, eqn.res)
  eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
  for j=1:length(eqn.res_vec) eqn.res_vec[j] = eqn.Minv[j]*eqn.res_vec[j] end
  return calcNorm(eqn, eqn.res_vec)
end


