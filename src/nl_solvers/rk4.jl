# rk4.jl
# Runge Kutta 4th order solver for ODEs
# Anthony Ashley

# base RK4 method:
# dxdt = f(t,x)

# Inputs:
#   f:      function, that accepts input: (scalar t, vector x_old, vector x_new)
#   h:      delta t
#   x_ic:   initial condition for x
#   t_max:  length of time to step through
# Outputs:
#   x:      solved x at t_max

function rk4(f, h, x_new, x_ic, t_max, extra_args)

  t = 0.0
  t_steps = int(floor(t_max/h))
  println("t_steps: ",t_steps)

  (m,) = size(x_ic)
#   x = Array(Float64,3,t_steps+2)
#  x = Array(Float64,m,t_steps+1)

  iter = 1

#  x[:,1] = x_ic

  x_old = x_ic
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
#     println("iter: ",iter)
    iter += 1
#    println("in rk4, iter = ", iter)
#    println("in rk4, t = ", t)

#    x_old = x[:,iter-1]

    f(t,x_old, k1, extra_args)
    x2[:] = x_old + (h/2)*k1

    f(t + h/2, x2, k2, extra_args)
    x3[:] = x_old + (h/2)*k2

    f(t + h/2, x3, k3, extra_args)
    x4[:] = x_old + h*k3

    f(t + h, x4, k4, extra_args)
    
    x_old = x_old + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
    fill!(k1, 0.0)
    fill!(k2, 0.0)
    fill!(k3, 0.0)
    fill!(k4, 0.0)

#    x[:,iter] = x_old + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
#    println("==== RK4 ==== i: ",i)
#     println("x[:,iter]: ",x[:,iter])
#    println("k1: ",k1)
    t = t + h

  end

#  writedlm("rk4_output.dat",x,",")
#   writecsv("rk4_output.dat",x," ")
#  return x[:, t_steps+1], x
  return x_old

end


# base RK4 method:
# dxdt = f(t,x)

# Inputs:
#   f:      function, that accepts input: (scalar t, vector x_old, vector x_new)
#   h:      delta t
#   x_ic:   initial condition for x
#   t_max:  length of time to step through
# Outputs:
#   x:      solved x at t_max


function rk4(f::Function, h::FloatingPoint, t_max::FloatingPoint, mesh::AbstractMesh, sbp::SBPOperator, eqn::AbstractEquation; res_tol = -1.0) 
#function rk4(f, h, x_new, x_ic, t_max, extra_args)

# res_tol is alternative stopping criteria
  SL0 = eqn.SL0
  SL = eqn.SL
#  extra_args = (mesh, sbp, eqn)

  t = 0.0
  t_steps = round(Int, t_max/h)
  println("t_steps: ",t_steps)

  (m,) = size(SL0)
#   x = Array(Float64,3,t_steps+2)
#  x = Array(Float64,m,t_steps+1)

  iter = 1

  f1 = open("convergence.dat", "w")

#  x[:,1] = x_ic

#  x_old = SL0
  x_old = zeros(SL0)
  x_old[:] = SL0
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
     println("iter: ",i)
    iter += 1
#    println("in rk4, iter = ", iter)
#    println("in rk4, t = ", t)

#    x_old = x[:,iter-1]

#    println("eqn.SL0 = ", eqn.SL0)

 #   eqn.SL0 = x_old
    eqn.SL[:] = 0.0
    f( mesh, sbp, eqn, eqn.SL0, eqn.SL, t)

    k1[:] = eqn.SL
    x2[:] = x_old + (h/2)*k1

    sol_norm = norm(eqn.SL)
    if (sol_norm < res_tol)
      println("breaking due to res_tol")
      break
    end


    eqn.SL0[:] = x2
    eqn.SL[:] = k2
    f( mesh, sbp, eqn,  eqn.SL0, eqn.SL, t + h/2)

    k2[:] = eqn.SL
    x3[:] = x_old + (h/2)*k2

    eqn.SL0[:] = x3
    eqn.SL[:] = k3
    f( mesh, sbp, eqn,  eqn.SL0, eqn.SL, t + h/2)

    k3[:] = eqn.SL
    x4[:] = x_old + h*k3

    eqn.SL0[:] = x4
    eqn.SL[:] = k4
    f( mesh, sbp, eqn,  eqn.SL0, eqn.SL, t + h)
    k4 = eqn.SL[:]

    x_old[:] = x_old + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
    eqn.SL0[:] = x_old

    write(f1, string(i, "   ", norm(eqn.SL)/mesh.numDof, "\n"))

    fill!(k1, 0.0)
    fill!(k2, 0.0)
    fill!(k3, 0.0)
    fill!(k4, 0.0)


#    x[:,iter] = x_old + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
#    println("==== RK4 ==== i: ",i)
#     println("x[:,iter]: ",x[:,iter])
#    println("k1: ",k1)
    t = t + h

  end

  close(f1)
#=
  # final result needs to be returned in a different variable for AD
  println("coping x_old to SL")
  println("x_old = ", x_old)
  for i = 1:length(x_old)
    SL[i] = x_old[i]
  end
=#
#  println("eqn.SL = ", eqn.SL)
#  println("SL = ", SL)

#  writedlm("rk4_output.dat",x,",")
#   writecsv("rk4_output.dat",x," ")
#  return x[:, t_steps+1], x
  return nothing

end
