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
    * f  : function to call
    * h  : time step size
    * t_max : time value to stop time stepping (time starts at 0)
    * mesh : AbstractMesh
    * sbp : SBPOperator 
    * eqn : AbstractSolutionData
    * opts : options dictionary
    * res_tol : keyword arg, residual topping tolerance
"""->
function rk4(f, h::FloatingPoint, t_max::FloatingPoint, mesh, sbp, eqn, opts; res_tol = -1.0) 
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
   if iter % 100 == 0
     println("iter: ",i)
  end


    iter += 1
#    println("in rk4, iter = ", iter)
#    println("in rk4, t = ", t)

#    x_old = x[:,iter-1]

#    println("eqn.SL0 = ", eqn.SL0)

 #   eqn.SL0 = x_old
    eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.SL0)
    f( mesh, sbp, eqn, opts, t)

    eqn.SL[:] = 0.0
    eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.SL)
    k1[:] = eqn.SL
    x2[:] = x_old + (h/2)*k1

    sol_norm = norm(eqn.SL)/mesh.numDof

    if iter % 100 == 0
      println("writing to convergence.dat")
      write(f1, string(i, "   ", sol_norm, "\n"))
    end

    if iter % 1000 == 0
      println("flushing convergence.dat to disk")
#      close(f1)
#      f1 = open("convergence.dat", "a+")
      flush(f1)
    end

    if iter % 1000 == 0

      saveSolutionToMesh(mesh, SL0)
      writeVisFiles(mesh, "solution_rk$iter")
    end



    if (sol_norm < res_tol)
      println("breaking due to res_tol")
     # put solution into SL0
#     fill!(eqn.SL0, 0.0)
#     eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.SL0)

      # put residual into eqn.SL
#      eqn.SL[:] = res_0
 
      break
    end


    eqn.SL0[:] = x2
    eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.SL0)
    f( mesh, sbp, eqn, opts, t + h/2)

    fill!(eqn.SL, 0.0)
    eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.SL)
    k2[:] = eqn.SL
    x3[:] = x_old + (h/2)*k2

    eqn.SL0[:] = x3
    eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.SL0)
    f( mesh, sbp, eqn, opts, t + h/2)

    fill!(eqn.SL, 0.0)
    eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.SL)

 
    k3[:] = eqn.SL
    x4[:] = x_old + h*k3

    eqn.SL0[:] = x4

    eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.SL0)
    f( mesh, sbp, eqn, opts, t + h)

    fill!(eqn.SL, 0.0)
    eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.SL)
 

    k4 = eqn.SL[:]

    x_old[:] = x_old + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
    eqn.SL0[:] = x_old

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

  # put solution into SL0a
#  fill!(eqn.SL0, 0.0)
#  eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.SL0)

    # evaluate residual at final q value
    eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.SL0)
    f( mesh, sbp, eqn, opts, t)

    eqn.SL[:] = 0.0
    eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.SL)
 
  # put residual into eqn.SL
#  eqn.SL[:] = k1
 
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
