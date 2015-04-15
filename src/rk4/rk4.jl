# rk4.jl
# Runge Kutta 4th order solver for ODEs
# Anthony Ashley

# base RK4 method:
# dxdt = f(t,x)

# Inputs:
#   f:      function, that accepts input: (scalar t, vector x)
#   h:      delta t
#   x_ic:   initial condition for x
#   t_max:  length of time to step through
# Outputs:
#   x:      solved x at t_max

function rk4(f, h, x_ic, t_max)

  t = 0
  t_steps = int(floor(t_max/h))
  println("t_steps: ",t_steps)

#   x = Array(Float64,3,t_steps+2)
  x = Array(Float64,3,t_steps+1)

  iter = 1

  x[:,1] = x_ic


  while t < t_max

#     println("iter: ",iter)
    iter += 1

    x_old = x[:,iter-1]

    k1 = f(t,x_old)
    k2 = f(t + h/2, x_old + (h/2)*k1)
    k3 = f(t + h/2, x_old + (h/2)*k2)
    k4 = f(t + h, x_old + h*k3)

    x[:,iter] = x_old + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
#     println(x[:,iter])
    t = t + h

  end

  writedlm("rk4_output.dat",x,",")
#   writecsv("rk4_output.dat",x," ")
  return x

end
