
# rk4(

function lorenz(t, x_vect)

  sigma = 10
  beta = 8/3
  rho = 28

  x = x_vect[1]
  y = x_vect[2]
  z = x_vect[3]

  f = Array{Float64}(3)



  f[1] = sigma*(y - x)
  f[2] = x*(rho - z) - y
  f[3] = x*y - beta*z

  return f

end
