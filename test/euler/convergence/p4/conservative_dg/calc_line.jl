# calculate the least squares line of the errors
function calc_line()
  dirs = [1, 2]
  npts = length(dirs)
  println("maximum(dirs) = ", maximum(dirs))
  err_vals = zeros(maximum(dirs))
  h_vals = zeros(maximum(dirs))
  for i in dirs
    println("i = ", i)
    fname = "m$i/error_calc.dat"
    fname2 = "m$i/counts_0.txt"
    val_i = readdlm(fname)
    err_vals[i] = val_i[1]

    counts_i = readdlm(fname2)
    h_vals[i] = 1/(sqrt(counts_i[3]/2))

  end

  min_idx = minimum(dirs)
  max_idx = maximum(dirs)
  err_vals = err_vals[min_idx:max_idx]
  h_vals = h_vals[min_idx:max_idx]
  npts = max_idx - min_idx + 1


  println("h values = ", h_vals)
  println("error values = ", err_vals)


  err_vals_log = zeros(npts)
  h_vals_log = zeros(npts)
  for i=1:npts
    err_vals_log[i] = log(err_vals[i])
    h_vals_log[i] = log(h_vals[i])
  end



  # now form system 

  A = ones(npts, 2)
  A[:, 2] = h_vals_log


  # solve system
  AAt = A.'*A
  rhs = A.'*err_vals_log

  x = AAt\rhs

  println("x = ", x)

  for i=1:npts
    h_i = h_vals_log[i]
    err_exp = x[1] + h_i*x[2]
    err_act = err_vals_log[i]
    err_pct = 100*(err_exp - err_act)/err_act
    println("point $i percent error of best fit line = ", err_pct, "%")
  end

  # calculate slope from last two points
  slope2 = (err_vals_log[npts - 1] - err_vals_log[npts])/(h_vals_log[npts - 1] - h_vals_log[npts])

  println("2 point slope: ", slope2)
  println("least squares slope: ", x[2])


  # write data to file for plotting
  data_arr = hcat(h_vals, err_vals)
  writedlm("err_data.dat", data_arr)

  return x[2]  # return least squares slope
end
