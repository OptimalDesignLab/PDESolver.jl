# calculate the least squares line of the errors

dirs = Array(1:3)  # do directories 1 through 3.  This range does not have
                   # to start at 1, but it must be contiguous
npts = length(dirs)
err_vals = zeros(maximum(dirs))
h_vals = zeros(maximum(dirs))
for i in dirs
  println("i = ", i)
  fname = "m$i/error_calc.dat"  # this file is written by the solver
  val_i = readdlm(fname)
  err_vals[i] = val_i[1]  # get error value
  h_vals[i] = val_i[2]  # get nominal mesh spacing h

end

min_idx = minimum(dirs)
max_idx = maximum(dirs)
err_vals = err_vals[min_idx:max_idx]
h_vals = h_vals[min_idx:max_idx]
npts = max_idx - min_idx + 1

# print the raw data
println("h  |   error")
println(hcat(h_vals, err_vals))

# tag the log of the data points
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

# compute error of each point compared to best fit line
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
writedlm("err_data.dat", data_arr, ' ')

