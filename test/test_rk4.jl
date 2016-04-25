# this file tests the rk4 function

#=
using ODLCommonTools
nl_solvers_path = joinpath(Pkg.dir("PDESolver", "src/NonlinearSolvers"))
println("nl_solvers_path = ", nl_solvers_path)
push!(LOAD_PATH, nl_solvers_path)
using NonlinearSolvers
=#

# this uses rk4 to integrate the function t^4 + t^3 + t^2 +t + 1
# thus the function supplied to rk4 is d/dt of that function
# and q_vec_rk[1] is initially function(t=0) = 1
function true_f(t)
#   return t + 1
  return t^4 + t^3 + t^2 + t + 1
end

function test_f(q_vec, res_vec, opts, t)
#   return res_vec[1] = 1
  res_vec[1] = 4*t^3 + 3*t^2 + 2*t + 1
end

opts_rk = Dict{Any, Any}("output_freq" => 1, "write_vis" => false, "use_itermax" => false)

q_vec_rk = [1.0]
res_vec_rk = [0.0]
ctx = (q_vec_rk, res_vec_rk)
pre_func = (a...) -> 3
post_func = (q_vec, res_vec, opts) -> norm(res_vec)
delta_t = 0.1
t_max = 1.0

t = rk4(test_f, delta_t, t_max, q_vec_rk, res_vec_rk, pre_func, post_func, ctx, opts_rk, real_time=true)

f_approx = q_vec_rk[1]
f_exact = true_f(t)
err = abs(f_approx - f_exact)
println("returned t = ", t)
println("at end of rk4, q_vec_rk[1] = ", f_approx)
println("true_f(t) = ", f_exact)
println("error = ", err)

@fact err --> less_than(1e-14)  # machine precision
