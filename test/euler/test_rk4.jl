# this file tests the rk4 function
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

pre_func = (a...) -> 3
function post_func(q_vec, res_vec, opts; calc_norm=true) 
  return norm(res_vec)
end

"""
  Use rk4 to solve an 4th order polynoial ODE, for which it is exact
"""
function test_rk4()
    facts("----- testing rk4 -----") do
    opts_rk = Dict{Any, Any}("run_type" => 1, "smb_name" => "", "output_freq" => 1, "write_vis" => false, "use_itermax" => false, "numBC" => 0)
    read_input(opts_rk)

    q_vec_rk = [1.0]
    res_vec_rk = [0.0]
    ctx = (q_vec_rk, res_vec_rk)

    delta_t = 0.1
    t_max = 1.0

    t = rk4(test_f, delta_t, t_max, q_vec_rk, res_vec_rk, pre_func, post_func, ctx, opts_rk, real_time=true)

    f_approx = q_vec_rk[1]
    f_exact = true_f(t)
    err = abs(f_approx - f_exact)

    @fact t --> roughly(t_max, atol=1e-14)
    @fact err --> less_than(1e-14)  # machine precision
  end  # end facts block

  facts("----- testing lserk54 -----") do
    opts_rk = Dict{Any, Any}("run_type" => 30, "smb_name" => "", "output_freq" => 1, "write_vis" => false, "use_itermax" => false, "numBC" => 0)
    read_input(opts_rk)

    q_vec_rk = [1.0]
    res_vec_rk = [0.0]
    ctx = (q_vec_rk, res_vec_rk)

    delta_t = 0.1
    t_max = 1.0

    t = lserk54(test_f, delta_t, t_max, q_vec_rk, res_vec_rk, pre_func, post_func, ctx, opts_rk, real_time=true)

    f_approx = q_vec_rk[1]
    f_exact = true_f(t)
    err = abs(f_approx - f_exact)

    @fact err --> less_than(1e-14)  # machine precision

  end

  return nothing
end

#test_rk4()
add_func1!(EulerTests, test_rk4, [TAG_NLSOLVERS, TAG_SHORTTEST])
