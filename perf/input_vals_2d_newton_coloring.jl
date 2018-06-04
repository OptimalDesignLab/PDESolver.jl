# run for 60 Newton iterations
# force the use of coloring for the Jacobian

arg_dict = Dict{Any, Any}(
"physics" => "Euler",
"run_type" => 5,
"jac_type" => 2,
"jac_method" => 2,
"order" => 1,
"use_DG" => true,
"Flux_name" => "RoeFlux",
"IC_name" => "ICIsentropicVortex",
"numBC" => 1,
"BC1" => [ 0, 1, 2, 3],
"BC1_name" => "isentropicVortexBC",
"delta_t" => 0.001,
"t_max" => 50.000,
"smb_name" => "SRCMESHES/square_benchmarksmall.smb",
"dmg_name" => ".null",
"res_abstol" => -1e-12,  # do itermax steps
"res_reltol" => -1e-10,
"step_tol" => -1e-10,
"itermax" => 40,
"output_freq" => 10,
"write_timing" => true,
"calc_jac_explicit" => false,
)
