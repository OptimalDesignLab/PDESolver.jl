# this user supplied file creates a dictionary of arguments
# if a key is repeated, the last use of the key is used
# it is a little bit dangerous letting the user run arbitrary code
# as part of the solver
# now that this file is read inside a function, it is better encapsulated

arg_dict = Dict{String, Any}(
"physics" => "Euler",
"run_type" => 5,
"jac_type" => 2,
"jac_method" => 2,
"operator_type" => "SBPDiagonalE",
"dimensions" => 2,
"order" => 1,
"IC_name" => "ICIsentropicVortex",
"numBC" => 1,
"BC1" => [0,1,2,3],
"BC1_name" => "isentropicVortexBC",
"delta_t" => 0.005,
"t_max" => 500.000,
"smb_name" => "SRCMESHES/vortex_curve.smb",
"dmg_name" => ".null",
"res_abstol" => 1e-12,
"res_reltol" => -1.0,
"step_tol" => 1e-9,
"itermax" => 20,
"solve" => false,
"use_lps" => true,
# make this DG
"Flux_name" => "RoeFlux",
"use_DG" => true, 
"p_free" => 2.0, #2.0,
"T_free" => 50.0, #50,
"calc_jac_explicit" => false,
"force_mesh_complex" => true,
"need_adjoint" => true,
)
