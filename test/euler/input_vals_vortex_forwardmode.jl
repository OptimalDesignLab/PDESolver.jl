# this user supplied file creates a dictionary of arguments
# if a key is repeated, the last use of the key is used
# it is a little bit dangerous letting the user run arbitrary code
# as part of the solver
# now that this file is read inside a function, it is better encapsulated
arg_dict = Dict{String, Any}(
"physics" => "Euler",
"run_type" => 11,
"jac_method" => 2,
"order" => 1,
"real_time" => true,
"variable_type" => :conservative,
"IC_name" => "ICIsentropicVortex",
"numBC" => 4,
"BC1" => [0],
"BC1_name" => "isentropicVortexBC",
"BC2" => [1],
"BC2_name" => "isentropicVortexBC",
"BC3" => [2],
"BC3_name" => "isentropicVortexBC",
"BC4" => [3],
"BC4_name" => "noPenetrationBC",
"aoa" => 0.0,
"smb_name" => "SRCMESHES/gvortex1.smb",
"dmg_name" => ".null",
"use_DG" => true,
"Flux_name" => "RoeFlux",
"t_max" => 5.0,
"res_abstol" => 1e-10,
"itermax" => 20,
"res_abstol" => 1e-10,
"res_reltol" => 1e-9,
"step_tol" => 1e-10,
"itermax" => 30, # 30

"num_functionals" => 1,
"functional_name1" => "lift",
"functional_bcs1" => [4],
# "analytical_functional_val" => -1/1.4,

# Objective Functional computation keys
"objective_function" => "drag",
"objective_bcs" => [4],

# Adjoint computation keys
"need_adjoint" => true,
"calc_adjoint" => false,
"writeq" => false,
"write_edge_vertnums" => false,
"write_face_vertnums" => false,
"write_qic" => false,
"writeboundary" => false,
"write_res" => false,
"print_cond" => false,
"write_counts" => false,
"write_vis" => false,
"output_freq" => 1,
"solve" => true, # This SHOULD be true, false for debugging
"do_postproc" => true,
"exact_soln_func" => "ICIsentropicVortex"
)
