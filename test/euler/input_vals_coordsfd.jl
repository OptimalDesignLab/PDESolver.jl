# this user supplied file creates a dictionary of arguments
# if a key is repeated, the last use of the key is used
# it is a little bit dangerous letting the user run arbitrary code
# as part of the solver
# now that this file is read inside a function, it is better encapsulated

arg_dict = Dict{Any, Any}(
"physics" => "Euler",
"operator_type" => "SBPOmega",
"dimensions" => 2,
"run_type" => 5,
"jac_method" => 2,
"jac_type" => 2,
"order" => 1,
"IC_name" => "ICIsentropicVortex",
"use_DG" => true,
"Flux_name" => "RoeFlux",
"numBC" => 3,
"BC1" => [0],
"BC1_name" => "isentropicVortexBC",  # outlet
"BC2" => [2],
"BC2_name" => "isentropicVortexBC", # inlet
"BC3" => [1, 3],
"BC3_name" => "noPenetrationBC",  # was noPenetrationBC
"aoa" => 0.0,
#"SRCname" => "SRCp3",
#"delta_t" => 0.005,
"smb_name" => "SRCMESHES/vortex_3x3_.smb",
#"smb_name" => "small.smb",
"dmg_name" => ".null",
"itermax" => 20,
#"print_cond" => true,
#"write_eigs" => true,
"res_abstol" => 1e-9,
"res_reltol" => 1e-9,
#"step_tol" => 1e-10,
"write_vis" => true,
"write_adjoint_vis" => true,
"output_freq" => 1,
"do_postproc" => true,
"exact_soln_func" => "ICIsentropicVortex",
"force_solution_complex" => true,
"force_mesh_complex" => true,
"need_adjoint" => true,
)
