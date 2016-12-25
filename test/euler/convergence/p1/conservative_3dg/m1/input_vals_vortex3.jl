# this user supplied file creates a dictionary of arguments
# if a key is repeated, the last use of the key is used
# it is a little bit dangerous letting the user run arbitrary code
# as part of the solver
# now that this file is read inside a function, it is better encapsulated

arg_dict = Dict{Any, Any} (
"var1" => 1,
"var2" => "a",
"var3" => 3.5,
"var4" => [1,2,3],
"var3" => 4,
"run_type" => 5,
"jac_method" => 2,
"jac_type" => 2,
"order" => 1,
"use_DG" => true,
"dimensions" => 3,
"Flux_name" => "RoeFlux",
"IC_name" => "ICExp",
"variable_type" => :conservative,
"SRCname" => "SRCExp",
#"SRCname" => "SRCExp",
#"IC_name" => "ICFile",
#"ICfname" => "start.dat",
#"Relfunc_name" => "ICRho1E2U3",
"numBC" => 1,
"BC1" => [ 0, 1, 2, 3, 4, 5],
"BC1_name" => "ExpBC",
#"BC2" => [4, 10],
#"BC2_name" => "noPenetrationBC",
#"BC2_name" => "ExpBC",
"delta_t" => 0.005,
"t_max" => 500.000,
"smb_name" => "SRCMESHES/cube_small.smb",
"dmg_name" => ".null",
"res_abstol" => 1e-8,
"res_reltol" => 1e-10,
"step_tol" => 1e-10,
"itermax" => 10,
"do_postproc" => true,
"exact_soln_func" => "ICExp",
"solve" => true,
)
