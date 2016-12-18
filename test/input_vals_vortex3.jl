# this user supplied file creates a dictionary of arguments
# if a key is repeated, the last use of the key is used
# it is a little bit dangerous letting the user run arbitrary code
# as part of the solver
# now that this file is read inside a function, it is better encapsulated

arg_dict = Dict{ASCIIString, Any}(
"var1" => 1,
"var2" => "a",
"var3" => 3.5,
"var4" => [1,2,3],
"var3" => 4,
"run_type" => 4,
"jac_type" => 2,
"order" => 1,
"IC_name" => "ICIsentropicVortex",
"numBC" => 1,
"BC1" => [0,1,2,3],
"BC1_name" => "isentropicVortexBC",
#"BC2_name" => "isentropicVortexBC",
"delta_t" => 0.005,
"t_max" => 10.000,
"smb_name" => "SRCMESHES/squarevortex_small.smb",
"res_abstol" => 1e-9,
"res_reltol" => -1.0,
"step_tol" => 1e-10,
"itermax" => 20,
"write_sparsity" => true,
"write_jac" => true,
"write_edge_vertnums" => true,
)
