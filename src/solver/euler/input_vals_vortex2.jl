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
"var3" => 7,
"run_type" => 1,
"order" => 1,
"IC_name" => "ICIsentropicVortex",
#"IC_name" => "ICFile",
#"ICfname" => "SL01.dat_old",
"numBC" => 1,
#"BC1" => [4, 10, 7, 13],
"BC1" => [ 0],
"BC1_name" => "isentropicVortexBC",
#"BC2" => [7, 13],
#"BC2_name" => "noPenetrationBC",
#"BC2_name" => "isentropicVortexBC",
"delta_t" => 0.005,
"t_max" => 5000.000,
"smb_name" => "../../mesh_files/quarter_vortex3l.smb",
"dmg_name" => ".null",
"res_tol" => 1e-5,
"step_tol" => 1e-9,
"itermax" => 20,
)
