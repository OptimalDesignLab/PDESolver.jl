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
"run_type" => 1,
"order" => 1,
# "IC_name" => "ICIsentropicVortex",
"numBC" => 1,
#"BC1" => [ 7, 13],
"BC1" => [0],
"BC1_name" => "flux1",
#"BC2" => [4, 10],
#"BC2_name" => "noPenetrationBC",
"delta_t" => 0.005,
"t_max" => 5.000,
"smb_name" => "src/mesh_files/Test1el.smb",
#"dmg_name" => "../../mesh_files/vortex.dmg",
"dmg_name" => ".null",
"res_tol" => 1e-10,
"step_tol" => 1e-6,
"itermax" => 20,
"writeCounts" => true,
)
