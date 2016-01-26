# this user supplied file creates a dictionary of arguments
# if a key is repeated, the last use of the key is used
# it is a little bit dangerous letting the user run arbitrary code
# as part of the solver
# now that this file is read inside a function, it is better encapsulated

arg_dict = Dict{Any, Any} (
"run_type" => 1,
"order" => 1,
"IC_name" => "ICx5plusy5",
"numBC" => 1,
"BC1" => [0],
"BC1_name" => "x5plusy5BC",
#"BC2" => [4, 10],
#"BC2_name" => "noPenetrationBC",
"delta_t" => 0.005,
"t_max" => 5.000,
"smb_name" => "src/mesh_files/channel.smb",
#"dmg_name" => "../../mesh_files/vortex.dmg",
"dmg_name" => ".null",
"res_tol" => 1e-10,
# "step_tol" => 1e-6,
"itermax" => 20,
"res_abstol" => 1e-10,
"res_reltol" => 1e-9,
"step_tol" => 1e-10,
"itermax" => 30,
"writeq" => false,
"itermax" => 30,
"writeq" => false,
"write_edge_vertnums" => false,
"write_face_vertnums" => false,
"write_qic" => false,
"writeboundary" => false,
"write_res" => false,
"print_cond" => false,
"write_counts" => false,
"write_vis" => false,
"solve" => true,
)
