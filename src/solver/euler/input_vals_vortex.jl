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
"run_type" => 4,
"jac_type" => 1,
"order" => 1,
"IC_name" => "ICIsentropicVortex",
#"numBC" => 2,
"numBC" => 1,
#"BC1" => [ 7, 13],
"BC1" => [0],
"BC1_name" => "isentropicVortexBC",
#"BC2" => [4, 10],
#"BC2_name" => "noPenetrationBC",
#"BC2_name" => "isentropicVortexBC",
"delta_t" => 0.001, # default 0.005
"t_max" => 500000.000,
"smb_name" => "src/mesh_files/vortex3_1.smb",
#"dmg_name" => "src/mesh_files/vortex.dmg",
"dmg_name" => ".null",
"res_abstol" => 1e-10,
"res_tol" => 1e-10,  # Only For RK method
"res_reltol" => 1e-9,
"step_tol" => 1e-10,
"itermax" => 30,
"edgestab_gamma" => 0.0, # default = -0.9,
"writeq" => false,
#"perturb_ic" => true,
#"perturb_mag" => 0.001,
#"write_sparsity" => true,
#"write_jac" => true,
"write_edge_vertnums" => false,
"write_face_vertnums" => false,
"write_qic" => false,
"writeboundary" => false,
"write_res" => false,
"print_cond" => true,
#"write_counts" => true,
"solve" => true,
)
