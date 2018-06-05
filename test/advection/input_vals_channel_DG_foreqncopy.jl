# this user supplied file creates a dictionary of arguments
# if a key is repeated, the last use of the key is used
# it is a little bit dangerous letting the user run arbitrary code
# as part of the solver
# now that this file is read inside a function, it is better encapsulated

#---------------------------------------------------------------------
# NOTE: Only for eqn_copy testing!!!
#   the eqn_deepcopy tests randomly change values in the eqn object,
#   so this input file cannot be used for any other tests
#---------------------------------------------------------------------

arg_dict = Dict{String, Any}(
"physics" => "Advection",
"run_type" => 1,
"order" => 1,
"use_DG" => true,
"Flux_name" => "LFFlux",
"IC_name" => "ICx5plusy5",
#"IC_name" => "ICFile",
#"ICfname" => "q_vec1.dat_old",
"numBC" => 1,
#"BC1" => [4, 10, 7, 13],
"BC1" => [ 0 ],
"BC1_name" => "x5plusy5BC",
#"BC2" => [7, 13],
#"BC2_name" => "noPenetrationBC",
#"BC2_name" => "Rho1E2U3BC",
"delta_t" => 0.005,
"t_max" => 500.000,
#"smb_name" => "../../mesh_files/channel.smb",
#"dmg_name" => "../../mesh_files/channel.dmg",
"smb_name" => "SRCMESHES/channel_large.smb",
#"dmg_name" => ".null",
"res_abstol" => 1e-12,
"res_reltol" => -1.0,
"step_tol" => 1e-9,
"itermax" => 20,
"write_face_vertnums" => true,
"write_edge_vertnums" => true,
"solve" => false
)
