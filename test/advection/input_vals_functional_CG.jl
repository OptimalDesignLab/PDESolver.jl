# this user supplied file creates a dictionary of arguments
# if a key is repeated, the last use of the key is used
# it is a little bit dangerous letting the user run arbitrary code
# as part of the solver
# now that this file is read inside a function, it is better encapsulated

arg_dict = Dict{ASCIIString, Any}(
"physics" => "Advection",
"run_type" => 5,
"jac_method" => 2,
"order" => 1,
"real_time" => true,
"IC_name" => "ICexp_xy",
"numBC" => 4,
"BC1" => [0],
"BC1_name" => "exp_xyBC",
"BC2" => [1],
"BC2_name" => "exp_xyBC",
"BC3" => [2],
"BC3_name" => "exp_xyBC",
"BC4" => [3],
"BC4_name" => "exp_xyBC",
"use_src_term" => true,
"SRCname" => "SRCexp_xy",
#"BC2" => [4, 10],
#"BC2_name" => "noPenetrationBC",
"delta_t" => 0.002,
"t_max" => 2.0,
"smb_name" => "SRCMESHES/gsquare2.smb",
# "smb_name" => "SRCMESHES/tri8l.smb",
#"smb_name" => "SRCMESHES/channel_verylarge.smb",
#"dmg_name" => "../../mesh_files/vortex.dmg",
"dmg_name" => ".null",
"res_abstol" => 1e-10,
# "step_tol" => 1e-6,
"itermax" => 20,
"res_abstol" => 1e-10,
"res_reltol" => 1e-9,
"step_tol" => 1e-10,
"itermax" => 30,
"writeq" => false,
"itermax" => 30,
"calc_functional" => true,
"num_functionals" => 1,
"functional_name1" => "qflux",
"analytical_functional_val" => 2*(exp(1) - 1),
"geom_edges_functional1" => [2,3],
"writeq" => false,
"write_edge_vertnums" => false,
"write_qic" => false,
"writeboundary" => false,
"write_res" => false,
"print_cond" => false,
"write_counts" => false,
"write_vis" => true,
"output_freq" => 100,
"solve" => true,
"do_postproc" => true,
"exact_soln_func" => "ICexp_xy",
"write_face_vertnums" => false
)