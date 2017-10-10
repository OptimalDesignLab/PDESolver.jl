# this user supplied file creates a dictionary of arguments
# if a key is repeated, the last use of the key is used
# it is a little bit dangerous letting the user run arbitrary code
# as part of the solver
# now that this file is read inside a function, it is better encapsulated

arg_dict = Dict{Any, Any}(
"physics" => "Advection",
"run_type" => 5,
# "jac_type" => 1,
"jac_method" => 2,
"order" => 1,
"real_time" => true,
"IC_name" => "ICxplusy",
"numBC" => 4,
"BC1" => [0],
"BC1_name" => "xplusyBC",
"BC2" => [1],
"BC2_name" => "xplusyBC",
"BC3" => [2],
"BC3_name" => "xplusyBC",
"BC4" => [3],
"BC4_name" => "xplusyBC",
"use_src_term" => true,
"SRCname" => "SRCxplusy",  # "SRC2",
"delta_t" => 0.002,
"t_max" => 2.0,
"smb_name" => "SRCMESHES/gsquare2.smb", # For the latest version of adjoint_capability
"dmg_name" => ".null",
"use_DG" => true,
"Flux_name" => "LFFlux",
"res_tol" => 1e-10,
# "step_tol" => 1e-6,
"itermax" => 20,
"res_abstol" => 1e-10,
"res_reltol" => 1e-9,
"step_tol" => 1e-10,
"itermax" => 30,
"writeq" => false,
"itermax" => 30,

# Compute functional value
# "calc_functional" => true,
"num_functionals" => 1,
"functional_name1" => "qflux",
# "functional_error" => true,
"geom_faces_functional1" => [1,2], # 0 based indexing for mesh edges
# "analytical_functional_val" => 3.0,

# Compute the objective function
"objective_function" => "qflux",
"geom_faces_objective" => [1,2],

# Compute adjoint vector
"calc_adjoint" => false,
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
"exact_soln_func" => "ICxplusy",
"write_face_vertnums" => false
)
