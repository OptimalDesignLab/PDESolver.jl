# this user supplied file creates a dictionary of arguments
# if a key is repeated, the last use of the key is used
# it is a little bit dangerous letting the user run arbitrary code
# as part of the solver
# now that this file is read inside a function, it is better encapsulated

arg_dict = Dict{String, Any}(
"physics" => "Euler",
"jac_method" => 2, # Default = 1,
"run_type" => 11,# 5,
"jac_type" => 2,
"order" => 1,
"use_DG" => true,
"dimensions" => 3,

"var1" => 1,
"var2" => "a",
"var3" => 3.5,
"var4" => [1,2,3],
"var3" => 4,
"Flux_name" => "RoeFlux",
"IC_name" => "ICExp",
"variable_type" => :conservative,
"numBC" => 3,
"BC1" => [ 0, 1, 2, 3],
"BC1_name" => "ExpBC",
"BC2" => [4],
"BC2_name" => "noPenetrationBC",
"BC3" => [10],
"BC3_name" => "noPenetrationBC",
#"BC2_name" => "isentropicVortexBC",
"delta_t" => 0.005,
"t_max" => 500.000,
"smb_name" => "SRCMESHES/tet8cube.smb",
"dmg_name" => ".null",

# Objective Functional computation keys
"objective_function" => "drag",
"geom_faces_objective" => [4],

# Adjoint computation keys
"need_adjoint" => true,

# Output options
"res_abstol" => 1e-12,
"res_reltol" => 1e-10,
"step_tol" => 1e-10,
"itermax" => 10000,
"use_edgestab" => false,
"use_GLS2" => false,
"tau_type" => 3,
"edgestab_gamma" => -0.9,
"use_filter" => false,
"use_dissipation" => false,
"dissipation_name" => "damp1",
"dissipation_const" => 12.00,
"write_finalsolution" => true,
"write_finalresidual" => true,
"write_counts" => true,
"write_rhs" => false,
"do_postproc" => true,
"exact_soln_func" => "ICIsentropicVortex",
"solve" => false,
)
