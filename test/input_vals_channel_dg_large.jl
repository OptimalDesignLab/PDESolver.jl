# this user supplied file creates a dictionary of arguments
# if a key is repeated, the last use of the key is used
# it is a little bit dangerous letting the user run arbitrary code
# as part of the solver
# now that this file is read inside a function, it is better encapsulated

arg_dict = Dict{Any, Any}(
"var1" => 1,
"var2" => "a",
"var3" => 3.5,
"var4" => [1,2,3],
"var3" => 7,
"run_type" => 4,
"order" => 1,
"IC_name" => "ICRho1E2U3",
"use_DG" => true,
"Flux_name" => "RoeFlux",
"numBC" => 0,
#"numBC" => 1,
#"BC1" => [ 0 ],
#"BC1_name" => "Rho1E2U3BC",
"delta_t" => 0.005,
"t_max" => 500.000,
"smb_name" => "SRCMESHES/tri8l.smb",
#"smb_name" => "/users/creanj/fasttmp/entropy_test/channel/abc.smb",
"res_abstol" => 1e-12,
"res_reltol" => -1.0,
"step_tol" => 1e-9,
"itermax" => 20,
"write_face_vertnums" => true,
"solve" => false
)
