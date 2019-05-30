# this user supplied file creates a dictionary of arguments
# if a key is repeated, the last use of the key is used
# it is a little bit dangerous letting the user run arbitrary code
# as part of the solver
# now that this file is read inside a function, it is better encapsulated

arg_dict = Dict{Any, Any}(
# specify physcs, SBP operator, and solve type (steady Newton's method)
"physics" => "Euler",
"operator_type" => "SBPOmega",
"dimensions" => 2,
"run_type" => 5,
"jac_method" => 2,
"jac_type" => 2,
"order" => 1,
"use_DG" => true,
"Flux_name" => "RoeFlux",

# use implicit Euler globalization for Newton's method
"newton_globalize_euler" => true,
#"euler_tau" => 1e-1,  # the pseudo-time step can be specified here, but
#                      # the default is fine for this case
"real_time" => false,

# set up initial and boundary conditions
"IC_name" => "ICFreeStream",
"numBC" => 2,
"BC1" => [8],
"BC1_name" => "FreeStreamBC",
"BC2" => [5],
"BC2_name" => "noPenetrationBC",

# problem-specific parameters
"aoa" => 2.0,
"Ma" => 0.5,

# mesh and geometry
"smb_name" => "SRCMESHES/2D_airfoil.smb",
"dmg_name" => "SRCMESHES/2D_airfoil.dmg",

# solve options
"res_abstol" => 1e-10,
"res_reltol" => 1e-9,
"itermax" => 50,
"output_freq" => 1,
"do_postproc" => false,
)
