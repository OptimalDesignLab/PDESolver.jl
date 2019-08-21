# this user supplied file creates a dictionary of arguments
# if a key is repeated, the last use of the key is used
# it is a little bit dangerous letting the user run arbitrary code
# as part of the solver
# now that this file is read inside a function, it is better encapsulated

arg_dict = Dict{Any, Any}(
# specify physcs, SBP operator, and solve type (homotopy continuation)
"physics" => "Euler",
"operator_type" => "SBPDiagonalE",
"dimensions" => 2,
"run_type" => 40,  # use homotopy to solve the nonlinear problem
"jac_method" => 2,
"jac_type" => 2,
"order" => 1,
"use_DG" => true,
"use_lps" => true,  # use LPS stabilization
"Flux_name" => "HLLFlux",  # HLL Riemann solver used for interface flux

"real_time" => false,
"volume_integral_type" => 2,
"Volume_flux_name" => "IRFlux",  # entropy stable volume integrals

# set up initial and boundary conditions
"IC_name" => "ICFreeStream",
"numBC" => 2,
"BC1" => [36],
"BC1_name" => "FreeStreamBC",
"BC2" => [5],
"BC2_name" => "noPenetrationESBC",  # entropy stable Euler wall BC

# problem-specific parameters
"aoa" => 1.25,
"Ma" => 0.8,

# mesh and geometry
"smb_name" => "../meshes/airfoil40_mesh1_curve_.smb",
"dmg_name" => ".null",  # if using Simmetrix tools, can be airfoil40.smd
"addShockCapturing" => true,  # enable shock capturing
"shock_capturing_name" => "Volume",  # volume dissipation term
"shock_sensor_name" => "SensorPP",  # Persons and Peraire's shock sensor
"homotopy_addBoundaryIntegrals" => true,  # include boundary dissipation in
                                          # homotopy function
# solve options
"res_abstol" => 1e-10,
"res_reltol" => 1e-10,
"itermax" => 500,
"output_freq" => 1,
"do_postproc" => false,
)
