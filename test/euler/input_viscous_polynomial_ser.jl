arg_dict   = Dict{Any, Any}(

  # RK4 settings
  "run_type"           => 1,
  "delta_t"            => 0.01,
  # "t_max"              => 40.0,
  # "itermax"            => 10000,
  "t_max"              => 40.0,
  "itermax"            => 10,

  # because RK4 is selected and not steady Newton, we have to explicitly set parallel_data
  "parallel_data"      => "element",

  "operator_type"      => "SBPGamma",
  # "operator_type"     => "SBPOmega",
  # "operator_type"     => "SBDiagE",
  "use_DG"             => true,
  # "delta_t"            => 5.0e-4,
  # "t_max"              => 500.0,
  "res_abstol"         => 1e-13,
  "res_reltol"         => 1e-13,
  "step_tol"           => 1e-13,
  "krylov_itermax"     => 30,
  "order"              => 1,
  "do_postproc"        => true,
  "use_staggered_grid" => false,
  "variable_type"      => :conservative,
  "dimensions"         => 2,
  "use_stagger_grid"   => false,
  "dmg_name"           => ".null",

  # "smb_name" => "/users/yanj4/Downloads/meshfiles/channel.smb",
  # "smb_name" => "/users/yanj4/Downloads/meshfiles/channel_12x12.smb",
  # "smb_name" => "/users/yanj4/Downloads/meshfiles/square128x128.smb",
  # "smb_name" => "/users/yanj4/Downloads/meshfiles/square81x81.smb",
  # "smb_name" => "/users/yanj4/Downloads/meshfiles/square64x64.smb",
  # "smb_name" => "/users/yanj4/Downloads/meshfiles/square32x32.smb",
  # "smb_name" => "/users/yanj4/Downloads/meshfiles/square8x8.smb",
  # "smb_name" => "/users/yanj4/Downloads/meshfiles/square16x16.smb",
  # "smb_name" => "/fasttmp/ashlea/jf_meshes/square16x16.smb",
  # "smb_name" => "/fasttmp/ashlea/nssolver_tests/0aa_meshes/square2_aa_par2.smb",
  # "smb_name" => "/fasttmp/ashlea/nssolver_tests/0aa_meshes/square2_aa_ser.smb",
  "smb_name" => "/fasttmp/ashlea/nssolver_tests/0aa_meshes/square_pert/square02_pert_ser.smb",
 
  #
  # Channel MMS: nonperiodic
  #
  "isViscous"          => true,
  "use_src_term"       => true,
  "numBC"              => 2,
  "BC1"                => [0, 2],
  "BC2"                => [1, 3],
  "BC1_name"           => "FreeStreamBC",
  "BC2_name"           => "FreeStreamBC",
  # "IC_name"            => "ICFreeStream",
  "IC_name"            => "ICPolynomial",
  # "IC_name"            => "ICPolynomial_fluxdbg",
  "perturb_ic"         => false,
  "perturb_mag"        => 0.01,
  "SRCname"            => "SRCPolynomial",
  "exact_soln_func"    => "ICPolynomial",

 
  "Flux_name"          => "RoeFlux",
  "SAT_type"           => "Hartman",
  "Cip"                => 1.0,
  "Ma"                 => 0.1,
  "Re"                 => 1.0e3,
  "aoa"                => 0.0,
  "sideslip_angle"     => 0.0,
  "CFL"                => 0.1,
  "jac_type"           => 3,
  "jac_method"         => 2,
  "solve"              => true,
  "write_eigs"         => false,
  "write_jac"          => false,
  "real_time"          => false,

  "write_timing"       => true,
  "physics"            => "Euler",

  "output_freq"        => 1,

  "write_finalresidual" => true,
  "write_finalsolution" => true,
)
