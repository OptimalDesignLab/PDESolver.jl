arg_dict = Dict{Any, Any}(

  "solve"               => true,

  "dimensions"          => 2,
  "physics"             => "Euler",

  #--- perturbing Ma
  "perturb_Ma"          => false,
  # "perturb_Ma_magnitude" => 1e-100,
  # "stabilize_v"         => false,
  # "stabilization_method"  => "clipJacFast",

  "Ma"                  => 0.25 + 1e-8,
  "aoa"                 => 2.0,

  "operator_type"       => "SBPOmega",

  "run_type"            => 101,     # rk4_ds
  "jac_method"          => 2,
  "jac_type"            => 3,
  "real_time"           => true,

  "t_max"               => 2000.0,
  "delta_t"             => 0.00039,       # Must hard code delta_t, not use CFL, if testing FD Ma. Ma changes wave speed.
  "itermax"             => 3,

  "order"               => 1,
  "use_DG"              => true,
  "Flux_name"           => "RoeFlux",

  "variable_type"       => :conservative,
  "smb_name"            => "SRCMESHES/naca0012_3997el_ser.smb",

  "IC_name"             => "ICFreeStream",

  "numBC"               => 2,
  "BC1"                 => [5],
  "BC1_name"            => "noPenetrationBC",
  "BC2"                 => [8],
  "BC2_name"            => "FreeStreamBC",

  "exact_soln_func"     => "ICFreeStream",

  "res_tol"             => 1e-10,
  "res_abstol"          => 1e-10,
  "res_reltol"          => 1e-9,
  "step_tol"            => 1e-10,

  #--- data writing
  "write_timing"        => true,
  "write_finalsolution" => true,
  "write_finalresidual" => true,

  "writeq"              => false,
  "write_edge_vertnums" => false,
  "write_face_vertnums" => false,
  "write_qic"           => false,
  "writeboundary"       => false,
  "write_res"           => false,
  "print_cond"          => false,
  "write_counts"        => false,

  "do_postproc"         => true,
  "force_solution_complex" => true,

  "num_functionals" => 1,
  "functional_name1" => "drag",
  "functional_bcs1" => [1],

  #--- drag & L2vnorm writing
  "write_drag"          => true,
  "write_L2vnorm"       => false,

  "write_vis"           => true,
  "output_freq"         => 1,

  "use_checkpointing"   => true,
  "checkpoint_freq"     => 5000,

)
