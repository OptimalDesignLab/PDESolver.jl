arg_dict = Dict{Any, Any}(
  # "operator_type"   => "SBPGamma",
  "operator_type"   => "SBPOmega",
  # "operator_type"   => "SBPDiagE",
  "physics"         => "Elliptic",
  "use_DG"          => true,
  "delta_t"         => 1.0e-3,
  "t_max"           => 5.0,
  "res_tol"         => 1e-12,
  "step_tol"        => 1e-12,
  "itermax"         => 1,
  "order"           => 2,
  "Relfunc_name"    => true,
  "smb_name"        => "SRCMESHES/square16x16.smb",
  "dmg_name"        => ".null",
  "numBC"           => 1,
  "BC1"             => [0, 1, 2, 3],
  "exactSolution"   =>"ExactExpTrig",
  "BC1_name"        => "DirichletTrig",
  "Diffusion"       => "poly0th", 
  "SRC_name"        => "SrcExpTrigPoly0thDiffn",
  "Functional"      => "volumeAverage",
  "exactFunctional" => 1.846230857168755189823348538541e-02,
  "IC_name"         => "ICRandom",

  # "exactSolution"   =>"ExactTrig",
  # "BC1_name"        => "DirichletTrig",
  # "Diffusion"       => "poly2nd", # lambda = [x^2+1 & xy \\ xy & y^2+1]
  # "SRC_name"        => "SrcTrigPoly2ndDiffn",

  # "exactSolution"   =>"ExactTrig",
  # "BC1_name"        => "DirichletTrig",
  # "Diffusion"       => "poly6th",
  # "SRC_name"        => "SrcTrigPoly6thDiffn",

  # "exactSolution"   =>"ExactPoly2nd",
  # "Diffusion"       => "DiffnPoly0th"
  # "BC1_name"        => "DirichletPolynial2nd",
  # "SRC_name"        => "SrcPoly2nd",

  # "SRC_name"        => "SRC0",		# for energy stability test

  # "Flux_name"       => "SAT",			# SAT-BR2
  "Flux_name"       => "SAT0",		# SAT-SIPG
  # "Flux_name"       => "Shahbazi",	# Shahbazi
  "Cip"             => -3.0,
  "run_type"        => 5,
  "jac_type"        => 2,
  "jac_method"      => 2,
  "write_eigs"      => false,
  "write_jac"       => false,
  "write_conditioning" => false,
  "write_energy"    => true,
  "solve"           => true,
  "real_time"       => false
)


