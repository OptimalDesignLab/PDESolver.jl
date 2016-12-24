arg_dict = Dict{Any, Any}(
"dimensions" => 3,
"use_DG" => true,
"run_type" => 1,
"real_time" => true,
"order" => 1,
"Flux_name" => "avgFlux",
"IC_name" => "ICsinwavexy",
"numBC" => 0,
"delta_t" => 0.0005,
"t_max" => 1.375,
"smb_name" => "SRCMESHES/tet5_periodic_p4.smb",
"dmg_name" => ".null",
"output_freq" => 100,
"write_vis" => false,
"do_postproc" => true,
"exact_soln_func" => "ICsinwavexy",
"solve" => true
)


