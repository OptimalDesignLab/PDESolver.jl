# options checking specific to the Euler Equations

"""
  This function checks the options in the options dictionary after default
  values have been supplied and throws exceptions if unsupported/incompatable
  options are specified

  Inputs:
    opts: options dictionary

  Outputs:
    none
"""
function checkOptions(opts)

  if opts["physics"] != PhysicsName
    error("physics not specified as $PhysicsName, are you lost?")
  end


  if opts["use_edge_res"]
    error("use_edge_res does not work correctly, do not use")
  end

  if opts["volume_integral_type"] != 1
    integral_type = opts["volume_integral_type"]
    error("SimpleODE equation supports volume integral type 1 only, type $integral_type not supported")
  end

  if opts["face_integral_type"] != 1
    integral_type = opts["face_integral_type"]
    error("SimpleODE equation supports face integral type 1 only, type $integral_type not supported")
  end

  if opts["use_edgestab"]
    error("SimpleODE does not support edge stabilization")
  end

  if opts["use_filter"]
    error("SimpleODE does not support filter stabilization")
  end

  if opts["use_dissipation"] || opts["use_dissipation_prec"]
    error("SimpleODE does not support dissipation stabilization")
  end

  if opts["numBC"] != 0
    error("SimpleODE does not support boundary conditions")
  end


  get!(opts, "calc_jac_explicit", false)
  get!(opts, "preallocate_jacobian_coloring", true)

  return nothing
end
