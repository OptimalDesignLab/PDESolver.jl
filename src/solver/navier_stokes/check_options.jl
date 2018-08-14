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

  if opts["use_staggered_grid"]
    error("staggered grid not supported for $PhysicsName")
  end

  get!(opts, "calc_jac_explicit", false)

  if opts["calc_jac_explicit"]
    get!(opts, "preallocate_jacobian_coloring", false)
  else
    get!(opts, "preallocate_jacobian_coloring", true)
  end

  get!(opts, "Cip", 1.0)
  get!(opts, "isViscous", true)

  if !opts["isViscous"]
    error("cannot use $PhysicsName with isViscous == false")
  end

  if !(opts["precompute_q_face"] && opts["precompute_face_flux"] &&
       opts["precompute_q_bndry"] && opts["precompute_boundary_flux"])

    error("precomputing face and boundary quantities required for $PhysicsName")
  end

  return nothing
end
