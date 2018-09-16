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

  # this is a problem for nested physics
#  if opts["physics"] != PhysicsName
#    error("physics not specified as $PhysicsName, are you lost?")
#  end

  println("doing physics specific options checking")

  if opts["use_edge_res"]
    error("use_edge_res does not work correctly, do not use")
  end

  if opts["use_staggered_grid"]
    if opts["volume_integral_type"] != 2
      error("cannot use staggered grids with non entropy stable volume integrals")
    end
#=
    if opts["face_integral_type"] != 2
      error("cannot use staggered grids with non entropy stable face integrals")
    end
=#
    if opts["operator_type2"] == "SBPNone"
      error("must specify a real SBP operator for staggered grid")
    end
  end


  run_type = opts["run_type"]
  println("about to set calc_jac_explicit")
  #TODO: what to do for jac_type 4?
  if opts["use_DG"] && opts["face_integral_type"] == 1 && 
     opts["Flux_name"] in keys(FluxDict_diff) && 
     (run_type == 5 || run_type == 40 || run_type == 20) && opts["jac_type"] >= 2 #&& 
#     !opts["use_src_term"]
    println("setting calc_jac_explicit to true")
    get!(opts, "calc_jac_explicit", true)
  else
    println("setting calc_jac_explicit to false")
    get!(opts, "calc_jac_explicit", false)
    if opts["calc_jac_explicit"]
      error("cannot calculate jacobian explicitly with given options")
    end
  end

  if opts["calc_jac_explicit"]
    get!(opts, "preallocate_jacobian_coloring", false)
  else
    get!(opts, "preallocate_jacobian_coloring", true)
  end


  println("opts[calc_jac_explicit] = ", opts["calc_jac_explicit"])

  return nothing
end
