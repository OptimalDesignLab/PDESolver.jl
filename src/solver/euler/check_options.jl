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
    if opts["volume_integral_type"] != 2
      error("cannot use staggered grids with non entropy stable volume integrals")
    end

    if opts["face_integral_type"] != 2
      error("cannot use staggered grids with non entropy stable face integrals")
    end

    if opts["operator_type2"] == "SBPNone"
      error("must specify a real SBP operator for staggered grid")
    end
  end



  return nothing
end
