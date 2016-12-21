# run tests related to parallel functionality that are not actually parallel
"""
  Run serial rk4 and newton cases, and generate input files for the
  equivalent parallel cases
"""
function test_parallel()
  ARGS[1] = "input_vals_vortex3.jl"

  cd("./rk4/serial")
  include(STARTUP_PATH)
  cd("../parallel")

  arg_dict["smb_name"] = "SRCMESHES/psquare2.smb"
  make_input(arg_dict, "input_vals_parallel")

  cd("../../newton/serial")
  ARGS[1] = "input_vals_vortex3.jl"
  include(STARTUP_PATH)


  cd("../parallel")
  arg_dict["smb_name"] = "SRCMESHES/psquare2.smb"
  make_input(arg_dict, "input_vals_parallel")

  return nothing
end

#test_parallel()
add_func1!(EulerTests, test_parallel)
