# run tests related to parallel functionality that are not actually parallel
"""
  Run serial rk4 and newton cases, and generate input files for the
  equivalent parallel cases
"""
function test_parallel()
  start_dir = pwd()
  ARGS[1] = "input_vals_vortex3.jl"

  # rk4
  cd("./rk4/serial")
  mesh, sbp, eqn, opts = solvePDE(ARGS[1])
  cd("../parallel")

  opts["smb_name"] = "SRCMESHES/psquare2.smb"
  make_input(opts, "input_vals_parallel")

  # staggered grid
  cd("../staggered_serial")
  mesh, sbp, eqn, opts = solvePDE(ARGS[1])

  cd("../staggered_parallel")
  opts["smb_name"] = "SRCMESHES/psquare2.smb"
  make_input(opts, "input_vals_parallel")

  cd("../../")
  # lserk
  if !isdir("./lserk")
    mkdir("./lserk")
  end
  if !isdir("./lserk/serial")
    mkdir("./lserk/serial")
  end
  if !isdir("lserk/parallel")
    mkdir("./lserk/parallel")
  end

  opts = Input.read_input_file("./rk4/serial/input_vals_vortex3.jl")
  cd("./lserk/serial")
  opts["run_type"] = 30
  make_input(opts, "input_vals_vortex3")
  mesh, sbp, eqn,  opts = solvePDE(ARGS[1])

  cd("../parallel")
  opts["smb_name"] = "SRCMESHES/psquare2.smb"
  make_input(opts, "input_vals_parallel")


  # newton
  cd("../../newton/serial")
  ARGS[1] = "input_vals_vortex3.jl"
  mesh, sbp, eqn, opts = solvePDE(ARGS[1])


  cd("../parallel")
  opts["smb_name"] = "SRCMESHES/psquare2.smb"
  make_input(opts, "input_vals_parallel")

  cd(start_dir)

  return nothing
end

#test_parallel()
add_func1!(EulerTests, test_parallel, [TAG_PARALLEL, TAG_SHORTTEST])
