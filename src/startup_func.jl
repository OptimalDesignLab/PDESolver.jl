# provide a way to invoke any solver

"""
  This function provides a way to invoke any physics solver based on the
  specification of the physics in the input file.
  This requires loading the input file twice, once to figure out the physics,
  and a second time when the physics-specific startup function is called
"""
function run_solver(input_file::AbstractString)

  if !MPI.Initialized()
    MPI.Init()
  end

  # load the input file in order to figure out what physics to solve
  opts = read_input(input_file)

  physics_name = opts["physics"]

  # get the module and startup function
  mod, func = retrieve_physics(physics_name)

  # call the function
  mesh, sbp, eqn, opts = func(input_file)

  return mesh, sbp, eqn, opts
end


