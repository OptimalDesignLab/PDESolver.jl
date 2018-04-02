# provide a way to invoke any solver

"""
  This function provides a way to invoke any physics solver based on the
  specification of the physics in the input file.
  This requires loading the input file twice, once to figure out the physics,
  and a second time when the physics-specific startup function is called

  The physics module must have already been registered using [`register_physics`](@ref)

  Inputs:

    input_file: an AbstractString specifying the path to the input file

  Outputs:

    mesh: the AbstractMesh object used during the solve
    sbp: the SBP operator used by the solver
    eqn: the AbstractSolutionData object during the solve.  At exit,
         eqn.q_vec should have the final solution in it
    opts: the options dictionary
"""
function run_solver(input_file::AbstractString)

  #TODO; `make physcs` a required key, don't load input file twice
  # load the input file in order to figure out what physics to solve
  opts = read_input(input_file)

  physics_name = opts["physics"]

  # get the module and startup function
  mod, func = retrieve_physics(physics_name)

  println("Physics function: ", func)

  # call the function
  mesh, sbp, eqn, opts = func(input_file)

  return mesh, sbp, eqn, opts
end


