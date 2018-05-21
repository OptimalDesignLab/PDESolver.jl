# provide a way to invoke any solver

"""
  This function provides a way to invoke any physics solver based on the
  specification of the physics in the input file.

  The physics module must have already been registered using [`register_physics`](@ref)

  This function is equivalent to [`solvePDE`](@ref), it is maintained here
  for legacy purposes only.

  **Inputs**

   * input_file: an AbstractString specifying the path to the input file

  **Outputs**

   * mesh: the AbstractMesh object used during the solve
   * sbp: the SBP operator used by the solver
   * eqn: the AbstractSolutionData object during the solve.  At exit,
          eqn.q_vec should have the final solution in it
   * opts: the options dictionary
"""
function run_solver(input_file::AbstractString)

  return solvePDE(input_file)

end


