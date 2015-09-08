
facts("--- Testing Sparse/Dense Jacobian ---") do

  resize!(ARGS, 1)
  ARGS[1] = "input_vals_vortex.jl"

  include("../src/solver/euler/startup.jl")

  @fact calcNorm(eqn, eqn.SL) => less_than(1e-9)

  resize!(ARGS, 1)
  ARGS[1] = "input_vals_vortex2.jl"

  include("../src/solver/euler/startup.jl")

  @fact calcNorm(eqn, eqn.SL) => less_than(1e-9)

  resize!(ARGS, 1)
  ARGS[1] = "input_vals_vortex3.jl"

  include("../src/solver/euler/startup.jl")

  @fact calcNorm(eqn, eqn.SL) => less_than(1e-9)

  resize!(ARGS, 1)
  ARGS[1] = "input_vals_vortex4.jl"

  include("../src/solver/euler/startup.jl")

  @fact calcNorm(eqn, eqn.SL) => less_than(1e-9)

end
