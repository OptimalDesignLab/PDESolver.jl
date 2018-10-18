module Input

  import MPI
  using PETSc2
  using Utils
  export read_input, make_input, read_input_file, registerOptionsChecker,
         printOpts

  include("physics_specific.jl")
  include("read_input.jl")
  include("make_input.jl")
  include("known_keys.jl")
  include("utils.jl")

end  # end module
