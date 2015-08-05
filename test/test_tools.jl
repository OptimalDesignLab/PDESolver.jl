# test the tools
include(joinpath(Pkg.dir("PDESolver"), "src/tools/misc.jl"))

facts("--- Testing tools/misc.jl ---") do

  A = rand(3,3)
  x = rand(3)
  b = A*x
  b2 = zeros(3)
  smallmatvec!(A, x, b2)

  @fact b => roughly(b2, atol=1e-14)

end
