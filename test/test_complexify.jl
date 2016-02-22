using FactCheck
include(joinpath(Pkg.dir("PDESolver"), "src/solver/euler/complexify.jl"))

facts("----- Testing Complexify.jl ------") do

  a = 1.0
  b = -1.0
  @fact absvalue(a) --> roughly(abs(a), atol=1e-15)
  @fact absvalue(b) --> roughly(abs(b), atol=1e-15)

  c = complex(1.0, 1.0)
  d = complex(1.0, -1.0)
  f = complex(-1.0, 1.0)
  g = complex(-1.0, -1.0)

  @fact absvalue(c) --> roughly(complex(1.0, 1.0), atol=1e-15)
  @fact absvalue(d) --> roughly(complex(1.0, -1.0), atol=1e-15)
  @fact absvalue(f) --> roughly(complex(1.0, -1.0), atol=1e-15)
  @fact absvalue(g) --> roughly(complex(1.0, 1.0), atol=1e-15)

end
