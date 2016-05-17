push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))
using Utils
using FactCheck

facts("---- Testing IO -----") do
  fname = "iotest.dat"
  f = open(fname, "w")
  f2 = BufferedIO(f)

  println(f2, 1)
  println(f2, 2)

  fstat = stat(fname)
  @fact fstat.size --> 0
  flush(f2)
#  close(f)

  data = readdlm(fname)
  @fact data[1] --> 1
  @fact data[2] --> 2

  close(f)
end
