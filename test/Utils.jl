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

facts("----- Testing Timings -----") do
  t = Timings()

  t.t_volume = 1.0
  t.t_face = 2.0
  t.t_source = 3.0
  t.t_sharedface = 4.0
  t.t_bndry = 5.0
  t.t_send = 6.0
  t.t_wait = 7.0
  t.t_allreduce = 8.0
  t.t_jacobian = 9.0
  t.t_solve = 10.0
  t.t_barrier = 11.0
  t.t_barrier2 = 12.0
  t.t_barrier3 = 13.0
  t.t_barriers[1:7] = 14.0:20.0

  write_timings(t, "timing_breakdown")
  vals = readdlm("timing_breakdown.dat")
  for i=1:20
    @fact vals[i] --> roughly(Float64(i))
  end

end
