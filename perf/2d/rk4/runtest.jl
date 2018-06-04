using ODLCommonTools
using PDESolver
using EulerEquationMod
using AdvectionEquationMod
using Utils

function runtest(;fname="input_vals_vortex3.jl", outname="perf_hist.txt")

  startup_path = joinpath(Pkg.dir("PDESolver"), "src/solver/euler/startup.jl")
  resize!(ARGS, 1)
  ARGS[1] = fname
  println("----- warming up -----")
#  val, telapsed, tgc, gc_bytes = @time_all include(startup_path)
  val, telapsed, tgc, gc_bytes = @time_all run_solver(fname)

  println("----- Final run -----")
  val, telapsed, tgc, gc_bytes = @time_all run_solver(fname)

  f = open(outname, "a+")
  tstr = getTimeString()
  branchname = getBranchName()
  println(f, tstr, " ", branchname, " ", telapsed, " ", gc_bytes)
  close(f)

end

function runtest_profile(;fname="input_vals_vortex3.jl")

  println("----- warming up -----")
  @time include(fname)

  Profile.clear_malloc_data()
  println("----- Final run -----")
  @time include(fname)
end

#runtestprofile()
runtest()

