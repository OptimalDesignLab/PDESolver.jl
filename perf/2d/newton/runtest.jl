using ODLCommonTools
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))
using Utils

function runtest_rk(;fname="input_vals_vortex3.jl", outname="perf_hist.txt")

  startup_path = joinpath(Pkg.dir("PDESolver"), "src/solver/euler/startup.jl")
  resize!(ARGS, 1)
  ARGS[1] = fname
  println("----- warming up -----")
  val, telapsed, tgc, gc_bytes = @time_all include(startup_path)

  println("----- Final run -----")
  val, telapsed, tgc, gc_bytes = @time_all include(startup_path)

  f = open(outname, "a+")
  tstr = getTimeString()
  branchname = getBranchName()
  println(f, tstr, " ", branchname, " ", telapsed, " ", gc_bytes)
  close(f)


end

runtest_rk()

