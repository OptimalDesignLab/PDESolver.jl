using ODLCommonTools
#push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))
using PDESolver
using EulerEquationMod
using Utils

function runtest_rk(;fname="input_vals_vortex3.jl", outname="perf_hist.txt")

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

function runtest_rkprofile(;fname="input_vals_vortex3.jl")

  startup_path = joinpath(Pkg.dir("PDESolver"), "src/solver/euler/startup.jl")
  resize!(ARGS, 1)
  ARGS[1] = fname
  println("----- warming up -----")
  @time include(startup_path)

  Profile.clear_malloc_data()
  println("----- Final run -----")
  @time include(startup_path)
end

# advection version
function runtestadvection_rkprofile(;fname="input_vals1.jl")

  startup_path = joinpath(Pkg.dir("PDESolver"), "src/solver/advection/startup_advection.jl")
  resize!(ARGS, 1)
  ARGS[1] = fname
  println("----- warming up -----")
  @time include(startup_path)

  Profile.clear_malloc_data()
  println("----- Final run -----")
  @time include(startup_path)
end


#runtest_rkprofile()
#runtestadvection_rkprofile()
runtest_rk()

