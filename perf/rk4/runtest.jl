using Coverage

function runtest_rk(;fname="input_vals_vortex3.jl")

  startup_path = joinpath(Pkg.dir("PDESolver"), "src/solver/euler/startup.jl")
  resize!(ARGS, 1)
  ARGS[1] = fname
  println("----- warming up -----")
  @time include(startup_path)

  println("----- Final run -----")
  @time include(startup_path)
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

runtest_rkprofile()
#runtest_rk()

