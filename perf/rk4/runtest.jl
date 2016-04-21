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
runtestadvection_rkprofile()
#runtest_rk()

