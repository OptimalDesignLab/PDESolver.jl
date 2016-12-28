TestFinalizeMPI = false
cd("./advection")
include(joinpath(pwd(), "runtests_parallel4.jl"))
cd("../euler")
resize!(ARGS, 0)  # reset to default
TestFinalizeMPI = true
include(joinpath(pwd(), "runtests_parallel4.jl"))
cd("../")
