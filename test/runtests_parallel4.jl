# run all 4 processor tests in a single session
args_orig = copy(ARGS)

TestFinalizeMPI = false
cd("./advection")
include(joinpath(pwd(), "runtests_parallel4.jl"))

cd("../euler")
resize!(ARGS, length(args_orig))  # reset to default
copy!(ARGS, args_orig)
TestFinalizeMPI = true
include(joinpath(pwd(), "runtests_parallel4.jl"))

cd("../")


FactCheck.exitstatus()
