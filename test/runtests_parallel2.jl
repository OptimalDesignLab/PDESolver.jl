# run all 2 processor tests in single session

args_orig = copy(ARGS)

TestFinalizeMPI = false
cd("./advection")
include(joinpath(pwd(), "runtests_parallel.jl"))
println("after advection tests, pwd = ", pwd())

cd("../euler")
resize!(ARGS, length(args_orig))  # reset to default
copy!(ARGS, args_orig)
TestFinalizeMPI = true
include(joinpath(pwd(), "runtests_parallel.jl"))

cd("../")

FactCheck.exitstatus()
