# run all 2 processor tests in single session

args_orig = copy(ARGS)

cd("./advection")
include(joinpath(pwd(), "runtests_parallel2.jl"))

cd("../euler")
resize!(ARGS, length(args_orig))  # reset to default
copy!(ARGS, args_orig)
include(joinpath(pwd(), "runtests_parallel2.jl"))

cd("../")
