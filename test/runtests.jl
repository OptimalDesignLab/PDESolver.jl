# run all serial tests in a single session

args_orig = copy(ARGS)


cd("./advection")
include(joinpath(pwd(), "runtests.jl"))
println("after advection tests, pwd = ", pwd())

cd("../euler")
resize!(ARGS, length(args_orig))  # reset to default
copy!(ARGS, args_orig)
include(joinpath(pwd(), "runtests.jl"))
println("after Euler tests, pwd = ", pwd())

cd("../simpleODE")
resize!(ARGS, length(args_orig))  # reset to default
copy!(ARGS, args_orig)
include(joinpath(pwd(), "runtests.jl"))

cd("../elliptic")
resize!(ARGS, length(args_orig))
copy!(ARGS, args_orig)
include(joinpath(pwd(), "runtests.jl"))

cd("../")
