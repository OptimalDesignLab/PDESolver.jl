#!/users/pandak/julia/julia

# script for running a test case on Prof. Hicken's euler code
# to run this:
#   nohup ./run_test.jl | tee output.txt &
# (may only work on diplomacy, mastermind, & balderdash)

push!(LOAD_PATH, "/users/ashlea/.julia/v0.4/PUMI")
using PumiInterface
using PdePumiInterface
using SummationByParts
include("test_solver.jl")

println("==== done loading everything ====")
println(" ")
println("=================================")
println("Start time: ")
timestr = Libc.strftime(time())
println(timestr)
println("=================================")

println("TestSolver.runtest(degree=3, Nx=20, Ny=20, cfl=0.005, output=true, energyhist=true, matrix=false)")
TestSolver.runtest(degree=3, Nx=20, Ny=20, cfl=0.005, output=true, energyhist=true, matrix=false)

# println("TestSolver.runtest(degree=1, Nx=10, Ny=10, cfl=0.1, output=true, energyhist=true, matrix=false)")
# TestSolver.runtest(degree=1, Nx=10, Ny=10, cfl=0.1, output=true, energyhist=true, matrix=false)

# println("TestSolver.runtest(degree=1, Nx=60, Ny=60, cfl=0.03, output=true, energyhist=true, matrix=false)")
# TestSolver.runtest(degree=1, Nx=60, Ny=60, cfl=0.03, output=true, energyhist=true, matrix=false)

# println("TestSolver.runtest(degree=1, Nx=40, Ny=40, cfl=0.05, output=true, energyhist=true, matrix=false)")
# TestSolver.runtest(degree=1, Nx=40, Ny=40, cfl=0.05, output=true, energyhist=true, matrix=false)

# println("TestSolver.runtest(degree=1, Nx=5, Ny=5, cfl=0.5, output=true, energyhist=true, matrix=false)")
# TestSolver.runtest(degree=1, Nx=5, Ny=5, cfl=0.5, output=true, energyhist=true, matrix=false)

println("=================================")
println("End time: ")
timestr = Libc.strftime(time())
println(timestr)
println("=================================")
