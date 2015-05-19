How to run Prof. Hicken's euler code
====================================

julia environment:

push!(LOAD_PATH, "/users/ashlea/.julia/v0.4/PUMI")
using PumiInterface
using PdePumiInterface
using SummationByParts
include("test_solver.jl")
TestSolver.runtest(degree=1, Nx=5, Ny=5, cfl=0.5, output=true, energyhist=true, matrix=false)



