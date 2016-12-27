__precompile__(false)
module PDESolver

# from registration.jl
export register_physics, retrieve_physics

# load paths for all the components of PDESolver
push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/euler"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))

# add physics modules to load path (but don't load them, because that would
# create a circular dependency)
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/advection"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/euler"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/simpleODE"))

include("registration.jl")

# load the modules
using ODLCommonTools
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
using ForwardDiff
using NonlinearSolvers   # non-linear solvers
using ArrayViews
using Utils
import ODLCommonTools.sview
using MPI




# package code goes here

end # module
