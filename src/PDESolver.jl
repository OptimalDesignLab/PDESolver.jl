__precompile__(false)
module PDESolver

# from registration.jl
export register_physics, retrieve_physics

# from interface.jl
export evalResidual

# from initialization.jl
export createMeshAndOperator, call_nlsolver

# from startup_func
export run_solver

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

include("registration.jl")  # registering physics modules
include("interface.jl")  # functions all physics modules need to implement
include("initialization.jl")  # startup related functions
include(joinpath(Pkg.dir("PDESolver"), "src/input/read_input.jl"))
include("startup_func.jl")  # unified solver invokation

# package code goes here

end # module
