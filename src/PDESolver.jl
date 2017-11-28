__precompile__(false)
module PDESolver

# from registration.jl
export register_physics, retrieve_physics, registerIC, registerBC

# from interface.jl
export evalResidual, evalHomotopy

# from initialization.jl
export createMeshAndOperator, loadRestartState, call_nlsolver

# from startup_func
export run_solver

# from interative.jl
export printICNames, printBCNames

# load paths for all the components of PDESolver
push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/euler"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/linearsolvers"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/input"))

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
using LinearSolvers
using NonlinearSolvers   # non-linear solvers
using ArrayViews
using Utils
import ODLCommonTools.sview
using MPI
using Input
using PETSc

if !MPI.Initialized()
  MPI.Init()
  mpi_inited = true
else
  mpi_inited =false
end

if PetscInitialized() == 0
  PetscInitialize()
  petsc_inited = true
else
  petsc_inited = false
end

function finalizePetsc()
  if PetscInitialized() == 0
    PetscFinalize()
  end
end

function finalizeMPI()
  if MPI.Initialized()
    MPI.Finalize()
  end
end

if petsc_inited
  atexit( () -> finalizePetsc() )
end

if mpi_inited
  atexit( () -> finalizeMPI() )
end

include("registration.jl")  # registering physics modules
include("interface.jl")  # functions all physics modules need to implement
include("initialization.jl")  # startup related functions
include("startup_func.jl")  # unified solver invokation
include("interactive.jl")
# package code goes here

end # module
