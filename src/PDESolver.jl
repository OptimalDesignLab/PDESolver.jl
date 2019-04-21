__precompile__(false)
module PDESolver

# defs.jl
export AssembleElementData, AbstractShockSensor

# from registration.jl
export register_physics, retrieve_physics, registerIC, registerBC

# from interface.jl
export evalResidual, evalJacobian, evalHomotopy, evalHomotopyJacobian,
       evalJacobianStrong, createFunctional,
       evalFunctional, evalFunctionalDeriv_q, evalFunctionalDeriv_m,
       updateMetricDependents,
       solvePDE, evalResidual_revm, evalResidual_revq, evaldRdqProduct,
       _getSparsityPattern, getShockSensor, setShockSensor, setShockSensorAlpha,
       createShockSensor

# from interface2.jl
export createObjects, createLinearSolver, getSparsityPattern

# from startup_func
export run_solver

# from interative.jl
export printICNames, printBCNames

# load paths for all the components of PDESolver
push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/linearsolvers"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/jacobian"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Debugging"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/input"))

# add physics modules to load path (but don't load them, because that would
# create a circular dependency)
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/advection"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/euler"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/elliptic"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/simpleODE"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/elliptic"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/navier_stokes"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/optimization"))

# load the modules
using MPI
using ODLCommonTools
import ODLCommonTools: evalFunctional, _evalFunctional,
                       evalFunctionalDeriv_m, _evalFunctionalDeriv_m,
                       evalFunctionalDeriv_q, _evalFunctionalDeriv_q

function finalizeMPI()
  if MPI.Initialized()
    MPI.Finalize()
  end
end

if !MPI.Initialized()
  MPI.Init()
  atexit(finalizeMPI)
end

using PETSc2

function finalizePetsc()
  if PetscInitialized()
    PetscFinalize()
  end
end

if !PetscInitialized()
  PetscInitialize()
  atexit(finalizePetsc)
end


using Input
using ODLCommonTools
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
#using LinearSolvers
#using NonlinearSolvers   # non-linear solvers
using ArrayViews
import ArrayViews.view
using Utils
import ODLCommonTools.sview
#using Input

include("defs.jl")  # common definitions
include("registration.jl")  # registering physics modules
include("interface.jl")  # functions all physics modules need to implement
include("interface2.jl") # functions physics modules don't have to implement
include("startup_func.jl")  # unified solver invokation
include("interactive.jl")
# package code goes here

end # module
