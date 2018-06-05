module EllipticEquationMod

using PDESolver
using ArrayViews
import ArrayViews.view
using SolverCommon
using ODLCommonTools
using SummationByParts
using PdePumiInterface
using NonlinearSolvers
using LinearSolvers
using Input
# using ForwardDiff
using Utils
using MPI
using Input  # input file processing
using PETSc2

import ODLCommonTools.sview

export AbstractEllipticData, EllipticData, EllipticData_, run_elliptic



abstract type AbstractEllipticData{Tsol, Tres} <: AbstractSolutionData{Tsol, Tres} end
abstract type EllipticData{Tsol, Tres, Tdim} <: AbstractEllipticData{Tsol, Tres} end

# now that EllipticData is declared, include other files that use it
include(joinpath(Pkg.dir("PDESolver"), "src/solver/debug.jl"))  # debug macro
include("elliptic.jl")
# include("time_advance.jl")
include("bc.jl")
include("ic.jl")
include("flux.jl")
include("check_options.jl")
include("source.jl")
#include("output.jl")
include("exactSolution.jl")
include("diffusion.jl")
include("functional.jl")
include("types.jl")
include("startup_func.jl")  # function for invoking the solver

global const PhysicsName = "Elliptic"
register_physics(PhysicsName, EllipticEquationMod, createObjects, checkOptions)
end # end of module
