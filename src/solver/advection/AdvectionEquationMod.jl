module AdvectionEquationMod

using PDESolver  # setupf LOAD_PATH to find PDESolver components
using SolverCommon
using ArrayViews
using ODLCommonTools
using SummationByParts
using PdePumiInterface
using NonlinearSolvers
using LinearSolvers
using MPI
using Utils
#using PETSc

import ODLCommonTools.sview
using Input


export AdvectionData, AdvectionData_, run_advection #getMass, array3DTo1D, array1DTo3D

# include("advection.jl")
# include("getMass.jl")

"""
  Direct subtype of [`AbstractSolutionData`](@ref), inheriting `Tsol` and
  `Tres` as static parameter
"""
abstract AbstractAdvectionData{Tsol, Tres} <: AbstractSolutionData{Tsol, Tres}

"""
  Subtype of [`AbstractAdvectionData`](@ref), inheriting its static parameters
  and adding `Tdim`.
"""
abstract AdvectionData{Tsol, Tres, Tdim} <: AbstractAdvectionData{Tsol, Tres}

include("types.jl")
include("functionals.jl")
include(joinpath(Pkg.dir("PDESolver"), "src/solver/debug.jl"))  # debug macro
include("advection.jl")
include("common_funcs.jl")
include("bc.jl")
include("bc_solvers.jl")
include("ic.jl")
include("GLS.jl")
include("GLS2.jl")
include("boundary_functional.jl")
include("functional_deriv.jl")
include("source.jl")
include("flux.jl")
include("check_options.jl")
include("eqn_deepcopy.jl")
include("startup_func.jl")  # function to invoke the solver

# register this physics module
"""
  This physics is named `Advection`
"""
global const PhysicsName = "Advection"
register_physics(PhysicsName, AdvectionEquationMod, createObjects, checkOptions)

end # end module
