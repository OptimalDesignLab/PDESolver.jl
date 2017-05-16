module AdvectionEquationMod

using PDESolver  # setupf LOAD_PATH to find PDESolver components
using ArrayViews
using ODLCommonTools
using SummationByParts
using PdePumiInterface
using ForwardDiff
using NonlinearSolvers
using MPI
using Utils
import ODLCommonTools.sview
using Input


export AdvectionData, AdvectionData_, run_advection #getMass, assembleSolution, disassembleSolution

# include("advection.jl")
# include("getMass.jl")


abstract AbstractAdvectionData{Tsol, Tres} <: AbstractSolutionData{Tsol, Tres}
abstract AdvectionData{Tsol, Tres, Tdim} <: AbstractAdvectionData{Tsol, Tres}

include("types.jl")
include(joinpath(Pkg.dir("PDESolver"), "src/solver/debug.jl"))  # debug macro
include("advection.jl")
include("common_funcs.jl")
include("bc.jl")
include("bc_solvers.jl")
include("ic.jl")
include("GLS.jl")
include("GLS2.jl")
include("boundary_functional.jl")
include("adjoint.jl")
include("source.jl")
include("flux.jl")
include("check_options.jl")
include("eqn_copy.jl")
include("startup_func.jl")  # function to invoke the solver

# register this physics module
global const PhysicsName = "Advection"
register_physics(PhysicsName, AdvectionEquationMod, run_advection)

end # end module
