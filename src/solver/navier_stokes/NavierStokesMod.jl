module NavierStokesMod

using PDESolver  # setup LOAD_PATH to find all PDESolver components
using SolverCommon
#using ArrayViews
using ODLCommonTools  # abstract type definitions + common functions
using SummationByParts
using PdePumiInterface  # common mesh interface implementation - pumi
#using NonlinearSolvers
using Jacobian
using LinearSolvers
using Utils
import ODLCommonTools.sview
using MPI
using Input  # input file processing
using PETSc2
using EulerEquationMod

"""
  Abstract supertype of all NavierStokesData implementations.  Functions should
  use this type for the `eqn` argument unless they rely on the behavior of a
  specific implementation
"""
abstract type NSData{Tsol, Tres, Tdim} <: AbstractSolutionData{Tsol, Tres} end

include("types.jl")            # type definitions
include("startup_func.jl")     # startup functions
include("check_options.jl")    # physics specific options checkgin
include("navier_stokes.jl")    # main residual evaluation
include("util.jl")     # utilities
include("bc.jl")       # viscous boundary condition
include("flux.jl")     # viscous face flux
include("viscous_func.jl")     # more viscous stuff
include("source.jl")   # viscous source term
include("ic.jl")       # initial conditions


"""
  This physics is named `NavierStokes`
"""
global const PhysicsName = "NavierStokes"

register_physics(PhysicsName, NavierStokesMod, createObjects, checkOptions)


end # end module
