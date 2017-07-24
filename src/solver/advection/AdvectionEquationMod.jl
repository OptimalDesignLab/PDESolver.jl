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
using PETSc

import ODLCommonTools.sview
using Input


export AdvectionData, AdvectionData_, run_advection #getMass, assembleSolution, disassembleSolution

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
include("startup_func.jl")  # function to invoke the solver

# register this physics module
"""
  This physics is named `Advection`
"""
global const PhysicsName = "Advection"
register_physics(PhysicsName, AdvectionEquationMod, run_advection)

@doc """
###AdvectionEquationMod.createObjectiveFunctionalData

Function for creating an object for functional and adjoint computation where the
functional is an objective function in an optimization.

**Arguments**

* `mesh` : Abstract PUMI mesh
* `sbp`  : Summation-by-parts operator
* `eqn`  : Advection equation object
* `opts` : Options dictionary

"""->

function createObjectiveFunctionalData{Tsol}(mesh::AbstractMesh, sbp::AbstractSBP,
                                             eqn::AdvectionData{Tsol}, opts)

  functional_faces = opts["geom_faces_objective"]
  if opts["objective_function"] == "qflux"
    objective = QfluxData{Tsol}(mesh, sbp, eqn, opts, functional_faces)
    objective.is_objective_fn = true
  end

  return objective
end

@doc """
###AdvectionEquationMod.createFunctionalData

Creates an object for functional computation. This function needs to be called
the same number of times as the number of functionals EXCLUDING the objective
function are being computed

**Arguments**

* `mesh` : Abstract PUMI mesh
* `sbp`  : Summation-by-parts operator
* `eqn`  : Advection equation object
* `opts` : Options dictionary
* `functional_number` : Which functional object is being generated. Default = 1

"""->

function createFunctionalData{Tsol}(mesh::AbstractMesh, sbp::AbstractSBP,
                                    eqn::AdvectionData{Tsol}, opts,
                                    functional_number::Int=1)

  dict_val = string("functional_name", functional_number)
  key = string("geom_faces_functional", functional_number)
  functional_faces = opts[key]

  if opts[dict_val] == "qflux"
    functional = QfluxData{Tsol}(mesh, sbp, eqn, opts, functional_faces)
  end

  return functional
end

end # end module
