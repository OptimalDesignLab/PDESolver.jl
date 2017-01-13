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
export OptimizationData

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
include("startup_func.jl")  # function to invoke the solver

# register this physics module
global const PhysicsName = "Advection"
register_physics(PhysicsName, AdvectionEquationMod, run_advection)

@doc """
###AdvectionEquationMod.qfluxData

Data type for storing relevant information pertaining to an a functional or an
objective function.

**Members**

*  `target_qFlux` : Target value for the functional qFlux

"""->

type QfluxData{Topt} <: AbstractOptimizationData
  target_qflux::Topt
  function QfluxData()
    target_qflux = zero(Topt)
    return new(target_qflux)
  end
end

@doc """
###AdvectionEquationMod.OptimizationData

"""->

type OptimizationData{Topt} <: AbstractOptimizationData

  is_objective_fn::Bool
  val::Topt
  qflux_obj::QfluxData{Topt}

  function OptimizationData(mesh::AbstractMesh, sbp::AbstractSBP, opts)

    functional = new()
    functional.val = zero(Topt)

    for i = 1:opts["num_functionals"]
      dict_val = string("functional_name",i)
      
      if opts[dict_val] == "qflux" || opts["objective_function"] == "qflux"
        if opts["objective_function"] == "qflux"
          functional.is_objective_fn = true
        else
          functional.is_objective_fn = false
        end
        functional.qflux_obj = QfluxData{Topt}()
      end # End if opts[dict_val]

    end  # End  for i = 1:opts["num_functionals"]
    
    return functional
  end # End inner constructor
end # End type OptimizationData

end # end module
