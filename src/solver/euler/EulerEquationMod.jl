# Description:
#   Declare the Euler data object
#   Includes all the files for the Euler module

module EulerEquationMod

using PDESolver  # setup LOAD_PATH to find all PDESolver components
using ArrayViews
using ODLCommonTools  # abstract type definitions + common functions
using SummationByParts
using PdePumiInterface  # common mesh interface implementation - pumi
using NonlinearSolvers
using ForwardDiff
using Utils
import ODLCommonTools.sview
using MPI
using Input  # input file processing
using PETSc

# the AbstractEquation type is declared in ODLCommonTools
# every equation will have to declare a new type that is a subtype of AbstractEquation

export AbstractEulerData, EulerData, EulerData_, run_euler

# add a layer of abstraction
@doc """
### EulerEquationMod.AbstractEulerData{Tsol, Tres}

  This abstract type should be the supertype of *all* solution data objects
  that are related to the Euler equations.

  It should be used for specify the type of a function argument only when
  the function does no operations on the solution data object itself, it just
  passes it onto other functions that do the work (thus AbstractEulerData
  should be used for only the highest level functions).

  Another way of saying say it is that this type should only be used when
  the function only needs to ensure that it is solving the Euler equations,
  but does not care even a little bit about how.
"""->
abstract AbstractEulerData{Tsol, Tres} <: AbstractSolutionData{Tsol, Tres}
#=
Use this type to leverage multiple dispatch.
This type holds any data, including large arrays of solution values,
  that are specific to the equation for the Euler equations.
  This includes the solutoin variables q, the fluxes in the xi and eta
  direction, and the result of the calculation the inverse mass matrix is
  stored here as well (not sure if that fits better here or in the mesh object)
  things like the coordinate field, the jacobian etc. are stored
  in the mesh object

The aux_vars array holds all auxiliary variables that are stored over
  the entire mesh.
Although it can be accessed directly, the preferred method is to use the macros
  defined in the euler_macros.jl file.
These macros return either a scalar or an ArrayView of the specified indices,
  depending on if the quantity requested is scalar or vector.
Every time a new variable is added to the array, the size must be updated
  and a new macro must be created.

The uses of aux_vars should mirror that of eqn.q, in that entire columns
  should be passed to low level functions and the low level functions
  use the macros to access individual variables.
The advantages of macros vs functions for access to variables remains unclear
  if aux_vars is a fixed size.
If it is variable sized then macros give the advantage of doing location lookup
  at compile time
=#

#=
  Initial Condition:
    All initial condition functions use conservative variables, and are
    converted to the variables used to solve the equation.
=#


#=
  Boundary conditions:
    All boundary condition calculations are done in conservative variables.
    If using some other variables, they are converted before being passed
    to the boundary condition functions.

    The functors that compute boundary conditions are gathered from BCDict
    during init() and stored in mesh.bndry_funcs.  The functors compute the
    flux from the boundary condition at a node.  The boundary fluxes for all
    the nodes on the boundary are stored in mesh.bndryflux, to be integrated
    later using boundaryintegrate!.
=#

#=
  Variable Conversion:
    If solving the equations is some set of variables other than conservative,
    (lets call them the v variables, and the conservative variables the q
     variables), the following functions must be defined:

    # convert to conservative variables at a node
    convertToConservative(params::ParamType, qe::AbstractArray{Tsol, 1},
                          qc::AbstractArray{Tsol, 1}

    # convert from conservative to the other variables at a node
    convertToEntropy(params::ParamType, qe::AbstractArray{Tsol, 1},
                     qc::AbstractArray{Tsol, 1})

    # convert from conservative variables to the variables the equation
    # is being solved in
    convertFromConservativeToWorkingVars(params::ParamType{Tdim, :var_type,
                        qe::AbstractArray{Tsol, 1}, qc::AbstractArray{Tsol, 1})

    # calculate the coefficient matrix of the time derivate (dq/dv) at a node
    calcA0(params::ParamType, q::AbstractArray{Tsol, 1},
           A0::AbstractArray{Tsol, 2})
    # calculate inv(dq/dv)
    calcA0inv(params::ParamType, q::AbstractArray{Tsol, 1},
              A0::AbstractArray{Tsol, 2})

    # multiply a 3D array by inv(A0)
    matVecA0Inv(mesh, sbp, eqn, opts, arr::AbstractArray{Tsol, 3})

    # calculate the Euler flux at a node
    calcEulerFlux(params::ParamType, q::AbstractArray{Tsol, 1},
                  aux_vars::AbstractArray{Tsol, 1}, dir::AbstractArray{Tmsh},
                  F::AbstractArray{Tsol, 1})

=#

#=
  Calculating Other Quantities:
    To calculate physical quantities from the solution variables (for which
    the calculation is different for the different sets of variables) at a node,
    a separate function should be used, usnig the params object to dispatch
    to the right function for the variables being used.

    Some examples:

    calcPressure(params::ParamType, q::AbstractArray{Tsol, 1})
    calcEntropy(params::ParamType, q::AbstractArray{Tsol, 1})
    calcSpeedOfSound(params::ParamType, q::AbstractArray{Tsol, 1})
=#

@doc """
### EulerEquationMod.EulerData

  This type, although abstract, is the type functions should use for their
  input arguments if they do any operations on the solution data object.
  It stores all data used in evaluting the Euler Equations.

  It is paramaterized on the types Tsol, the type of the
  conservative variables q, and Tdim, the dimension of the equation

  It should have the following fields:
   * res_type : datatype of residual (depreciated)
    * q  : 3D array holding conservative variables
    * q_vec  : vector to assemble q into
    * aux_vars : 3D array holding auxiliary variables
    * flux_parametric : 4D array [ndof per node, nnodes per element, nelements, Tdim]
             holding the Euler flux in the xi and eta directions
    * res  : 3D array holding residual
    * res_vec   : vector form of res
    * edgestab_alpha : paramater used for edge stabilization, 4d array
    * bndryflux : 3D array holding boundary flux data
    * stabscale : 2D array holding edge stabilization scale factor
    * M : vector holding the mass matrix
    * Minv :  vector holding inverse mass matrix
    # Minv3D :  3D array holding inverse mass matrix for application to res (not res_vec)
"""->
abstract EulerData{Tsol, Tres, Tdim, var_type} <: AbstractEulerData{Tsol, Tres}

"""
  Functor type for faceElementIntegrals.  These integrals operate on a face,
  but require data from the entirety of the elements that make up the
  face, rather than data interpolated to the face
"""
abstract FaceElementIntegralType
# high level functions should take in an AbstractEulerData, remaining
# agnostic to the dimensionality of the equation
# Mid level function should take in an EulerData{Tsol, Tdim}, so they
# know the dimensionality of the equation, but have a single method that
# can handle all dimensions

# low level functions should take in EulerData{Tsol, 2} or EulerData{Tsol, 3}
# this allows them to have different methods for different dimension equations.

include("types.jl")  # type definitions
include(joinpath(Pkg.dir("PDESolver"), "src/solver/debug.jl"))  # debug macro
include("euler_macros.jl")
include("common_funcs.jl")
include("euler_funcs.jl")
include("conversion.jl")
include("euler.jl")
include("ic.jl")
include("bc.jl")
include("stabilization.jl")
include("faceElementIntegrals.jl")
include("flux.jl")
# include("artificialViscosity.jl")
# include("constant_diff.jl")
include("GLS2.jl")
include("boundary_functional.jl")
include("adjoint.jl")
include("source.jl")
include("PressureMod.jl")
include("entropy_flux.jl")
include("eigensystem.jl")
include("check_options.jl")
include("startup_func.jl")  # function for invoking the solver
include("./deriv/differentiateByMetrics.jl")

global const PhysicsName = "Euler"
register_physics(PhysicsName, EulerEquationMod, run_euler)

@doc """
### EulerEquationMod.PressureData

Subtype of AbstractOptimizationData. Stores all the information relevant to computing
an objective function pertaining to pressure coefficeint

**Members**

*  `targetCp_arr` : An array of arrays that stores the target coefficient of
                    pressure. length(targetCp_arr) = number of geometric edges
                    over which the functional is being computed. Each sub array
                    has dimensions (sbpface.numnodes, nfaces) *(from calcBoundarFlux
                    in bc.jl)*
*  `nodal_info` : 1D array of indices for one node needed to acces `targetCp_arr`
                  at a particular data point.
                  nodal_info[1] = geometric edge number
                  nodal_info[2] = sbpface node number
                  nodal_info[3] = element face number on the geometric edge

"""->
type PressureData{Tpress} <: AbstractOptimizationData
  targetCp_arr::Array{Array{Tpress,2},1}

  function PressureData(mesh::AbstractMesh, g_edges::AbstractArray{Int,1},
                        nface_arr::AbstractArray{Int,1})

    targetCp_arr = Array(Array{Tpress,2},length(g_edges))
    for i = 1:length(g_edges)
      targetCp_arr[i] = zeros(Tpress, mesh.sbpface.numnodes, nface_arr[i])
    end

    return new(targetCp_arr)

  end # End inner constructor
end   # End type PressureData

@doc """
### EulerEquationMod.LiftData

Subtype of AbstractOptimizationData. Stores all the information relevant to computing
an objective function pertaining to lift. Presently its an empty type
"""->

type LiftData{Topt} <: AbstractOptimizationData
  target_lift::Topt
  function LiftData()
    target_lift = zero(Topt)
    return new(target_lift)
  end
end

@doc """
### EulerEquationMod.DragData

Subtype of AbstractOptimizationData. Stores all the information relevant to
computing an objective function pertaining to drag. Presently its an empty type

"""->

type DragData{Topt} <: AbstractOptimizationData
  target_drag::Topt
  function DragData()
    target_drag = zero(Topt)
    return new(target_drag)
  end
end

@doc """
###EulerEquationMod.createObjectiveFunctionalData

Function for create an object for functional and adjoint computation where the
functional is an objective function in an optimization.

**Arguments**

* `mesh` : Abstract PUMI mesh
* `sbp`  : Summation-by-parts operator
* `eqn`  : Euler equation object
* `opts` : Options dictionary

"""->
function createObjectiveFunctionalData{Tsol}(mesh::AbstractMesh, sbp::AbstractSBP,
                                             eqn::EulerData{Tsol}, opts)

  functional_faces = opts["geom_faces_objective"]

  if opts["objective_function"] == "lift"
    objective = BoundaryForceData{Tsol, :lift}(mesh, sbp, eqn, opts, functional_faces)
    objective.is_objective_fn = true
  elseif opts["objective_function"] == "drag"
    objective = BoundaryForceData{Tsol, :drag}(mesh, sbp, eqn, opts, functional_faces)
    objective.is_objective_fn = true
  end # End if opts["objective_function"]

  return objective
end # End function createObjectiveFunctionalData(mesh, sbp, eqn, opts)

@doc """
###EulerEquationMod.createFunctionalData

Creates an object for functional computation. This function needs to be called
the same number of times as the number of functionals EXCLUDING the objective
function are being computed

**Arguments**

* `mesh` : Abstract PUMI mesh
* `sbp`  : Summation-by-parts operator
* `eqn`  : Euler equation object
* `opts` : Options dictionary
* `functional_number` : Which functional object is being generated. Default = 1

"""->

function createFunctionalData{Tsol}(mesh::AbstractMesh, sbp::AbstractSBP,
                                    eqn::EulerData{Tsol}, opts,
                                    functional_number::Int=1)

  dict_val = string("functional_name", functional_number)
  key = string("geom_faces_functional", functional_number)
  functional_faces = opts[key]

  if opts[dict_val] == "lift"
    functional = BoundaryForceData{Tsol, :lift}(mesh, sbp, eqn, opts, functional_faces)
  elseif opts[dict_val] == "drag"
    functional = BoundaryForceData{Tsol, :drag}(mesh, sbp, eqn, opts, functional_faces)
  end

  return functional
end

#=
@doc """
### EulerEquationMod.OptimizationData

A composite data type with data types pertaining to all possible objective
functions or boundary functionals. While dealing with an objective function,
make sure to use the bool

**Members**

* `pressCoeff_obj` : data object for functional corresponding to pressure
                     coefficients
* `lift_obj`       : data object for functional lift
* `drag_obj`       : data object for functional drag

**Constructor**

In order to perform optimization, a variable of type OptimizationData can be
constructed as

```
objective = OptimizationData(mesh, sbp, opts)
```

"""->

type OptimizationData{Topt} <: AbstractOptimizationData

  is_objective_fn::Bool
  ndof::Int
  val::Topt
  pressCoeff_obj::PressureData{Topt} # Objective function related to pressure coeff
  lift_obj::LiftData{Topt} # Objective function is lift
  drag_obj::DragData{Topt} # Objective function is drag
  force_obj::BoundaryForceData{Topt} # Objective function is boundaryForce

  function OptimizationData(mesh::AbstractMesh, sbp::AbstractSBP, opts)

    functional = new()
    functional.val = zero(Topt)

    for i = 1:opts["num_functionals"]
      dict_val = string("functional_name", i)

      if opts[dict_val] == "targetCp" || opts["objective_function"] == "targetCp"

        if opts["objective_function"] == "targetCp"
          functional.is_objective_fn = true
        else
          functional.is_objective_fn = false
        end
        g_edges = opts[string("geom_edges_functional", i)]
        nface_arr = zeros(Int, length(g_edges))
        for j = 1:length(nface_arr)
          nface_arr[j] = getnFaces(mesh, g_edges[j])
        end
        functional.pressCoeff_obj = PressureData{Topt}(mesh, g_edges, nface_arr)
        functional.ndof = 1

      elseif opts[dict_val] == "lift" || opts["objective_function"] == "lift"

        if opts["objective_function"] == "lift"
          functional.is_objective_fn = true
        else
          functional.is_objective_fn = false
        end
        functional.lift_obj = LiftData{Topt}()
        functional.dof = 1

      elseif opts[dict_val] == "drag" || opts["objective_function"] == "drag"

        if opts["objective_function"] == "drag"
          functional.is_objective_fn = true
        else
          functional.is_objective_fn = false
        end
        functional.drag_obj = DragData{Topt}()
        functional.ndof = 1


      end # End if
    end   # End for i = 1:opts["num_functionals"]

    return functional
  end  # End inner constructor
end    # End OptimizationData
=#

end # end module
