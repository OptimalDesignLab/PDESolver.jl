# Description:
#   Declare the Euler data object
#   Includes all the files for the Euler module

module EulerEquationMod

using PDESolver  # setup LOAD_PATH to find all PDESolver components
using SolverCommon
using ArrayViews
using ODLCommonTools  # abstract type definitions + common functions
using SummationByParts
using PdePumiInterface  # common mesh interface implementation - pumi
using NonlinearSolvers
using LinearSolvers
using ForwardDiff
using Utils
import ODLCommonTools.sview
using MPI
using Input  # input file processing
#using PETSc
# using FreeFormDeformation
# using MeshMovement

# the AbstractEquation type is declared in ODLCommonTools
# every equation will have to declare a new type that is a subtype of AbstractEquation

export AbstractEulerData, EulerData, EulerData_, run_euler, eval_dJdaoa

# add a layer of abstraction
@doc """
### EulerEquationMod.AbstractEulerData{Tsol, Tres}

  This abstract type should be the supertype of *all* solution data objects
  that are related to the Euler equations.

  It should be used for specifying the type of a function argument only when
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
  This includes the solution variables q, the fluxes in the xi and eta
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
  at compile time.
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
include("functionals.jl")
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
include("functional_deriv.jl")
include("source.jl")
include("PressureMod.jl")
include("entropy_flux.jl")
include("eigensystem.jl")
include("check_options.jl")
include("eqn_deepcopy.jl")
include("startup_func.jl")  # function for invoking the solver
include("dataprep_rev.jl")
include("evaldRdm.jl")
include("homotopy.jl")

# Jacobian calculation
include("evalJacobian.jl")
include("euler_funcs_diff.jl")
include("flux_diff.jl")
include("bc_solvers_diff.jl")
include("bc_diff.jl")


"""
  This physics is named `Euler`
"""
global const PhysicsName = "Euler"

register_physics(PhysicsName, EulerEquationMod, run_euler)


end # end module
