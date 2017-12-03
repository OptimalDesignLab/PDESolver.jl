# functional definitions

@doc """
###AdvectionEquationMod.QfluxData

Data type for storing relevant information pertaining to an a functional or an
objective function.

**Members**

*  `bcnums` : boundary condition groups on which the functional is to be
                             computed.
*  `val` : Computed value of the functional
*  `target_qFlux` : Target value for the functional qFlux

"""->
type QfluxData{Topt} <: AbstractIntegralFunctional{Topt}
  bcnums::Array{Int,1}
  val::Topt
  target_qflux::Topt
end # End type OfluxData

"""
  Constructor for QfluxData
"""
function QfluxDataConstructor{Topt}(::Type{Topt}, mesh, sbp, eqn, opts, bcnums)

  val = zero(Topt)
  target_qflux = zero(Topt)
  return QfluxData{Topt}(bcnums, val, target_qflux)
end

"""
  Functional that integrates the solution q over the specified boundary(/ies)

  **Fields**

   * bcnums: the boundary condition groups the functional is computed over
   * val: the value of the functional, initially 0.0

"""
type IntegralQData{Topt} <: AbstractIntegralFunctional{Topt}
  bcnums::Array{Int, 1}
  val::Topt
  bndries::Array{Boundary, 1}
end

"""
  Constructor for IntegralQData
"""
function IntegralQDataConstructor{Topt}(::Type{Topt}, mesh, sbp, eqn, opts,
                                        bcnums)

  val = 0.0
  bndries = getBoundaries(mesh, bcnums)
  return IntegralQData{Topt}(bcnums, val, bndries)
end


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

  functional_name = opts["objective_function"]
  functional_bcs = opts["objective_bcs"]
#  functional_regions = opts["geom_regions_objective"]
  func_constructor = FunctionalDict[functional_name]
  objective = func_constructor(Tsol, mesh, sbp, eqn, opts, functional_bcs)

  #=
  if opts["objective_function"] == "qflux"
    objective = QfluxData{Tsol}(mesh, sbp, eqn, opts, functional_faces)
  end
  =#
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

  dict_key = string("functional_name", functional_number)
  key = string("functional_bcs", functional_number)
  func_name = opts[dict_key]
  functional_faces = opts[key]

  func_constructor = FunctionalDict[func_name]
  functional = func_constructor(Tsol, mesh, sbp, eqn, opts, functional_faces)
#=
  if func_name == "qflux"
    functional = QfluxData{Tsol}(mesh, sbp, eqn, opts, functional_faces)
  elseif func_name "integralq"
    functional = IntegralQData{Tsol}(mesh, sbp, eqn, opts, functional_faces)
  end
=#

  return functional
end


"""
  Maps names of functionals to their outer constructor functions

  See `FunctionalDict` in the Euler module for descriptions of the constructors
  and arguments.
"""
global const FunctionalDict = Dict{ASCIIString, Function}(
"qflux" => QfluxDataConstructor,
"integralq" => IntegralQDataConstructor,
)


