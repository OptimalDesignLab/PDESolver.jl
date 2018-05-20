# functional definitions

import PDESolver.createFunctional

@doc """
###EulerEquationMod.BoundaryForceData

Composite data type for storing data pertaining to the boundaryForce. It holds
lift and drag values

"""->

type BoundaryForceData{Topt, fname} <: AbstractIntegralFunctional
  bcnums::Array{Int,1}  #TODO: make this non-abstract
  ndof::Int
  bndry_force::Array{Topt,1}
  isLift::Bool
  lift_val::Topt
  drag_val::Topt
  dLiftdaoa::Topt # Partial derivative of lift w.r.t. angle of attack
  dDragdaoa::Topt # Partial derivative of drag w.r.t. angle of attack

end

"""
  Constructor for BoundaryForceData{Topt, :lift}
"""
function LiftForceDataConstructor{Topt}(::Type{Topt}, mesh, sbp, eqn, opts, bcnums)

  ndof = mesh.dim
  bndry_force = zeros(Topt, mesh.dim)
  isLift = true
  lift_val = 0.0
  drag_val = 0.0
  dLiftdaoa = 0.0
  dDragdaoa = 0.0

  return BoundaryForceData{Topt, :lift}(bcnums, ndof,
                           bndry_force, isLift, lift_val, drag_val, dLiftdaoa,
                           dDragdaoa)
end

"""
  Constructor for BoundaryForceData{Topt, :drag}
"""
function DragForceDataConstructor{Topt}(::Type{Topt}, mesh, sbp, eqn, opts,
                             bcnums)

  ndof = mesh.dim
  bndry_force = zeros(Topt, mesh.dim)
  isLift = false
  lift_val = 0.0
  drag_val = 0.0
  dLiftdaoa = 0.0
  dDragdaoa = 0.0

  return BoundaryForceData{Topt, :drag}(bcnums, ndof,
                           bndry_force, isLift, lift_val, drag_val, dLiftdaoa,
                           dDragdaoa)
end

"""
  Functional for computing lift coefficient.  Uses the lift functional to
  compute the force and then divides by the (non-dimensional) dynamic pressure
  0.5*rho_free*Ma^2.  Note that this assumes the chord length (in 2d) is 1
"""
type LiftCoefficient{Topt} <: AbstractIntegralFunctional
  lift::BoundaryForceData{Topt, :lift}
  val::Topt
  bcnums::Array{Int, 1}
end

"""
  Constructor for LiftCoefficient functional
"""
function LiftCoefficientConstructor{Topt}(::Type{Topt}, mesh, sbp, eqn, opts,
                                          bcnums)

  lift = LiftForceDataConstructor(Topt, mesh, sbp, eqn, opts, bcnums)
  val = 0.0

  return LiftCoefficient{Topt}(lift, val, bcnums)
end




"""
  Type for computing the mass flow rate over a boundary (integral rho*u dot n
  dGamma)
"""
type MassFlowData{Topt} <: AbstractIntegralFunctional
  bcnums::Array{Int, 1}
  ndof::Int
  val::Topt
end

"""
  Constructor for MassFlowData.  This needs to have a different name from
  the type so it can be put in a dictionary
"""
function MassFlowDataConstructor{Topt}(::Type{Topt}, mesh, sbp, eqn, opts, 
                            bcnums)
  return MassFlowData{Topt}(bcnums, 1, 0.0)
end


"""
  Creates a functional object.

**Arguments**

 * `mesh` : Abstract PUMI mesh
 * `sbp`  : Summation-by-parts operator
 * `eqn`  : Euler equation object
 * `opts` : Options dictionary
 * `functional_name`: the name of the functional (in [`FunctionalDict`](@ref)
 * `functional_bcs`: the boundary condition numbers the functional is
                     computed on.
"""
function createFunctional{Tsol, I<:Integer}(mesh::AbstractMesh, sbp::AbstractSBP,
                                    eqn::EulerData{Tsol}, opts,
                                    functional_name::AbstractString,
                                    functional_bcs::Vector{I})

  func_constructor = FunctionalDict[functional_name]
  objective = func_constructor(Tsol, mesh, sbp, eqn, opts, functional_bcs)

  return objective
end

"""
  Maps functional names to their outer constructors.

  All outer constructors must have the signature

  MyTypeNameConstructor{Topt}(::Type{Topt}, mesh, sbp, eqn, opts, bcnums)

  where MyTypeName is the name of the type, bcnums are the
  boundary conditions that the functional is
  defined on (Array{Int, 1}), 

  Currently only boundary functionals are supported.

  For non-boundary functionals
  the region numbers associated with the faces would also be needed.
  Consider:

                | -> edge 1
     face 1     |         face2   
                |
                |

  The mesh edges lying on geometric edge 1 have two possible parent elements,
  one of face 1 and one on face 2.  `geom_regions` picks between them.  This
  effectively determines the direction of the normal vector.

  Note that the outer constructor name is the type name with the suffix "Constructor"
"""
global const FunctionalDict = Dict{ASCIIString, Function}(
"lift" => LiftForceDataConstructor,
"drag" => DragForceDataConstructor,
"massflow" => MassFlowDataConstructor,
"liftCoefficient" => LiftCoefficientConstructor
)


