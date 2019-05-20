# functional definitions

import PDESolver: createFunctional, _evalFunctional, _evalFunctionalDeriv_m,
                 _evalFunctionalDeriv_q
import ODLCommonTools: getParallelData, setupFunctional


#------------------------------------------------------------------------------
# Lift and drag

@doc """
###EulerEquationMod.BoundaryForceData

Composite data type for storing data pertaining to the boundaryForce. It holds
lift and drag values

"""->

mutable struct BoundaryForceData{Topt, fname} <: AbstractBoundaryFunctional{Topt}
  bcnums::Array{Int,1}

  # factors to multiply x, y, z momenta by, determines if lift or drag is
  # calculated
  facx::Topt
  facy::Topt
  facz::Topt

  # things needed for the calculation
  qg::Vector{Topt}
  euler_flux::Vector{Topt}
end

"""
  Constructor for BoundaryForceData{Topt, :lift}
"""
function LiftForceDataConstructor(::Type{Topt}, mesh, sbp, eqn, opts, bcnums) where Topt

  facx = 0
  facy = 0
  facz = 0

  qg = zeros(Topt, mesh.numDofPerNode)
  euler_flux = zeros(Topt, mesh.numDofPerNode)

  return BoundaryForceData{Topt, :lift}(bcnums, facx, facy, facz, qg,
                                        euler_flux)
end

"""
  Constructor for BoundaryForceData{Topt, :drag}
"""
function DragForceDataConstructor(::Type{Topt}, mesh, sbp, eqn, opts,
                             bcnums) where Topt

  facx = 0
  facy = 0
  facz = 0

  qg = zeros(Topt, mesh.numDofPerNode)
  euler_flux = zeros(Topt, mesh.numDofPerNode)

  return BoundaryForceData{Topt, :drag}(bcnums, facx, facy, facz, qg,
                                        euler_flux)
end


function setupFunctional(mesh::AbstractMesh, sbp, eqn::AbstractSolutionData,
                         opts::Dict, func::BoundaryForceData{Topt, :lift}) where {Topt}

  if mesh.dim == 2
    func.facx = -sin(eqn.params.aoa)
    func.facy =  cos(eqn.params.aoa)
  else
    func.facx = -sin(eqn.params.aoa)
    func.facy =  0
    func.facz =  cos(eqn.params.aoa)
  end

  return nothing
end

function setupFunctional(mesh::AbstractMesh, sbp, eqn::AbstractSolutionData,
                         opts::Dict, func::BoundaryForceData{Topt, :drag}) where {Topt}

  if mesh.dim == 2
    func.facx = cos(eqn.params.aoa)
    func.facy = sin(eqn.params.aoa)
  else
    func.facx = cos(eqn.params.aoa)
    func.facy = 0
    func.facz = sin(eqn.params.aoa)
  end

  return nothing
end


"""
  Functional for computing lift coefficient.  Uses the lift functional to
  compute the force and then divides by the (non-dimensional) dynamic pressure
  0.5*rho_free*Ma^2.  Note that this assumes the chord length (in 2d) is 1
"""
mutable struct LiftCoefficient{Topt} <: AbstractBoundaryFunctional{Topt}
  func::BoundaryForceData{Topt, :lift}
  bcnums::Array{Int, 1}
end

"""
  Constructor for LiftCoefficient functional
"""
function LiftCoefficientConstructor(::Type{Topt}, mesh, sbp, eqn, opts,
                                    bcnums) where Topt

  lift = LiftForceDataConstructor(Topt, mesh, sbp, eqn, opts, bcnums)

  return LiftCoefficient{Topt}(lift, bcnums)
end



"""
  Functional for computing drag coefficient.  See [`LiftCoefficient`](@ref)
  for the non-dimensionalization.
"""
mutable struct DragCoefficient{Topt} <: AbstractBoundaryFunctional{Topt}
  func::BoundaryForceData{Topt, :drag}
  bcnums::Array{Int, 1}
end

"""
  Constructor for DragCoefficient functional
"""
function DragCoefficientConstructor(::Type{Topt}, mesh, sbp, eqn, opts,
                                    bcnums) where Topt

  drag = DragForceDataConstructor(Topt, mesh, sbp, eqn, opts, bcnums)

  return DragCoefficient{Topt}(drag, bcnums)
end


"""
  Typedef for functionals that compute aerodynamic coefficients (lift, drag)
"""
const AeroCoefficients{Topt} = Union{LiftCoefficient{Topt}, DragCoefficient{Topt}}

function setupFunctional(mesh::AbstractMesh, sbp, eqn::AbstractSolutionData,
                         opts::Dict, func::AeroCoefficients)

  setupFunctional(mesh, sbp, eqn, opts, func.func)

end


#------------------------------------------------------------------------------
# other functionals


"""
  Type for computing the mass flow rate over a boundary (integral rho*u dot n
  dGamma)
"""
mutable struct MassFlowData{Topt} <: AbstractBoundaryFunctional{Topt}
  bcnums::Array{Int, 1}
end

"""
  Constructor for MassFlowData.  This needs to have a different name from
  the type so it can be put in a dictionary
"""
function MassFlowDataConstructor(::Type{Topt}, mesh, sbp, eqn, opts, 
                            bcnums) where Topt
  return MassFlowData{Topt}(bcnums)
end

"""
  Type for computing the entropy flux rate over a boundary (integral S * u_i dot n_i
  dGamma), where S is the entropy function (from the IR entropy variables).
"""
mutable struct EntropyFluxData{Topt} <: AbstractBoundaryFunctional{Topt}
  bcnums::Array{Int, 1}
end

"""
  Constructor for EntropyFlux.  This needs to have a different name from
  the type so it can be put in a dictionary
"""
function EntropyFluxConstructor(::Type{Topt}, mesh, sbp, eqn, opts, 
                            bcnums) where Topt
  return EntropyFluxData{Topt}(bcnums)
end

"""
  Functional that computes function `w^T d(u)`, where `w` are the entropy
  variables and `d(u)` is the entropy stable dissipation computed by
  [`ELFPenaltyFaceIntegral`](@ref).  Note that this is an integral
  over interior faces and not boundary faces
"""
mutable struct EntropyDissipationData{Topt} <: EntropyPenaltyFunctional{Topt}
  func::ELFPenaltyFaceIntegral
  func_sparseface::FluxType
  func_sparseface_revq::FluxType_revq
  func_sparseface_revm::FluxType_revm
end

"""
  Constructor for [`EntropyDissipationData`](@ref)

  This function takes `bcnums` as an argument for consistency with 
  the boundary functional constructors, but doesn't use it.
"""
function EntropyDissipationConstructor(::Type{Topt}, mesh, sbp, eqn, opts,
                                       bcnums) where Topt

  func = ELFPenaltyFaceIntegral(mesh, eqn)
  func_sparseface = LFPenalty()
  func_sparseface_revq = LFPenalty_revq()
  func_sparseface_revm = LFPenalty_revm()

  return EntropyDissipationData{Topt}(func, func_sparseface, func_sparseface_revq, func_sparseface_revm)
end


"""
  Returns the negative of [`EntropyDissipationData`](@ref).  That functional
  is always negative, so this one is always positive.
"""
mutable struct NegEntropyDissipationData{Topt} <: EntropyPenaltyFunctional{Topt}
  func::EntropyDissipationData{Topt}
end

"""
  Constructor for [`NegEntropyDissipationData`](@ref).  `bcnums` argument is
  unused.
"""
function NegEntropyDissipationConstructor(::Type{Topt}, mesh, sbp, eqn, opts,
                                       bcnums) where Topt

  func = EntropyDissipationConstructor(Topt, mesh, sbp, eqn, opts, bcnums)
  return NegEntropyDissipationData{Topt}(func)
end


"""
  Functional that computes function `w^T d(u)`, where `w` are the entropy
  variables and `d(u)` is the entropy stable dissipation computed by
  [`ELFPenaltyFaceIntegral`](@ref).  Note that this is an integral
  over interior faces and not boundary faces
"""
mutable struct EntropyJumpData{Topt} <: EntropyPenaltyFunctional{Topt}
  func::EntropyJumpPenaltyFaceIntegral
end

"""
  Constructor for [`EntropyDissipationData`](@ref)

  This function takes `bcnums` as an argument for consistency with 
  the boundary functional constructors, but doesn't use it.
"""
function EntropyJumpConstructor(::Type{Topt}, mesh, sbp, eqn, opts,
                                       bcnums) where Topt

  func = EntropyJumpPenaltyFaceIntegral(mesh, eqn)

  return EntropyJumpData{Topt}(func)
end


"""
  Computes entropy dissipation by computing
  
  -W^T * R_face + \int_\Gamma_interior (psi_L - psi_R)

  This works for flux functions that do not have the EC + dissipation form,
  and should be < 0 for any entropy stable scheme.
"""
mutable struct EntropyDissipation2Data{Topt} <: EntropyPenaltyFunctional{Topt}
  flux_functional::EntropyDissipationData{Topt}
  potential_flux::PotentialFlux
  potential_flux_revq::PotentialFlux_revq
  potential_flux_revm::PotentialFlux_revm
end

function EntropyDissipation2Constructor(::Type{Topt}, mesh, sbp, eqn, opts,
                                       bcnums) where Topt

  flux_functional = EntropyDissipationConstructor(Topt, mesh, sbp, eqn, opts,
                                                  bcnums)

  # configure the functional to compute (wL - wR)^T F, where F is the flux
  # currently used by the discretization
  flux_functional.func_sparseface = eqn.flux_func
  flux_functional.func_sparseface_revq = eqn.flux_func_revq
  flux_functional.func_sparseface_revm = eqn.flux_func_revm

  potential_flux = FluxDict["PotentialFlux"]
  potential_flux_revq = FluxDict_revq["PotentialFlux"]
  potential_flux_revm = FluxDict_revm["PotentialFlux"]

  return EntropyDissipation2Data{Topt}(flux_functional, potential_flux,
                                       potential_flux_revq, potential_flux_revm)
end

"""
  Computes negative of EntropyDissipation2
"""
mutable struct NegEntropyDissipation2Data{Topt} <: EntropyPenaltyFunctional{Topt}
  func::EntropyDissipation2Data{Topt}
end

function NegEntropyDissipation2Constructor(::Type{Topt}, mesh, sbp, eqn, opts,
                                           bcnums) where Topt
  func = EntropyDissipation2Constructor(Topt, mesh, sbp, eqn, opts, bcnums)
  return NegEntropyDissipation2Data{Topt}(func)
end


function getParallelData(obj::EntropyPenaltyFunctional)
  return PARALLEL_DATA_ELEMENT
end

const NegEntropyDissipations = Union{NegEntropyDissipationData,
                                     NegEntropyDissipation2Data}


"""
  Functional that computes -delta_w^T lambda_max Y^T Y delta_w for
  boundary faces.  delta_w is the difference between the (entropy variable)
  numerical solution and the boundary state.  This only works for boundary
  conditions that have [`getDirichletState`](@ref) defined.
"""
mutable struct BoundaryEntropyDiss{Topt} <: AbstractBoundaryFunctional{Topt}
  bcnums::Array{Int, 1}
end

function BoundaryEntropyDissConstructor(::Type{Topt}, mesh, sbp, eqn, opts, 
                            bcnums) where Topt
  return BoundaryEntropyDiss{Topt}(bcnums)
end

"""
  Negative of [`BoundaryEntropyDiss`](@ref).  That one is always negative
  so this one is always positive.
"""
mutable struct NegBoundaryEntropyDiss{Topt} <: AbstractBoundaryFunctional{Topt}
  func::BoundaryEntropyDiss{Topt}
end

function NegBoundaryEntropyDissConstructor(::Type{Topt}, mesh, sbp, eqn, opts, 
                            bcnums) where Topt
  obj = BoundaryEntropyDissConstructor(Topt, mesh, sbp, eqn, opts, bcnums)
  return NegBoundaryEntropyDiss{Topt}(obj)
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
function createFunctional(mesh::AbstractMesh, sbp::AbstractOperator,
                  eqn::EulerData{Tsol}, opts,
                  functional_name::AbstractString,
                  functional_bcs::Vector{I}) where {Tsol, I<:Integer}

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
global const FunctionalDict = Dict{String, Function}(
"lift" => LiftForceDataConstructor,
"drag" => DragForceDataConstructor,
"liftCoefficient" => LiftCoefficientConstructor,
"dragCoefficient" => DragCoefficientConstructor,
"massflow" => MassFlowDataConstructor,
"entropyflux" => EntropyFluxConstructor,
"entropydissipation" => EntropyDissipationConstructor,
"entropydissipation2" => EntropyDissipation2Constructor,
"negentropydissipation" => NegEntropyDissipationConstructor,
"negentropydissipation2" => NegEntropyDissipation2Constructor,
"entropyjump" => EntropyJumpConstructor,
"boundaryentropydiss" => BoundaryEntropyDissConstructor,
"negboundaryentropydiss" => NegBoundaryEntropyDissConstructor,
)


