# types needed for shock capturing

#------------------------------------------------------------------------------
# Persson and Peraire stuff

"""
  Struct holding data required to transform back and form from a nodal to
  a modal solution representation.
"""
struct VandermondeData
  degree::Int
  Vp::Matrix{Float64}  # Vandermonde matrix for degree p
  Vpinv::Matrix{Float64}  # pseudo-inverse of Vp
  filt::Matrix{Float64}  # V*Vp (u_modal = Vp*u, u_filtered = V*u_modal)
  filtT::Matrix{Float64} # transpose of above
end

function VandermondeData(sbp::AbstractOperator, degree::Int)

  coords = calcnodes(sbp)
  Vp = SummationByParts.calcvandermondproriol(coords.', degree)
  Vpinv = pinv(Vp)
  filt = Vp*Vpinv
  filtT = filt.'

  return VandermondeData(degree, Vp, Vpinv, filt, filtT)
end

"""
  Shock sensor from Persson and Peraire's method, "Sub-Cell Shock Caputirng for
  Discontinuous Galerkin Methods", AIAA 2006.
"""
mutable struct ShockSensorPP{Tsol, Tres} <: AbstractShockSensor
  Vp::VandermondeData  # Vandermonde matrix for degree p and pseudo-invers
  Vp1::VandermondeData  # degree p-1

  # constants (this struct is mutable so these can be changed at runtime)
  s0::Float64
  kappa::Float64
  e0::Float64

  # storage
  up::Vector{Tsol}
  up_tilde::Vector{Tsol}
  up1_tilde::Vector{Tsol}
  
  num_dot::Vector{Tres}
  den_dot::Vector{Tres}

  function ShockSensorPP{Tsol, Tres}(mesh::AbstractMesh, sbp::AbstractOperator, opts) where {Tsol, Tres}

    Vp = VandermondeData(sbp, sbp.degree)
    Vp1 = VandermondeData(sbp, sbp.degree-1)

    # constants from Barter's thesis
    s0 = -(4 + 4.25*log10(sbp.degree))
    kappa = 0.5
    e0 = 10  # this is a bit weird, because PP says it should be O(h/p)
    
    up = zeros(Tsol, sbp.numnodes)
    up_tilde = zeros(Tsol, sbp.numnodes)
    up1_tilde = zeros(Tsol, sbp.numnodes)

    num_dot = zeros(Tres, sbp.numnodes)
    den_dot = zeros(Tres, sbp.numnodes)


    return new(Vp, Vp1, s0, kappa, e0, up, up_tilde, up1_tilde,
               num_dot, den_dot)
  end
end


"""
  Shock sensor that errors if called.  This is used when shock capturing is
  not supposed to be added.
"""
struct ShockSensorNone{Tsol, Tres} <: AbstractShockSensor

  function ShockSensorNone{Tsol, Tres}(mesh, sbp, opts) where {Tsol, Tres}
    return new()
  end
end



#------------------------------------------------------------------------------
# Projection-based shock capturing

mutable struct ProjectionShockCapturing{Tsol, Tres} <: AbstractVolumeShockCapturing
  filt::Matrix{Float64}  # the filter operator

  # storage
  w::Matrix{Tsol}
  t1::Matrix{Tres}
  t2::Matrix{Tres}
  Se_jac::Matrix{Tres}
  ee_jac::Matrix{Tres}
  A0inv::Matrix{Tsol}

  function ProjectionShockCapturing{Tsol, Tres}(mesh::AbstractMesh, sbp::AbstractOperator, opts) where {Tsol, Tres}

    numDofPerNode = mesh.numDofPerNode
    filt = zeros(Float64, sbp.numnodes, sbp.numnodes)
    getFilterOperator!(sbp, filt)
    println("max entry in filt = ", maximum(abs.(filt)))

    w = zeros(Tsol, numDofPerNode, sbp.numnodes)
    t1 = zeros(Tres, numDofPerNode, sbp.numnodes)
    t2 = zeros(Tres, numDofPerNode, sbp.numnodes)

    Se_jac = zeros(Tres, numDofPerNode, sbp.numnodes)
    ee_jac = zeros(Tres, numDofPerNode, sbp.numnodes)
    A0inv = zeros(Tsol, numDofPerNode, numDofPerNode)

    return new(filt, w, t1, t2, Se_jac, ee_jac, A0inv)
  end
end


#------------------------------------------------------------------------------
# for DG-type shock capturing, need to know interfaces

"""
  Interface that also stores the index that interface corresponds to in the
  parent mesh.
"""
struct RelativeInterface
  iface::Interface  # needed for SBP
  idx_orig::Int
  # add peer?
end


#TODO: docs
mutable struct ShockedElements{Tres}
  elnums_shock::Vector{Int}  # elements (global numbering)  where shock
                             # indicator is non-zero
  elnums_neighbor::Vector{Int}  # elements (global numbering) that neighbor
                                # elnums_shock but are in elnums_shock
  elnums_all::Vector{Int}  # union of elnums_shock and elnums_neighbor
  elnums_mesh::Vector{Int}  # temporary array, length mesh.numEl
                            # contains the indices of the selected elements
                            # in elnums_all
                            # contains the indices of elements in elnums_all
                            # TODO: can this be a BitArray?
  ee::Vector{Tres}  # the numerical viscoscity of each element in
                    # neighbor_elnums
  ifaces::Vector{RelativeInterface}

  # current indices in elnums_shock and elnums_neighbor
  idx_shock::Int

  numShock::Int  # number of elements with shocks in them.  Because we try
                 # to pre-allocate elnums_shock, its length may be greater than
                 # the number of elements with shocks
  numNeighbor::Int  # number of elements that neighbor an element with a shock
                    # in it, but don't have a shock in it
  numInterfaces::Int  # number of interfaces in ifaces
  numEl::Int  # total number of elements
  function ShockedElements{Tres}(mesh::AbstractMesh) where {Tres}

    # try to guess initial size
    size_guess = max(div(mesh.numEl, 10), 1)
    elnums_shock = zeros(Int, size_guess)
    elnums_neighbor = zeros(Int, size_guess)  # the size here is approximate
    elnums_all = Array{Int}(0)  # do this later
    elnums_mesh = zeros(mesh.numEl)
    ee = Array{Tres}(size_guess)
    ifaces = Array{Interface}(0)

    idx_shock = 1

    numShock = 0
    numNeighbor = 0
    numInterfaces = 0
    numEl = 0

    return new(elnums_shock, elnums_neighbor, elnums_all, elnums_mesh, ee,
               ifaces, idx_shock, numShock, numNeighbor, numInterfaces, numEl)
  end
end

"""
  Shock capturing using the entropy-stable varient of the Local Discontinuous
  Galerkin method.  [`allocateArrays`](@ref) must be called after the
  shock mesh is known before this type is usable.

  The sizes of these arrays may be larger than required.

  **Fields**

   * w_el: `numDofPerNode` x `numNodesPerElement` x `shockmesh.numEl` array
           for storing entropy variables
   * q_j: `numDofPerNode` x `numNodesPerElement` x `dim` x `shockmesh.numEl`
          array for storing theta_j and q_j
   * convert_entropy: a function for converting conservative variables to
                      entropy variables, must have signature
                      `convert_entropy(params::ParamType, q::AbstractVector,
                                       w::AbstractVector)`, where `q` and `w`
                      are of length `numDofPerNode`.  `w` should be overwritten
   * `flux`: an [`AbstractLDGFlux`](@ref)
   * diffusion: an [`AbstractDiffusion`](@ref).
"""
mutable struct LDGShockCapturing{Tsol, Tres} <: AbstractFaceShockCapturing
  # Note: the variable names are from Chen and Shu's "Entropy Stable High Order
  #       DG Methods" paper
  w_el::Array{Tsol, 3}  # entropy variables for all elements in elnums_all
  q_j::Array{Tres, 4}  # auxiliary equation solution for elements in
                       # elnums_all (this gets used for both theta and q)
  convert_entropy::Any  # convert to entropy variables
  flux::AbstractLDGFlux
  diffusion::AbstractDiffusion
  function LDGShockCapturing{Tsol, Tres}() where {Tsol, Tres}
    # we don't know the right sizes yet, so just make them zero size
    w_el = Array{Tsol}(0, 0, 0)
    q_j = Array{Tsol}(0, 0, 0, 0)

    # default values
    convert_entropy = convertToIR_
    flux = LDG_ESFlux()
    diffusion = ShockDiffusion{Tres}()

    return new(w_el, q_j, convert_entropy, flux, diffusion)
  end

  function LDGShockCapturing{Tsol, Tres}(mesh::AbstractMesh, sbp::AbstractOperator, opts) where {Tsol, Tres}
    return LDGShockCapturing{Tsol, Tres}()
  end

end


"""
  This function allocates the arrays of [`LDGShockCapturing`](@ref) that
  depend on the shock mesh.  The arrays are only re-allocated if they are
  too small.

  **Inputs**

   * capture: [`LDGShockCapturing`](@ref)
   * mesh
   * shockmesh: a `ShockedElements` object, fully initialized
"""
function allocateArrays(capture::LDGShockCapturing{Tsol, Tres}, mesh::AbstractMesh,
                        shockmesh::ShockedElements) where {Tsol, Tres}

  # can't resize multi-dimension arrays, so reallocate
  if size(capture.q_j, 4) < shockmesh.numEl
    capture.q_j = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement,
                              mesh.dim, shockmesh.numEl)
    # the final dimension must be numEl, not numShock, because q_j is zero
    # for the neighboring elements later
  end

  if size(capture.w_el, 3) < shockmesh.numEl
    capture.w_el = Array{Tsol}(mesh.numDofPerNode, mesh.numNodesPerElement,
                               shockmesh.numEl)
  end

  setDiffusionArray(capture.diffusion, shockmesh.ee)

  return nothing
end


#------------------------------------------------------------------------------
# Viscoscity and fluxes for LDG
"""
  Diagonal viscoscity (constant for each element), used for shock capturing
"""
mutable struct ShockDiffusion{T} <: AbstractDiffusion
  ee::Vector{T}
  function ShockDiffusion{T}() where {T}
    ee = Vector{T}(0)
    return new(ee)
  end

  function ShockDiffusion{T}(ee::AbstractVector{T}) where {T}
    return new(ee)
  end
end

function ShockDiffusion(ee::AbstractVector{T}) where {T}
  return ShockDiffusion{T}(ee)
end

"""
  Function to set the elementwise diffusion coefficient
"""
function setDiffusionArray(obj::ShockDiffusion, vals::AbstractArray)
  obj.ee = vals
end


"""
  Entropy stable LDG flux
"""
struct LDG_ESFlux  <: AbstractLDGFlux
  alpha::Float64
  beta::Float64  # Chen and Shu say this can be an arbitrary vector, but my
                 # analysis says all entries must the the same, so it acts
                 # like a scalar
  function LDG_ESFlux()
    alpha = 1
    beta = 1

    return new(alpha, beta)
  end
end


"""
  Shock capturing type that errors out.  Used when shock capturing is not
  supposed to be added.
"""
struct ErrorShockCapturing{Tsol, Tres} <: AbstractShockCapturing

  function ErrorShockCapturing{Tsol, Tres}(mesh::AbstractMesh, sbp::AbstractSBP,
                                           opts) where {Tsol, Tres}
    return new()
  end
end


#------------------------------------------------------------------------------
# BR2 Shock Capturing

"""
  Abstract type for the different diffusion penalties that can be expressed
  in the framework of Yan et.al. "Interior Penalties for Summation-by-Parts
  Discretization os Linear Second-Order Differential Equations"
"""
abstract type AbstractDiffusionPenalty end

struct BR2Penalty{Tsol, Tres} <: AbstractDiffusionPenalty
end

"""
  BR2 shock capturing

  **Fields**

   * grad_w: stores Lambda * D * q at volume nodes, numDofPerNode x
             numNodesPerElement x dim x shockmesh.numEl
"""
mutable struct SBPParabolicSC{Tsol, Tres} <: AbstractShockCapturing
  w_el::Array{Tsol, 3}
  grad_w::Array{Tres, 4}
  convert_entropy::Any  # convert to entropy variables
  flux::AbstractLDGFlux
  diffusion::AbstractDiffusion
  penalty::AbstractDiffusionPenalty

  function SBPParabolicSC{Tsol, Tres}() where {Tsol, Tres}
    # we don't know the right sizes yet, so just make them zero size
    w_el = Array{Tsol}(0, 0, 0)
    grad_w = Array{Tres}(0, 0, 0, 0)

    # default values
    convert_entropy = convertToIR_
    flux = LDG_ESFlux()
    diffusion = ShockDiffusion{Tres}()
    penalty = BR2Penalty{Tsol, Tres}()

    return new(w_el, grad_w, convert_entropy, flux, diffusion, penalty)
  end

  function SBPParabolicSC{Tsol, Tres}(mesh::AbstractMesh, sbp::AbstractOperator, opts) where {Tsol, Tres}
    return SBPParabolicSC{Tsol, Tres}()
  end


end

function allocateArrays(capture::SBPParabolicSC{Tsol, Tres}, mesh::AbstractMesh,
                        shockmesh::ShockedElements) where {Tsol, Tres}

  # can't resize multi-dimension arrays, so reallocate
  if size(capture.grad_w, 4) < shockmesh.numEl
    capture.grad_w = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement,
                                 mesh.dim, shockmesh.numEl)
  end

  if size(capture.w_el, 3) < shockmesh.numEl
    capture.w_el = Array{Tsol}(mesh.numDofPerNode, mesh.numNodesPerElement,
                               shockmesh.numEl)
  end

  setDiffusionArray(capture.diffusion, shockmesh.ee)

  return nothing
end



#------------------------------------------------------------------------------
# Creating shock sensors

global const ShockSensorDict = Dict{String, Type{T} where T <: AbstractShockSensor}(
"SensorNone" => ShockSensorNone,
"SensorPP" => ShockSensorPP,
)

function getShockSensor(mesh, sbp, eqn::EulerData{Tsol, Tres}, opts) where {Tsol, Tres}

  name = opts["shock_sensor_name"]
  obj = ShockSensorDict[name]{Tsol, Tres}(mesh, sbp, opts)
  eqn.shock_sensor = obj

  return nothing
end


#------------------------------------------------------------------------------
# Creating shock capturing

global const ShockCapturingDict = Dict{String, Type{T} where T <: AbstractShockCapturing}(
"ShockCapturingNone" => ErrorShockCapturing,
"ElementProjection" => ProjectionShockCapturing,
"LDG" => LDGShockCapturing
)


"""
  This function populates the `eqn.shock_capturing` field of the mesh with the
  shock capturing object for the scheme.  The shock capturing scheme is
  determined by opts["shock_capturing_name"]
"""
function getShockCapturing(mesh, sbp, eqn::EulerData{Tsol, Tres}, opts) where {Tsol, Tres}

  name = opts["shock_capturing_name"]
  obj = ShockCapturingDict[name]{Tsol, Tres}(mesh, sbp, opts)
  eqn.shock_capturing = obj

  return nothing
end
