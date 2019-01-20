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

  function ShockSensorPP{Tsol, Tres}(sbp::AbstractOperator) where {Tsol, Tres}

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


#------------------------------------------------------------------------------
# Projection-based shock capturing

mutable struct ProjectionShockCapturing{Tsol, Tres} <: AbstractShockCapturing
  filt::Matrix{Float64}  # the filter operator

  # storage
  w::Matrix{Tsol}
  t1::Matrix{Tres}
  t2::Matrix{Tres}
  Se_jac::Matrix{Tres}
  ee_jac::Matrix{Tres}
  A0inv::Matrix{Tsol}

  function ProjectionShockCapturing{Tsol, Tres}(sbp::AbstractOperator, numDofPerNode::Integer) where {Tsol, Tres}

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

# This is needed because of Julia issue 24909
mutable struct PushVector{A <: AbstractVector{T}} <: AbstractVector{T}
  data::A
  len::Int
end

function PushVector(A::T) where {T<:AbstractVector}

  return PushVector{T}(A, length(A))
end

function PushVector{T}(len::Integer=0)
  data = Array{T}(len)
  return PushVector(data)
end

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
  function ShockedElement{Tres}(mesh::AbstractMesh) where {Tres}

    # try to guess initial size
    size_guess = max(div(mesh.numEl, 10), 1)
    elnums_shock = zeros(Int, size_guess)
    elnums_neighbor = zeros(Int, size_guess)  # the size here is approximate
    elnums_all = Array{Int}(0)  # do this later
    elnums_mesh = zeros(mesh.numEl)
    ee = Array{Tsol}(size_guess)
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
  Galerkin method.
"""
mutable struct LDGShockCapturing{Tsol, Tres} <: AbstractShockCapturing
  # Note: the variable names are from Chen and Shu's "Entropy Stable High Order
  #       DG Methods" paper
  w_el::Array{Tsol, 3}  # entropy variables for all elements in elnums_all
  q_j::Array{Tres, 4}  # auxiliary equation solution for elements in
                       # elnums_all (this gets used for both theta and q)
  function LDGShockCapturing{Tsol, Tres}()
    # we don't know the right sizes yet, so just make them zero size
    w_el = Array{Tsol}(0, 0, 0)
    q_j = Array{Tsol}(0, 0, 0, 0)

    return new(w_el, q_j)
  end
end




"""
  Diagonal viscoscity (constant for each element), used for shock capturing
"""
struct ShockDiffusion{T} <: AbstractDiffusion
  ee::Vector{T}
end

function ShockDiffusion(ee::AbstractVector{T}) where {T}
  return ShockDiffusion{T}(ee)
end


"""
  Entropy stable LDG flux
"""
struct LDG_ESFlux
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


