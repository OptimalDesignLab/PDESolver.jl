# types needed for shock capturing

#------------------------------------------------------------------------------
# Persson and Peraire stuff

struct VandermondeData
  degree::Int
  Vp::Matrix{Float64}  # Vandermonde matrix for degree p
  Vpinv::Matrix{Float64}  # pseudo-inverse of Vp
  filt::Matrix{Float64}  # V*Vp (u_modal = Vp*u, u_filtered = V*u_modal)
end

function VandermondeData(sbp::AbstractOperator, degree::Int)

  coords = calcnodes(sbp)
  Vp = calcvandermondproriol(coords, degree)
  Vpinv = pinv(Vp)
  filt = Vp*Vpinv

  return VandermondeData(degree, Vp, Vpinv, filt)
end

"""
  Data for doing shock capturing using Persson and Peraire's method.
"""
mutable struct ShockPPData{Tsol}
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

  function ShockPPData{Tsol}(sbp::AbstractOperator) where {Tsol}

    Vp = VandermondeData(sbp, sbp.degree)
    Vp1 = VandermondeData(sbp, sbp.degree-1)

    # constants from Barter's thesis
    s0 = -(4 + 4.25*log10(sbp.degree))
    kappa = 0.5
    e0 = 1  # this is a bit weird, because PP says it should be O(h/p)
    
    up = zeros(Tsol, sbp.numnodes)
    up_tilde = zeros(Tsol, sbp.numnodes)
    up1_tilde = zeros(Tsol, sbp.numnodes)

    return new(Vp, Vp1, s0, kappa, e0, up, up_tilde, up1_tilde)
  end
end


#------------------------------------------------------------------------------
# Projection-based shock capturing

mutable struct ProjectionShockCapturing{Tsol, Tres} <: AbstractShockCapturing
  shock_sensor::ShockPPData{Tsol}
  filt::Matrix{Float64}  # the filter operator

  t1::Matrix{Tres}
  t2::Matrix{Tres}
  function ProjectionShockCapturing{Tsol, Tres}(sbp::AbstractOperator, numDofPerNode::Integer) where {Tsol, Tres}

    shock_sensor = ShockPPData{Tsol}(sbp)
    filt = zeros(Float64, sbp.numnodes, sbp.numnodes)
    getFilterOperator(sbp, filt)

    t1 = zeros(Tres, numDofPerNode, sbp.numnodes)
    t2 = zeros(Tres, numDofPerNode, sbp.numnodes)

    return new(shock_sensor, filt, t1, t2)
  end
end


