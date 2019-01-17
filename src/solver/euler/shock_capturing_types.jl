# types needed for shock capturing

#------------------------------------------------------------------------------
# Persson and Peraire stuff

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
  Data for doing shock capturing using Persson and Peraire's method.
"""
mutable struct ShockPPData{Tsol, Tres}
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

  function ShockPPData{Tsol, Tres}(sbp::AbstractOperator) where {Tsol, Tres}

    Vp = VandermondeData(sbp, sbp.degree)
    Vp1 = VandermondeData(sbp, sbp.degree-1)

    # constants from Barter's thesis
    s0 = -(4 + 4.25*log10(sbp.degree))
    kappa = 0.5
    e0 = 1  # this is a bit weird, because PP says it should be O(h/p)
    
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
  shock_sensor::ShockPPData{Tsol, Tres}
  filt::Matrix{Float64}  # the filter operator

  w::Matrix{Tsol}
  t1::Matrix{Tres}
  t2::Matrix{Tres}
  Se_jac::Matrix{Tres}
  ee_jac::Matrix{Tres}
  A0inv::Matrix{Tsol}

  function ProjectionShockCapturing{Tsol, Tres}(sbp::AbstractOperator, numDofPerNode::Integer) where {Tsol, Tres}

    shock_sensor = ShockPPData{Tsol, Tres}(sbp)
    filt = zeros(Float64, sbp.numnodes, sbp.numnodes)
    getFilterOperator!(sbp, filt)

    w = zeros(Tsol, numDofPerNode, sbp.numnodes)
    t1 = zeros(Tres, numDofPerNode, sbp.numnodes)
    t2 = zeros(Tres, numDofPerNode, sbp.numnodes)

    Se_jac = zeros(Tres, numDofPerNode, sbp.numnodes)
    ee_jac = zeros(Tres, numDofPerNode, sbp.numnodes)
    A0inv = zeros(Tsol, numDofPerNode, numDofPerNode)

    return new(shock_sensor, filt, w, t1, t2, Se_jac, ee_jac, A0inv)
  end
end


