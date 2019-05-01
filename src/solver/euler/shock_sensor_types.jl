# different shock sensors

#------------------------------------------------------------------------------
# Persson and Peraire stuff

"""
  Struct holding data required to transform back and form from a nodal to
  a modal solution representation.
"""
struct VandermondeData
  degree::Int
  nmodes::Int # number of modes
  nmodes1::Int  # number of modes in degree - 1 basis
  Vp::Matrix{Float64}  # Vandermonde matrix for degree p
  Vpinv::Matrix{Float64}  # pseudo-inverse of Vp
  filt::Matrix{Float64}  # Vp*Vpinv (u_modal = Vpinv*u, u_filtered = Vp*u_modal)
  filtT::Matrix{Float64} # transpose of above
  filt1::Matrix{Float64} # filter operator that only projects the modes from a
                         # degree-1 operator back to the nodal basis
  filt1T::Matrix{Float64} # transpose of above

  function VandermondeData(sbp::AbstractOperator, degree::Integer)

    coords = calcnodes(sbp)
    Vp = SummationByParts.calcvandermondproriol(coords.', degree)
    nmodes = size(Vp, 2)
    Vpinv = pinv(Vp)
    filt = Vp*Vpinv
    filtT = filt.'

    # this is a rather inefficient way to compute the number of modes
    Vp1 = SummationByParts.calcvandermondproriol(coords.', degree - 1)
    nmodes1 = size(Vp1, 2)
    nmodes_diff = nmodes - nmodes1

    filt1 = Vp[:, 1:(end-nmodes_diff)]*Vpinv[1:(end-nmodes_diff), :]
    filt1T = filt1.'

    return new(degree, nmodes, nmodes1, Vp, Vpinv, filt, filtT,
                           filt1, filt1T)
  end


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
  _e0::Float64  # original e0 value
  e0::Float64  # e0 scaled by alpha

  # storage
  up::Vector{Tsol}
  up_tilde::Vector{Tsol}
  up1_tilde::Vector{Tsol}
  
  num_dot::Vector{Tres}
  den_dot::Vector{Tres}
  ee_dot::Vector{Tres}
  lambda_max_dot::Matrix{Tres}

  up_tilde_bar::Vector{Tsol}
  up1_tilde_bar::Vector{Tsol}
  up_bar::Vector{Tsol}

  function ShockSensorPP{Tsol, Tres}(mesh::AbstractMesh, sbp::AbstractOperator, opts) where {Tsol, Tres}

    Vp = VandermondeData(sbp, sbp.degree)
    Vp1 = VandermondeData(sbp, sbp.degree-1)  #TODO: unneded?

    # constants from Barter's thesis
    s0 = -(4 + 4.25*log10(sbp.degree))  # was -(4 + 4.25*log10(sbp.degree))
    kappa = 0.5  # was 0.5
    _e0 = 0.01
    e0 = _e0
    
    up = zeros(Tsol, sbp.numnodes)
    up_tilde = zeros(Tsol, sbp.numnodes)
    up1_tilde = zeros(Tsol, sbp.numnodes)

    num_dot = zeros(Tres, sbp.numnodes)
    den_dot = zeros(Tres, sbp.numnodes)
    ee_dot = zeros(Tres, mesh.numNodesPerElement)
    lambda_max_dot = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)

    up_tilde_bar = zeros(Tsol, mesh.numNodesPerElement)
    up1_tilde_bar = zeros(Tsol, mesh.numNodesPerElement)
    up_bar = zeros(Tsol, mesh.numNodesPerElement)

    return new(Vp, Vp1, s0, kappa, _e0, e0, up, up_tilde, up1_tilde,
               num_dot, den_dot, ee_dot, lambda_max_dot,
               up_tilde_bar, up1_tilde_bar, up_bar)
  end
end


function setAlpha(obj::ShockSensorPP, alpha::Number)

  obj.e0 = obj._e0*alpha
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


function setAlpha(obj::ShockSensorNone, alpha::Number)

end

#------------------------------------------------------------------------------
# Sensor for testing only: there is a shock everywhere

"""
  Shock sensor that always says there is a shock and returns a viscoscity 
  of 1
"""
mutable struct ShockSensorEverywhere{Tsol, Tres} <: AbstractShockSensor
  alpha::Float64

  function ShockSensorEverywhere{Tsol, Tres}(mesh::AbstractMesh, sbp::AbstractSBP, opts) where {Tsol, Tres}
    alpha = 1.0

    return new(alpha)
  end
end

function setAlpha(obj::ShockSensorEverywhere, alpha::Number)

  obj.alpha = alpha
end




#------------------------------------------------------------------------------
# ShockSensoryVelocity

"""
  Shock sensor where ee[d, i] = q[d+1, i] (ie. the viscoscity in direction
  `d` is the momentum in direction `d`.  This is useful for testing a non
  constant viscoscity with non-zero derivative wrt q.
"""
mutable struct ShockSensorVelocity{Tsol, Tres} <: AbstractShockSensor
  alpha::Float64

  function ShockSensorVelocity{Tsol, Tres}(mesh::AbstractMesh, sbp::AbstractSBP, opts) where {Tsol, Tres}
    alpha = 1.0

    return new(alpha)
  end
end

function setAlpha(obj::ShockSensorVelocity, alpha::Number)

  obj.alpha = alpha
end




#------------------------------------------------------------------------------
# Hartmanns Isotropic sensor

"""
  Struct for computing the norm of the strong form residual
"""
struct StrongFormData{Tsol, Tres}

  flux::Array{Tres, 3}
  nrm::Vector{Tres}
  aux_vars::Vector{Tres}
  work::Array{Tres, 3}

  flux_jac::Array{Tres, 5}
  Dx::Array{Float64, 3}

  flux_bar::Array{Tres, 3}

  function StrongFormData{Tsol, Tres}(mesh, sbp, opts) where {Tsol, Tres}

    numDofPerNode = mesh.numDofPerNode
    numNodesPerElement = mesh.numNodesPerElement
    dim = mesh.dim

    flux = zeros(Tres, numDofPerNode, numNodesPerElement, dim)
    nrm = zeros(Tres, dim)
    aux_vars = Tres[]
    work = zeros(Tres, numDofPerNode, numNodesPerElement, dim)

    # this only needs to be 4D, but applyOperatorJac doesn't have a method
    # for that
    flux_jac = zeros(Tres, numDofPerNode, numDofPerNode, dim,
                   numNodesPerElement, numNodesPerElement)
    Dx = zeros(numNodesPerElement, numNodesPerElement, dim)

    flux_bar = zeros(Tres, numDofPerNode, numNodesPerElement, dim)

    return new(flux, nrm, aux_vars, work,
               flux_jac, Dx,
               flux_bar)
  end
end





"""
  Sensor from Hartmann's "Adaptive Discontinuous Galerkin Finite Element
  Methods for the Compressible Euler Equations".  This sensor is suitable for
  isotropic meshes only.  See the paper "for the Compressible Navier-Stokes
  Equations" for the anisotropic one.
"""
mutable struct ShockSensorHIso{Tsol, Tres} <: AbstractShockSensor
  _C_eps::Float64  # original value
  C_eps::Float64  # value scaled by alpha
  beta::Float64

  strongdata::StrongFormData{Tsol, Tres}
  res::Matrix{Tres}

  res_jac::Array{Tres, 4}

  function ShockSensorHIso{Tsol, Tres}(mesh::AbstractMesh, sbp::AbstractSBP,
                                       opts) where {Tsol, Tres}
    numDofPerNode = mesh.numDofPerNode
    numNodesPerElement = mesh.numNodesPerElement

    _C_eps = 1/25  # was 1/25
    C_eps = _C_eps
    beta = 1/10
    strongdata = StrongFormData{Tsol, Tres}(mesh, sbp, opts)
    res = zeros(Tres, numDofPerNode, numNodesPerElement)
    res_jac = zeros(Tres, numDofPerNode, numDofPerNode, numNodesPerElement,
                        numNodesPerElement)

    return new(_C_eps, C_eps, beta, strongdata, res, res_jac)
  end
end


function setAlpha(obj::ShockSensorHIso, alpha::Number)

  obj.C_eps = obj._C_eps*alpha
end



#------------------------------------------------------------------------------
# Odens gradient based sensor

"""
  Shock sensor from Baumann and Oden's "An Adaptive-Order Discontinuous
  Galerkin Method for the Solution of the Euler Equations of Gas Dynamics".

  Basically a gradient sensor.
"""
mutable struct ShockSensorBO{Tsol, Tres} <: AbstractShockSensor
  _alpha::Float64  # arbitrary coeffcieint
  alpha::Float64  # scaled by user-defined alpha
  function ShockSensorBO{Tsol, Tres}(mesh::AbstractMesh, sbp::AbstractSBP,
                                       opts) where {Tsol, Tres}
    _alpha = 0.01  # was 1.0
    alpha = _alpha

    return new(_alpha, alpha)
  end
end


function setAlpha(obj::ShockSensorBO, alpha::Number)

  obj.alpha = obj._alpha*alpha
end



#------------------------------------------------------------------------------
# Hartmann's high order sensor

"""
  This is an approximate version of the sensor from Hartmann's paper:
  "Higher-order and adaptive discontinuous Galerkin methods with
   shock-capturing applied to transonic turbulent delta wing flow"
"""
mutable struct ShockSensorHHO{Tsol, Tres} <: AbstractShockSensor
  _C_eps::Float64  # original value
  C_eps::Float64  # value scaled by alpha

  strongdata::StrongFormData{Tsol, Tres}
  h_k_tilde::Matrix{Tres}
  h_k_tilde_bar::Matrix{Tres}

  p_dot::Vector{Tsol}
  press_el::Matrix{Tsol}
  press_dx::Array{Tres, 3}
  work::Array{Tres, 3}
  res::Array{Tres, 2}
  Rp::Vector{Tres}
  fp::Vector{Tres}


  p_jac::Array{Tsol, 3}
  res_jac::Array{Tres, 4}
  p_hess::Matrix{Tsol}
  Dx::Array{Tres, 3}
  px_jac::Array{Tres, 5}
  val_dot::Matrix{Tres}
  Rp_jac::Array{Tres, 3}
  fp_jac::Array{Tres, 3}

  fp_bar::Vector{Tres}
  Rp_bar::Vector{Tres}
  press_el_bar::Matrix{Tres}
  press_dx_bar::Array{Tres, 3}
  p_dot_bar::Vector{Tres}
  res_bar::Matrix{Tres}

  function ShockSensorHHO{Tsol, Tres}(mesh::AbstractMesh, sbp::AbstractSBP,
                                       opts) where {Tsol, Tres}

    @unpack mesh numDofPerNode numNodesPerElement dim

    _C_eps = 1/5  # was 1/5
    C_eps = _C_eps
    strongdata = StrongFormData{Tsol, Tres}(mesh, sbp, opts)
    h_k_tilde = Array{Tres}(mesh.dim, mesh.numEl)
    h_k_tilde_bar = Array{Tres}(mesh.dim, mesh.numEl)
    calcAnisoFactors(mesh, sbp, opts, h_k_tilde)

    p_dot = zeros(Tsol, numDofPerNode)
    press_el = zeros(Tsol, 1, numNodesPerElement)
    press_dx = zeros(Tres, 1, numNodesPerElement, dim)
    work = zeros(Tres, 1, numNodesPerElement, dim)
    res = zeros(Tres, numDofPerNode, numNodesPerElement)
    Rp = zeros(Tres, numNodesPerElement)
    fp = zeros(Tres, numNodesPerElement)


    p_jac = zeros(Tsol, 1, numDofPerNode, numNodesPerElement)
    res_jac = zeros(Tres, numDofPerNode, numDofPerNode, numNodesPerElement,
                          numNodesPerElement)
    p_hess = zeros(Tsol, numDofPerNode, numDofPerNode)
    Dx = zeros(numNodesPerElement, numNodesPerElement, dim)
    px_jac = zeros(Tres, 1, numDofPerNode, dim, numNodesPerElement,
                         numNodesPerElement)
    val_dot = zeros(Tres, numDofPerNode, numNodesPerElement)

    Rp_jac = zeros(Tres, numNodesPerElement, numDofPerNode, numNodesPerElement)
    fp_jac = zeros(Tres, numDofPerNode, numNodesPerElement, numNodesPerElement)

    fp_bar = zeros(Tres, numNodesPerElement)
    Rp_bar = zeros(Tres, numNodesPerElement)
    press_el_bar = zeros(Tres, 1, numNodesPerElement)
    press_dx_bar = zeros(Tres, 1, numNodesPerElement, dim)
    p_dot_bar = zeros(Tres, numDofPerNode)
    res_bar = zeros(Tres, numDofPerNode, numNodesPerElement)


    return new(_C_eps, C_eps, strongdata, h_k_tilde, h_k_tilde_bar,
               p_dot, press_el, press_dx, work, res, Rp, fp,
               p_jac, res_jac, p_hess, Dx, px_jac, val_dot, Rp_jac, fp_jac,
               fp_bar, Rp_bar, press_el_bar, press_dx_bar, p_dot_bar, res_bar)
  end
end



#------------------------------------------------------------------------------
# Harten's High Order sensor, with elementwise constant viscoscity

mutable struct ShockSensorHHOConst{Tsol, Tres} <: AbstractShockSensor
  _C_eps::Float64  # original value
  C_eps::Float64  # value scaled by alpha

  strongdata::StrongFormData{Tsol, Tres}
  h_k_tilde::Matrix{Tres}
  h_k_tilde_bar::Matrix{Tres}

  p_dot::Vector{Tsol}
  press_el::Matrix{Tsol}
  press_dx::Array{Tres, 3}
  work::Array{Tres, 3}
  res::Array{Tres, 2}
  Rp::Vector{Tres}
  fp::Vector{Tres}
  epsilon::Matrix{Tres}


  p_jac::Array{Tsol, 3}
  res_jac::Array{Tres, 4}
  p_hess::Matrix{Tsol}
  Dx::Array{Tres, 3}
  px_jac::Array{Tres, 5}
  val_dot::Matrix{Tres}
  Rp_jac::Array{Tres, 3}
  fp_jac::Array{Tres, 3}
  epsilon_dot::Array{Tres, 4}
  val2_dot::Matrix{Tres}

  fp_bar::Vector{Tres}
  Rp_bar::Vector{Tres}
  press_el_bar::Matrix{Tres}
  press_dx_bar::Array{Tres, 3}
  p_dot_bar::Vector{Tres}
  res_bar::Matrix{Tres}

  function ShockSensorHHOConst{Tsol, Tres}(mesh::AbstractMesh, sbp::AbstractSBP,
                                       opts) where {Tsol, Tres}

    @unpack mesh numDofPerNode numNodesPerElement dim

    _C_eps = 1/5  # was 1/5
    C_eps = _C_eps
    strongdata = StrongFormData{Tsol, Tres}(mesh, sbp, opts)
    h_k_tilde = Array{Tres}(mesh.dim, mesh.numEl)
    h_k_tilde_bar = Array{Tres}(mesh.dim, mesh.numEl)
    calcAnisoFactors(mesh, sbp, opts, h_k_tilde)

    p_dot = zeros(Tsol, numDofPerNode)
    press_el = zeros(Tsol, 1, numNodesPerElement)
    press_dx = zeros(Tres, 1, numNodesPerElement, dim)
    work = zeros(Tres, 1, numNodesPerElement, dim)
    res = zeros(Tres, numDofPerNode, numNodesPerElement)
    Rp = zeros(Tres, numNodesPerElement)
    fp = zeros(Tres, numNodesPerElement)
    epsilon = zeros(Tres, 1, numNodesPerElement)


    p_jac = zeros(Tsol, 1, numDofPerNode, numNodesPerElement)
    res_jac = zeros(Tres, numDofPerNode, numDofPerNode, numNodesPerElement,
                          numNodesPerElement)
    p_hess = zeros(Tsol, numDofPerNode, numDofPerNode)
    Dx = zeros(numNodesPerElement, numNodesPerElement, dim)
    px_jac = zeros(Tres, 1, numDofPerNode, dim, numNodesPerElement,
                         numNodesPerElement)
    val_dot = zeros(Tres, numDofPerNode, numNodesPerElement)

    Rp_jac = zeros(Tres, numNodesPerElement, numDofPerNode, numNodesPerElement)
    fp_jac = zeros(Tres, numDofPerNode, numNodesPerElement, numNodesPerElement)
    epsilon_dot = zeros(Tres, 1, numDofPerNode, numNodesPerElement,
                              numNodesPerElement)
    val2_dot = zeros(Tres, numDofPerNode, numNodesPerElement)

    fp_bar = zeros(Tres, numNodesPerElement)
    Rp_bar = zeros(Tres, numNodesPerElement)
    press_el_bar = zeros(Tres, 1, numNodesPerElement)
    press_dx_bar = zeros(Tres, 1, numNodesPerElement, dim)
    p_dot_bar = zeros(Tres, numDofPerNode)
    res_bar = zeros(Tres, numDofPerNode, numNodesPerElement)


    return new(_C_eps, C_eps, strongdata, h_k_tilde, h_k_tilde_bar,
               p_dot, press_el, press_dx, work, res, Rp, fp, epsilon,
               p_jac, res_jac, p_hess, Dx, px_jac, val_dot, Rp_jac, fp_jac,
               epsilon_dot, val2_dot,
               fp_bar, Rp_bar, press_el_bar, press_dx_bar, p_dot_bar, res_bar)
  end
end


const HHOSensors{Tsol, Tres} = Union{ShockSensorHHO{Tsol, Tres}, ShockSensorHHOConst{Tsol, Tres}}

function setAlpha(obj::HHOSensors, alpha::Number)

  obj.C_eps = obj._C_eps*alpha
end


function updateMetrics(mesh, sbp, opts, sensor::HHOSensors)

  calcAnisoFactors(mesh, sbp, opts, sensor.h_k_tilde)

  return nothing
end

function initForRevm(sensor::HHOSensors)

  fill!(sensor.h_k_tilde_bar, 0)

  return nothing
end

function finishRevm(mesh::AbstractMesh, sbp::AbstractOperator, eqn::EulerData,
                    opts, sensor::HHOSensors)

  calcAnisoFactors_revm(mesh, sbp, opts, sensor.h_k_tilde, sensor.h_k_tilde_bar)

  return nothing
end


#------------------------------------------------------------------------------
# Creating shock sensors

global const ShockSensorDict = Dict{String, Type{T} where T <: AbstractShockSensor}(
"SensorNone" => ShockSensorNone,
"SensorEverywhere" => ShockSensorEverywhere,
"SensorVelocity" => ShockSensorVelocity,
"SensorPP" => ShockSensorPP,
"SensorHIso" => ShockSensorHIso,
"SensorBO" => ShockSensorBO,
"SensorHHO" => ShockSensorHHO,
"SensorHHOConst" => ShockSensorHHOConst,
)


import PDESolver: createShockSensor
"""
  Constructs and returns the shock sensor specified in the options dictionary.
"""
function createShockSensor(mesh::AbstractMesh, sbp::AbstractOperator,
                        eqn::EulerData{Tsol, Tres}, opts,
                        name=opts["shock_sensor_name"]) where {Tsol, Tres}

  obj = ShockSensorDict[name]{Tsol, Tres}(mesh, sbp, opts)

  return obj
end


