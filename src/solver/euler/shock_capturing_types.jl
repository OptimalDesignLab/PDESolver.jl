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
# Sensor for testing only: there is a shock everywhere

"""
  Shock sensor that always says there is a shock and returns a viscoscity 
  of 1
"""
mutable struct ShockSensorEverywhere{Tsol, Tres} <: AbstractShockSensor

  function ShockSensorEverywhere{Tsol, Tres}(mesh::AbstractMesh, sbp::AbstractSBP, opts) where {Tsol, Tres}
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

  function ProjectionShockCapturing{Tsol, Tres}(mesh::AbstractMesh, sbp::AbstractOperator, eqn, opts) where {Tsol, Tres}

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

  function LDGShockCapturing{Tsol, Tres}(mesh::AbstractMesh, sbp::AbstractOperator, eqn, opts) where {Tsol, Tres}
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
                                           eqn, opts) where {Tsol, Tres}
    return new()
  end
end


#------------------------------------------------------------------------------
# BR2 Shock Capturing

"""
  Penalty for BR2
"""
struct BR2Penalty{Tsol, Tres} <: AbstractDiffusionPenalty

  delta_w_n::Matrix{Tsol}
  qL::Array{Tres, 3}
  qR::Array{Tres, 3}
  t1L::Array{Tres, 3}
  t1R::Array{Tres, 3}
  t2L::Array{Tres, 3}
  t2R::Array{Tres, 3}

  function BR2Penalty{Tsol, Tres}(mesh, sbp, opts) where {Tsol, Tres}

    delta_w_n = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
    qL = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
    qR = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
    t1L = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
    t1R = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
    t2L = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace, mesh.dim)
    t2R = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace, mesh.dim)

    return new(delta_w_n, qL, qR, t1L, t1R, t2L, t2R)
  end

end

"""
  Parallel communications for alpha
"""
struct AlphaComm
  tag::Cint
  send_req::Vector{MPI.Request}
  recv_req::Vector{MPI.Request}
  recv_waited::Vector{Bool}
  send_bufs::Vector{Vector{UInt8}}  # make the buffers UInt8 so (hopefully)
  recv_bufs::Vector{Vector{UInt8}}  # the send will be done eagerly
  comm::MPI.Comm
  peer_parts::Vector{Int}  # the list of all processes that might be
                           # communicated with (same order as mesh.peer_parts)
  peer_parts_red::Vector{Cint}  # list of peer parts that will actually be
                                # communicated with
  function AlphaComm(comm::MPI.Comm, peer_parts::Vector{Int})

    tag = getNextTag(TagManager)
    send_req = Vector{MPI.Request}(0)
    recv_req = Vector{MPI.Request}(0)
    recv_waited = Vector{Bool}(0)
    send_bufs = Vector{Vector{UInt8}}(0)
    recv_bufs = Vector{Vector{UInt8}}(0)
    peer_parts_red = Vector{Cint}(0)

    return new(tag, send_req, recv_req, recv_waited, send_bufs, recv_bufs,
               comm, peer_parts, peer_parts_red)
  end
end

function allocateArrays(comm::AlphaComm, shockmesh::ShockedElements)

  # make sure previous communications have finished before reallocating arrays
  MPI.Waitall!(comm.send_req)
  assertReceivesWaited(comm)

  resize!(comm.send_bufs, shockmesh.npeers)
  resize!(comm.recv_bufs, shockmesh.npeers)
  resize!(comm.send_req, shockmesh.npeers)
  resize!(comm.recv_req, shockmesh.npeers)
  resize!(comm.recv_waited, shockmesh.npeers)
  resize!(comm.peer_parts_red, shockmesh.npeers)

  # extract needed MPI ranks
  for i=1:shockmesh.npeers
    comm.peer_parts_red[i] = comm.peer_parts[shockmesh.peer_indices[i]]
    comm.recv_bufs[i] = Vector{UInt8}(shockmesh.numSharedInterfaces[i])
    comm.send_bufs[i] = Vector{UInt8}(shockmesh.numSharedInterfaces[i])
  end

  return nothing
end


import MPI: Waitall!, Waitany!, Isend, Irecv!

"""
  Wrapper around MPI.Waitall! for the receive requests
"""
function Waitall!(comm::AlphaComm)
  
  MPI.Waitall!(comm.recv_req)
end

"""
  Wrapper around MPI.Waitany! for the receive requests.
  Returns only the index, not the Status object.
"""
function Waitany!(comm::AlphaComm)

  idx, stat = MPI.Waitany!(comm.recv_req)
  comm.recv_req[idx] = MPI.REQUEST_NULL
  comm.recv_waited[idx] = true

  return idx
end


import Utils: assertReceivesWaited

"""
  Asserts all receives have been waited on
"""
function assertReceivesWaited(comm::AlphaComm)

  if length(comm.recv_waited) > 0
    for i=1:length(comm.recv_waited)
      @assert comm.recv_waited[i]
    end
  end

  return nothing
end


"""
  Starts communication for a given peer index (note: this is the index in
  the range 1:shockmesh.npeer)

  **Inputs**

   * comm: `AlphaComm` object
   * peeridx: peer index
"""
function Isend(comm::AlphaComm, peeridx::Integer)

  req = MPI.Isend(comm.send_bufs[peeridx], comm.peer_parts_red[peeridx],
                  comm.tag, comm.comm)
  comm.send_req[peeridx] = req

  return nothing
end

"""
  Counterpart to Isend, starts the receive for a given peer index.  Same
  arguments as Isend.
"""
function Irecv!(comm::AlphaComm, peeridx::Integer)

  req = MPI.Irecv!(comm.recv_bufs[peeridx], comm.peer_parts_red[peeridx],
                   comm.tag, comm.comm)
  comm.recv_waited[peeridx] = false
  comm.recv_req[peeridx] = req

  return nothing
end


"""
  Computes any diffusion scheme from the SBPParabolic paper.

  Note that the sizes for the arrays are the *minimum* sizes.  In practice
  they may be larger.  Do not call `size()` on these arrays, use the appropriate
  integer field of `mesh` or `shockmesh`.

  **Fields**

   * w_el: entropy variables at each volume node of each element,
           `numDofPerNode` x `numNodesPerElement` x `shockmesh.numEl`
   * grad_w: stores Lambda * D * q at volume nodes, numDofPerNode x
             numNodesPerElement x dim x shockmesh.numEl.  Should be zero
             for all neighboring elements (see `ShockedElement`) because
             the viscoscity is zero there, but storing values in this array
             significantly simplifies the code
   * convert_entropy: a function that converts from the conservative to
                      entropy variables, must have signature:
                      `convert_entropy(params::ParamType, q_c, q_e)`
                      where `q_c` contains the conservative variables ata
                      a node and `q_e` is overwritten with the entropy
                      variables.  Defaults to `convertToIR_`.
  * diffusion: an [`AbstractDiffusion`](@ref) object that specifies the
               diffusion tensor. Defaults to `ShockDiffusion`
  * penalty: an [`AbstractDiffusionPenalty`](@ref) object that specifies which
             scheme to use.  Defaults to opts["DiffusionPenalty"]
  * alpha: a 2 x `shockmesh.numInterfaces` array containing alpha_gk and
           alpha_gn for each interface.

  Note that `convert_entropy`, `diffusion` and `penalty` are abstract types
  and should only be accessed through a function barrier in performance
  critical functions.  The user is allowed to change these fields at any
  time to customize the behavior of the scheme.

  Note that this function is not fully constructed by the constructor,
  see the `allocateArrays` function.

"""
mutable struct SBPParabolicSC{Tsol, Tres} <: AbstractFaceShockCapturing
  w_el::Array{Tsol, 3}
  grad_w::Array{Tres, 4}
  convert_entropy::Any  # convert to entropy variables
  diffusion::AbstractDiffusion
  penalty::AbstractDiffusionPenalty
  alpha::Array{Float64, 2}  # 2 x numInterfaces
  alpha_parallel::Array{Array{Float64, 2}, 1}  # one array for each peer

  #------------------
  # temporary arrays
  
  # getFaceVariables
  w_faceL::Matrix{Tsol}
  w_faceR::Matrix{Tsol}
  grad_faceL::Matrix{Tres}
  grad_faceR::Matrix{Tres}

  # applyDgkTranspose
  temp1L::Matrix{Tres}
  temp1R::Matrix{Tres}
  temp2L::Array{Tres, 3}
  temp2R::Array{Tres, 3}
  temp3L::Array{Tres, 3}
  temp3R::Array{Tres, 3}
  work::Array{Tres, 3}

  alpha_comm::AlphaComm # MPI communications for alpha

  function SBPParabolicSC{Tsol, Tres}(mesh::AbstractMesh, sbp::AbstractOperator, eqn, opts) where {Tsol, Tres}
    # we don't know the right sizes yet, so make them zero size
    w_el = Array{Tsol}(0, 0, 0)
    grad_w = Array{Tres}(0, 0, 0, 0)

    # default values
    convert_entropy = convertToIR_
    diffusion = ShockDiffusion{Tres}()
    penalty = getDiffusionPenalty(mesh, sbp, eqn, opts)
    alpha = zeros(Float64, 0, 0)
    alpha_parallel = Array{Array{Float64, 2}}(0)


    # temporary arrays
    w_faceL = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
    w_faceR = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
    grad_faceL = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
    grad_faceR = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)

    temp1L = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
    temp1R = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
    temp2L = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
    temp2R = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
    temp3L = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
    temp3R = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
    work = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)

    alpha_comm = AlphaComm(mesh.comm, mesh.peer_parts)

    return new(w_el, grad_w, convert_entropy, diffusion, penalty, alpha,
               alpha_parallel, w_faceL, w_faceR, grad_faceL, grad_faceR,
              temp1L, temp1R, temp2L, temp2R, temp3L, temp3R, work, alpha_comm)
  end
end


"""
  Sets up the shock capturing object for use, using the fully initialized
  `shockmesh`.  Also does some other misc. setup stuff that needs `shockmesh`.
"""
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
  allocateArrays(capture.alpha_comm, shockmesh)
  computeAlpha(capture, mesh, shockmesh)

  return nothing
end



#------------------------------------------------------------------------------
# Creating shock sensors

global const ShockSensorDict = Dict{String, Type{T} where T <: AbstractShockSensor}(
"SensorNone" => ShockSensorNone,
"SensorEverywhere" => ShockSensorEverywhere,
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
"LDG" => LDGShockCapturing,
"SBPParabolic" => SBPParabolicSC,
)


"""
  This function populates the `eqn.shock_capturing` field of the mesh with the
  shock capturing object for the scheme.  The shock capturing scheme is
  determined by opts["shock_capturing_name"]
"""
function getShockCapturing(mesh, sbp, eqn::EulerData{Tsol, Tres}, opts) where {Tsol, Tres}

  name = opts["shock_capturing_name"]
  obj = ShockCapturingDict[name]{Tsol, Tres}(mesh, sbp, eqn, opts)
  eqn.shock_capturing = obj

  return nothing
end


#------------------------------------------------------------------------------
# Diffusion penalties for SBPParabolic

global const DiffusionPenaltyDict = Dict{String, Type{T} where T <: AbstractDiffusionPenalty}(
"BR2" => BR2Penalty,
)

"""
  Returns a fully constructed [`AbstractDiffusionPenalty`](@ref).

  **Inputs**
  
   * mesh
   * sbp
   * eqn
   * opts
   * name: name of penalty to get, defaults to opts["DiffusionPenalty"]

  **Outputs**

   * obj: `AbstractDiffusionPenalty` object
"""
function getDiffusionPenalty(mesh, sbp, eqn::EulerData{Tsol, Tres}, opts,
                             name::String=opts["DiffusionPenalty"]
                            ) where {Tsol, Tres}

  obj = DiffusionPenaltyDict[name]{Tsol, Tres}(mesh, sbp, opts)
  assertArraysUnique(obj)

  return obj
end

  
