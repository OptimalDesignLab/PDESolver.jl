# types needed for shock capturing

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
  e0::Float64

  # storage
  up::Vector{Tsol}
  up_tilde::Vector{Tsol}
  up1_tilde::Vector{Tsol}
  
  num_dot::Vector{Tres}
  den_dot::Vector{Tres}
  ee_dot::Vector{Tres}
  lambda_max_dot::Matrix{Tres}

  function ShockSensorPP{Tsol, Tres}(mesh::AbstractMesh, sbp::AbstractOperator, opts) where {Tsol, Tres}

    Vp = VandermondeData(sbp, sbp.degree)
    Vp1 = VandermondeData(sbp, sbp.degree-1)  #TODO: unneded?

    # constants from Barter's thesis
    s0 = -(4 + 4.25*log10(sbp.degree))  # was -(4 + 4.25*log10(sbp.degree))
    kappa = 0.5
    e0 = 1  # this is a bit weird, because PP says it should be O(h/p)
    
    up = zeros(Tsol, sbp.numnodes)
    up_tilde = zeros(Tsol, sbp.numnodes)
    up1_tilde = zeros(Tsol, sbp.numnodes)

    num_dot = zeros(Tres, sbp.numnodes)
    den_dot = zeros(Tres, sbp.numnodes)
    ee_dot = zeros(Tres, mesh.numNodesPerElement)
    lambda_max_dot = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)

    return new(Vp, Vp1, s0, kappa, e0, up, up_tilde, up1_tilde,
               num_dot, den_dot, ee_dot, lambda_max_dot)
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
   * entropy_vars: an [`AbstractVariables`](@ref) specifying which entropy
                   variables to use
   * `flux`: an [`AbstractLDGFlux`](@ref)
   * diffusion: an [`AbstractDiffusion`](@ref).
"""
mutable struct LDGShockCapturing{Tsol, Tres} <: AbstractFaceShockCapturing
  # Note: the variable names are from Chen and Shu's "Entropy Stable High Order
  #       DG Methods" paper
  w_el::Array{Tsol, 3}  # entropy variables for all elements in elnums_all
  q_j::Array{Tres, 4}  # auxiliary equation solution for elements in
                       # elnums_all (this gets used for both theta and q)
  entropy_vars::AbstractVariables  # convert to entropy variables
  flux::AbstractLDGFlux
  diffusion::AbstractDiffusion
  function LDGShockCapturing{Tsol, Tres}() where {Tsol, Tres}
    # we don't know the right sizes yet, so just make them zero size
    w_el = Array{Tsol}(0, 0, 0)
    q_j = Array{Tsol}(0, 0, 0, 0)

    # default values
    entropy_vars = IRVariables()
    flux = LDG_ESFlux()
    diffusion = ShockDiffusion{Tres}()

    return new(w_el, q_j, entropy_vars, flux, diffusion)
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

#------------------------------------------------------------------------------
# Shock Viscoscity Model
"""
  Diagonal viscoscity (constant for each element), used for shock capturing
"""
mutable struct ShockDiffusion{T} <: AbstractDiffusion
  ee::Vector{T}
  ee_dot::Array{T, 3}  # derivative of ee wrt solution
  is_nonlinear::BitArray{1}  # true if ee depends on the solution for each
                             # element
  function ShockDiffusion{T}() where {T}
    ee = Vector{T}(0)
    ee_dot = Array{T}(0, 0, 0)
    is_nonlinear = BitArray(0)
    return new(ee, ee_dot, is_nonlinear)
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
  Function to set the elementwise diffusion coefficient and its derivative
"""
function setDiffusionArray_diff(obj::ShockDiffusion, vals::AbstractArray,
                                vals_dot::AbstractArray{T, 3},
                                is_nonlinear::BitArray{1}) where {T}

  obj.ee =vals
  obj.ee_dot = vals_dot
  obj.is_nonlinear = is_nonlinear
end



#------------------------------------------------------------------------------

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
  * entropy_vars: an [`AbstractVariables`](@ref) specifying which entropy
                   variables to use.
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
  entropy_vars::AbstractVariables
  diffusion::AbstractDiffusion
  penalty::AbstractDiffusionPenalty
  alpha::Array{Float64, 2}  # 2 x numInterfaces
  alpha_b::Array{Float64, 1}  # numBoundaryFaces
  alpha_parallel::Array{Array{Float64, 2}, 1}  # one array for each peer

  #------------------
  # temporary arrays
  
  # getFaceVariables
  w_faceL::Matrix{Tsol}
  w_faceR::Matrix{Tsol}
  grad_faceL::Matrix{Tres}
  grad_faceR::Matrix{Tres}

  # getFaceVariables_diff
  wL_dot::Array{Tsol, 3}
  wR_dot::Array{Tsol, 3}
  Dx::Array{Float64, 3}
  t1::Array{Tres, 3}
  t1_dot::Array{Tres, 5}
  t2L_dot::Array{Tres, 5}
  t2R_dot::Array{Tres, 5}
  t3L_dot::Array{Tres, 5}
  t3R_dot::Array{Tres, 5}

  # applyDgkTranspose
  temp1L::Matrix{Tres}
  temp1R::Matrix{Tres}
  temp2L::Array{Tres, 3}
  temp2R::Array{Tres, 3}
  temp3L::Array{Tres, 3}
  temp3R::Array{Tres, 3}
  work::Array{Tres, 3}

  # applyDgkTranspose_diff
  t3L_dotL::Array{Tres, 5}
  t3L_dotR::Array{Tres, 5}
  t3R_dotL::Array{Tres, 5}
  t3R_dotR::Array{Tres, 5}

  t4L_dotL::Array{Tres, 5}
  t4L_dotR::Array{Tres, 5}
  t4R_dotL::Array{Tres, 5}
  t4R_dotR::Array{Tres, 5}

  t5L_dotL::Array{Tres, 5}
  t5L_dotR::Array{Tres, 5}
  t5R_dotL::Array{Tres, 5}
  t5R_dotR::Array{Tres, 5}

  # computeBoundaryTerm_diff
  w_dot::Array{Tsol, 3}
  t2_dot::Array{Tres, 5}
  t3_dot::Array{Tres, 5}
  t4_dot::Array{Tres, 4}

  alpha_comm::AlphaComm # MPI communications for alpha

  function SBPParabolicSC{Tsol, Tres}(mesh::AbstractMesh, sbp::AbstractOperator, eqn, opts) where {Tsol, Tres}

    @unpack mesh numDofPerNode numNodesPerFace dim numNodesPerElement

    # we don't know the right sizes yet, so make them zero size
    w_el = Array{Tsol}(0, 0, 0)
    grad_w = Array{Tres}(0, 0, 0, 0)

    # default values
#    entropy_vars = IRVariables()
    entropy_vars = ConservativeVariables()
    diffusion = ShockDiffusion{Tres}()
    penalty = getDiffusionPenalty(mesh, sbp, eqn, opts)
    alpha = zeros(Float64, 0, 0)
    alpha_b = zeros(Float64, 0)
    alpha_parallel = Array{Array{Float64, 2}}(0)


    # temporary arrays

    # getFaceVariables
    w_faceL = zeros(Tsol, numDofPerNode, numNodesPerFace)
    w_faceR = zeros(Tsol, numDofPerNode, numNodesPerFace)
    grad_faceL = zeros(Tres, numDofPerNode, numNodesPerFace)
    grad_faceR = zeros(Tres, numDofPerNode, numNodesPerFace)

    # getFaceVariables_diff
    wL_dot = zeros(Tsol, numDofPerNode, numDofPerNode, numNodesPerElement)
    wR_dot = zeros(Tsol, numDofPerNode, numDofPerNode, numNodesPerElement)
    Dx = zeros(numNodesPerElement, numNodesPerElement, dim)
    t1 = zeros(Tres, numDofPerNode, numNodesPerElement, dim)
    t1_dot = zeros(Tres, numDofPerNode, numDofPerNode, dim, numNodesPerElement,
                         numNodesPerElement)
    t2L_dot = zeros(Tres, numDofPerNode, numDofPerNode, dim, numNodesPerElement,
                         numNodesPerElement)
    t2R_dot = zeros(Tres, numDofPerNode, numDofPerNode, dim, numNodesPerElement,
                         numNodesPerElement)
    t3L_dot = zeros(Tres, numDofPerNode, numDofPerNode, dim, numNodesPerFace,
                          numNodesPerElement)
    t3R_dot = zeros(Tres, numDofPerNode, numDofPerNode, dim, numNodesPerFace,
                          numNodesPerElement)


    # applyDgk transpose
    temp1L = zeros(Tres, numDofPerNode, numNodesPerFace)
    temp1R = zeros(Tres, numDofPerNode, numNodesPerFace)
    temp2L = zeros(Tres, numDofPerNode, numNodesPerElement, mesh.dim)
    temp2R = zeros(Tres, numDofPerNode, numNodesPerElement, mesh.dim)
    temp3L = zeros(Tres, numDofPerNode, numNodesPerElement, mesh.dim)
    temp3R = zeros(Tres, numDofPerNode, numNodesPerElement, mesh.dim)
    work = zeros(Tres, numDofPerNode, numNodesPerElement, mesh.dim)

    # applyDgk_transpose_diff
    t3L_dotL = zeros(Tres, numDofPerNode, numDofPerNode, dim, numNodesPerFace,
                           numNodesPerElement)
    t3L_dotR = zeros(Tres, numDofPerNode, numDofPerNode, dim, numNodesPerFace,
                           numNodesPerElement)
    t3R_dotL = zeros(Tres, numDofPerNode, numDofPerNode, dim, numNodesPerFace,
                           numNodesPerElement)
    t3R_dotR = zeros(Tres, numDofPerNode, numDofPerNode, dim, numNodesPerFace,
                           numNodesPerElement)

    t4L_dotL = zeros(Tres, numDofPerNode, numDofPerNode, dim, numNodesPerElement,
                           numNodesPerElement)
    t4L_dotR = zeros(Tres, numDofPerNode, numDofPerNode, dim, numNodesPerElement,
                           numNodesPerElement)
    t4R_dotL = zeros(Tres, numDofPerNode, numDofPerNode, dim, numNodesPerElement,
                           numNodesPerElement)
    t4R_dotR = zeros(Tres, numDofPerNode, numDofPerNode, dim, numNodesPerElement,
                           numNodesPerElement)

    t5L_dotL = zeros(Tres, numDofPerNode, numDofPerNode, dim, numNodesPerElement,
                           numNodesPerElement)
    t5L_dotR = zeros(Tres, numDofPerNode, numDofPerNode, dim, numNodesPerElement,
                           numNodesPerElement)
    t5R_dotL = zeros(Tres, numDofPerNode, numDofPerNode, dim, numNodesPerElement,
                           numNodesPerElement)
    t5R_dotR = zeros(Tres, numDofPerNode, numDofPerNode, dim, numNodesPerElement,
                           numNodesPerElement)


    # computeBoundaryTerm_diff
    w_dot = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode,
                        mesh.numNodesPerElement)
    #t1 = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
    #t1_dot = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode, mesh.dim,
    #                     mesh.numNodesPerElement, mesh.numNodesPerElement)
    t2_dot = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode, mesh.dim,
                         mesh.numNodesPerElement, mesh.numNodesPerElement)

    t3_dot = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode, mesh.dim,
                         mesh.numNodesPerFace, mesh.numNodesPerElement)

    t4_dot = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode,
                         mesh.numNodesPerFace, mesh.numNodesPerElement)


    alpha_comm = AlphaComm(mesh.comm, mesh.peer_parts)

    return new(w_el, grad_w, entropy_vars, diffusion, penalty, alpha,
               alpha_b, alpha_parallel, w_faceL, w_faceR, grad_faceL,
               grad_faceR,
               wL_dot, wR_dot, Dx, t1, t1_dot, t2L_dot, t2R_dot, t3L_dot,
               t3R_dot,
               temp1L, temp1R, temp2L, temp2R, temp3L, temp3R, work,
               t3L_dotL, t3L_dotR, t3R_dotL, t3R_dotR,
               t4L_dotL, t4L_dotR, t4R_dotL, t4R_dotR,
               t5L_dotL, t5L_dotR, t5R_dotL, t5R_dotR,
               w_dot, t2_dot, t3_dot, t4_dot,
               alpha_comm)
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
  assertArraysUnique(obj)
  eqn.shock_capturing = obj

  return nothing
end


 
