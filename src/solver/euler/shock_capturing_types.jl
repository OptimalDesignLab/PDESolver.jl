# types needed for different shock capturing schemes (basically different
# diffusion terms)

#------------------------------------------------------------------------------
# Projection-based shock capturing

mutable struct ProjectionShockCapturing{Tsol, Tres} <: AbstractVolumeShockCapturing
  filt::Matrix{Float64}  # the filter operator

  # storage
  w::Matrix{Tsol}
  t1::Matrix{Tres}
  t2::Matrix{Tres}
  Se_jac::Array{Tres, 4}
  ee_jac::Array{Tres, 4}
  A0inv::Matrix{Tsol}

  sensor::AbstractShockSensor

  function ProjectionShockCapturing{Tsol, Tres}(mesh::AbstractMesh, sbp::AbstractOperator, eqn, opts, sensor::AbstractShockSensor) where {Tsol, Tres}

    numDofPerNode = mesh.numDofPerNode
    filt = zeros(Float64, sbp.numnodes, sbp.numnodes)
    getFilterOperator!(sbp, filt)
    println("max entry in filt = ", maximum(abs.(filt)))

    w = zeros(Tsol, numDofPerNode, sbp.numnodes)
    t1 = zeros(Tres, numDofPerNode, sbp.numnodes)
    t2 = zeros(Tres, numDofPerNode, sbp.numnodes)

    Se_jac = zeros(Tres, mesh.dim, numDofPerNode, mesh.numNodesPerElement,
                          mesh.numNodesPerElement)
    ee_jac = zeros(Tres, mesh.dim, numDofPerNode, mesh.numNodesPerElement,
                          mesh.numNodesPerElement)

    A0inv = zeros(Tsol, numDofPerNode, numDofPerNode)

    return new(filt, w, t1, t2, Se_jac, ee_jac, A0inv, sensor)
  end
end

#------------------------------------------------------------------------------
# LDG shock capturing



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
  sensor::AbstractShockSensor
  function LDGShockCapturing{Tsol, Tres}(mesh::AbstractMesh, sbp,
                      eqn, opts, sensor::AbstractShockSensor) where {Tsol, Tres}
    # we don't know the right sizes yet, so just make them zero size
    w_el = Array{Tsol}(0, 0, 0)
    q_j = Array{Tsol}(0, 0, 0, 0)

    # default values
    entropy_vars = getSCVariables(opts)
    flux = LDG_ESFlux()
    diffusion = ShockDiffusion{Tres}(mesh, sensor)

    return new(w_el, q_j, entropy_vars, flux, diffusion, sensor)
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

"""
  Shock capturing type that errors out.  Used when shock capturing is not
  supposed to be added.
"""
mutable struct ErrorShockCapturing{Tsol, Tres} <: AbstractShockCapturing

  sensor::AbstractShockSensor
  diffusion::AbstractDiffusion
  function ErrorShockCapturing{Tsol, Tres}(mesh::AbstractMesh, sbp::AbstractSBP,
                                           eqn, opts,
                                           sensor::AbstractShockSensor) where {Tsol, Tres}

    diffusion = ShockDiffusion{Tres}(mesh, sensor)
    return new(sensor, diffusion)
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
  sensor::AbstractShockSensor
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

  function SBPParabolicSC{Tsol, Tres}(mesh::AbstractMesh, sbp::AbstractOperator,
                                      eqn, opts, sensor::AbstractShockSensor) where {Tsol, Tres}

    @unpack mesh numDofPerNode numNodesPerFace dim numNodesPerElement

    # we don't know the right sizes yet, so make them zero size
    w_el = Array{Tsol}(0, 0, 0)
    grad_w = Array{Tres}(0, 0, 0, 0)

    # default values
    entropy_vars = getSCVariables(opts)
    diffusion = ShockDiffusion{Tres}(mesh, sensor)
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

    return new(w_el, grad_w, entropy_vars, diffusion, penalty, sensor, alpha,
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

  allocateArrays(capture.alpha_comm, shockmesh)
  computeAlpha(capture, mesh, shockmesh)

  return nothing
end




#------------------------------------------------------------------------------
# SBPParabolicReducedSC

"""
  Computes a reduced form of the scheme from the SBPParabolic paper.

  Note that the sizes for the arrays are the *minimum* sizes.  In practice
  they may be larger.  Do not call `size()` on these arrays, use the appropriate
  integer field of `mesh` or `shockmesh`.

  **Fields**

   * w_el: entropy variables at each volume node of each element,
           `numDofPerNode` x `numNodesPerElement` x `shockmesh.numEl`
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
mutable struct SBPParabolicReducedSC{Tsol, Tres} <: AbstractFaceShockCapturing
  w_el::Array{Tsol, 3}
  entropy_vars::AbstractVariables
  diffusion::AbstractDiffusion
  penalty::AbstractDiffusionPenalty
  sensor::AbstractShockSensor
  sensor_const::ShockSensorHApprox{Tsol, Tres}  # shock sensor used for
                                                   # face terms
  alpha::Array{Float64, 2}  # 2 x numInterfaces
  alpha_b::Array{Float64, 1}  # numBoundaryFaces
  alpha_parallel::Array{Array{Float64, 2}, 1}  # one array for each peer

  #------------------
  # temporary arrays
  
  # getFaceVariables
  w_faceL::Matrix{Tsol}
  w_faceR::Matrix{Tsol}

  # getFaceVariables_diff
  wL_dot::Array{Tsol, 3}
  wR_dot::Array{Tsol, 3}

  # computeBoundaryTerm_diff
  w_dot::Array{Tsol, 3}
  t2_dot::Array{Tres, 5}
  t3_dot::Array{Tres, 5}
  t4_dot::Array{Tres, 4}

  alpha_comm::AlphaComm # MPI communications for alpha

  function SBPParabolicReducedSC{Tsol, Tres}(mesh::AbstractMesh, sbp::AbstractOperator,
                                      eqn, opts, sensor::AbstractShockSensor) where {Tsol, Tres}

    @unpack mesh numDofPerNode numNodesPerFace dim numNodesPerElement

    # we don't know the right sizes yet, so make them zero size
    w_el = Array{Tsol}(0, 0, 0)

    # default values
    entropy_vars = getSCVariables(opts)
    diffusion = ShockDiffusion{Tres}(mesh, sensor)
    penalty = getDiffusionPenalty(mesh, sbp, eqn, opts)
    sensor_const = ShockSensorHApprox{Tsol, Tres}(mesh, sbp, opts)
    alpha = zeros(Float64, 0, 0)
    alpha_b = zeros(Float64, 0)
    alpha_parallel = Array{Array{Float64, 2}}(0)


    # temporary arrays

    # getFaceVariables
    w_faceL = zeros(Tsol, numDofPerNode, numNodesPerFace)
    w_faceR = zeros(Tsol, numDofPerNode, numNodesPerFace)

    # getFaceVariables_diff
    wL_dot = zeros(Tsol, numDofPerNode, numDofPerNode, numNodesPerElement)
    wR_dot = zeros(Tsol, numDofPerNode, numDofPerNode, numNodesPerElement)

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

    return new(w_el, entropy_vars, diffusion, penalty, sensor, sensor_const,
               alpha,
               alpha_b, alpha_parallel, w_faceL, w_faceR,
               wL_dot, wR_dot,
               w_dot, t2_dot, t3_dot, t4_dot,
               alpha_comm)
  end
end


"""
  Sets up the shock capturing object for use, using the fully initialized
  `shockmesh`.  Also does some other misc. setup stuff that needs `shockmesh`.
"""
function allocateArrays(capture::SBPParabolicReducedSC{Tsol, Tres},
                        mesh::AbstractMesh,
                        shockmesh::ShockedElements) where {Tsol, Tres}

  if size(capture.w_el, 3) < shockmesh.numEl
    capture.w_el = Array{Tsol}(mesh.numDofPerNode, mesh.numNodesPerElement,
                               shockmesh.numEl)
  end

  allocateArrays(capture.alpha_comm, shockmesh)
  computeAlpha(capture, mesh, shockmesh)

  return nothing
end

function getSparsityPattern(capture::SBPParabolicReducedSC)

  return INVISCID
end



#------------------------------------------------------------------------------
# VolumeShockCapturing


mutable struct VolumeShockCapturing{Tsol, Tres} <: AbstractVolumeShockCapturing
  entropy_vars::AbstractVariables
  diffusion::AbstractDiffusion
  sensor::AbstractShockSensor

  function VolumeShockCapturing{Tsol, Tres}(mesh::AbstractMesh, 
                                sbp::AbstractOperator, eqn,  opts,
                                sensor::AbstractShockSensor) where {Tsol, Tres}
    entropy_vars = getSCVariables(opts)
    diffusion = ShockDiffusion{Tres}(mesh, sensor)
    return new(entropy_vars, diffusion, sensor)
  end
end


"""
  Get the entropy variables specified by the options dictionary for use in the
  shock capturing dissipation term.
"""
function getSCVariables(opts)

  name = opts["shock_capturing_variables"]
  if name == "IR"
    return IRVariables()
  elseif name == "conservative"
    return Conservative()
  else
    error("unrecognized variables for shock capturing dissipation: $name")
  end
end


#------------------------------------------------------------------------------
# Creating shock capturing

global const ShockCapturingDict = Dict{String, Type{T} where T <: AbstractShockCapturing}(
"ShockCapturingNone" => ErrorShockCapturing,
"ElementProjection" => ProjectionShockCapturing,
"LDG" => LDGShockCapturing,
"SBPParabolic" => SBPParabolicSC,
"SBPParabolicReduced" => SBPParabolicReducedSC,
"Volume" => VolumeShockCapturing,
)


"""
  This function populates the `eqn.shock_capturing` field of the mesh with the
  shock capturing object for the scheme.  The shock capturing scheme is
  determined by opts["shock_capturing_name"]
"""
function getShockCapturing(mesh, sbp, eqn::EulerData{Tsol, Tres}, opts,
                           sensor::AbstractShockSensor) where {Tsol, Tres}

  name = opts["shock_capturing_name"]
  obj = ShockCapturingDict[name]{Tsol, Tres}(mesh, sbp, eqn, opts, sensor)
  assertArraysUnique(obj)
  eqn.shock_capturing = obj

  return nothing
end


 
