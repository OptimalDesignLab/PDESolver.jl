# declare some datatypes used for parallel communication

"""
  Default MPI tag manager
"""
global const TagManager = MPITagManager()

# allow sending of arbitrary face/element-based data
global const PARALLEL_DATA_FACE = 1001
global const PARALLEL_DATA_ELEMENT = 1002

"""
  Given the string describing what data is shared in parallel, returns the
  enum value.  An exception is thrown if an unrecognized value is supplied.

  **Inputs**

   * pdata: a string describing what data is shared in parallel

  **Outputs**

   * an Int enum
"""
function getParallelDataEnum(pdata::String)
  if pdata == "face"
    return PARALLEL_DATA_FACE
  elseif pdata == "element"
    return PARALLEL_DATA_ELEMENT
  else
    error("unrecognized parallel data string: $pdata")
  end
end

"""
  Inverse of [`getParallelDataEnum`](@ref), takes the enum and returns the
  string.

  **Inputs**

   * pdata: Int enum

  **Outputs**

   * a String describing the parallel data
"""
function getParallelDataEnum(pdata::Int)
  if pdata == PARALLEL_DATA_FACE
    return "face"
  elseif pdata == PARALLEL_DATA_ELEMENT
    return "element"
  else
    error("unrecognized parallel data enum: $pdata")
  end
end


"""
  This type holds all the data necessary to perform MPI communication with
  a given peer process that shared mesh edges (2D) or faces (3D) with the
  current process.

  **Fields**

   * peernum: the MPI rank of the peer process
   * peeridx: the index of this peer in mesh.peer_parts
   * myrank: MPI rank of the current process
   * comm: MPI communicator used to define the above
   * pdata: enum describing if face or element data is sent in parallel

   * q_send: the send buffer, a 3D array of n x m x d.  While these dimensions
            are arbitrary, there are two commonly used case.  If
            `PARALLEL_DATA_FACE, then m is mesh.numNodesPerFace and
            d is the number of faces shared with peernum.
            If `PARALLEL_DATA_ELEMENT, then 
            m = mesh.numNodesPerElement and d is the number of elements that
            share faces with peernum.
   * q_recv: the receive buffer.  Similar to q_send, except the size needs to
            to be the number of entities on the *remote* process.

   * send_waited: has someone called MPI.Wait() on send_req yet?  Some MPI
                 implementations complain if Wait() is called on a Request
                 more than once, so use this field to avoid doing so.
   * send_req: the MPI.Request object for the Send/Isend/whatever other type of
              Send
   * send_status: the MPI.Status object returned by calling Wait() on send_req

   * recv_waited: like send_waited, but for the receive
   * recv_req: like send_req, but for the receive
   * recv_status: like send_status, but for the receive

   * bndries_local: Vector of Boundaries describing the faces from the local
                   side of the interface
   * bndries_remote: Vector of Boundaries describing the facaes from the remote
                    side (see the documentation for PdePumiInterface before
                    using this field)
   * interfaces: Vector of Interfaces describing the faces from both sides (see
                the documentation for PdePumiInterfaces, particularly the
                mesh.shared_interfaces field, before using this field

"""
mutable struct SharedFaceData{T} <: AbstractSharedFaceData{T}
  peernum::Int
  peeridx::Int
  myrank::Int
  comm::MPI.Comm
  pdata::Int
  tag::Cint
  q_send::Array{T, 3}  # send buffer
  q_recv::Array{T, 3}  # receive buffer
  send_waited::Bool
  send_req::MPI.Request
  send_status::MPI.Status
  recv_waited::Bool
  recv_req::MPI.Request
  recv_status::MPI.Status
  # keep these fields or not?  They are not strictly required, but are good
  # examples
  bndries_local::Array{Boundary, 1}
  bndries_remote::Array{Boundary, 1}
  interfaces::Array{Interface, 1}

  # default inner constructor (for now)
end

#------------------------------------------------------------------------------
# constructors

"""
  Outer constructor for SharedFaceData.

  **Inputs**

   * mesh: a mesh object
   * peeridx: the index of a peer in mesh.peer_parts
   * pdata: either integer enum or string describing what data is shared in
            parallel
   * tag: MPI tag to use for all communications
   * q_send: the send buffer
   * q_recv: the receive buffer

"""
function SharedFaceData(mesh::AbstractMesh, peeridx::Integer, pdata::Integer,
                     tag::Integer,
                     q_send::Array{T, 3}, q_recv::Array{T, 3}) where T
# create a SharedFaceData for a given peer

  peernum = mesh.peer_parts[peeridx]
  myrank = mesh.myrank
  comm = mesh.comm

  send_waited = true  # don't wait on the request
  send_req = MPI.REQUEST_NULL
  send_status= MPI.Wait!(send_req)

  recv_waited = true
  recv_req = MPI.REQUEST_NULL
  recv_status = MPI.Wait!(recv_req)

  bndries_local = mesh.bndries_local[peeridx]
  bndries_remote = mesh.bndries_remote[peeridx]
  interfaces = mesh.shared_interfaces[peeridx]

  # it would be nice to have a finalizer to free the MPI tag, but finalizers
  # don't run in a defined order, so its possible the MPI manager could be
  # destroyed before the finalizer for the SharedFaceData is run.

  return SharedFaceData{T}(peernum, peeridx, myrank, comm, pdata, tag,
                           q_send, q_recv,
                           send_waited, send_req, send_status,
                           recv_waited, recv_req, recv_status,
                           bndries_local, bndries_remote, interfaces)

end


function SharedFaceData(mesh::AbstractMesh, peeridx::Integer, pdata::String,
                     tag::Integer,
                     q_send::Array{T, 3}, q_recv::Array{T, 3}) where T

  pdata_enum = getParallelDataEnum(pdata)
  return SharedFaceData(mesh, peeridx, pdata_enum, tag, q_send, q_recv)
end

import Base.copy
"""
  Copy function for SharedFaceData.  Note that this does *not* retain the
  send_req/send_status (and similarly for the recceive) state
  of the original object.  Instead, they are initialized the same as the
  constructor.

  This function may only be called after receiving is complete,
  otherwise an exception is thrown.
"""
function copy(data::SharedFaceData{T}) where T


  @assert data.recv_waited 

  peernum = data.peernum
  peeridx = data.peeridx
  myrank = data.myrank
  comm = data.comm
  pdata = data.pdata

  q_send = copy(data.q_send)
  q_recv = copy(data.q_recv)

  # don't copy the state of the MPI operations
  send_waited = true  # don't wait on the request
  send_req = MPI.REQUEST_NULL
  send_status= MPI.Wait!(send_req)

  recv_waited = true
  recv_req = MPI.REQUEST_NULL
  recv_status = MPI.Wait!(recv_req)

  bndries_local = data.bndries_local
  bndries_remote = data.bndries_remote
  interfaces = data.interfaces

  return SharedFaceData{T}(peernum, peeridx, myrank, comm, pdata,
                           q_send, q_recv,
                           send_waited, send_req, send_status,
                           recv_waited, recv_req, recv_status,
                           bndries_local, bndries_remote, interfaces)
end

import Base.copy!
"""
  In-place copy for SharedFaceData.  This copies the buffers, but does not
  retain the state of the Request and Status fields.  Instead they are
  initialized the same as the constructor.

  This function may only be called after receiving is complete,
  otherwise an exception is thrown.
"""
function copy!(dest::SharedFaceData{T}, src::SharedFaceData{T}) where T

  # if a communication is in progress, copying would leave the buffers in
  # an undefind state
  @assert src.recv_waited

  # we assume all the bookkeeping fields don't need to be copy
  copy!(dest.q_send, src.q_send)
  dest.send_waited = true
  dest.send_req = MPI.REQUEST_NULL
  dest.send_status = MPI.Wait!(dest.send_req)
  
  copy!(dest.q_recv, src.q_recv)
  dest.recv_waited = true
  dest.recv_req = MPI.REQUEST_NULL
  dest.recv_status = MPI.Wait!(dest.recv_req)

  return nothing
end

"""
  This function returns a vector of SharedFaceData objects, one for each
  peer processes the current process shares mesh edge (2D) or face (3D) with.
  This function is intended to be used by the AbstractSolutionData constructors,
  although it can be used to create additional vectors of SharedFaceData
  objects.

  if `pdata` == "face", then the send and receive buffers
  are numDofPerNode x numNodesPerFace x number of shared faces.

  if `pdata` == "element", the send and receive buffers are
    numDofPerNode x numNodesPerElement x number of elements that share the
    faces.  Note that the number of elements that share the faces can be
    different for the send and receive buffers.

  **Inputs**

   * Tsol: element type of the arrays
   * mesh: an AbstractMesh object
   * sbp: an SBP operator
   * opts: the options dictonary
   * pdata: string describing what data is shared in parallel

  **Keyword Arguments**

   * tag: integer specifying the MPI tag that will be used for all
          communications by the returned objects.  If negative (default),
          a new unused tag will be determined.  If the user supplies a
          non-negative tag value, it must not be in use according to the
          `TagManager`.
    
  **Outputs**

   * data_vec: Vector{SharedFaceData}.  data_vec[i] corresponds to 
               mesh.peer_parts[i]
"""
function getSharedFaceData(::Type{Tsol}, mesh::AbstractMesh, sbp::AbstractSBP, opts, pdata::String; tag::Integer=-1) where Tsol
# return the vector of SharedFaceData used by the equation object constructor

  @assert mesh.isDG

  pdata_enum = getParallelDataEnum(pdata)
  data_vec = Array{SharedFaceData{Tsol}}(mesh.npeers)

  if tag < 0
    tag = getNextTag(TagManager)
  else
    markTagUsed(TagManager, tag)
  end

  if pdata_enum == PARALLEL_DATA_FACE
    dim2 =  mesh.numNodesPerFace
    dim3_send = mesh.peer_face_counts
    dim3_recv = mesh.peer_face_counts
  elseif pdata_enum == PARALLEL_DATA_ELEMENT
    dim2 = mesh.numNodesPerElement
    dim3_send = mesh.local_element_counts
    dim3_recv = mesh.remote_element_counts
  else
    throw(ErrorException("Unsupported parallel type requested: $pdata"))
  end

  for i=1:mesh.npeers
    qsend = Array{Tsol}(mesh.numDofPerNode, dim2,dim3_send[i])
    qrecv = Array{Tsol}(mesh.numDofPerNode, dim2, dim3_recv[i])
    data_vec[i] = SharedFaceData(mesh, i, pdata_enum, tag, qsend, qrecv)
  end

  return data_vec
end

#------------------------------------------------------------------------------
# API for SharedFaceData
"""
  This function verifies all the receives have been waited on for the 
  supplied SharedFaceData objects
"""
function assertReceivesWaited(shared_data::Vector{SharedFaceData{T}}) where T

  for i=1:length(shared_data)
    data_i = shared_data[i]
    if !shared_data[i].recv_waited
      throw(ErrorException("Process $(data_i.myrank) has not yet completed receive from proccess $(data_i.peernum)"))
    end
  end

  return nothing
end

"""
  Verify either all or none of the sends have been waited on.  Throw an
  exception otherwise.

  Inputs:

    shared_data: Vector of SharedFaceData objects

  Output:

    val: number of receives that have been waited on
"""
function assertSendsConsistent(shared_data::Vector{SharedFaceData{T}}) where T

  nvals = length(shared_data)
  val = 0
  for i=1:nvals
    if shared_data[i].send_waited
      val += 1
    end
  end

  if val != nvals && val != 0  # either all have been waited on or none
    throw(ErrorException("send requests partially waited on: $val of $nvals"))
  end

  return val
end

"""
  Like assertSendsConsistent, but for the receives
"""
function assertReceivesConsistent(shared_data::Vector{SharedFaceData{T}}) where T

  nvals = length(shared_data)
  val = 0
  for i=1:nvals
    if shared_data[i].recv_waited
      val += 1
    end
  end

  if val != nvals && val != 0  # either all have been waited on or none
    throw(ErrorException("receive requests partially waited on: $val of $nvals"))
  end

  return val
end



"""
  This function is like MPI.Waitall, operating on the sends of a vector of 
  SharedFaceData objects
"""
function waitAllSends(shared_data::Vector{SharedFaceData{T}}) where T

  # MPI.Requests are (currently) mutable, and therefore have reference semantics
  npeers = length(shared_data)
  send_reqs = Array{MPI.Request}(npeers)
  for i=1:npeers
    send_reqs[i] = shared_data[i].send_req
  end

  # this might be better than calling Wait() repeatedly (especially if MPICH
  # busywaits)
  stats = MPI.Waitall!(send_reqs)

  for i=1:npeers
    data_i = shared_data[i]
    data_i.send_req = send_reqs[i]  # just in case Request has value semantics
    data_i.send_status = stats[i]
    data_i.send_waited = true
  end

  return nothing
end

"""
  This function is like MPI.Waitall, operating on the recvs of a vector of 
  SharedFaceData objects
"""
function waitAllReceives(shared_data::Vector{SharedFaceData{T}}) where T

  # MPI.Requests are (currently) mutable, and therefore have reference semantics
  npeers = length(shared_data)
  recv_reqs = Array{MPI.Request}(npeers)
  for i=1:npeers
    recv_reqs[i] = shared_data[i].recv_req
  end

  # this might be better than calling Wait() repeatedly (especially if MPICH
  # busywaits)
  stats = MPI.Waitall!(recv_reqs)

  for i=1:npeers
    data_i = shared_data[i]
    data_i.recv_req = recv_reqs[i]  # just in case Request has value semantics
    data_i.recv_status = stats[i]
    data_i.recv_waited = true
  end

  return nothing
end


# static buffer of Request objects
global const _waitany_reqs = Array{MPI.Request}(50)
"""
  Like MPI.WaitAny, but operates on the receives of  a vector of SharedFaceData.
  Only the index of the Request that was waited on is returned, 
  the Status and recv_waited fields of hte SharedFaceData are updated internally
"""
function waitAnyReceive(shared_data::Vector{SharedFaceData{T}}) where T

  npeers = length(shared_data)
  if length(_waitany_reqs) != npeers
    resize!(_waitany_reqs, npeers)
  end

  for i=1:npeers
    _waitany_reqs[i] = shared_data[i].recv_req
  end

  idx, stat = MPI.Waitany!(_waitany_reqs)
 
  shared_data[idx].recv_req = MPI.REQUEST_NULL  # make sure this request is not
                                                # used again
  shared_data[idx].recv_status = stat
  shared_data[idx].recv_waited = true

  return idx
end
