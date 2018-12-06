# parallel communication primatives

include("tag_manager.jl")  # managment of MPI tags
include("parallel_types.jl")  # datatype
include("parallel_rev1.jl")
include("parallel_rev2.jl")


"""
  This function is a thin wrapper around exchangeData().  It is used for the
  common case of sending and receiving the solution variables to other processes.
  It uses eqn.shared_data to do the parallel communication.
  eqn.shared_data *must* be passed into the corresponding finishDataExchange
  call.

  It is safe to call this function in the serial case.

  **Inputs**

   * mesh: an AbstractMesh
   * sbp: an SBP operator
   * eqn: an AbstractSolutionData
   * opts: options dictionary

  **Keyword arguments**

   * wait: wait for sends and receives to finish before exiting
"""
function startSolutionExchange(mesh::AbstractMesh, sbp::AbstractOperator,
                                  eqn::AbstractSolutionData, opts;
                                  wait=false)

  if mesh.npeers == 0
    return nothing
  end

  pdata = eqn.shared_data[1].pdata  # assume all are same
  if pdata == PARALLEL_DATA_FACE
    populate_buffer = getSendDataFace
  elseif pdata == PARALLEL_DATA_ELEMENT
    populate_buffer = getSendDataElement
  else
    throw(ErrorException("unsupported parallel_type = $(getParallelDataString(pdata))"))
  end

  exchangeData(mesh, sbp, eqn, opts, eqn.shared_data, populate_buffer, wait=wait)

  return nothing
end


"""
  This function posts the MPI sends and receives for a vector of SharedFaceData.  It works for both `PARALLEL_DATA_FACE or `PARALLEL_DATA_ELEMENT.  The only
  difference between these two cases is the populate_buffer() function.

  The previous receives using these SharedFaceData objects should have
  completed by the time this function is called.  An exception is throw
  if this is not the case.

  The previous sends are likely to have completed by the time this function
  is called, but they are waited on if not.  This function might not perform
  well if the previous sends have not completed.
  #TODO: fix this using WaitAny

  Inputs:
    mesh: an AbstractMesh
    sbp: an SBPOperator
    eqn: an AbstractSolutionData
    opts: the options dictionary
    populate_buffer: function with the signature:
                     populate_buffer(mesh, sbp, eqn, opts, data::SharedFaceData)
                     that populates data.q_send
  Inputs/Outputs:
    shared_data: vector of SharedFaceData objects representing the parallel
                 communication to be done

  Keyword Arguments:

    wait: wait for the sends and receives to finish before returning.  This
          is a debugging option only.  It will kill parallel performance.
"""
function exchangeData(mesh::AbstractMesh, sbp::AbstractOperator,
                      eqn::AbstractSolutionData, opts,
                      shared_data::Vector{SharedFaceData{T}},
                      populate_buffer::Function;
                      wait=false) where T

  npeers = length(shared_data)

  # bail out early if there is no communication to do
  # not sure if the rest of this function runs correctly if npeers == 0
  if npeers == 0
    return nothing
  end

  # this should already have happened.  If it hasn't something else has
  # gone wrong in the solver.  Throw an exception
  assertReceivesWaited(shared_data)

  # post the receives first
  for i=1:npeers
    Irecv!(shared_data[i])
  end

  # verify the sends are consistent
  assertSendsConsistent(shared_data)

  for i=1:npeers
    # wait for these in order because doing the waitany trick doesn't work
    # these should have completed long ago, so it shouldn't be a performance
    # problem

    # the waitany trick doesn't work because this loop posts new sends, reusing
    # the same array of MPI_Requests.

    # TODO: use 2 arrays for the Requests: old and new, so the WaitAny trick
    #       works

    idx = i
    data_i = shared_data[idx]
    tag = data_i.tag

    # wait on the previous send if it hasn't been waited on yet
    if !data_i.send_waited
      waitSend(data_i)
    end

    # move data to send buffer
    populate_buffer(mesh, sbp, eqn, opts, data_i)

    # post the send
    Isend(data_i)
  end

  if wait
    waitAllSends(shared_data)
    waitAllReceives(shared_data)
  end

  return nothing
end




"""
  This is the counterpart of exchangeData.  This function finishes the
  receives started in exchangeData.

  This function (efficiently) waits for a receive to finish and calls
  a function to do calculations for on that data. If `PARALLEL_DATA_FACE`
  == PARALLEL_DATA_FACE, it also permutes the data in the receive buffers to agree
  with the ordering of elementL.  For `PARALLEL_DATA_ELEMENT,
  users should call SummationByParts.interiorFaceInterpolate to interpolate
  the data to the face while ensuring proper permutation.

  Calling this function repeatedly without calling [`exchangeData`](@ref)
  in between is supported.  In this case, it only waits for the MPI
  communication to finish the first time.

  **Inputs**

   * mesh: an AbstractMesh
   * sbp: an SBPOperator
   * eqn: an AbstractSolutionData
   * opts: the options dictonary
   * calc_func: function that does calculations for a set of shared faces
               described by a single SharedFaceData.  It must have the signature
               calc_func(mesh, sbp, eqn, opts, data::SharedFaceData)

  **Inputs/Outputs**

   * shared_data: vector of SharedFaceData, one for each peer process that
                 needs to be communicated with.  By the time calc_func is
                 called, the SharedFaceData passed to it has its q_recv field
                 populated.  See note above about data permutation.
"""
function finishExchangeData(mesh, sbp, eqn, opts,
                            shared_data::Vector{SharedFaceData{T}},
                            calc_func::Function) where T

  npeers = length(shared_data)

  if mesh.npeers == 0
    return nothing
  end

  val = assertReceivesConsistent(shared_data)
  
  for i=1:npeers
    if val == 0  # request have not been waited on previously
      eqn.params.time.t_wait += @elapsed idx = waitAnyReceive(shared_data)
    else
      idx = i
    end

    data_idx = shared_data[idx]
    if data_idx.pdata == PARALLEL_DATA_FACE && val == 0
      # permute the received nodes to be in the elementR orientation
      permuteinterface!(mesh.sbpface, data_idx.interfaces, data_idx.q_recv)
    end

    calc_func(mesh, sbp, eqn, opts, data_idx)

  end

  return nothing
end


 
@doc """
### Utils.verifyCommunication

  This function checks the data provided by the Status object to verify a 
  communication completed successfully.  The sender's rank and the number of
  elements is checked agains the expected sender and the buffer size

  Inputs:
    data: a SharedFaceData
"""->
function verifyReceiveCommunication(data::SharedFaceData{T}) where T
# verify a communication occured correctly by checking the fields of the 
# Status object
# if the Status came from a send, then peer should be comm_rank ?

  sender = MPI.Get_source(data.recv_status)
  @assert sender == data.peernum

  ndata = MPI.Get_count(data.recv_status, T)
  @assert ndata == length(data.q_recv)

  return nothing
end

"""
  This function populates the send buffer from `eqn.q` for 
  `PARALLEL_DATA_FACE`

  Inputs:
    mesh: a mesh
    sbp: an SBP operator
    eqn: an AbstractSolutionData
    opts: options dictonary

  Inputs/Outputs:
    data: a SharedFaceData.  data.q_send will be overwritten
"""
function getSendDataFace(mesh::AbstractMesh, sbp::AbstractOperator,
                         eqn::AbstractSolutionData, opts, data::SharedFaceData)


  idx = data.peeridx
  bndryfaces = data.bndries_local
  boundaryinterpolate!(mesh.sbpface, bndryfaces, eqn.q, data.q_send)

  return nothing
end


"""
  This function populates the send buffer from eqn.q for 
  `PARALLEL_DATA_ELEMENT`

  Inputs:

    mesh: a mesh
    sbp: an SBP operator
    eqn: an AbstractSolutionData
    opts: options dictonary

  Inputs/Outputs:

    data: a SharedFaceData.  data.q_send will be overwritten
"""
function getSendDataElement(mesh::AbstractMesh, sbp::AbstractOperator,
                         eqn::AbstractSolutionData, opts, data::SharedFaceData)

  # copy data into send buffer
  idx = data.peeridx
  local_els = mesh.local_element_lists[idx]
  send_buff = data.q_send
  for j=1:length(local_els)
    el_j = local_els[j]
    for k=1:size(eqn.q, 2)
      for p=1:size(eqn.q, 1)
        send_buff[p, k, j] = eqn.q[p, k, el_j]
      end
    end
  end

  return nothing
end




@doc """
### Utils.mpi_master

  This macro introduces an if statement that causes the expression to be 
  executed only if the variable myrank is equal to zero.  myrank must exist
  in the scope of the caller

"""->
macro mpi_master(ex)
  return quote
#    println("myrank = ", esc(myrank))
    if $(esc(:(myrank == 0)))
      $(esc(ex))
    end
  end
end

@doc """
### Utils.time_all 

  This macro returns the value produced by the expression as well as 
  the execution time, the GC time, and the amount of memory allocated
"""->
macro time_all(ex)
  quote
    local stats = Base.gc_num()
    local elapsedtime = time_ns()
    local val = $(esc(ex))
    elapsedtime = time_ns() - elapsedtime
    local diff = Base.GC_Diff(Base.gc_num(), stats)
    (val, elapsedtime/1e9, diff.total_time/1e9, diff.allocd)
  end
end

function print_time_all(f, t_elapsed, t_gc, alloc)
    println(f, t_elapsed, " seconds, ", t_gc, " GC seconds, ", alloc, " bytes allocated")
end

#------------------------------------------------------------------------------
# Debugging functions

"""
  Verify the contents of the receive buffers are current by doing a parallel
  communication and comparing the new contents of the buffers with the old
  ones.

  Throws an error if contents are not consistent.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * f: IO, defaults to STDOUT
"""
function checkBufferConsistency(mesh, sbp, eqn::AbstractSolutionData{Tsol}, opts, f::IO=STDOUT) where Tsol

  println(f, "\nChecking buffer consistency")

  # make sure communication is finished
  assertReceivesWaited(eqn.shared_data)

  # copy buffers
  nbufs = length(eqn.shared_data)
  old_bufs = Array{Array{Tsol, 3}}(nbufs)
  for i=1:nbufs
    old_bufs[i] = copy(eqn.shared_data[i].q_recv)
    println(f, "initially, norm of buffer ", i, " = ", vecnorm(old_bufs[i]))
  end

  # do the parallel communications
  startSolutionExchange(mesh, sbp, eqn, opts)

  # finish the communication
  finishExchangeData(mesh, sbp, eqn, opts, eqn.shared_data, tmpCalcFunc)

  # compare buffers
  for i=1:nbufs
    buf1 = old_bufs[i]
    buf2 = eqn.shared_data[i].q_recv
    nrm_i = vecnorm(buf1 - buf2)
    println(f, "buffer ", i, " diffnorm = ", nrm_i)
    println(f, "buffer ", i, " diffnorm real real = ", vecnorm(real(buf1) - real(buf2)))
    println(f, "buffer ", i, " diffnorm real complex = ", vecnorm(real(buf1) - buf2))
    println(f, "buffer ", i, " diffnorm complex real = ", vecnorm(buf1 - real(buf2)))
#    @assert abs(nrm_i) < 1e-13
  end

  return nothing
end

# used for debugging
function tmpCalcFunc(mesh, sbp, eqn, opts, data::SharedFaceData)

  return nothing
end
