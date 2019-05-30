# functions for doing parallel computations that are precisely the reverse mode
# of the regular parallel communication
# This can be suboptimal because it requires doing 2 parallel communications
# one after another (one for q and another one to return q_bar), rather than
# doing one big communication.
# The idea behind the functions in this file is:
#
# Primal: q_local + q_remote -> res_local
# Reverse: q_local + q_remote + res_bar_local -> q_bar_local + q_bar_remote
#
# The fact that q_bar_remote is computed necessitates the second communication

"""
  This function combines the functionality of finishExchangeData for the
  solution and exchangeData for `q_bar`.  This function does 4 things:

    1. waits to receive the solution from the remote process (typically started
       with [`startSolutionExchange`](@ref)
    2. calls `calc_func` to populate `send_data_bar.q_recv` with `q_bar`
    3. sends `q_bar` to remote process
    4. posts receives for `q_bar` from remote processes

  Users should call [`finishSolutionExchange_rev`](@ref) to complete
  the receives for `q_bar`.

  It is safe to call this function even if `shared_data` has already been
  waited on.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * shared_data: the vector of `SharedFaceData` objects for the solution
   * shared_data_bar: the vector of `SharedFaceData` objects for `q_bar`
   * calc_func: function that takes an element of `shared_data` and populates
                the *receive* buffer of the corresponding element of
                `shared_data_bar`.  See below

  **Keyword Arguments**

   * wait: if true, waits for the `shared_data_bar` sends and receives to finish
           before exiting.  Default false

  `calc_func` should have the signature

  ```
    calc_func(mesh, sbp, eqn, opts, data::SharedFaceData,
              data_bar::SharedFaceData)
  ```

  where `data` is one element of `shared_data` and `data_bar` is the
  corresponding element of `shared_data_bar`.  The receive buffer of
  `data` will have been populated before `calc_func` is called.  `calc_func`
  must populate the *receive* buffer of `data_bar` with `q_bar` values.
  Note that this is opposite the usual usage of `SharedFaceData` objects.
  After `calc_func` returns, the data in `data_bar.q_recv` will be sent
  to the remote MPI process.  The receive buffer will be zeroed out
  before being passed to `calc_func`.


  Currently `shared_data_bar` must have a `parallel_data` setting of PARALLEL_DATA_ELEMENT
"""
function exchangeData_rev(mesh::AbstractMesh, sbp::AbstractOperator,
                      eqn::AbstractSolutionData, opts,
                      shared_data::Vector{SharedFaceData{T}},
                      shared_data_bar::Vector{SharedFaceData{T}},
                      calc_func::Function;
                      wait=false) where T

  npeers = length(shared_data)

  # bail out early if there is no communication to do
  # not sure if the rest of this function runs correctly if npeers == 0
  if npeers == 0
    return nothing
  end

  # this should already have happened.  If it hasn't something else has
  # gone wrong in the solver.  Throw an exception
  assertSendsWaited(shared_data_bar)

  # post the receives first
  for i=1:npeers
    Irecv_rev!(shared_data_bar[i])
  end

  # verify the sends are consistent
  val = assertReceivesConsistent(shared_data)
  assertReceivesConsistent(shared_data_bar)


  for i=1:npeers
    # wait on a primal communication to finish
    if val == 0  # requests not previously waited on
      eqn.params.time.t_wait += @elapsed idx = waitAnyReceive(shared_data)
    else
      idx = i
    end

    data_i = shared_data[idx]
    data_bar_i = shared_data_bar[idx]
    tag = data_bar_i.tag

    if data_i.pdata == PARALLEL_DATA_FACE && val == 0
      # permute the received nodes to be in the elementR orientation
      permuteinterface!(mesh.sbpface, data_i.interfaces, data_i.q_recv)
    end

    # wait on the previous bar send if it hasn't been waited on yet
    # this should have completed long ago
    if !data_bar_i.recv_waited
      waitReceive(data_bar_i)
    end
    fill!(data_bar_i.q_recv, 0)

    # call calc_func to populate bar q_recv
    calc_func(mesh, sbp, eqn, opts, data_i, data_bar_i)

    # post the bar send
    Isend_rev(data_bar_i)
  end

  if wait
    waitAllSends(shared_data)
    waitAllSends(shared_data_bar)
    waitAllReceives(shared_data)
    waitAllReceives(shared_data_bar)
  end

  return nothing
end


"""
  This function completes the receive of the [`shared_data_bar`](@ref)
  started by [`exchangeData_rev`](@ref).

  It is safe to call this function even if the `shared_data_bar` objects
  have been waited on already, however calling this function more than once
  will result in `populate_buffer_rev` being called multipled times, which is
  likely not what the user wants.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * shared_data_bar: vector of [`SharedFaceData`](@ref) objects.  Must be
                      the same vector passed to the `shared_data_bar` argument
                      of `exchangeData_rev`.
   * populate_buffer_rev: function that will be called after every receive.
                          This function is often the reverse mode of
                          the `populate_buffer` argument to `exchangeData`

  `populate_buffer_rev` must have the signature

  ```
    populate_buffer_rev(mesh, sbp, eqn, opts, data_bar::SharedFaceData)
  ```

  And should take the data in the *send* buffer of the `SharedFaceData` object
  and do some computation with it, typically updating `eqn.q_bar`.  The reason
  the send buffer should be used is that this is a reverse mode function, which
  reverses the typical rules of send and receive.

  The parallel_data setting of the `shared_data_bar` objects must be
  `PARALLEL_DATA_ELEMENT`.
"""
function finishExchangeData_rev(mesh, sbp, eqn, opts,
                            shared_data_bar::Vector{SharedFaceData{T}},
                            populate_buffer_rev::Function) where T


  npeers = length(shared_data_bar)
  val = assertSendsConsistent(shared_data_bar)

  for i=1:npeers
    if val == 0  # request have not been waited on previously
      eqn.params.time.t_wait += @elapsed idx = waitAnySend(shared_data_bar)
    else
      idx = i
    end

    data_idx = shared_data_bar[idx]

    # to support PARALLEL_DATA_FACE: add permutation function here

    populate_buffer_rev(mesh, sbp, eqn, opts, data_idx)
  end

  # wait on all receives to complete.  This is not strictly necessary, but 
  # avoids a problem with exchangeData_rev2 where it expects the receives
  # to already have been waited on
  waitAllReceives(shared_data_bar)


  return nothing
end


"""
  This is a thin wrapper around [`finishExchangeData_rev`](@ref) for the
  common case of receiving `q_bar` from other processes and adding it to
  `eqn.q_bar.

  **Inputs**

   * mesh
   * sbp
   * eqn: communication for `eqn.shared_data_bar` completed
   * opts
"""
function finishSolutionBarExchange(mesh::AbstractMesh, sbp::AbstractOperator,
                                   eqn::AbstractSolutionData, opts)

  if mesh.npeers == 0
    return nothing
  end

  pdata = eqn.shared_data_bar[1].pdata  # assume all are same
  if pdata == PARALLEL_DATA_FACE
    error("PARALLEL_DATA_FACE not supported")
  elseif pdata == PARALLEL_DATA_ELEMENT
    populate_buffer_rev = getSendDataElement_rev
  else
    throw(ErrorException("unsupported parallel_type = $(getParallelDataString(pdata))"))
  end

  finishExchangeData_rev(mesh, sbp, eqn, opts, eqn.shared_data_bar, populate_buffer_rev)

  return nothing
end


"""
  Reverse mode of [`getSendDataElement`](@ref).  Takes data from the send
  buffer of `data` and adds it to `eqn.q_bar`
"""
function getSendDataElement_rev(mesh::AbstractMesh, sbp::AbstractOperator,
                         eqn::AbstractSolutionData, opts, data::SharedFaceData)

  # copy data into send buffer
  idx = data.peeridx
  local_els = mesh.local_element_lists[idx]
  send_buff = data.q_send
  for j=1:length(local_els)
    el_j = local_els[j]
    for k=1:size(eqn.q, 2)
      for p=1:size(eqn.q, 1)
        eqn.q_bar[p, k, el_j] += send_buff[p, k, j]
      end
    end
  end

  return nothing
end


