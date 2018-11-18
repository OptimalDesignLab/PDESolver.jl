# functions for doing the parallel communication required for reverse mode
# calculations
# this is different than the communication in parallel_rev2.jl because it
# sends both q and res_bar at the beginning, and therefore does not require
# a parallel communication at the end to return q_bar to the process that
# ownes it.  Unfortunately, this makes the code the physics modules must
# implement a little bit different in reverse mode.  The idea is:
#
# Primal: q_local + q_remote -> res_local
# Reverse: q_local + q_remote + res_bar_local + res_bar_remote -> q_bar_local
#
# Notice that only q_bar_local needs to be computed because q_bar_remote will
# be computed by the process that ownes it.  It is able to do this because
# it has access to res_bar_local and res_bar_remote



"""
  Similar to [`startSolutionExchange`](@ref), this function starts
  communciation of the solution `eqn.res_bar` and possibly `eqn.q`.
  See [`finishExchangeData_rev`](@ref) for a way to finish communication
  for reverse-mode calculations. `eqn.shared_data` and `eqn.shared_data_bar`
  are used for the communications.

  **Inputs**

   * mesh: an AbstractMesh
   * sbp: an SBP operator
   * eqn: an AbstractSolutionData
   * opts: options dictionary

  **Keyword arguments**

   * send_q: if true, start communication for the solution `eqn.q`.
             Default true.  Callers should carefully consider if this can
             be set to false, which may improve performance.
   * wait: wait for sends and receives to finish before exiting, default false


  **Implementation Notes**

  This function is slightly sub-optimial in that it posts all the sends and
  receives for the solution and then posts the sends and receives for
  `res_bar`.  It would be slightly better to post all receives first and then
  all sends.
"""
function startSolutionExchange_rev2(mesh::AbstractMesh, sbp::AbstractSBP,
                                  eqn::AbstractSolutionData, opts;
                                  send_q=true, wait=false)

  if mesh.npeers == 0
    return nothing
  end
  
  if send_q
    startSolutionExchange(mesh, sbp, eqn, opts, wait=wait)
  end

  pdata = eqn.shared_data_bar[1].pdata  # assume all are same
  if pdata == PARALLEL_DATA_FACE
    populate_buffer = getSendDataFace_resbar
  elseif pdata == PARALLEL_DATA_ELEMENT
    populate_buffer = getSendDataElement_resbar
  else
    throw(ErrorException("unsupported parallel_type = $(getParallelDataString(pdata))"))
  end

  exchangeData(mesh, sbp, eqn, opts, eqn.shared_data_bar, populate_buffer,
               wait=wait)

  return nothing
end





"""
  Similar to [`finishExchangeData`](@ref), but accepts two vectors of
  `SharedFaceData` objects, one for the solution data and the other for some
  `_bar` quantity.  The `calc_func` must have the signature:

  ```
    calc_func(mesh, sbp, eqn, opts, data::SharedFaceData,
              data_bar::SharedFaceData)
  ```

  **Inputs**

   * mesh: an AbstractMesh
   * sbp: an SBPOperator
   * eqn: an AbstractSolutionData
   * opts: the options dictonary
   * calc_func: function that does calculations for a set of shared faces
               described by a single SharedFaceData.  It must have the signature
               calc_func(mesh, sbp, eqn, opts, data::SharedFaceData)

  **Inputs/Outputs**

   * shared_data: vector of `SharedFaceData` for solution variables
   * shared_data_bar: vector of `SharedFaceData` for `_bar` variables

  In reality, whether or not `shared_data` contains the solution and
  `shared_data_bar` contains `_bar` variables of some kind does not matter.
  This function is simply the extensions of the `finishExchangeData` where
  two vector of `SharedFaceData` are provided.

  **Implementation Notes**

  The implementation uses MPI's Waitany function on `shared_data_bar`.  This
  is because the more likely case for reverse-mode calculations is that the
  solution was sent previously and only `shared_face_bar` needs to be
  waited on.

"""
function finishExchangeData_rev2(mesh, sbp, eqn, opts,
                                shared_data::Vector{SharedFaceData{T}},
                                shared_data_bar::Vector{SharedFaceData{T2}},
                                calc_func::Function) where {T, T2}

  
  npeers = length(shared_data)
  @assert npeers == length(shared_data_bar)
  val = assertReceivesConsistent(shared_data)
  val_bar = assertReceivesConsistent(shared_data_bar)

  println(eqn.params.f, "entered finishExchangeData_rev2")
  println(eqn.params.f, "val = ", val, ", val_bar = ", val_bar)

  for i=1:npeers

    #TODO: there is a slightly better algorithm for waiting on two sets of
    #      Requestions: Waitany on the first, check if the corresponding
    #      entry in the second set has finished.  If not, Waitany on the 
    #      second set.  Continue ping-ponging back and forth unitil a matched
    #      set is ready
    if val_bar == 0  # request have not been waited on previously
      eqn.params.time.t_wait += @elapsed idx = waitAnyReceive(shared_data_bar)
    else
      idx = i
    end

    if val == 0
      waitReceive(shared_data, idx)
    end


    data_bar_idx = shared_data_bar[idx]
    if data_bar_idx.pdata == PARALLEL_DATA_FACE && val_bar == 0
      # permute the received nodes to be in the elementR orientation
      permuteinterface!(mesh.sbpface, data_bar_idx.interfaces, data_bar_idx.q_recv)
    end

    data_idx = shared_data[idx]
    if data_idx.pdata == PARALLEL_DATA_FACE && val == 0
      # permute the received nodes to be in the elementR orientation
      permuteinterface!(mesh.sbpface, data_idx.interfaces, data_idx.q_recv)
    end

    calc_func(mesh, sbp, eqn, opts, data_idx, data_bar_idx)

  end

  return nothing
end



"""
  Like [`getSendDataFace`](@ref), but takes the data from `eqn.res_bar`
  and puts it in the send buffer.
"""
function getSendDataFace_resbar(mesh::AbstractMesh, sbp::AbstractSBP,
                         eqn::AbstractSolutionData, opts, data::SharedFaceData)


  idx = data.peeridx
  bndryfaces = data.bndries_local
  boundaryinterpolate!(mesh.sbpface, bndryfaces, eqn.res_bar, data.q_send)

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
function getSendDataElement_resbar(mesh::AbstractMesh, sbp::AbstractSBP,
                         eqn::AbstractSolutionData, opts, data::SharedFaceData)

  # copy data into send buffer
  idx = data.peeridx
  local_els = mesh.local_element_lists[idx]
  send_buff = data.q_send
  for j=1:length(local_els)
    el_j = local_els[j]
    for k=1:size(eqn.q, 2)
      for p=1:size(eqn.q, 1)
        send_buff[p, k, j] = eqn.res_bar[p, k, el_j]
      end
    end
  end

  return nothing
end


