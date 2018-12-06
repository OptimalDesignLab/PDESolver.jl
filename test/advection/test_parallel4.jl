# run tests in parallel (np=4)

"""
  Test the parallel communication primatives in the Utils module
"""
function test_parallel2_comm()

  fname = "input_vals_parallel2.jl"
  mesh, sbp, eqn, opts = solvePDE(fname)

  myrank = mesh.myrank
  commsize = mesh.commsize

  AdvectionEquationMod.ICDict["ICp1"](mesh, sbp, eqn, opts, eqn.q_vec)
  array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)
  # add offset so we can tell where data came from
  for i=1:length(eqn.q)
    eqn.q[i] += myrank
  end


#=
  function wrap(i, commsize)
    if i > (commsize-1)
      i = 0
    elseif i < 0
      i = commsize - 1
    end

    return i
  end

  function populate_buffer(mesh, sbp, eqn, opts, data_i)

    return nothing
  end
=#
  @testset "----- Testing Parallel Communication -----" begin


    # parallel_data == PARALLEL_DATA_FACE
    populate_buffer = Utils.getSendDataFace
    shared_data = getSharedFaceData(Float64, mesh, sbp, opts, opts["parallel_data"])
    shared_data_bar = getSharedFaceData(Float64, mesh, sbp, opts, opts["parallel_data"])

    exchangeData(mesh, sbp, eqn, opts, shared_data, populate_buffer)
    calc_func = (mesh, sbp, eqn, opts, data) -> return nothing
    finishExchangeData(mesh, sbp, eqn, opts, shared_data, calc_func)

    testSharedDataFace(mesh, sbp, eqn, opts, shared_data)

    # test exchanging element data
    setParallelData(shared_data, PARALLEL_DATA_ELEMENT)
    populate_buffer = Utils.getSendDataElement

    exchangeData(mesh, sbp, eqn, opts, shared_data, populate_buffer)
    calc_func = (mesh, sbp, eqn, opts, data) -> return nothing
    finishExchangeData(mesh, sbp, eqn, opts, shared_data, calc_func)


    testSharedDataElement(mesh, sbp, eqn, opts, shared_data)


    # test finishExchangeData_rev2

    fill!(eqn.q_bar, 0)
    for i=1:mesh.npeers
      rand!(eqn.shared_data_bar[i].q_recv)
    end
    #TODO: fill recv buffers with random values to test zeroing out
    startSolutionExchange(mesh, sbp, eqn, opts)
    exchangeData_rev(mesh, sbp, eqn, opts, eqn.shared_data, eqn.shared_data_bar, calc_func_rev)
    finishSolutionBarExchange(mesh, sbp, eqn, opts)
    testSharedDataElement_rev(mesh, sbp, eqn, opts)


    
    calc_func2 = (mesh, sbp, eqn, opts, data, data2) -> return nothing
    setParallelData(shared_data, PARALLEL_DATA_FACE)
    setParallelData(shared_data_bar, PARALLEL_DATA_FACE)
    copy!(eqn.res_bar, eqn.q)
    for i=1:length(eqn.res_bar)
      eqn.res_bar[i] += 1   # add offset to distinguish q and res_bar
    end

    populate_buffer = Utils.getSendDataFace
    populate_buffer_rev = Utils.getSendDataFace_resbar

    #TODO: repalce with startSolutionExchange_rev2
    exchangeData(mesh, sbp, eqn, opts, shared_data, populate_buffer)
    exchangeData(mesh, sbp, eqn, opts, shared_data_bar, populate_buffer_rev)

    finishExchangeData_rev2(mesh, sbp, eqn, opts, shared_data, shared_data_bar,
                           calc_func2)
    testSharedDataFace(mesh, sbp, eqn, opts, shared_data)
    testSharedDataFace(mesh, sbp, eqn, opts, shared_data_bar, 1)


    # make parallel_data different
    setParallelData(shared_data, PARALLEL_DATA_ELEMENT)
    setParallelData(shared_data_bar, PARALLEL_DATA_FACE)

    populate_buffer = Utils.getSendDataElement
    populate_buffer_rev = Utils.getSendDataFace_resbar

    exchangeData(mesh, sbp, eqn, opts, shared_data, populate_buffer)
    exchangeData(mesh, sbp, eqn, opts, shared_data_bar, populate_buffer_rev)

    finishExchangeData_rev2(mesh, sbp, eqn, opts, shared_data, shared_data_bar,
                           calc_func2)

    testSharedDataElement(mesh, sbp, eqn, opts, shared_data)
    testSharedDataFace(mesh, sbp, eqn, opts, shared_data_bar, 1)

  end  # end facts block

  return nothing
end


function calc_func_rev(mesh, sbp, eqn, opts, data::SharedFaceData, data_bar::SharedFaceData)

  # use Remote metrics to add the polynomial + comm_rank to q_recv
  idx = data_bar.peeridx
  coords = mesh.remote_metrics[idx].coords
  myrank = mesh.myrank

  
  for i=1:size(coords, 3)
    for j=1:size(coords, 2)
      coords_j = sview(coords, :, j, i)
      q_exact = AdvectionEquationMod.calc_p1(eqn.params, coords_j, 0.0) + myrank
      data_bar.q_recv[1, j, i] += q_exact
    end
  end

  return nothing
end


function testSharedDataElement_rev(mesh, sbp, eqn, opts)

  # test that each element of eqn.q_bar has the sum of the contributions of
  # each remote process

  q_bar2 = zeros(eqn.q_bar)
  for i=1:mesh.npeers
    peer = mesh.peer_parts[i]
    for el in mesh.local_element_lists[i]
      for j=1:mesh.numNodesPerElement
        coords_j = sview(mesh.coords, :, j, el)
        q_exact = AdvectionEquationMod.calc_p1(eqn.params, coords_j, 0.0) + peer
        q_bar2[1, j, el] += q_exact
      end
    end
  end

  @test maximum(abs.(q_bar2 - eqn.q_bar)) < 1e-13

  return nothing
end



"""
  Test that face parallel data was sent correctly.  This requires the solution
  on each process is linear and has a constant added to it equal to the
  MPI rank of the process.
"""
function testSharedDataFace(mesh, sbp, eqn, opts, shared_data, offset=0)

  # test the interpolation was exact + the data came from the right place
  # using the peer rank offset
  for i=1:mesh.npeers
    peer_rank = mesh.peer_parts[i]
    coords_peer = mesh.coords_sharedface[i]
    data_i = shared_data[i]
    for j=1:length(data_i.bndries_local)
      for k=1:mesh.numNodesPerFace
        coords_k = sview(coords_peer, :, k, j)
        q_exact = AdvectionEquationMod.calc_p1(eqn.params, coords_k, 0.0) + peer_rank + offset
        @test abs(data_i.q_recv[1, k, j] - q_exact) < 1e-13
      end
    end
  end

  return nothing
end


"""
  Tests element parallel data was sent correctly.  Same requirements as the
  above function
"""
function testSharedDataElement(mesh, sbp, eqn, opts, shared_data)

  q_faceL = zeros(eltype(eqn.q), mesh.numDofPerNode, mesh.numNodesPerFace)
  q_faceR = zeros(q_faceL)
  for i=1:mesh.npeers
    peer_rank = mesh.peer_parts[i]
    data_i = shared_data[i]
    el_offset = mesh.shared_element_offsets[i]
    for j=1:length(data_i.interfaces)
      iface_i = data_i.interfaces[j]
    
      qL = sview(eqn.q, :, :, iface_i.elementL)
      elR = iface_i.elementR - el_offset + 1
      qR = sview(data_i.q_recv, :, :, elR)

      interiorFaceInterpolate!(mesh.sbpface, iface_i, qL, qR, q_faceL, q_faceR)

      for k=1:mesh.numNodesPerFace
        @test abs(q_faceR[1, k] - peer_rank - (q_faceL[1, k] - mesh.myrank)) < 1e-13
      end
    end
  end


  return nothing
end




add_func1!(AdvectionTests, test_parallel2_comm, [TAG_SHORTTEST, TAG_TMP])
