# run tests in parallel (np=4)
"""
  Test the parallel communication primatives in the Utils module
"""
function test_parallel2_comm()

  ARGS[1] = "input_vals_parallel2.jl"
  mesh, sbp, eqn, opts = run_advection(ARGS[1])

  myrank = mesh.myrank
  commsize = mesh.commsize
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

  facts("----- Testing Parallel Communication -----") do
    peer_up = wrap(myrank+1, commsize)
    peer_down = wrap(myrank-1, commsize)

    mesh.npeers = 2
    mesh.peer_parts = [peer_down, peer_up]
    #=
    mesh.send_reqs = Array(MPI.Request, mesh.npeers)
    mesh.recv_reqs = Array(MPI.Request, mesh.npeers)
    mesh.recv_waited= Array(Bool, mesh.npeers)
    mesh.send_waited = Array(Bool, mesh.npeers)

    initMPIStructures(mesh, opts)
    =#

    shared_data = getSharedFaceData(Float64, mesh, sbp, opts)

#    send_data = Array(Array{Float64, 1}, mesh.npeers)
#    recv_data = Array(Array{Float64, 1}, mesh.npeers)
    for i=1:mesh.npeers
      shared_data[i].q_send = reshape(Float64[myrank + i, myrank + i + 1], 2, 1, 1)
      shared_data[i].q_recv = reshape(Array(Float64, mesh.npeers), 2, 1, 1)
    end

    exchangeData(mesh, sbp, eqn, opts, shared_data, populate_buffer, wait=true)

    # peer down: the sent to its peer up
    data = shared_data[1]
    @fact data.q_recv[1] --> peer_down + 2
    @fact data.q_recv[2] --> peer_down + 3

    # peer up: sent to its peer down
    data = shared_data[2]
    @fact data.q_recv[1] --> peer_up + 1
    @fact data.q_recv[2] --> peer_up + 2


    # test exchangeElementData
  #  fill!(eqn.q, 42)
    for i=1:mesh.npeers
      mesh.local_element_lists[i] = [i]
      data_i = shared_data[i]
      data_i.q_send = zeros(Float64, mesh.numDofPerNode, mesh.numNodesPerElement, 1)
      data_i.q_recv = zeros(Float64, mesh.numDofPerNode, mesh.numNodesPerElement, 1)
      
      eqn.q[:,:, i] = i + myrank
    end
    fill!(mesh.recv_waited, true)
    fill!(mesh.send_waited, true)


    exchangeData(mesh, sbp, eqn, opts, shared_data, Utils.getSendDataElement, wait=true)

    data = shared_data[1].q_recv
    for j in data
      @fact j --> peer_down + 2
    end
    data = shared_data[2].q_recv
    for j in data
      @fact j --> peer_up + 1
    end


  end  # end facts block

  return nothing
end

add_func1!(AdvectionTests, test_parallel2_comm, [TAG_SHORTTEST])

#=
function test_adjoint_parallel()

  facts("--- Testing Functional Computation on a Geometric Boundary ---") do
    resize!(ARGS, 1)
    ARGS[1] = "input_vals_functional_DG_parallel.jl"
    include(STARTUP_PATH)

    @fact mesh.isDG --> true
    @fact opts["functional_name1"] --> "qflux"
    @fact opts["functional_error"] --> true
    @fact opts["smb_name"] --> "src/mesh_files/gsquare2np2.smb"
    @fact opts["analytical_functional_val"] --> roughly(2*(exp(1) - 1), atol=1e-12)
    @fact opts["geom_edges_functional1"] --> [1,2]

    fname = "./functional_error1.dat"
    error = readdlm(fname)

    @fact error[1] --> roughly(0.00681567877682826, atol=1e-6)


  end # End facts("--- Testing Functional Computation on a Geometric Boundary ---")

  return nothing
end

add_func1(AdvectionTests, test_adjoint_parallel, [TAG_ADJOINT])
=#
