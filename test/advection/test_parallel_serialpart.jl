# serial part of parallel tests

"""
  Test the parallel communication constructs.
"""
function test_parallel_mpi()
  fname = "input_vals_parallel.jl"
  mesh, sbp, eqn, opts = solvePDE(fname)

  # test the Utils parallel functions
  @testset "----- Testing Parallel Functions -----" begin

    test_tagmanager()


    Tsol = eltype(eqn.q)
    mesh.npeers = 1
    nfaces = length(mesh.bndryfaces)
    shared_data = Array{SharedFaceData{Tsol}}(mesh.npeers)
    resize!(mesh.bndries_local, mesh.npeers)
    resize!(mesh.bndries_remote, mesh.npeers)
    resize!(mesh.shared_interfaces, mesh.npeers)
    resize!(mesh.peer_parts, mesh.npeers)
    resize!(mesh.peer_face_counts, mesh.npeers)
    resize!(mesh.local_element_counts, mesh.npeers)
    resize!(mesh.remote_element_counts, mesh.npeers)
    fill!(mesh.peer_parts, mesh.myrank)
    fill!(mesh.peer_face_counts, nfaces)
    fill!(mesh.local_element_counts, nfaces)
    fill!(mesh.remote_element_counts, nfaces)
    for i=1:mesh.npeers
      println("creating SharedFaceData for peer ", i, ", Tsol = ", Tsol)
      mesh.bndries_local[i] = mesh.bndryfaces
      mesh.bndries_remote[i] = mesh.bndryfaces
      mesh.shared_interfaces[i] = Array{Interface}(0)
      shared_data[i] = SharedFaceData(Tsol, mesh, i, opts["parallel_data"], 100)
    end

    for i=1:mesh.npeers
      @test ( shared_data[i].send_req )== MPI.REQUEST_NULL
      @test ( shared_data[i].recv_req )== MPI.REQUEST_NULL
    end

    data_i = shared_data[1]
    Utils.getSendDataFace(mesh, sbp, eqn, opts, data_i)
    buff = data_i.q_send

    # verify that the inteperpolation was exact
    eqn.params.alpha_x = 1
    eqn.params.alpha_y = 1
    for i=1:length(mesh.bndryfaces)
      for j=1:mesh.numNodesPerFace
        coords = mesh.coords_bndry[:, j, i]
        val_exp = AdvectionEquationMod.calc_p4(eqn.params, coords, 0.0)
        @test isapprox( buff[1, j, i], val_exp) atol=1e-13
      end
    end

    rand!(data_i.q_send)
    data_i.send_req = MPI.Isend(data_i.q_send, data_i.myrank, 1, data_i.comm)
    data_i.recv_req = MPI.Irecv!(data_i.q_recv, data_i.myrank, 1, data_i.comm)

    data_i.recv_status = MPI.Wait!(data_i.recv_req)
    # if the assertiosn do not trigger, then the test passes
    verifyReceiveCommunication(data_i)
    @test isapprox( data_i.q_send, data_i.q_recv) atol=1e-13

    # test calcSharedFaceIntegrals
    mesh.npeers = 1
    mesh.peer_face_counts = [2]

    # do a dummy send
    a = [1]
    ar = [1]
    data_i.send_req = MPI.Isend(a, mesh.myrank, 1, mesh.comm)
    data_i.recv_req = MPI.Irecv!(ar, mesh.myrank, 1, mesh.comm)
    data_i.q_send = Array{Float64}(mesh.numDofPerNode, mesh.numNodesPerFace, 2)
    data_i.q_recv = Array{Float64}(mesh.numDofPerNode, mesh.numNodesPerFace, 2)
    mesh.nrm_sharedface = Array{Array{Float64, 3}}(1)
    mesh.nrm_sharedface[1] = zeros(2, mesh.numNodesPerFace, 2)
    nrm_arr = mesh.nrm_sharedface[1]

    interiorfaceinterpolate!(mesh.sbpface, mesh.interfaces, eqn.q, eqn.q_face)
    for i=1:mesh.peer_face_counts[1]
      iface = mesh.interfaces[i]
      for k=1:mesh.numNodesPerFace
        nrm_arr[:, k, i] = mesh.nrm_face[ :, k, i]
        for j=1:mesh.numDofPerNode
          data_i.q_send[j, k, i] = eqn.q_face[j, 1, k, i]
          data_i.q_recv[j, k, i] = eqn.q_face[j, 2, k, i]
        end
      end
    end

    # get the face flux using calcFaceFlux
    ifaces_orig = mesh.interfaces
    mesh.interfaces = ifaces_orig[1:2]
    fill!(eqn.res, 0.0)
    eqn.flux_face = zeros(mesh.numDofPerNode, mesh.numNodesPerFace, 2)
    AdvectionEquationMod.calcFaceFlux(mesh, sbp, eqn, eqn.flux_func, mesh.interfaces, eqn.flux_face)

    # now get it using Boundary integrate
    eqn.flux_sharedface = Array{Array{Float64, 3}}(1)
    eqn.flux_sharedface[1] = Array{Float64}(mesh.numDofPerNode, mesh.numNodesPerFace, 2)

    data_i.interfaces = ifaces_orig[1:2]
    # create the boundary array
    data_i.bndries_local = Array{Boundary}(2)
    for i=1:2
      data_i.bndries_local[i] = Boundary(ifaces_orig[i].elementL, ifaces_orig[i].faceL)
    end
    # use the identity permutation because the values were already permuted
    for i=1:mesh.numNodesPerFace
      mesh.sbpface.nbrperm[i, 1] = i
    end
    AdvectionEquationMod.calcSharedFaceIntegrals(mesh, sbp, eqn, opts, data_i)

    for i=1:length(eqn.flux_sharedface[1])
      @test isapprox( eqn.flux_sharedface[1][i], eqn.flux_face[i]) atol=1e-13
    end

    # test changing buffers
    @test getParallelData(shared_data) == PARALLEL_DATA_FACE
    @test size(data_i.q_send, 1) == mesh.numDofPerNode
    @test size(data_i.q_send, 2) == mesh.numNodesPerFace
    @test size(data_i.q_send, 3) == mesh.peer_face_counts[1]

    @test size(data_i.q_recv, 1) == mesh.numDofPerNode
    @test size(data_i.q_recv, 2) == mesh.numNodesPerFace
    @test size(data_i.q_recv, 3) == mesh.peer_face_counts[1]

    setParallelData(shared_data, PARALLEL_DATA_ELEMENT)
    @test getParallelData(shared_data) == PARALLEL_DATA_ELEMENT
    @test size(data_i.q_send, 1) == mesh.numDofPerNode
    @test size(data_i.q_send, 2) == mesh.numNodesPerElement
    @test size(data_i.q_send, 3) == mesh.local_element_counts[1]

    @test size(data_i.q_recv, 1) == mesh.numDofPerNode
    @test size(data_i.q_recv, 2) == mesh.numNodesPerElement
    @test size(data_i.q_recv, 3) == mesh.remote_element_counts[1]

    # test copy
    data_i2 = copy(data_i)
    data_i2.q_send[1] = 10
    data_i2.q_recv[1] = 10
    @test data_i2.q_send[1] != data_i.q_send[1]
    @test data_i2.q_recv[1] != data_i.q_recv[1]
    @test pointer(data_i2.q_send) != pointer(data_i.q_send)
    @test pointer(data_i2.q_recv) != pointer(data_i.q_recv)


    copy!(data_i2, data_i)
    @test data_i2.q_send[1] == data_i.q_send[1]
    @test data_i2.q_recv[1] == data_i.q_recv[1]

    old_tag = shared_data[1].tag
    setNewTag(shared_data)
    @test old_tag != shared_data[1].tag

  end  # end facts block

  return nothing
end

function test_tagmanager()
  @testset "MPITagManager" begin
    mgr = MPITagManager()

    for i=1:10
      @test getNextTag(mgr) == i
    end

    # test using tags
    markTagUsed(mgr, 12)
    @test getNextTag(mgr) == 11
    @test getNextTag(mgr) == 13

    freeTag(mgr, 3)
    freeTag(mgr, 4)
    @test getNextTag(mgr) == 3
    @test getNextTag(mgr) == 4
    @test getNextTag(mgr) == 14

    freeTag(mgr, 5)
    freeTag(mgr, 4)
    freeTag(mgr, 6)
    @test getNextTag(mgr) == 4
    @test getNextTag(mgr) == 5
    @test getNextTag(mgr) == 6

    freeTag(mgr, 3)
    @test_throws Exception freeTag(mgr, 3)

    @test_throws Exception markTagUsed(mgr, 4)


    # test starting tag = 5
    mgr = MPITagManager(5)

    for i=1:10
      @test getNextTag(mgr) == i + 5 - 1
    end

    markTagUsed(mgr, 16)
    @test getNextTag(mgr) == 15
    @test getNextTag(mgr) == 17

    freeTag(mgr, 9)
    freeTag(mgr, 7)
    freeTag(mgr, 8)
    @test getNextTag(mgr) == 7
    @test getNextTag(mgr) == 8
    @test getNextTag(mgr) == 9
    @test getNextTag(mgr) == 18

    @test_throws Exception freeTag(mgr, 2)
    @test_throws Exception markTagUsed(mgr, 2)

    freeTag(mgr, 7)
    @test_throws Exception freeTag(mgr, 7)


  end

  return nothing
end


#test_parallel_mpi()
add_func1!(AdvectionTests, test_parallel_mpi, [TAG_SHORTTEST])

function test_parallel_serialpart()
  fname = "input_vals_parallel.jl"
  mesh, sbp, eqn, opts = solvePDE(fname)

  # do a serial rk4 run to compare against later
  start_dir = pwd()
  cd("./rk4/serial")
  opts["order"] = 1
  opts["solve"] = true
  fname = "input_vals_parallel_run"
  make_input(opts, fname)

  fname2 = string(fname, ".jl")
  mesh, sbp, eqn, opts = solvePDE(fname2)
  cd(start_dir)

  # make the parallel version
  opts["smb_name"] = "SRCMESHES/psquare2.smb"
  opts["dmg_name"] = "SRCMESHES/psquare2.dmg"
  make_input( opts, string("./rk4/parallel/", fname, "p"))

  # serial newton
  start_dir = pwd()
  cd("./newton/serial")
  fname = "input_vals_serial.jl"
  mesh, sbp, eqn, opts = solvePDE(fname)
  cd(start_dir)
  opts["smb_name"] = "SRCMESHES/psquare2.smb"
  opts["dmg_name"] = "SRCMESHES/psquare2.dmg"
  make_input(opts, string("./newton/parallel/", "input_vals_parallel"))


  # 3D rk4
  start_dir = pwd()
  cd("./rk4_3d/serial")
  fname = "input_vals_rk4_3d.jl"
  mesh, sbp, eqn, opts = solvePDE(fname)

  # make the parallel version
  cd(start_dir)
  opts["smb_name"] = "SRCMESHES/ptet8cube.smb"
  make_input(opts, string("./rk4_3d/parallel/", "input_vals_parallel"))

  # 3D Newton's method
  opts["run_type"] = 5
  opts["jac_method"] = 2
  opts["jac_type"] = 3
  opts["parallel_type"] = 2
  opts["parallel_data"] = PARALLEL_DATA_ELEMENT
  make_input(opts, string("./newton_3d/parallel/", "input_vals_parallel"))

  opts["smb_name"] = "SRCMESHES/tet8cube.smb"
  make_input(opts, string("./newton_3d/serial/", "input_vals_serial"))

  cd("./newton_3d/serial")
  fname = "input_vals_serial.jl"
  mesh, sbp, eqn, opts = solvePDE(fname)

  cd(start_dir)
end

#test_parallel_serialpart()
add_func1!(AdvectionTests, test_parallel_serialpart, [TAG_SHORTTEST])
