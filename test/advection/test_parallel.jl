# serial part of parallel tests

"""
  Test the parallel communication constructs.
"""
function test_parallel_mpi()
  ARGS[1] = "input_vals_parallel.jl"
  mesh, sbp, eqn, opts = run_advection(ARGS[1])

  # test the Utils parallel functions
  facts("----- Testing Parallel Functions -----") do
    mesh.npeers = 1
    nfaces = length(mesh.bndryfaces)
    shared_data = Array(SharedFaceData, mesh.npeers)
    resize!(mesh.bndries_local, mesh.npeers)
    resize!(mesh.bndries_remote, mesh.npeers)
    resize!(mesh.shared_interfaces, mesh.npeers)
    resize!(mesh.peer_parts, mesh.npeers)
    fill!(mesh.peer_parts, mesh.myrank)
    for i=1:mesh.npeers
      mesh.bndries_local[i] = mesh.bndryfaces
      mesh.bndries_remote[i] = mesh.bndryfaces
      mesh.shared_interfaces[i] = Array(Interface, 0)
      q_send = zeros(mesh.numDofPerNode, mesh.numNodesPerFace, nfaces)
      q_recv = zeros(q_send)
      shared_data[i] = SharedFaceData(mesh, i, q_send, q_recv)
    end

    for i=1:mesh.npeers
      @fact shared_data[i].send_req --> MPI.REQUEST_NULL
      @fact shared_data[i].recv_req --> MPI.REQUEST_NULL
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
        @fact buff[1, j, i] --> roughly(val_exp, atol=1e-13)
      end
    end

    rand!(data_i.q_send)
    data_i.send_req = MPI.Isend(data_i.q_send, data_i.myrank, 1, data_i.comm)
    data_i.recv_req = MPI.Irecv!(data_i.q_recv, data_i.myrank, 1, data_i.comm)

    data_i.recv_status = MPI.Wait!(data_i.recv_req)
    # if the assertiosn do not trigger, then the test passes
    verifyReceiveCommunication(data_i)
    @fact data_i.q_send --> roughly(data_i.q_recv, atol=1e-13)

    # test calcSharedFaceIntegrals
    mesh.npeers = 1
    mesh.peer_face_counts = [2]

    # do a dummy send
    a = [1]
    ar = [1]
    data_i.send_req = MPI.Isend(a, mesh.myrank, 1, mesh.comm)
    data_i.recv_req = MPI.Irecv!(ar, mesh.myrank, 1, mesh.comm)
    data_i.q_send = Array(Float64, mesh.numDofPerNode, mesh.numNodesPerFace, 2)
    data_i.q_recv = Array(Float64, mesh.numDofPerNode, mesh.numNodesPerFace, 2)
    mesh.nrm_sharedface = Array(Array{Float64, 3}, 1)
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
    eqn.flux_sharedface = Array(Array{Float64, 3}, 1)
    eqn.flux_sharedface[1] = Array(Float64, mesh.numDofPerNode, mesh.numNodesPerFace, 2)

    data_i.interfaces = ifaces_orig[1:2]
    # create the boundary array
    data_i.bndries_local = Array(Boundary, 2)
    for i=1:2
      data_i.bndries_local[i] = Boundary(ifaces_orig[i].elementL, ifaces_orig[i].faceL)
    end
    # use the identity permutation because the values were already permuted
    for i=1:mesh.numNodesPerFace
      mesh.sbpface.nbrperm[i, 1] = i
    end
    AdvectionEquationMod.calcSharedFaceIntegrals(mesh, sbp, eqn, opts, data_i)

    for i=1:length(eqn.flux_sharedface[1])
      @fact eqn.flux_sharedface[1][i] --> roughly(eqn.flux_face[i], atol=1e-13)
    end
  end  # end facts block

  return nothing
end

#test_parallel_mpi()
add_func1!(AdvectionTests, test_parallel_mpi, [TAG_SHORTTEST])

function test_parallel_serialpart()
  ARGS[1] = "input_vals_parallel.jl"
  mesh, sbp, eqn, opts = run_advection(ARGS[1])

  # do a serial rk4 run to compare against later
  start_dir = pwd()
  cd("./rk4/serial")
  opts["order"] = 1
  opts["solve"] = true
  fname = "input_vals_parallel_run"
  make_input(opts, fname)

  ARGS[1] = string(fname, ".jl")
  mesh, sbp, eqn, opts = run_advection(ARGS[1])
  cd(start_dir)

  # make the parallel version
  opts["smb_name"] = "SRCMESHES/psquare2.smb"
  opts["dmg_name"] = "SRCMESHES/psquare2.dmg"
  make_input( opts, string("./rk4/parallel/", fname, "p"))

  # serial newton
  start_dir = pwd()
  cd("./newton/serial")
  ARGS[1] = "input_vals_serial.jl"
  mesh, sbp, eqn, opts = run_advection(ARGS[1])
  cd(start_dir)
  opts["smb_name"] = "SRCMESHES/psquare2.smb"
  opts["dmg_name"] = "SRCMESHES/psquare2.dmg"
  make_input(opts, string("./newton/parallel/", "input_vals_parallel"))


  # 3D rk4
  start_dir = pwd()
  cd("./rk4_3d/serial")
  ARGS[1] = "input_vals_rk4_3d.jl"
  mesh, sbp, eqn, opts = run_advection(ARGS[1])

  # make the parallel version
  cd(start_dir)
  opts["smb_name"] = "SRCMESHES/ptet8cube.smb"
  make_input(opts, string("./rk4_3d/parallel/", "input_vals_parallel"))

  # 3D Newton's method
  opts["run_type"] = 5
  opts["jac_method"] = 2
  opts["jac_type"] = 3
  opts["parallel_type"] = 2
  opts["parallel_data"] = "element"
  make_input(opts, string("./newton_3d/parallel/", "input_vals_parallel"))

  opts["smb_name"] = "SRCMESHES/tet8cube.smb"
  make_input(opts, string("./newton_3d/serial/", "input_vals_serial"))

  cd("./newton_3d/serial")
  ARGS[1] = "input_vals_serial.jl"
  mesh, sbp, eqn, opts = run_advection(ARGS[1])

  cd(start_dir)
end

#test_parallel_serialpart()
add_func1!(AdvectionTests, test_parallel_serialpart, [TAG_SHORTTEST])
