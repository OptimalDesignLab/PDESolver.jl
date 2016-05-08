#=
push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/advection"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))
include(joinpath(Pkg.dir("PDESolver"), "src/input/make_input.jl"))

using PDESolver
#using Base.Test
using FactCheck
using ODLCommonTools
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
using AdvectionEquationMod
using ForwardDiff
using NonlinearSolvers   # non-linear solvers
using ArrayViews

global const STARTUP_PATH = joinpath(Pkg.dir("PDESolver"), "src/solver/advection/startup_advection.jl")
=#
resize!(ARGS, 1)

ARGS[1] = "input_vals_parallel.jl"
include(STARTUP_PATH)

# test the Utils parallel functions
facts("----- Testing Parallel Functions -----") do
  mesh.npeers = 2
  mesh.send_reqs = Array(MPI.Request, mesh.npeers)
  mesh.recv_reqs = Array(MPI.Request, mesh.npeers)
  mesh.send_stats = Array(MPI.Status, mesh.npeers)
  mesh.recv_stats = Array(MPI.Status, mesh.npeers)
  mesh.send_waited = Array(Bool, mesh.npeers)
  mesh.recv_waited = Array(Bool, mesh.npeers)

  initMPIStructures(mesh, opts)

  for i=1:mesh.npeers
    @fact mesh.send_reqs[i] --> MPI.REQUEST_NULL
    @fact mesh.recv_reqs[i] --> MPI.REQUEST_NULL
  end

  buff = zeros(Float64, mesh.numDofPerNode, mesh.numNodesPerFace, length(mesh.bndryfaces))
  getSendData(mesh, opts, eqn.q, mesh.bndryfaces, buff, MPI.REQUEST_NULL, true)

  # verify that the inteperpolation was exact
  for i=1:length(mesh.bndryfaces)
    for j=1:mesh.numNodesPerFace
      coords = mesh.coords_bndry[:, j, i]
      val_exp = AdvectionEquationMod.calc_p4(coords, 1, 1, 0.0)
      @fact buff[1, j, i] --> roughly(val_exp, atol=1e-13)
    end
  end

  nvals = 5
  vals = rand(nvals)
  vals_recv = zeros(Float64, nvals)
  req = MPI.Isend(vals, mesh.myrank, 1, mesh.comm)
  req2 = MPI.Irecv!(vals_recv, mesh.myrank, 1, mesh.comm)

  stat = MPI.Wait!(req2)
  # if the assertiosn do not trigger, then the test passes
  verifyCommunication(mesh, opts, vals_recv, mesh.myrank, stat)
  @fact vals_recv --> roughly(vals, atol=1e-13)

  req = MPI.Isend(vals, mesh.myrank, 1, mesh.comm)
  req2 = MPI.Irecv!(vals_recv, mesh.myrank, 1, mesh.comm)

  stat = MPI.Wait!(req2)
  @fact_throws AssertionError verifyCommunication(mesh, opts, vals_recv[1:2], mesh.myrank, stat)


  # test calcSharedFaceIntegrals
  mesh.npeers = 1
  mesh.peer_face_counts = [2]
  mesh.send_reqs = Array(MPI.Request, mesh.npeers)
  mesh.recv_reqs = Array(MPI.Request, mesh.npeers)
  mesh.send_stats = Array(MPI.Status, mesh.npeers)
  mesh.recv_stats = Array(MPI.Status, mesh.npeers)
  mesh.send_waited = Array(Bool, mesh.npeers)
  mesh.recv_waited = Array(Bool, mesh.npeers)


  initMPIStructures(mesh, opts)

  # do a dummy send
  a = [1]
  ar = [1]
  mesh.send_reqs[1] = MPI.Isend(a, mesh.myrank, 1, mesh.comm)
  mesh.recv_reqs[1] = MPI.Irecv!(ar, mesh.myrank, 1, mesh.comm)
  eqn.q_face_send = Array(Array{Float64, 3}, 1)
  eqn.q_face_recv = Array(Array{Float64, 3}, 1)
  eqn.q_face_send[1] = Array(Float64, mesh.numDofPerNode, mesh.numNodesPerFace, 2)
  eqn.q_face_recv[1] = Array(Float64, mesh.numDofPerNode, mesh.numNodesPerFace, 2)
  mesh.dxidx_sharedface = Array(Array{Float64, 4}, 1)
  mesh.dxidx_sharedface[1] = zeros(2,2, mesh.numNodesPerFace, 2)
  dxidx_arr = mesh.dxidx_sharedface[1]

  interiorfaceinterpolate!(mesh.sbpface, mesh.interfaces, eqn.q, eqn.q_face)
  for i=1:mesh.peer_face_counts[1]
    iface = mesh.interfaces[i]
    for k=1:mesh.numNodesPerFace
      dxidx_arr[:, :, k, i] = mesh.dxidx[:, :, 1, iface.elementL]
      for j=1:mesh.numDofPerNode
        eqn.q_face_send[1][j, k, i] = eqn.q_face[j, 1, k, i]
        eqn.q_face_recv[1][j, k, i] = eqn.q_face[j, 2, k, i]
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

  mesh.shared_interfaces = Array(Array{Interface, 1}, 1)
  mesh.shared_interfaces[1] = ifaces_orig[1:2]
  # create the boundary array
  mesh.bndries_local = Array(Array{Boundary, 1}, 1)
  mesh.bndries_local[1] = Array(Boundary, 2)
  for i=1:2
    mesh.bndries_local[1][i] = Boundary(ifaces_orig[i].elementL, ifaces_orig[i].faceL)
  end
  # use the identity permutation because the values were already permuted
  for i=1:mesh.numNodesPerFace
    mesh.sbpface.nbrperm[i, 1] = i
  end
  AdvectionEquationMod.calcSharedFaceIntegrals(mesh, sbp, eqn, opts, eqn.flux_func)

  for i=1:length(eqn.flux_sharedface[1])
    @fact eqn.flux_sharedface[1][i] --> roughly(eqn.flux_face[i], atol=1e-13)
  end


  # to a serial rk4 run to compare against later
  start_dir = pwd()
  cd ("./rk4/serial")
  opts["order"] = 1
  opts["solve"] = true
  fname = "input_vals_parallel_run"
  make_input(opts, fname)

  ARGS[1] = string(fname, ".jl")
  include(STARTUP_PATH)
  cd(start_dir)

  # make the parallel version
  opts["smb_name"] = "src/mesh_files/psquare2.smb"
  opts["dmg_name"] = "src/mesh_files/psquare2.dmg"
  make_input( opts, string("./rk4/parallel/", fname, "p"))

end
