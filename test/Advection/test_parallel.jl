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

  initMPIStructures(mesh, opts)

  for i=1:mesh.npeers
    @fact mesh.send_reqs[i] --> MPI.REQUEST_NULL
    @fact mesh.recv_reqs[i] --> MPI.REQUEST_NULL
  end

  buff = zeros(Float64, mesh.numDofPerNode, mesh.numNodesPerFace, length(mesh.bndryfaces))
  getSendData(mesh, opts, eqn.q, mesh.bndryfaces, buff, MPI.REQUEST_NULL)

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



end
