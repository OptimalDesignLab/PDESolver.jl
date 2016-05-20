# run tests in parallel


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
ARGS[1] = "input_vals_parallel2.jl"
include(STARTUP_PATH)

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

facts("----- Testing Parallel Communication -----") do
  peer_up = wrap(myrank+1, commsize)
  peer_down = wrap(myrank-1, commsize)

  mesh.npeers = 2
  mesh.peer_parts = [peer_down, peer_up]
  initMPIStructures(mesh, opts)

  send_data = Array(Array{Float64, 1}, mesh.npeers)
  recv_data = Array(Array{Float64, 1}, mesh.npeers)
  for i=1:mesh.npeers
    send_data[i] = Float64[myrank + i, myrank + i + 1]
    recv_data[i] = Array(Float64, mesh.npeers)
  end

  exchangeFaceData(mesh, opts, send_data, recv_data, wait=true)

  # peer down: the sent to its peer up
  data = recv_data[1]
  @fact data[1] --> peer_down + 2
  @fact data[2] --> peer_down + 3

  # peer up: sent to its peer down
  data = recv_data[2]
  @fact data[1] --> peer_up + 1
  @fact data[2] --> peer_up + 2

end


