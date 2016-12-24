push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/advection"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))

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
# insert a command line argument
resize!(ARGS, 1)
ARGS[1] = "input_vals_channel.jl"
include("../../src/solver/advection/startup_advection.jl")  # initialization and construction
include("../../src/solver/advection/DG_advection.jl")
include("../../src/solver/advection/bc_solvers.jl")
include("../../src/solver/euler/complexify.jl")

fill!(eqn.res_vec, 0.0)
println("eqn.q_vec = \n", eqn.q_vec)
println("eqn.q = \n", eqn.q)

facts("--- Testing evalInteriorFlux ---") do

  for (findex, face) in enumerate(mesh.interfaces)
    @fact sbp.numfacenodes --> 2
    
    nbrnodeindex = Array(sbp.numfacenodes:-1:1)
    @fact nbrnodeindex --> [2,1]

    iL = Array(Int, 2)
    iR = Array(Int, 2)
    iL[1] = sbp.facenodes[1, face.faceL]
    iL[2] = sbp.facenodes[2, face.faceL]
    iR[1] = sbp.facenodes[nbrnodeindex[1], face.faceR]
    iR[2] = sbp.facenodes[nbrnodeindex[2], face.faceR]
    
    @fact iL[1] --> 2
    @fact iL[2] --> 3
    @fact iR[1] --> 1
    @fact iR[2] --> 3

    evalInteriorFlux(mesh, sbp, eqn, opts)
    println("eqn.res = \n", eqn.res)
    #=
  	for j = 1:sbp.numfacenodes
	  il = sbp.facenodes[j, face.faceL]::Int
	  println("il = ", il)
	  ir = sbp.facenodes[nbrnodeindex[j], face.faceR]::Int
	  println("ir = ", ir)
  	end
  	=#
  end
  

end