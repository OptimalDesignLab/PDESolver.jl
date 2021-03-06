using PDESolver
#using Base.Test
using FactCheck
using ODLCommonTools
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
using Utils
using AdvectionEquationMod
using NonlinearSolvers   # non-linear solvers
using ArrayViews

global const STARTUP_PATH = joinpath(Pkg.dir("PDESolver"), "src/solver/advection/startup.jl")
# insert a command line argument
resize!(ARGS, 1)
fname = "input_vals_channel.jl"
mesh, sbp, eqn, opts = solvePDE(fname)

fill!(eqn.res_vec, 0.0)

@testset "--- Testing evalInteriorFlux ---" begin

  for (findex, face) in enumerate(mesh.interfaces)
    
    nbrnodeindex = Array(sbp.numfacenodes:-1:1)
    @test ( nbrnodeindex )== [2,1]

    iL = Array{Int}(2)
    iR = Array{Int}(2)
    iL[1] = mesh.facenodes[1, face.faceL]
    iL[2] = mesh.facenodes[2, face.faceL]
    iR[1] = mesh.facenodes[nbrnodeindex[1], face.faceR]
    iR[2] = mesh.facenodes[nbrnodeindex[2], face.faceR]
    
    @test ( iL[1] )== 2
    @test ( iL[2] )== 3
    @test ( iR[1] )== 1
    @test ( iR[2] )== 3

    evalInteriorFlux(mesh, sbp, eqn, opts)
    println("eqn.res = \n", eqn.res)
    #=
    for j = 1:sbp.numfacenodes
    il = mesh.facenodes[j, face.faceL]::Int
    println("il = ", il)
    ir = mesh.facenodes[nbrnodeindex[j], face.faceR]::Int
    println("ir = ", ir)
    end
    =#
  end
  

end
