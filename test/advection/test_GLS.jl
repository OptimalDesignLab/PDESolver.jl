# test_GLS.jl

# Write tests corresponding to GLS here.

using PDESolver
using Base.Test
using ODLCommonTools
using ArrayViews

include(joinpath(Pkg.dir("PDESolver"), "src/input/read_input.jl"))

# insert a command line argument
resize!(ARGS, 1)
ARGS[1] = "input_vals_GLS.jl"
opts = read_input(ARGS[1])
@testset "--- Check input arguments applied correctly ---" begin
  @test ( opts["use_GLS"] )== true
  @test ( opts["Q_transpose"] )== true
end


include("../../src/solver/advection/startup.jl")
include("../../src/solver/advection/GLS.jl")

@testset "--- Check functions in ../src/solver/advection/GLS.jl ---" begin

  @testset "Checking shape function derivatives" begin
    sbp2 = getTriSBPGamma{Float64}()
    shapefuncderiv = zeros(Tsol, sbp2.numnodes, sbp2.numnodes, Tdim)
    calcShapefuncDeriv(sbp2, shapefuncderiv)
    @test isapprox(shapefuncderiv[:,:,1] -  [-0.5  0.5  0.0
                                             -0.5  0.5  0.0
                                             -0.5  0.5  0.0]) atol=1e-13
    @test isapprox(shapefuncderiv[:,:,2] -  [-0.5  0.0  0.5
                                             -0.5  0.0  0.5
                                             -0.5  0.0  0.5]) atol=1e-13
  end

  @testset "Check calcAxidxi" begin
    Tdim = 2
    alpha_x = ones(Tsol, mesh.numNodesPerElement)
    alpha_y = ones(alpha_x)
    dxidx = sview(mesh.dxidx,:,:,:,1)
    shapefuncderiv = zeros(Tsol, sbp.numnodes, sbp.numnodes, Tdim)
    calcShapefuncDeriv(sbp, shapefuncderiv)
    AxiDxi = calcAxiDxi(mesh, dxidx, alpha_x, alpha_y, shapefuncderiv)
    @test isapprox(AxiDxi -[-0.5 0.25 0.25
                            -0.5 0.25 0.25
                            -0.5 0.25 0.25]) atol=1e-13
  end

  @testset "Check calcTau" begin
    fill!(eqn.params.alpha_x, 1.0)
    fill!(eqn.params.alpha_y, 1.0)
    tau = zeros(Tsol,mesh.numNodesPerElement, mesh.numEl)
    calcTau(mesh, sbp, eqn, tau)
    @test ( mesh.numEl )== 1
    @test isapprox( tau[:,1], [1.0,1.0,1.0]) 
  end

  @testset "Check GLS" begin
    fill!(eqn.params.alpha_x, 1.0)
    fill!(eqn.params.alpha_y, 1.0)
    fill!(eqn.res, 0.0)
    eqn.q[1,:,1] = eqn.q_vec
    GLS(mesh, sbp, eqn)
    @test isapprox( eqn.res[1,:,1], [8.0 -4.0 -4.0]) 

  end

end  # end facts
