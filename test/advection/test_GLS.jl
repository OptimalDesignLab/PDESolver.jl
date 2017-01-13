# test_GLS.jl

# Write tests corresponding to GLS here.

using PDESolver
using FactCheck
using ODLCommonTools
using ArrayViews

include(joinpath(Pkg.dir("PDESolver"), "src/input/read_input.jl"))

# insert a command line argument
resize!(ARGS, 1)
ARGS[1] = "input_vals_GLS.jl"
opts = read_input(ARGS[1])
facts("--- Check input arguments applied correctly ---") do
  @fact opts["use_GLS"] --> true
  @fact opts["Q_transpose"] --> true
end


include("../../src/solver/advection/startup.jl")
include("../../src/solver/advection/GLS.jl")

facts("--- Check functions in ../src/solver/advection/GLS.jl ---") do

  context("Checking shape function derivatives") do
  	sbp2 = TriSBP{Float64}()
    shapefuncderiv = zeros(Tsol, sbp2.numnodes, sbp2.numnodes, Tdim)
    calcShapefuncDeriv(sbp2, shapefuncderiv)
    @fact shapefuncderiv[:,:,1] --> roughly([-0.5  0.5  0.0
                                             -0.5  0.5  0.0
                                             -0.5  0.5  0.0])
    @fact shapefuncderiv[:,:,2] --> roughly([-0.5  0.0  0.5
                                             -0.5  0.0  0.5
                                             -0.5  0.0  0.5])
  end

  context("Check calcAxidxi") do
  	Tdim = 2
  	alpha_x = ones(Tsol, mesh.numNodesPerElement)
  	alpha_y = ones(alpha_x)
  	dxidx = sview(mesh.dxidx,:,:,:,1)
    shapefuncderiv = zeros(Tsol, sbp.numnodes, sbp.numnodes, Tdim)
    calcShapefuncDeriv(sbp, shapefuncderiv)
    AxiDxi = calcAxiDxi(mesh, dxidx, alpha_x, alpha_y, shapefuncderiv)
    @fact AxiDxi --> roughly([-0.5 0.25 0.25
    	                      -0.5 0.25 0.25
    	                      -0.5 0.25 0.25])
  end

  context("Check calcTau") do
    fill!(eqn.params.alpha_x, 1.0)
    fill!(eqn.params.alpha_y, 1.0)
    tau = zeros(Tsol,mesh.numNodesPerElement, mesh.numEl)
    calcTau(mesh, sbp, eqn, tau)
    @fact mesh.numEl --> 1
    @fact tau[:,1] --> roughly([1.0,1.0,1.0])
  end

  context("Check GLS") do
    fill!(eqn.params.alpha_x, 1.0)
    fill!(eqn.params.alpha_y, 1.0)
    fill!(eqn.res, 0.0)
    eqn.q[1,:,1] = eqn.q_vec
    GLS(mesh, sbp, eqn)
    @fact eqn.res[1,:,1] --> roughly([8.0 -4.0 -4.0])

  end

end  # end facts
