push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/euler"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))


using PDESolver
#using Base.Test
using FactCheck
using ODLCommonTools
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
using EulerEquationMod
using ForwardDiff
using NonlinearSolvers   # non-linear solvers
using ArrayViews
include( joinpath(Pkg.dir("PDESolver"), "src/solver/euler/complexify.jl"))
include( joinpath(Pkg.dir("PDESolver"), "src/input/make_input.jl"))
global const STARTUP_PATH = joinpath(Pkg.dir("PDESolver"), "src/solver/euler/startup.jl")
# insert a command line argument
resize!(ARGS, 1)
ARGS[1] = "input_vals_3d.jl"
include(STARTUP_PATH)

# create entropy variable param type
params_e = EulerEquationMod.ParamType{3, :entropy, Float64, Float64, Float64}(mesh, sbp, opts, 1)

facts("----- Testing Conversion -----") do

  # test unsafe conversion kernels
  q = [1., 2, 3, 4, 15]
  q2 = zeros(q)
  EulerEquationMod.convertToEntropy_(eqn.params, q, q2)
  q3 = zeros(q)
  EulerEquationMod.convertToConservative_(eqn.params, q2, q3)
  @fact q3 --> roughly(q, atol=1e-13)

  # test in-place conversion
  q2 = copy(q)
  EulerEquationMod.convertToEntropy_(eqn.params, q2, q2)
  EulerEquationMod.convertToConservative_(eqn.params, q2, q2)
  @fact q2 --> roughly(q, atol=1e-13)

  # test safe interface
  EulerEquationMod.convertToEntropy(eqn.params, q, q2)
  EulerEquationMod.convertToEntropy_(eqn.params, q, q3)
  @fact q2 --> roughly(q3, atol=1e-13)

  EulerEquationMod.convertToConservative(eqn.params, q, q2)
  @fact q2 --> roughly(q, atol=1e-13)

  EulerEquationMod.convertToEntropy_(eqn.params, q, q2)
  EulerEquationMod.convertToEntropy(params_e, q2, q3)
  @fact q3 --> roughly(q2, atol=1e-13)

  EulerEquationMod.convertToConservative(params_e, q2, q3)
  @fact q3 --> roughly(q, atol=1e-13)

  
end

facts("----- Testing Flux Calculation -----") do
  q = [1., 2, 3, 4, 15]
  F1 = [2, 4.2, 6, 8, 30.4]
  F2 = [3, 6, 9.2, 12, 45.6]
  F3 = [4, 8, 12, 16.2, 60.8]
  q2 = zeros(q)
  EulerEquationMod.convertToEntropy(eqn.params, q, q2)

  p = EulerEquationMod.calcPressure(eqn.params, q)
  p2 = EulerEquationMod.calcPressure(params_e, q2)
  @fact p --> roughly(0.2, atol=1e-12)
  @fact p2 --> roughly(p, atol=1e-10)
  aux_vars = [p]

  F = zeros(Float64, 5)
  Fe = zeros(F)
  nrm = [1., 0, 0]
  EulerEquationMod.calcEulerFlux(eqn.params, q, aux_vars, nrm, F)
  EulerEquationMod.calcEulerFlux(params_e, q2, aux_vars, nrm, Fe)
  @fact F --> roughly(F1, atol=1e-12)
  @fact Fe --> roughly(F, atol=1e-10)

  nrm = [0., 1, 0]
  EulerEquationMod.calcEulerFlux(eqn.params, q, aux_vars, nrm, F)
  EulerEquationMod.calcEulerFlux(params_e, q2, aux_vars, nrm, Fe)
  @fact F --> roughly(F2, atol=1e-12)
  @fact Fe --> roughly(F, atol=1e-10)

  nrm = [0., 0, 1]
  EulerEquationMod.calcEulerFlux(eqn.params, q, aux_vars, nrm, F)
  EulerEquationMod.calcEulerFlux(params_e, q2, aux_vars, nrm, Fe)
  @fact F --> roughly(F3, atol=1e-12)
  @fact Fe --> roughly(F, atol=1e-10)
end

facts("----- Testing Misc. Functions -----") do

  q = [1., 2, 3, 4, 15]
  q2 = zeros(q)
  EulerEquationMod.convertToEntropy(eqn.params, q, q2)
  a = EulerEquationMod.calcSpeedofSound(eqn.params, q)
  @fact a --> roughly(sqrt(0.28), atol=1e-13)
  a2 = EulerEquationMod.calcSpeedofSound(params_e, q2)
  @fact a2 --> roughly(a, atol=1e-13)

  s = EulerEquationMod.calcEntropy(eqn.params, q)
  @fact s --> roughly(log(0.4*0.5), atol=1e-12)
  s2 = EulerEquationMod.calcEntropy(params_e, q2)
  @fact s2 --> roughly(s, atol=1e-12)

end

facts("----- Testing Coefficient Matrices calculations -----") do
  q = [1., 2, 3, 4, 15]
  q2 = zeros(q)
  EulerEquationMod.convertToEntropy(eqn.params, q, q2)

  A0 = [-4 -8 -12 -16 -60;
         -8 -16.8 -24 -32 -121.6;
         -12 -24 -36.8 -48 -182.4;
         -16 -32 -48 -64.8 -243.2;
         -60 -121.6 -182.4 -243.2 -923.6]
  fac = -0.625
  A0 = A0*fac

  A02 = zeros(A0)
  EulerEquationMod.calcA0(params_e, q2, A02)
  for i=1:length(A0)
    @fact A0[i] --> roughly(A02[i], atol=1e-10)
  end

  A0inv = inv(A0)
  A02inv = zeros(A02)
  EulerEquationMod.calcA0Inv(params_e, q2, A02inv)
  for i=1:length(A0)
    @fact A0inv[i] --> roughly(A02inv[i], atol=1e-8)
  end



end
