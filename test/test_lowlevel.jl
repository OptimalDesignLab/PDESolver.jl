include("../src/solver/euler/startup.jl")  # initialization and construction

facts("--- Testing Euler Low Level Functions --- ") do

 q = [1.0, 2.0, 3.0, 7.0]
 aux_vars = [0.0]
 dir = [1.0, 0.0]
 F = zeros(4)
 coords = [1.0,  0.0]
 context("--- Testing calc functions ---") do

   @fact EulerEquationMod.calcPressure(q, eqn.params) => roughly(0.2)
   EulerEquationMod.calcEulerFlux(eqn.params, q, aux_vars, dir, F)
   @fact F => roughly([2.0, 4.2, 6, 14.4])

 end

 context("--- Testing common functions ---") do

   fill!(F, 0.0)
   EulerEquationMod.calcRho1Energy2(coords, eqn.params, F)
   @fact F[1] => 1.0
   @fact F[4] => 2.0

   fill!(F, 0.0)
   EulerEquationMod.calcRho1Energy2U3(coords, eqn.params, F)
   @fact F[1] => roughly(1.0)
   @fact F[2] => roughly(0.35355)
   @fact F[3] => roughly(0.35355)
   @fact F[4] => roughly(2.0)

   fill!(F, 0.0)
   EulerEquationMod.calcIsentropicVortex(coords, eqn.params, F)
   @fact F[1] => roughly(2.000)
   @fact F[2] => roughly(0.000)
   @fact F[3] => roughly(-1.3435)
   @fact F[4] => roughly(2.236960)




 end






 end
