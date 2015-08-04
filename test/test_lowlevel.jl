include("../src/solver/euler/startup.jl")  # initialization and construction

facts("--- Testing Mesh --- ") do

  @fact mesh.numVert => 4
  @fact mesh.numEdge => 5
  @fact mesh.numEl => 2
  @fact mesh.order => 1
  @fact mesh.numDof => 16
  @fact mesh.numNodes => 4
  @fact mesh.numDofPerNode => 4
  @fact mesh.numBoundaryEdges => 4
  @fact mesh.numInterfaces => 1
  @fact mesh.numNodesPerElement => 3
  @fact mesh.numNodesPerType => [1, 0 , 0]

  @fact mesh.bndry_funcs[1] => EulerEquationMod.Rho1E2U3BC()
  @fact mesh.bndryfaces[1].element => 1
  @fact mesh.bndryfaces[1].face => 2
  @fact mesh.bndryfaces[2].element => 2
  @fact mesh.bndryfaces[2].face => 2
  @fact mesh.bndryfaces[3].element => 1
  @fact mesh.bndryfaces[3].face => 1
  @fact mesh.bndryfaces[4].element => 2
  @fact mesh.bndryfaces[4].face => 3

  @fact mesh.interfaces[1].elementL => 2
  @fact mesh.interfaces[1].elementR => 1
  @fact mesh.interfaces[1].faceL => 1
  @fact mesh.interfaces[1].faceR => 3

  @fact mesh.coords[:, :, 1] => roughly([-1.0 1 1; -1 -1 1])
  @fact mesh.coords[:, :, 2] => roughly([-1.0 1 -1; -1 1 1])

  @fact mesh.dxidx[:, :, 1, 1] => roughly([1.0 -1; 0 1], atol=1e-14)

  @fact mesh.dxidx[:, :, 1, 1] => roughly([1.0 -1; 0 1], atol=1e-14)
  @fact mesh.dxidx[:, :, 2, 1] => roughly([1.0 -1; 0 1], atol=1e-14)
  @fact mesh.dxidx[:, :, 3, 1] => roughly([1.0 -1; 0 1], atol=1e-14)

  @fact mesh.dxidx[:, :, 1, 2] => roughly([1.0 0; -1 1], atol=1e-14)
  @fact mesh.dxidx[:, :, 2, 2] => roughly([1.0 0; -1 1], atol=1e-14)
  @fact mesh.dxidx[:, :, 3, 2] => roughly([1.0 0; -1 1], atol=1e-14)

  @fact mesh.jac => roughly(ones(3,2))


end




facts("--- Testing Euler Low Level Functions --- ") do

 q = [1.0, 2.0, 3.0, 7.0]
 qg = deepcopy(q)
 aux_vars = [0.0]
 dxidx = mesh.dxidx[:, :, 1, 1]  # arbitrary
 dir = [1.0, 0.0]
 F = zeros(4)
 coords = [1.0,  0.0]


 context("--- Testing calc functions ---") do

   @fact EulerEquationMod.calcPressure(q, eqn.params) => roughly(0.2)
   EulerEquationMod.calcEulerFlux(eqn.params, q, aux_vars, dir, F)
   @fact F => roughly([2.0, 4.2, 6, 14.4], atol=1e-14)

 end

  context("--- Testing Boundary Function ---") do
  
   nx = dxidx[1,1]*dir[1] + dxidx[2,1]*dir[2]
   ny = dxidx[1,2]*dir[1] + dxidx[2,2]*dir[2]
   nrm = [nx, ny]

   EulerEquationMod.calcEulerFlux(eqn.params, q, aux_vars, nrm, F)

   F_roe = zeros(4)
   EulerEquationMod.RoeSolver(q, qg, aux_vars, dxidx, dir, F_roe, eqn.params)
   @fact F_roe => roughly(-F) 


   # test that roe flux = euler flux of BC functions
   EulerEquationMod.calcIsentropicVortex(coords, eqn.params, q)

   func1 = EulerEquationMod.isentropicVortexBC()
   func1(q, aux_vars, coords, dxidx, dir, F_roe, eqn.params)
   EulerEquationMod.calcEulerFlux(eqn.params, q, aux_vars, nrm, F)
   @fact F_roe => roughly(-F) 

   q[3] = 0  # make flow parallel to wall
   func1 = EulerEquationMod.noPenetrationBC()
   func1(q, aux_vars, coords, dxidx, dir, F_roe, eqn.params)
   EulerEquationMod.calcEulerFlux(eqn.params, q, aux_vars, nrm, F)
   @fact F_roe => roughly(-F) 

   EulerEquationMod.calcRho1Energy2U3(coords, eqn.params, q)
   func1 = EulerEquationMod.Rho1E2U3BC()
   func1(q, aux_vars, coords, dxidx, dir, F_roe, eqn.params)
   EulerEquationMod.calcEulerFlux(eqn.params, q, aux_vars, nrm, F)
   @fact F_roe => roughly(-F) 



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
