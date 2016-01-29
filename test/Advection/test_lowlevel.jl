include("../../src/solver/advection/startup_advection.jl")  # initialization and construction
fill!(eqn.res_vec, 0.0)
using ArrayViews
facts("--- Testing Mesh --- ") do

  @fact mesh.numVert => 4
  @fact mesh.numEdge => 5
  @fact mesh.numEl => 2
  @fact mesh.order => 1
  @fact mesh.numDof => 4
  @fact mesh.numNodes => 4
  @fact mesh.numDofPerNode => 1
  @fact mesh.numBoundaryEdges => 4
  @fact mesh.numInterfaces => 1
  @fact mesh.numNodesPerElement => 3
  @fact mesh.numNodesPerType => [1, 0 , 0]

  println("mesh.bndryfaces = ", mesh.bndryfaces)
  @fact mesh.bndry_funcs[1] => AdvectionEquationMod.x5plusy5BC()
  @fact mesh.bndryfaces[1].element => 1
  @fact mesh.bndryfaces[1].face => 1
  @fact mesh.bndryfaces[2].element => 1
  @fact mesh.bndryfaces[2].face => 2
  @fact mesh.bndryfaces[3].element => 2
  @fact mesh.bndryfaces[3].face => 1
  @fact mesh.bndryfaces[4].element => 2
  @fact mesh.bndryfaces[4].face => 3

  println("mesh.interfaces = ",  mesh.interfaces)
  @fact mesh.interfaces[1].elementL => 2
  @fact mesh.interfaces[1].elementR => 1
  @fact mesh.interfaces[1].faceL => 2
  @fact mesh.interfaces[1].faceR => 3

  jac_fac = 0.25
  fac = 2
  println("mesh.coords = ", mesh.coords)
  @fact mesh.coords[:, :, 1] => roughly([4 4 0; 0 4 4.0])
  @fact mesh.coords[:, :, 2] => roughly([0.0 4 0; 0 0 4])

  println("mesh.dxidx = \n", mesh.dxidx)

  @fact mesh.dxidx[:, :, 1, 1] => roughly(fac*[1 1; -1 0.0], atol=1e-14)
  @fact mesh.dxidx[:, :, 2, 1] => roughly(fac*[1 1; -1 0.0], atol=1e-14)
  @fact mesh.dxidx[:, :, 3, 1] => roughly(fac*[1 1; -1 0.0], atol=1e-14)

  @fact mesh.dxidx[:, :, 1, 2] => roughly(fac*[1 0; 0 1.0], atol=1e-14)
  @fact mesh.dxidx[:, :, 2, 2] => roughly(fac*[1 0; 0 1.0], atol=1e-14)
  @fact mesh.dxidx[:, :, 3, 2] => roughly(fac*[1 0; 0 1.0], atol=1e-14)

  @fact mesh.jac => roughly(jac_fac*ones(3,2))

end


facts("--- Testing Functions Within AdvectionData_--- ") do
  Tsol = Float64
  u_vec = Tsol[1,2,3,4]
  u  = zeros(Tsol, 1, 3, 2)

  # checking disassembleSolution
  eqn.disassembleSolution(mesh, sbp, eqn, opts, u, u_vec)
  @fact u[1,1,1] => roughly(1.0)
  @fact u[1,2,1] => roughly(2.0)
  @fact u[1,3,1] => roughly(3.0)
  @fact u[1,1,2] => roughly(4.0)
  @fact u[1,2,2] => roughly(1.0)
  @fact u[1,3,2] => roughly(3.0)

  println("mesh.dofs = \n", mesh.dofs)
  #checking assembleSolution
  fill!(u_vec, 0.0)
  eqn.assembleSolution(mesh, sbp, eqn, opts, u, u_vec)
  @fact u_vec => roughly([2.0,2.0,6.0,4.0])

  # check mass matrix
  # just for testing, make jac != 1
  println("sbp.w = ", sbp.w)
  println("mesh.jac = ", mesh.jac)
  w_val = sbp.w[1]  # all of sbp.w entries are the same
  jac_val = 0.25
  M1 = 2*w_val/jac_val
  M2 = w_val/jac_val
  M3 = 2*w_val/jac_val
  M4 = w_val/jac_val
  M_test = [M1, M2, M3, M4]
  M_code = AdvectionEquationMod.calcMassMatrix(mesh, sbp, eqn)
  @fact M_code => roughly(M_test, atol=1e-13)

  Minv_test = 1./M_test
  Minv_test = AdvectionEquationMod.calcMassMatrixInverse(mesh, sbp, eqn)
  @fact Minv_test => roughly(Minv_test, atol=1e-13)
end

  context("--- Testing Boundary Function ---") do

 end



 context("--- Testing common functions ---") do




 end


 context("--- Testing dataPrep ---") do
 
  end


  context("--- Testing evalVolumeIntegrals ---")  do



  end


  context("--- Testing evalBoundaryIntegrals ---") do
  end

  context("--- Testing evalEuler --- ")  do

end # end facts block
