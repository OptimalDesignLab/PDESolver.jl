include("../src/solver/euler/startup.jl")  # initialization and construction
fill!(eqn.SL, 0.0)
using ArrayViews
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
 aux_vars = [EulerEquationMod.calcPressure(q, eqn.params)]
 dxidx = mesh.dxidx[:, :, 1, 1]  # arbitrary
 dir = [1.0, 0.0]
 F = zeros(4)
 coords = [1.0,  0.0]

 F_xi = zeros(4,2)

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

   # calc Euler fluxs needed by Roe solver
   F_roe = zeros(4)

   nrm1 = [dxidx[1,1], dxidx[1,2]]
   EulerEquationMod.calcEulerFlux(eqn.params, q, aux_vars, nrm1, view(F_xi, :, 1))
   nrm2 = [dxidx[2,1], dxidx[2,2]]
   EulerEquationMod.calcEulerFlux(eqn.params, q, aux_vars, nrm2, view(F_xi, :, 2))

   EulerEquationMod.RoeSolver(q, qg, F_xi, aux_vars, dxidx, dir, F_roe, eqn.params)
   @fact F_roe => roughly(-F) 


   # test that roe flux = euler flux of BC functions
   EulerEquationMod.calcIsentropicVortex(coords, eqn.params, q)

   func1 = EulerEquationMod.isentropicVortexBC()
   func1(q, F_xi, aux_vars, coords, dxidx, dir, F_roe, eqn.params)
   EulerEquationMod.calcEulerFlux(eqn.params, q, aux_vars, nrm, F)
   @fact F_roe => roughly(-F) 

   q[3] = 0  # make flow parallel to wall
   func1 = EulerEquationMod.noPenetrationBC()
   func1(q, F_xi, aux_vars, coords, dxidx, dir, F_roe, eqn.params)
   EulerEquationMod.calcEulerFlux(eqn.params, q, aux_vars, nrm, F)
   @fact F_roe => roughly(-F) 

   EulerEquationMod.calcRho1Energy2U3(coords, eqn.params, q)
   func1 = EulerEquationMod.Rho1E2U3BC()
   func1(q, F_xi, aux_vars, coords, dxidx, dir, F_roe, eqn.params)
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
   @fact F[1] => roughly(1.0, atol=1e-4)
   @fact F[2] => roughly(0.35355, atol=1e-4)
   @fact F[3] => roughly(0.35355, atol=1e-4)
   @fact F[4] => roughly(2.0, atol=1e-4)

   fill!(F, 0.0)
   EulerEquationMod.calcIsentropicVortex(coords, eqn.params, F)
   @fact F[1] => roughly(2.000, atol=1e-4)
   @fact F[2] => roughly(0.000, atol=1e-4)
   @fact F[3] => roughly(-1.3435, atol=1e-4)
   @fact F[4] => roughly(2.236960, atol=1e-4)




 end


 context("--- Testing dataPrep ---") do
 
   EulerEquationMod.dataPrep(mesh, sbp, eqn, eqn.SL0, opts)


   # test disassembleSolutino
   for i=1:mesh.numEl
     for j=1:mesh.numNodesPerElement
       @fact eqn.q[:, j, i] => roughly([1.0, 0.35355, 0.35355, 2.0], atol=1e-5)
     end
   end

   #=
   for i=1:mesh.numEl
     println("i = ", i)
     for j=1:mesh.numNodesPerElement
       println("j = ", j)
       aux_vars_i = eqn.aux_vars[ :, j, i]
       println("aux_vars_i = ", aux_vars_i)
       p = EulerEquationMod.@getPressure(aux_vars_i)
       @fact p => roughly(0.750001, atol=1e-5)
     end
   end
   =#

   # test calcEulerFlux
   for i=1:mesh.numNodesPerElement
     println("eq.F_xi[:, $i, 1, 1] = ", eqn.F_xi[:, i, 1, 1])
     @fact eqn.F_xi[:, i, 1, 1] => roughly([0.0, 0.750001, -0.750001, 0.0], atol=1e-5)
   end

   for i=1:mesh.numNodesPerElement
     @fact eqn.F_xi[:, i, 2, 1] => roughly([0.35355, 0.874999, 0.12499, 0.972263], atol=1e-5)
   end

   for i=1:mesh.numNodesPerElement
     @fact eqn.F_xi[:, i, 1, 2] => roughly([0.35355, 0.124998, 0.874999, 0.972263], atol=1e-5)
   end

   for i=1:mesh.numNodesPerElement
     @fact eqn.F_xi[:, i, 2, 2] => roughly([0.0, -0.750001, 0.750001, 0.0], atol=1e-5)
   end


   # test getBCFluxes
     for j= 1:sbp.numfacenodes
       @fact eqn.bndryflux[:, j, 1] => roughly([-0.35355, -0.874999, -0.124998, -0.972263], atol=1e-5)
     end

     for j= 1:sbp.numfacenodes
       @fact eqn.bndryflux[:, j, 2] => roughly([-0.35355,  -0.124998, -0.874999, -0.972263], atol=1e-5)
     end

     for j= 1:sbp.numfacenodes
       @fact eqn.bndryflux[:, j, 3] => roughly([0.35355,  0.124998, 0.874999, 0.972263], atol=1e-5)
     end

     for j= 1:sbp.numfacenodes
       @fact eqn.bndryflux[:, j, 4] => roughly([0.35355, 0.874999, 0.124998, 0.972263], atol=1e-5)
     end



  end


  context("--- Testing evalVolumeIntegrals ---")  do

    EulerEquationMod.evalVolumeIntegrals(mesh, sbp, eqn)

    el1_res = [-0.35355  0  0.35355;
                -0.874999  0.750001  0.124998;
		-0.124998  -0.750001  0.874999;
		-0.972263  0  0.972263]
    el2_res = [-0.35355  0.35355 0;
                -0.124998  0.874999 -0.75001;
		-0.874999 0.124998 0.75001;
		-0.972263  0.972263 0]
 
    @fact eqn.res[:, :, 1] => roughly(el1_res, atol=1e-5)
    @fact eqn.res[:, :, 2] => roughly(el2_res, atol=1e-5)



  end


  context("--- Testing evalBoundaryIntegrals ---") do
    fill!(eqn.res, 0.0)

    EulerEquationMod.evalBoundaryIntegrals( mesh, sbp, eqn)

    el1_res = [0.35355 0 -0.35355;
               0.124998 -0.750001 -0.874999;
	       0.874999 0.750001 -0.124998;
	       0.972263  0  -0.972263]
    el2_res = [0.35355 -0.35355 0;
               0.874999 -0.124998 0.750001;
	       0.124998  -0.874999  -0.750001;
	       0.972263  -0.972263  0]

    @fact eqn.res[:, :, 1] => roughly(el1_res, atol=1e-5)
    @fact eqn.res[:, :, 2] => roughly(el2_res, atol=1e-5)

  end

  context("--- Testing evalEuler --- ")  do

    fill!(eqn.SL, 0.0)
    fill!(eqn.res, 0.0)
    EulerEquationMod.evalEuler(mesh, sbp, eqn, opts, eqn.SL0, eqn.SL)

    for i=1:mesh.numDof
      eqn.SL[i] => roughly(0.0, atol=1e-14)
    end

  end

end # end facts block
