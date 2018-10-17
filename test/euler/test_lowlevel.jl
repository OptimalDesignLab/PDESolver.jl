#using ArrayViews
# using input_vals_channel.jl

"""
  This function tests whether the mesh counted the number of mesh entities
  is correct.  It also checks the coordinates and mapping jacobian, which
  means it is dependent on the element ordering
"""
function test_lowlevel_mesh(mesh, sbp, eqn, opts)
  @testset "--- Testing Mesh --- " begin

    @test ( mesh.numVert )== 4
    @test ( mesh.numEdge )== 5
    @test ( mesh.numEl )== 2
    @test ( mesh.order )== 1
    @test ( mesh.numDof )== 16
    @test ( mesh.numNodes )== 4
    @test ( mesh.numDofPerNode )== 4
    @test ( mesh.numBoundaryFaces )== 4
    @test ( mesh.numInterfaces )== 1
    @test ( mesh.numNodesPerElement )== 3
    @test ( mesh.numNodesPerType )== [1, 0 , 0]

    @test ( mesh.bndry_funcs[1] )== EulerEquationMod.Rho1E2U3BC(mesh, eqn)
    @test ( mesh.bndryfaces[1].element )== 1
    @test ( mesh.bndryfaces[1].face )== 3
    @test ( mesh.bndryfaces[2].element )== 2
    @test ( mesh.bndryfaces[2].face )== 1
    @test ( mesh.bndryfaces[3].element )== 1
    @test ( mesh.bndryfaces[3].face )== 2
    @test ( mesh.bndryfaces[4].element )== 2
    @test ( mesh.bndryfaces[4].face )== 2

    @test ( mesh.interfaces[1].elementL )== 1
    @test ( mesh.interfaces[1].elementR )== 2
    @test ( mesh.interfaces[1].faceL )== 1
    @test ( mesh.interfaces[1].faceR )== 3


  #=
    @test ( mesh.bndryfaces[1].element )== 1
    @test ( mesh.bndryfaces[1].face )== 2
    @test ( mesh.bndryfaces[2].element )== 2
    @test ( mesh.bndryfaces[2].face )== 2
    @test ( mesh.bndryfaces[3].element )== 1
    @test ( mesh.bndryfaces[3].face )== 1
    @test ( mesh.bndryfaces[4].element )== 2
    @test ( mesh.bndryfaces[4].face )== 3

    @test ( mesh.interfaces[1].elementL )== 2
    @test ( mesh.interfaces[1].elementR )== 1
    @test ( mesh.interfaces[1].faceL )== 1
    @test ( mesh.interfaces[1].faceR )== 3
  =#
    @test isapprox( mesh.coords[:, :, 2], [-1.0 1 1; -1 -1 1]) 
    @test isapprox( mesh.coords[:, :, 1], [-1.0 1 -1; -1 1 1]) 

    @test isapprox( mesh.dxidx[:, :, 1, 2], [1.0 -1; 0 1]) atol=1e-14

    @test isapprox( mesh.dxidx[:, :, 1, 2], [1.0 -1; 0 1]) atol=1e-14
    @test isapprox( mesh.dxidx[:, :, 2, 2], [1.0 -1; 0 1]) atol=1e-14
    @test isapprox( mesh.dxidx[:, :, 3, 2], [1.0 -1; 0 1]) atol=1e-14

    @test isapprox( mesh.dxidx[:, :, 1, 1], [1.0 0; -1 1]) atol=1e-14
    @test isapprox( mesh.dxidx[:, :, 2, 1], [1.0 0; -1 1]) atol=1e-14
    @test isapprox( mesh.dxidx[:, :, 3, 1], [1.0 0; -1 1]) atol=1e-14

    @test isapprox( mesh.jac, ones(3,2)) 


  end
end # end function

#test_lowlevel_mesh(mesh, sbp, eqn, opts)
add_func2!(EulerTests, test_lowlevel_mesh, "input_vals_channel.jl", [TAG_SHORTTEST])

"""
  This function tests calculation of quantities needed for using the
  entropy variables: conversion, A0, A1

"""
function test_lowlevel_entropyvars(mesh, sbp, eqn, opts)

  Tsol = eltype(eqn.q)
  Tres = eltype(eqn.res)
  Tmsh = eltype(mesh.dxidx)
  @testset "--- Testing Euler Low Level Functions --- " begin
    opts["variable_type"] = :entropy
    eqn_e = EulerData_{Tsol, Tres, 2, Tmsh, opts["variable_type"]}(mesh, sbp, opts)

    e_params = eqn_e.params
    opts["variable_type"] = :conservative

    q = [1.0, 2.0, 3.0, 7.0]
    qg = deepcopy(q)
    aux_vars = [EulerEquationMod.calcPressure(eqn.params, q)]
    dxidx = mesh.dxidx[:, :, 1, 1]  # arbitrary
    dir = [1.0, 0.0]
    F = zeros(4)
    Fe = zeros(4)
    coords = [1.0,  0.0]

    flux_parametric = zeros(4,2)

    v = zeros(4)
    EulerEquationMod.convertToEntropy(eqn.params, q, v)
    v_analytic = [-2*4.99528104378295, 4., 6, -2*1]
    @test isapprox( v, v_analytic) 
    # test inplace operation
    q2 = copy(q)
    EulerEquationMod.convertToEntropy(eqn.params, q2, q2)
    @test ( q2 )== v_analytic
    q_ret = zeros(4)
    EulerEquationMod.convertToConservative(e_params, v, q_ret)
    @test isapprox( q_ret, q) 
    
    # test inplace operation
    v2 = copy(v)
    EulerEquationMod.convertToConservative(e_params, v2, v2)
    @test isapprox( v2, q) 

    v = zeros(4)
    vIR = zeros(4)
    EulerEquationMod.convertToEntropy(eqn.params, q, v)
    v_analytic = [-2*4.99528104378295, 4., 6, -2*1]
    @test isapprox( v, v_analytic) 
    # test inplace operation
    q2 = copy(q)
    EulerEquationMod.convertToEntropy(eqn.params, q2, q2)
    @test ( q2 )== v_analytic

    # test IR variables
    EulerEquationMod.convertToIR(eqn.params, q, vIR)
    diff = vIR - q2./eqn.params.gamma_1
    @test isapprox( norm(diff), 0.0) atol=1e-12

    # test inplace operation
    vIR2 = copy(q)
    EulerEquationMod.convertToIR(eqn.params, vIR2, vIR2)
    diff = vIR2 - q2./eqn.params.gamma_1
    @test isapprox( norm(diff), 0.0) atol=1e-12

    # convert back
    EulerEquationMod.convertToConservativeFromIR_(eqn.params, vIR2, vIR2)
    diff = q - vIR2
    @test isapprox( norm(diff), 0.0) atol=1e-12

    q_ret = zeros(4)
    EulerEquationMod.convertToConservative(e_params, v, q_ret)
    @test isapprox( q_ret, q) 
    
    # test inplace operation
    v2 = copy(v)
    EulerEquationMod.convertToConservative(e_params, v2, v2)
    @test isapprox( v2, q) 


    # ------------------------------------------------------------------------
    # test matrices

    # test inv(A0)
    A0inv = zeros(4,4)
    A0inv2 = [170.4 -52 -78 24; 
              -52 18 24 -8; 
              -78 24 38 -12; 
              24 -8  -12 4]
    EulerEquationMod.calcA0Inv(e_params, v, A0inv)

    @test isapprox( A0inv, A0inv2) 

    # test A0
    A0 = zeros(4,4)
    A02 = inv(A0inv)
    EulerEquationMod.calcA0(e_params, v, A0)
   
    for i=1:16
      @test isapprox( A0[i], A02[i]) atol=1e-10
    end


    A0inv_c = zeros(4,4)
    EulerEquationMod.calcA0(eqn.params, q, A0inv_c)
    @test ( A0inv_c )== eye(4)

    A0_c = zeros(4,4)
    EulerEquationMod.calcA0Inv(eqn.params, q, A0_c)
    @test ( A0_c )== eye(4)

   # test calcIRA0
   A03 = A0 *eqn.params.gamma_1
   A04_test = zeros(A02)
   A04_code = zeros(A02)
   EulerEquationMod.getIRA0(eqn.params, q, A04_code)
   q3 = zeros(Complex128, length(q))
   q4 = zeros(Complex128, length(q))

   q3[:] = q
   EulerEquationMod.convertToEntropy_(eqn.params, q3, q4)
   for i=1:length(q)
     h = 1e-20
     q4[i] += complex(0, h)
     EulerEquationMod.convertToConservative_(eqn.params, q4, q3)
     A04_test[:, i] = imag(q3)/h
     q4[i] -= complex(0, h)
   end
   scale!(A04_test, eqn.params.gamma_1)

   for i=1:length(A04_test)
     @test isapprox( A04_code[i], A04_test[i]) atol=1e-13
     @test isapprox( A04_code[i], A03[i]) atol=1e-10
   end

   test_IRA0inv(eqn.params, q)

      # test A1
    A1 = zeros(4,4)
    q_tmp = ones(Tsol, 4)
    EulerEquationMod.calcA1(eqn.params, q_tmp, A1)
    @test isapprox(A1,    [0.0 1.0 0.0  0.0
                          -0.6 1.6 -0.4 0.4
                          -1.0 1.0 1.0 0.0
                          -0.6 0.6 -0.4 1.4]) atol=1e-13


    EulerEquationMod.calcA1(e_params, v, A1)
    fac = 0.3125
    A1_analytic = fac*[16 33.6 48 115.2;
                            33.6 73.6 100.8 248.32; 
                            48 100.8 147.2 355.2;
                            115.2 248.32 355.2 4*218.32]

    A1_diff = A1 - A1_analytic
    for i=1:16
      @test isapprox( A1[i], A1_analytic[i]) atol=1e-10
    end

    A2 = zeros(4,4)
    A2 = zeros(4,4)
    EulerEquationMod.calcA2(eqn.params, q_tmp, A2)
    @test isapprox(A2,  [0.0 0.0 1.0 0.0
                        -1.0 1.0 1.0 0.0
                        -0.6 -0.4 1.6 0.4
                        -0.6 -0.4 0.6 1.4]) atol=1e-13



    EulerEquationMod.calcA2(e_params, v, A2)
    A2_analytic = fac*[24. 48 73.6 172.8;
                            48 100.8 147.2 355.2; 
                            73.6 147.2 230.4 544.32;
                            172.8 355.2 544.32 1309.92]
    A2_diff = A2 - A2_analytic

    for i=1:16
      @test isapprox( A2[i], A2_analytic[i]) atol=1e-10
    end


    # check that checkDensity and checkPresure work
    @test_throws Exception  EulerEquationMod.checkDensity(eqn, mesh)
    @test_throws Exception  EulerEquationMod.checkPressure(eqn, mesh)

  end  # end facts block

  return nothing
end  # end function

function test_IRA0inv(params, q)

  h = 1e-20
  pert = Complex128(0, h)

  # test calcIRA0inv
  q3 = zeros(Complex128, length(q))
  copy!(q3, real(q))
  
  numDofPerNode = length(q)

  qe = zeros(q3)
  A0inv = zeros(Complex128, numDofPerNode, numDofPerNode)
  A0invc = zeros(Complex128, numDofPerNode, numDofPerNode)

  for i=1:numDofPerNode
    q3[i] += pert
    EulerEquationMod.convertToIR_(params, q3, qe)
    q3[i] -= pert
    A0invc[:, i] = imag(qe)/h
  end

  EulerEquationMod.getIRA0inv(params, q3, A0inv)

  @test maximum(abs.(A0inv - A0invc)) < 1e-12

  # test against the inverse of A0
  A0 = zeros(A0inv)
  EulerEquationMod.getIRA0(params, q3, A0)
  A0inv2 = inv(A0)
  @assert maximum(abs.(A0inv - A0inv2)) < 1e-9

  return nothing
end



#test_lowlevel_entropyvars(mesh, sbp, eqn, opts)
add_func2!(EulerTests, test_lowlevel_entropyvars, "input_vals_channel.jl", [TAG_ENTROPYVARS, TAG_SHORTTEST])

"""
  This function tests functions that are used to convert entire vectors/arrays
  of values between conservative and entropy variables
"""
function test_lowlevel_convert(mesh, sbp, eqn, opts)
  Tsol = eltype(eqn.q)
  Tres = eltype(eqn.res)
  Tmsh = eltype(mesh.dxidx)
  opts["variable_type"] = :entropy
  eqn_e = EulerData_{Tsol, Tres, 2, Tmsh, opts["variable_type"]}(mesh, sbp, opts)

  e_params = eqn_e.params
  opts["variable_type"] = :conservative
  #=
  q = [1.0, 2.0, 3.0, 7.0]
  qg = deepcopy(q)
  aux_vars = [EulerEquationMod.calcPressure(eqn.params, q)]
  dxidx = mesh.dxidx[:, :, 1, 1]  # arbitrary
  dir = [1.0, 0.0]
  F = zeros(4)
  Fe = zeros(4)
  coords = [1.0,  0.0]
  =#

  A0 = zeros(4,4)
  A0inv = zeros(4,4)

  @testset "--- Testing convert Functions ---" begin
    # for the case, the solution is uniform flow
    array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)
    v_arr = copy(eqn.q)
    v2 = zeros(4)

    EulerEquationMod.convertToEntropy(eqn.params, eqn.q[:, 1, 1], v2)
    EulerEquationMod.convertToEntropy(mesh, sbp, eqn, opts, v_arr)
    # test conversion to entropy variables
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        @test ( v_arr[:, j, i] )== v2
      end
    end

    eqn_e.q = v_arr # attach entropy variables to eqn_e

    v_vec = copy(eqn.q_vec)
    EulerEquationMod.convertToEntropy(mesh, sbp, eqn, opts, v_vec)
    for i=1:4:mesh.numDof
      @test ( v_vec[i:(i+3)] )== v2
    end

    eqn_e.q_vec = v_vec
    v_arr2 = copy(v_arr)

    # test multiplying an entire array by A0inv
    EulerEquationMod.calcA0Inv(e_params, v_arr2[:, 1, 1], A0inv)
    v2 = A0inv*v_arr2[:, 1, 1]
    EulerEquationMod.matVecA0inv(mesh, sbp, eqn_e, opts, v_arr2)
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        @test isapprox( v_arr2[:, j, i], v2) 
      end
    end

    v_arr3 = copy(v_arr)

    EulerEquationMod.calcA0(e_params, v_arr3[:, 1, 1], A0)
    v3 = A0*v_arr3[:, 1, 1]  # store original values
    EulerEquationMod.matVecA0(mesh, sbp, eqn_e, opts, v_arr3)
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        @test isapprox( v_arr3[:, j, i], v3) 
      end
    end

    # now test converting back to conservative
    EulerEquationMod.convertToConservative(mesh, sbp, eqn_e, opts, v_arr)
    for i =1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        @test isapprox( v_arr[:, j, i], eqn.q[:, j, i]) 
      end
    end

    EulerEquationMod.convertToConservative(mesh, sbp, eqn_e, opts, v_vec)
    for i=1:mesh.numDof
      @test isapprox( v_vec[i], eqn.q_vec[i]) 
    end



  end  # end facts block

  return nothing
end
     
#test_lowlevel_convert(mesh, sbp, eqn, opts)
add_func2!(EulerTests, test_lowlevel_convert, "input_vals_channel.jl", [TAG_ENTROPYVARS, TAG_SHORTTEST])



function test_lowlevel_eigsystem(mesh, sbp, eqn, opts)


  @testset "--- Testing Eigensystem ---" begin
    # compute flux jacobian in x and y directions
    params = eqn.params
    q = [1.0, 2.0, 3.0, 7.0]
    qc = convert(Array{Complex128}, q)
    aux_vars = Array{Complex128}(1)
    qg = deepcopy(q)
    dxidx = mesh.dxidx[:, :, 1, 1]  # arbitrary
    F = zeros(Complex128, 4)

    Ax = zeros(4,4)
    Ax3 = zeros(Ax)
    Ay = zeros(4,4)
    Ay3 = zeros(Ay)

    EulerEquationMod.calcA(params, q, Ax3)
    EulerEquationMod.calcB(params, q, Ay3)
    h = 1e-20
    pert = Complex128(0, h)
    for i=1:4
      qc[i] += pert
      dir = [1.0, 0.0]
      p = EulerEquationMod.calcPressure(params, qc)
      aux_vars[1] = p
      EulerEquationMod.calcEulerFlux(params, qc, aux_vars, dir, F)
      Ax[:, i] = imag(F)/h

      dir = [0.0, 1.0]
      EulerEquationMod.calcEulerFlux(params, qc, aux_vars, dir, F)
      Ay[:, i] = imag(F)/h

      qc[i] -= pert
    end

    # now compute Ax and and Ay from their eigensystem and compare
    Yx = zeros(4,4)
    Yy = zeros(4,4)
    Lambdax = zeros(4)
    Lambday = zeros(4)

    EulerEquationMod.calcEvecsx(params, q, Yx)
    EulerEquationMod.calcEvecsy(params, q, Yy)
    EulerEquationMod.calcEvalsx(params, q, Lambdax)
    EulerEquationMod.calcEvalsy(params, q, Lambday)

    Ax2 = Yx*diagm(Lambdax)*inv(Yx)
    Ay2 = Yy*diagm(Lambday)*inv(Yy)

    @test isapprox( Ax2, Ax) atol=1e-12
    @test isapprox( Ay2, Ay) atol=1e-12
    @test isapprox( Ax3, Ax) atol=1e-12
    @test isapprox( Ay3, Ay) atol=1e-12

    # check A0 = Y*S2*Y.'
    A0 = zeros(4,4)
    Sx = zeros(4)
    Sy = zeros(4)

    EulerEquationMod.getIRA0(params, q, A0)
    EulerEquationMod.calcEScalingx(params, q, Sx)
    EulerEquationMod.calcEScalingy(params, q, Sy)
    A02 = Yx*diagm(Sx)*Yx.'
    A03 = Yy*diagm(Sy)*Yy.'

    @test isapprox( A0, A02) atol=1e-12
    @test isapprox( A0, A03) atol=1e-12

  end  # end facts block

  return nothing
end  # end function

#test_lowlevel_eigsystem()

add_func2!(EulerTests, test_lowlevel_eigsystem, "input_vals_channel.jl", [TAG_SHORTTEST])


"""
  This function tests calculation of dependent quantities: pressure, speed of
  sound, Euler flux.
"""
function test_lowlevel_calc(mesh, sbp, eqn, opts)
    Tsol = eltype(eqn.q)
    Tres = eltype(eqn.res)
    Tmsh = eltype(mesh.dxidx)

   opts["variable_type"] = :entropy
   eqn_e = EulerData_{Tsol, Tres, 2, Tmsh, opts["variable_type"]}(mesh, sbp, opts)

   e_params = eqn_e.params
   opts["variable_type"] = :conservative
   
   q = [1.0, 2.0, 3.0, 7.0]
   v = zeros(q)
   EulerEquationMod.convertToEntropy(eqn.params, q, v)


  @testset "--- Testing calc functions ---" begin

    @test isapprox( EulerEquationMod.calcPressure(eqn.params, q), 0.2) 
    @test isapprox( EulerEquationMod.calcPressure(e_params, v), 0.2) 
    a_cons = EulerEquationMod.calcSpeedofSound(eqn.params, q)
    a_ent = EulerEquationMod.calcSpeedofSound(e_params, v)
    @test isapprox( a_cons, a_ent) 
    EulerEquationMod.calcEulerFlux(eqn.params, q, aux_vars, dir, F)
    EulerEquationMod.calcEulerFlux(e_params, v, aux_vars, dir, Fe)
    @test isapprox( F, [2.0, 4.2, 6, 14.4]) atol=1e-14
    @test isapprox( Fe, F) 
  end

  return nothing
end

#test_lowlevel_convert(mesh, sbp, eqn, opts)
add_func2!(EulerTests, test_lowlevel_convert, "input_vals_channel.jl", [TAG_ENTROPYVARS, TAG_SHORTTEST])

"""
  Ths function tests calculation of boundary conditions and the Roe solver.
"""
function test_lowlevel_boundary(mesh, sbp, eqn, opts)

  q = [1.0, 2.0, 3.0, 7.0]
  qg = deepcopy(q)
  aux_vars = [EulerEquationMod.calcPressure(eqn.params, q)]
  dxidx = mesh.dxidx[:, :, 1, 1]  # arbitrary
  dir = [1.0, 0.0]
  F = zeros(4)
  Fe = zeros(4)
  coords = [1.0,  0.0]
  nrm_xy = zeros(mesh.dim)


  flux_parametric = zeros(4,2)

  v = zeros(4)
  EulerEquationMod.convertToEntropy(eqn.params, q, v)

  @testset "--- Testing Boundary Function ---" begin

    nx = dxidx[1,1]*dir[1] + dxidx[2,1]*dir[2]
    ny = dxidx[1,2]*dir[1] + dxidx[2,2]*dir[2]
    nrm = [nx, ny]

    EulerEquationMod.calcEulerFlux(eqn.params, q, aux_vars, nrm, F)

    # calc Euler fluxs needed by Roe solver
    F_roe = zeros(4)

    nrm1 = [dxidx[1,1], dxidx[1,2]]
    EulerEquationMod.calcEulerFlux(eqn.params, q, aux_vars, nrm1, sview(flux_parametric, :, 1))
    nrm2 = [dxidx[2,1], dxidx[2,2]]
    EulerEquationMod.calcEulerFlux(eqn.params, q, aux_vars, nrm2, sview(flux_parametric, :, 2))

    EulerEquationMod.RoeSolver(eqn.params, q, qg, aux_vars, nrm, F_roe)
    @test isapprox( F_roe, F) 


    # test that roe flux = euler flux of BC functions
    EulerEquationMod.calcIsentropicVortex(eqn.params, coords, q)

    nrm1 = [dxidx[1,1], dxidx[1,2]]
    EulerEquationMod.calcEulerFlux(eqn.params, q, aux_vars, nrm1, sview(flux_parametric, :, 1))
    nrm2 = [dxidx[2,1], dxidx[2,2]]
    EulerEquationMod.calcEulerFlux(eqn.params, q, aux_vars, nrm2, sview(flux_parametric, :, 2))


    func1 = EulerEquationMod.isentropicVortexBC(mesh, eqn)
    EulerEquationMod.calcEulerFlux(eqn.params, q, aux_vars, nrm, F)
    calcBCNormal(eqn.params, dxidx, dir, nrm_xy)
    func1(eqn.params, q, aux_vars, coords, nrm_xy, F_roe,)

    @test isapprox( F_roe, F) 

    q[3] = 0  # make flow parallel to wall
    func1 = EulerEquationMod.noPenetrationBC(mesh, eqn)
    nrm1 = [dxidx[1,1], dxidx[1,2]]
    EulerEquationMod.calcEulerFlux(eqn.params, q, aux_vars, nrm1, sview(flux_parametric, :, 1))
    nrm2 = [dxidx[2,1], dxidx[2,2]]
    EulerEquationMod.calcEulerFlux(eqn.params, q, aux_vars, nrm2, sview(flux_parametric, :, 2))
   

    EulerEquationMod.calcEulerFlux(eqn.params, q, aux_vars, nrm, F)
    func1(eqn.params, q, aux_vars, coords, nrm_xy, F_roe)

    @test isapprox( F_roe, F) 

    # test the entropy stable BC
    func2 = EulerEquationMod.noPenetrationESBC(mesh, eqn)
    func2(eqn.params, q, aux_vars, coords, nrm_xy, F_roe)

    @test isapprox( F_roe, F) 

    EulerEquationMod.calcRho1Energy2U3(eqn.params, coords, q)
    func1 = EulerEquationMod.Rho1E2U3BC(mesh, eqn)
    nrm1 = [dxidx[1,1], dxidx[1,2]]
    EulerEquationMod.calcEulerFlux(eqn.params, q, aux_vars, nrm1, sview(flux_parametric, :, 1))
    nrm2 = [dxidx[2,1], dxidx[2,2]]
    EulerEquationMod.calcEulerFlux(eqn.params, q, aux_vars, nrm2, sview(flux_parametric, :, 2))


    EulerEquationMod.calcEulerFlux(eqn.params, q, aux_vars, nrm, F)
    func1(eqn.params, q, aux_vars, coords, nrm_xy, F_roe)

    @test isapprox( F_roe, F) 
  end  # end facts block

  return nothing
end

#test_lowlevel_boundary(mesh, sbp, eqn, opts)
add_func2!(EulerTests, test_lowlevel_boundary, "input_vals_channel.jl", [TAG_BC, TAG_FLUX, TAG_SHORTTEST])


"""
  This function tests things in commonfuncs.jl
"""
function test_lowlevel_commonfuncs(mesh, sbp, eqn, opts)

  coords = [1.0,  0.0]

  @testset "--- Testing common functions ---" begin

    F = zeros(4)
    EulerEquationMod.calcRho1Energy2(eqn.params, coords, F)
    @test ( F[1] )== 1.0
    @test ( F[4] )== 2.0

    fill!(F, 0.0)
    EulerEquationMod.calcRho1Energy2U3(eqn.params, coords, F)
    @test isapprox( F[1], 1.0) atol=1e-4
    @test isapprox( F[2], 0.35355) atol=1e-4
    @test isapprox( F[3], 0.35355) atol=1e-4
    @test isapprox( F[4], 2.0) atol=1e-4

    fill!(F, 0.0)
    EulerEquationMod.calcIsentropicVortex(eqn.params, coords, F)
    @test isapprox( F[1], 2.000) atol=1e-4
    @test isapprox( F[2], 0.000) atol=1e-4
    @test isapprox( F[3], -1.3435) atol=1e-4
    @test isapprox( F[4], 2.236960) atol=1e-4


    level = EulerEquationMod.getPascalLevel(1)
    @test ( level )== 1

    for i=2:3
      level = EulerEquationMod.getPascalLevel(i)
      @test ( level )== 2
    end

    for i=4:6
      level = EulerEquationMod.getPascalLevel(i)
      @test ( level )== 3
    end

    for i=7:10
      level = EulerEquationMod.getPascalLevel(i)
      @test ( level )== 4
    end

    for i=11:15
      level = EulerEquationMod.getPascalLevel(i)
      @test ( level )== 5
    end
  end  # end facts block
  
  return nothing
end

#test_lowlevel_commonfuncs(mesh, sbp, eqn, opts)
add_func2!(EulerTests, test_lowlevel_commonfuncs, "input_vals_channel.jl", [TAG_SHORTTEST])

"""
  This functinon tests dataprep and some auxiliary functions
"""
function test_lowlevel_dataprep(mesh, sbp, eqn, opts)

  @testset "--- Testing dataPrep ---" begin
 
    # the input file loaded a uniform flow into q_vev
    EulerEquationMod.array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)
    EulerEquationMod.dataPrep(mesh, sbp, eqn, opts)
    # test array1DTo3D
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        @test isapprox( eqn.q[:, j, i], [1.0, 0.35355, 0.35355, 2.0]) atol=1e-5
      end
    end

    # testing arrToVecAssign
    q_vec_orig = copy(eqn.q_vec)
    EulerEquationMod.arrToVecAssign(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    for i = 1:mesh.numDof
      @test isapprox( eqn.q_vec[i], q_vec_orig[i]) atol=1e-5
    end


    #=
    for i=1:mesh.numEl
      println("i = ", i)
      for j=1:mesh.numNodesPerElement
        println("j = ", j)
        aux_vars_i = eqn.aux_vars[ :, j, i]
        println("aux_vars_i = ", aux_vars_i)
        p = EulerEquationMod.@getPressure(aux_vars_i)
        @test isapprox( p, 0.750001) atol=1e-5
      end
    end
    =#

    # test calcEulerFlux
    for i=1:mesh.numNodesPerElement
      @test isapprox( eqn.flux_parametric[:, i, 1, 2], [0.0, -0.750001, 0.750001, 0.0]) atol=1e-5
    end

    for i=1:mesh.numNodesPerElement
      @test isapprox( eqn.flux_parametric[:, i, 2, 2], [0.35355,  0.12499, 0.874999 ,0.972263]) atol=1e-5
    end

    for i=1:mesh.numNodesPerElement
      @test isapprox( eqn.flux_parametric[:, i, 1, 1], [0.35355,  0.874999, 00.124998,.972263]) atol=1e-5
    end

    for i=1:mesh.numNodesPerElement
      @test isapprox( eqn.flux_parametric[:, i, 2, 1], [0.0, 0.750001, -0.750001, 0.0]) atol=1e-5
    end


    # test getBCFluxes
    for j= 1:mesh.numNodesPerFace
      @test isapprox( eqn.bndryflux[:, j, 1], [-0.35355, -0.874999, -0.124998, -0.972263]) atol=1e-5
    end

    for j= 1:mesh.numNodesPerFace
      @test isapprox( eqn.bndryflux[:, j, 2], [-0.35355,  -0.124998, -0.874999, -0.972263]) atol=1e-5
    end

    for j= 1:mesh.numNodesPerFace
      @test isapprox( eqn.bndryflux[:, j, 3], [0.35355,  0.124998, 0.874999, 0.972263]) atol=1e-5
    end

    for j= 1:mesh.numNodesPerFace
      @test isapprox( eqn.bndryflux[:, j, 4], [0.35355, 0.874999, 0.124998, 0.972263]) atol=1e-5
    end
  end  # end facts block

  return nothing
end

#test_lowlevel_dataprep(mesh, sbp, eqn, opts)
add_func2!(EulerTests, test_lowlevel_dataprep, "input_vals_channel.jl", [TAG_FLUX, TAG_SHORTTEST])

"""
  This function tests evaluating volume and boundary integrals over the
  entire mesh.
"""
function test_lowlevel_integrals(mesh, sbp, eqn, opts)

  EulerEquationMod.array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)
  EulerEquationMod.dataPrep(mesh, sbp, eqn, opts)

  @testset "--- Testing evalVolumeIntegrals ---" begin

    EulerEquationMod.evalVolumeIntegrals(mesh, sbp, eqn, opts)

    el1_res = [-0.35355  0  0.35355;
                -0.874999  0.750001  0.124998;
    -0.124998  -0.750001  0.874999;
    -0.972263  0  0.972263]
    el2_res = [-0.35355  0.35355 0;
                -0.124998  0.874999 -0.75001;
    -0.874999 0.124998 0.75001;
    -0.972263  0.972263 0]
 
    @test isapprox( eqn.res[:, :, 2], el1_res) atol=1e-4
    @test isapprox( eqn.res[:, :, 1], el2_res) atol=1e-4

  end  # end facts block


  @testset "--- Testing evalBoundaryIntegrals ---" begin
    fill!(eqn.res, 0.0)

    EulerEquationMod.evalBoundaryIntegrals(mesh, sbp, eqn, opts)

    el1_res = [0.35355 0 -0.35355;
               0.124998 -0.750001 -0.874999;
         0.874999 0.750001 -0.124998;
         0.972263  0  -0.972263]
    el2_res = [0.35355 -0.35355 0;
               0.874999 -0.124998 0.750001;
         0.124998  -0.874999  -0.750001;
         0.972263  -0.972263  0]

    @test isapprox( eqn.res[:, :, 2], el1_res) atol=1e-5
    @test isapprox( eqn.res[:, :, 1], el2_res) atol=1e-5

  end  # end facts block

  @testset "--- Testing evalResidual --- " begin

    fill!(eqn.res_vec, 0.0)
    fill!(eqn.res, 0.0)
    EulerEquationMod.evalResidual(mesh, sbp, eqn, opts)

    # for uniform flow, residual is zero
    for i=1:mesh.numDof
      @test isapprox( eqn.res_vec[i], 0.0) atol=1e-14
    end

  end  # end facts block


  return nothing
end

#test_lowlevel_integrals(mesh, sbp, eqn, opts)
add_func2!(EulerTests, test_lowlevel_integrals, "input_vals_channel.jl", [TAG_BC, TAG_FLUX, TAG_SHORTTEST])

