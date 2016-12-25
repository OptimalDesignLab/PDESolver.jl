
const test_3d_inputfile = "input_vals_3d.jl"  # input file used by all tests

"""
  Test converting back and forth between conservative and entropy variables
  in 3d
"""
function test_3d_conversion(mesh, sbp, eqn, opts)

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

    vIR = copy(q)
    EulerEquationMod.convertToIR(eqn.params, vIR, vIR)
    diff = vIR - q2./eqn.params.gamma_1
    @fact norm(diff) --> roughly(0.0, atol=1e-12)

    EulerEquationMod.convertToConservativeFromIR_(eqn.params, vIR, vIR)
    diff = vIR - q
    @fact norm(diff) --> roughly(0.0, atol=1e-12)
  

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
  end  # end facts block

  return nothing
end  # end function

#test_3d_conversion(mesh, sbp, eqn, opts)
add_func2!(EulerTests, test_3d_conversion,  test_3d_inputfile, [TAG_ENTROPYVARS])


function test_3d_eigensystem(mesh, sbp, eqn, opts)

  facts("----- Testing Eigensystem -----") do

    params = eqn.params
    q = [1., 2, 3, 4, 15]
    numDofPerNode = length(q)
    qc = convert(Array{Complex128}, q)
    aux_vars = Array(Complex128, 1)
    qg = deepcopy(q)
    dxidx = mesh.dxidx[:, :, 1, 1]  # arbitrary
    F = zeros(Complex128, numDofPerNode)

    Ax = zeros(numDofPerNode, numDofPerNode)
    Ax3 = zeros(Ax)
    Ay = zeros(Ax)
    Ay3 = zeros(Ay)
    Az = zeros(Ax)
    Az3 = zeros(Ax)

    EulerEquationMod.calcA(params, q, Ax3)
    EulerEquationMod.calcB(params, q, Ay3)
    EulerEquationMod.calcC(params, q, Az3)
    h = 1e-20
    pert = Complex128(0, h)
    for i=1:numDofPerNode
      qc[i] += pert
      dir = [1.0, 0.0, 0]
      p = EulerEquationMod.calcPressure(params, qc)
      aux_vars[1] = p
      EulerEquationMod.calcEulerFlux(params, qc, aux_vars, dir, F)
      Ax[:, i] = imag(F)/h

      dir = [0.0, 1.0, 0.0]
      EulerEquationMod.calcEulerFlux(params, qc, aux_vars, dir, F)
      Ay[:, i] = imag(F)/h

      dir = [0.0, 0.0, 1.0]
      EulerEquationMod.calcEulerFlux(params, qc, aux_vars, dir, F)
      Az[:, i] = imag(F)/h

      qc[i] -= pert
    end

    # now compute Ax and and Ay from their eigensystem and compare
    Yx = zeros(Ax)
    Yy = zeros(Ax)
    Yz = zeros(Ax)
    Lambdax = zeros(numDofPerNode)
    Lambday = zeros(numDofPerNode)
    Lambdaz = zeros(numDofPerNode)

    EulerEquationMod.calcEvecsx(params, q, Yx)
    EulerEquationMod.calcEvecsy(params, q, Yy)
    EulerEquationMod.calcEvecsz(params, q, Yz)
    EulerEquationMod.calcEvalsx(params, q, Lambdax)
    EulerEquationMod.calcEvalsy(params, q, Lambday)
    EulerEquationMod.calcEvalsz(params, q, Lambdaz)

    Ax2 = Yx*diagm(Lambdax)*inv(Yx)
    Ay2 = Yy*diagm(Lambday)*inv(Yy)
    Az2 = Yz*diagm(Lambdaz)*inv(Yz)

    @fact Ax2 --> roughly(Ax, atol=1e-12)
    @fact Ay2 --> roughly(Ay, atol=1e-12)
    @fact Az2 --> roughly(Az, atol=1e-12)
    @fact Ax3 --> roughly(Ax, atol=1e-12)
    @fact Ay3 --> roughly(Ay, atol=1e-12)
    @fact Az3 --> roughly(Az, atol=1e-12)

    # check A0 = Y*S2*Y.'
    A0 = zeros(Ax)
    Sx = zeros(numDofPerNode)
    Sy = zeros(numDofPerNode)
    Sz = zeros(numDofPerNode)

    EulerEquationMod.getIRA0(params, q, A0)
    EulerEquationMod.calcEScalingx(params, q, Sx)
    EulerEquationMod.calcEScalingy(params, q, Sy)
    EulerEquationMod.calcEScalingz(params, q, Sz)
    A02 = Yx*diagm(Sx)*Yx.'
    A03 = Yy*diagm(Sy)*Yy.'
    A04 = Yz*diagm(Sz)*Yz.'

    @fact A0 --> roughly(A02, atol=1e-12)
    @fact A0 --> roughly(A03, atol=1e-12) 
    @fact A0 --> roughly(A04, atol=1e-12) 

  end  # end facts block

  return nothing
end


add_func2!(EulerTests, test_3d_eigensystem,  test_3d_inputfile, [TAG_ENTROPYVARS])


"""
  Test calculation of Euler flux in 3D, for both conservative and entropy
  variables.
"""
function test_3d_flux(mesh, sbp, eqn, opts)
  params_e = EulerEquationMod.ParamType{3, :entropy, Float64, Float64, Float64}(mesh, sbp, opts, 1)

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
  end  # end facts block

  return nothing
end  # end function

#test_3d_flux(mesh, sbp, eqn, opts)
add_func2!(EulerTests, test_3d_flux,  test_3d_inputfile, [TAG_ENTROPYVARS, TAG_FLUX])

"""
  Test some auxiliary calculation functinos in 3D
"""
function test_3d_misc(mesh, sbp, eqn, opts)
  params_e = EulerEquationMod.ParamType{3, :entropy, Float64, Float64, Float64}(mesh, sbp, opts, 1)

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

  end  # end facts block

  return nothing
end  # end function

#test_3d_misc(mesh, sbp, eqn, opts)
add_func2!(EulerTests, test_3d_misc,  test_3d_inputfile, [TAG_ENTROPYVARS])

"""
  Test calculation of A0
"""
function test_3d_matrices(mesh, sbp, eqn, opts)
  params_e = EulerEquationMod.ParamType{3, :entropy, Float64, Float64, Float64}(mesh, sbp, opts, 1)

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

    # test IRA0 against A0 and a jacobian computed with complex step
    A03 = scale(A0, eqn.params.gamma_1)
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
      @fact A04_code[i] --> roughly(A04_test[i], atol=1e-13)
      @fact A04_code[i] --> roughly(A03[i], atol=1e-10)
    end

  end  # end facts block

  return nothing
end  # end function

#test_3d_matrices(mesh, sbp, eqn, opts)
add_func2!(EulerTests, test_3d_matrices,  test_3d_inputfile)

"""
  Test Roe solver in 3D
"""
function test_3d_bc(mesh, sbp, eqn, opts)
  facts("----- Testing BC Solvers -----") do

    q = [1., 2, 3, 4, 15]
    F = zeros(q)
    F2 = zeros(q)
    nrm = [1., 1, 1]
    dir = mesh.dxidx[:, :, 1, 1].'*nrm
    p = EulerEquationMod.calcPressure(eqn.params, q)
    aux_vars = [p]
    dxidx = mesh.dxidx[:, :, 1, 1]

    EulerEquationMod.RoeSolver(eqn.params, q, q, aux_vars, dxidx, nrm, F2)
    EulerEquationMod.calcEulerFlux(eqn.params, q, aux_vars, dir, F)

  #  println("F = \n", F)
  #  println("F2 = \n", F2)
    @fact F2 --> roughly(F, atol=1e-13)
  end  # end facts block

  return nothing
end  # end functiona

add_func2!(EulerTests, test_3d_bc,  test_3d_inputfile, [TAG_BC])
