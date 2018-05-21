
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
add_func2!(EulerTests, test_3d_conversion,  test_3d_inputfile, [TAG_ENTROPYVARS, TAG_SHORTTEST])


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


add_func2!(EulerTests, test_3d_eigensystem,  test_3d_inputfile, [TAG_ENTROPYVARS, TAG_SHORTTEST])


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
add_func2!(EulerTests, test_3d_flux,  test_3d_inputfile, [TAG_ENTROPYVARS, TAG_FLUX, TAG_SHORTTEST])

"""
  Test some auxiliary calculation functions in 3D
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
add_func2!(EulerTests, test_3d_misc,  test_3d_inputfile, [TAG_ENTROPYVARS, TAG_MISC, TAG_SHORTTEST])


function test_3d_secondary_quantities()

  fname = "input_vals_3d.jl"
  opts = read_input(fname)

  for p=1:4
    opts["order"] = p
    fname = "input_vals_3d_p$p"
    make_input(opts, fname)
    mesh, sbp, eqn, opts = solvePDE(fname*".jl")

    coords = mesh.coords[:, :, 1]
    dxidx = mesh.dxidx[:, :, :, 1]
    jac = mesh.jac[:, 1]
    q = zeros(mesh.numDofPerNode, mesh.numNodesPerElement)
    vorticity = zeros(3, mesh.numNodesPerElement)

    ufunc = (x, y, z) -> 2*x^p + y^p + z^p + 1
    vfunc = (x, y, z) -> x^p + 2*y^p + z^p + 1
    wfunc = (x, y, z) -> x^p + y^p + 2*z^p + 1

    dudx = (x, y, z) -> 2*p*x^(p-1)
    dudy = (x, y, z) -> 1*p*y^(p-1)
    dudz = (x, y, z) -> 1*p*z^(p-1)

    dvdx = (x, y, z) -> 1*p*x^(p-1)
    dvdy = (x, y, z) -> 2*p*y^(p-1)
    dvdz = (x, y, z) -> 1*p*z^(p-1)

    dwdx = (x, y, z) -> 1*p*x^(p-1)
    dwdy = (x, y, z) -> 1*p*y^(p-1)
    dwdz = (x, y, z) -> 2*p*z^(p-1)

    for itr=0:1
      # populate q
      for i=1:mesh.numNodesPerElement
        x = coords[1, i] + itr
        y = coords[2, i] + itr
        z = coords[3, i] + itr

        q[1, i] = 1
        q[2, i] = ufunc(x, y, z)
        q[3, i] = vfunc(x, y, z)
        q[4, i] = wfunc(x, y, z)
        q[5, i] = 10
      end
      EulerEquationMod.calcVorticity(eqn.params, sbp, q, dxidx, jac, vorticity)

      for i=1:mesh.numNodesPerElement
        x = coords[1, i] + itr
        y = coords[2, i] + itr
        z = coords[3, i] + itr

        vortx = dwdy(x, y, z) - dvdz(x, y, z)
        vorty = -dwdx(x, y, z) + dudz(x, y, z)
        vortz = dvdx(x, y, z) - dudy(x, y, z)

        @fact vortx --> roughly(vorticity[1, i], atol=1e-12)
        @fact vorty --> roughly(vorticity[2, i], atol=1e-12)
        @fact vortz --> roughly(vorticity[3, i], atol=1e-12)
      end
    end  # end loop itr
  end  # end loop over p

  # the enstrophy integral is O(2p + 1), so for integration to be exact, make
  # solution 2*p - 1
  for pprime = 3:4
    p = pprime - 2
#    p = 1
    opts["order"] = pprime
    fname = "input_vals_3d_p$p"
    make_input(opts, fname)
    mesh, sbp, eqn, opts = solvePDE(fname*".jl")

    ufunc = (x, y, z) -> 2*x^p + y^p + z^p + 1
    vfunc = (x, y, z) -> x^p + 2*y^p + z^p + 1
    wfunc = (x, y, z) -> x^p + y^p + 2*z^p + 1

    q = eqn.q
    for el = 1:mesh.numEl
      for i=1:mesh.numNodesPerElement
        x = mesh.coords[1, i, el]
        y = mesh.coords[2, i, el]
        z = mesh.coords[3, i, el]

        q[1, i, el] = 1
        q[2, i, el] = ufunc(x, y, z)
        q[3, i, el] = vfunc(x, y, z)
        q[4, i, el] = wfunc(x, y, z)
        q[5, i, el] = 10
      end
    end

    # ----- test enstrophy ------
    enstrophy_exact = (24*3^p - (24*p^2)/(2*p - 1) - 12*3^(2*p) + (8*3^(2*p)*p^2)/(2*p - 1) - 12)/(2*mesh.volume)


    enstrophy_numerical = EulerEquationMod.calcEnstrophy(mesh, sbp, eqn, opts, q)

    @fact enstrophy_exact --> roughly(enstrophy_numerical, atol=1e-12)

    #----- test Kinetic Energy -----
    q_vec = reshape(q, mesh.numDof)

    if p == 1
      ke_exact = 1992
    elseif p== 2
      ke_exact = 132712/15
    end

    ke_exact /= 2*mesh.volume

    ke_numerical = EulerEquationMod.calcKineticEnergy(mesh, sbp, eqn, opts, q_vec)

    @fact ke_exact --> roughly(ke_numerical, atol=1e-12)


    # ----- test kinetic energy dissipation rate ------
    t = 0.5
    ufunc = (x, y, z) -> 2*x^p + y^p + z^p + 1 + t^p
    vfunc = (x, y, z) -> x^p + 2*y^p + z^p + 1 + t^p
    wfunc = (x, y, z) -> x^p + y^p + 2*z^p + 1 + t^p
    rhofunc = (x, y, z) -> t^p

    dudt = (x, y, z) -> p*t^(p-1)
    dvdt = (x, y, z) -> p*t^(p-1)
    dwdt = (x, y, z) -> p*t^(p-1)
    drhodt = (x, y, z) -> p*t^(p-1)

    q = eqn.q
    res = zeros(q)
    for el = 1:mesh.numEl
      for i=1:mesh.numNodesPerElement
        x = mesh.coords[1, i, el]
        y = mesh.coords[2, i, el]
        z = mesh.coords[3, i, el]

        rho = rhofunc(x, y, z)
        q[1, i, el] = rho
        q[2, i, el] = rho*ufunc(x, y, z)
        q[3, i, el] = rho*vfunc(x, y, z)
        q[4, i, el] = rho*wfunc(x, y, z)
        q[5, i, el] = 10

        res[1, i, el] = drhodt(x, y, z)
        res[2, i, el] = rhofunc(x,  y, z)*dudt(x, y, z) + ufunc(x, y, z)*drhodt(x, y, z)
        res[3, i, el] = rhofunc(x,  y, z)*dvdt(x, y, z) + vfunc(x, y, z)*drhodt(x, y, z)

        res[4, i, el] = rhofunc(x,  y, z)*dwdt(x, y, z) + wfunc(x, y, z)*drhodt(x, y, z)

        res[5, i, el] = 0

      end
    end

    q_vec = reshape(q, mesh.numDof)
    res_vec = reshape(res, mesh.numDof)

    if p == 1
      dkedt_exact = 24*t + 216
    else
      dkedt_exact = 16*t*(3*t^2 + 55)
    end
    dkedt_exact *= rhofunc(1, 1, 1)/mesh.volume

    dkedt_numerical = EulerEquationMod.calcKineticEnergydt(mesh, sbp, eqn, opts, q_vec, res_vec)

    @fact dkedt_exact --> roughly(dkedt_numerical, atol=1e-12)


  end  # end loop over p



  return nothing

end 

add_func1!(EulerTests, test_3d_secondary_quantities, [TAG_MISC, TAG_SHORTTEST])



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
add_func2!(EulerTests, test_3d_matrices,  test_3d_inputfile, [TAG_SHORTTEST])

"""
  Test Roe solver in 3D
"""
function test_3d_bc(mesh, sbp, eqn, opts)
  facts("----- Testing BC Solvers -----") do

    params = eqn.params
    q = [1., 2, 3, 4, 15]
    F = zeros(q)
    F2 = zeros(q)
    nrm = [1., 1, 1]
    dir = mesh.dxidx[:, :, 1, 1].'*nrm
    p = EulerEquationMod.calcPressure(eqn.params, q)
    aux_vars = [p]
    dxidx = mesh.dxidx[:, :, 1, 1]
    nrm_xy = zeros(mesh.dim)

    calcBCNormal(params, dxidx, nrm, nrm_xy)
    EulerEquationMod.RoeSolver(params, q, q, aux_vars, nrm_xy, F2)
    EulerEquationMod.calcEulerFlux(params, q, aux_vars, dir, F)

    @fact F2 --> roughly(F, atol=1e-13)

    func1 = EulerEquationMod.noPenetrationESBC()
    # make velocity parallel to the boundary
    nrm = mesh.nrm_bndry[:, 1, 1]
    tngt = mesh.coords_bndry[:, 2, 1] - mesh.coords_bndry[:, 1, 1]
    tngt = tngt./norm(tngt)
    vval = 3
    q = [1., vval*tngt[1], vval*tngt[2], vval*tngt[3], 15]
    aux_vars = eqn.aux_vars[:, 1, 1]
    coords = mesh.coords_bndry[:, 1, 1]


    func1(params, q, aux_vars, coords, nrm, F2)
    EulerEquationMod.calcEulerFlux(params, q, aux_vars, nrm, F)

    @fact F2 --> roughly(F, atol=1e-13)

    func1 = EulerEquationMod.noPenetrationBC()
    func1(params, q, aux_vars, coords, nrm, F2)

    @fact F2 --> roughly(F, atol=1e-13)

  end  # end facts block

  return nothing
end  # end functiona

function test_3d_functional(mesh, sbp, eqn, opts)
  facts("----- Testing 3D functional -----") do

    q_bndry_orig = eqn.q_bndry
    q_bndry2  = zeros(eqn.q_bndry)
    eqn.q_bndry = q_bndry2
    s1 = sview(eqn.q_bndry, 1, :, :)
    fill!(s1, 1.0)
    s2 = sview(eqn.q_bndry, 2, :, :)
    fill!(s2, 0.355)
    s3 = sview(eqn.q_bndry, 3, :, :)
    fill!(s3, 0.355)
    s4 = sview(eqn.q_bndry, 4, :, :)
    fill!(s4, 0.355)
    s5 = sview(eqn.q_bndry, 5, :, :)
    fill!(s5, 3.0)


    moment_about = [2., 2, 2]
    mom = EulerEquationMod.calcMomentContribution!(mesh, eqn, [1], moment_about)
    println("moment = ", mom)
    for i=1:length(mom)
      @fact mom[i] --> roughly(0.0, atol=1e-13)
    end

    # test reverse mode
    moment_about = [0.0, 0, 0]
    nout = 3
    nin = mesh.dim*mesh.numNodesPerFace*mesh.numBoundaryFaces
    jac = zeros(nout, nin)

    mom_0 = EulerEquationMod.calcMomentContribution!(mesh, eqn, [1], moment_about)
    pert = 1e-6
    for i=1:nin
      mesh.nrm_bndry[i] += pert
      mom_i = EulerEquationMod.calcMomentContribution!(mesh, eqn, [1], moment_about)

      for j=1:nout
        jac[j, i] = (mom_i[j] - mom_0[j])/pert
      end

      mesh.nrm_bndry[i] -= pert
    end

    # reverse mode
    jac2 = zeros(jac)
    mom_bar = zeros(mom_0)
    for i=1:nout
      mom_bar[i] = 1
      fill!(mesh.nrm_bndry_bar, 0.0)
      EulerEquationMod.calcMomentContribution_revm!(mesh, eqn, [1], moment_about, mom_bar)

      for j=1:nin
        jac2[i, j] = mesh.nrm_bndry_bar[j]
      end

      mom_bar[i] = 0
    end

    #=
    for i=1:nin
      @fact norm(jac[:, i] - jac2[:, i])/size(jac, 1) --> roughly(0.0, atol=1e-5)
    end
    =#
    @fact norm(jac - jac2)/sqrt(length(jac)) --> roughly(0.0, atol=1e-5)


    # test other method
    println("checking second method")
    face_range = mesh.bndry_offsets[1]:(mesh.bndry_offsets[2] - 1)
    bndryfaces = mesh.bndryfaces[face_range]
    coords_faces = mesh.coords_bndry[:, :, face_range]
    nrm = mesh.nrm_bndry[:, :, face_range]

    dforce = EulerEquationMod.computeDForce(mesh, eqn, bndryfaces, nrm)
    mom_0 = EulerEquationMod.calcMomentContribution!(mesh.sbpface, coords_faces, dforce, moment_about)
    nin = length(dforce)
    nout = 3
    jac = zeros(nout, nin)
    for i=1:nin
      dforce[i] += pert
      mom_i = EulerEquationMod.calcMomentContribution!(mesh.sbpface, coords_faces, dforce, moment_about)

      for j=1:nout
        jac[j, i] = (mom_i[j] - mom_0[j])/pert
      end

      dforce[i] -= pert
    end

    jac2 = zeros(jac)
    dforce_bar = zeros(dforce)
    coords_bndry_bar = zeros(coords_faces)
    mom_bar = zeros(mom_0)
    for i=1:nout
      mom_bar[i] = 1
      fill!(dforce_bar, 0.0)
      fill!(coords_bndry_bar, 0.0)

      EulerEquationMod.calcMomentContribution_rev!(mesh.sbpface, coords_faces, coords_bndry_bar, dforce, dforce_bar, moment_about, mom_bar)
      for j=1:nin
        jac2[i, j] = dforce_bar[j]
      end

      mom_bar[i] = 0
    end

    #=
    for i=1:size(jac, 2)
      @fact norm(jac[:, i] - jac2[:, i]) --> roughly(0.0, atol=1e-5)
    end
    =#
    @fact norm(jac - jac2)/sqrt(length(jac)) --> roughly(0.0, atol=1e-5)


    eqn.q_bndry = q_bndry_orig
  end


  return nothing
end

add_func2!(EulerTests, test_3d_bc,  test_3d_inputfile, [TAG_BC, TAG_SHORTTEST])
add_func2!(EulerTests, test_3d_functional,  test_3d_inputfile, [TAG_REVERSEMODE, TAG_SHORTTEST])
