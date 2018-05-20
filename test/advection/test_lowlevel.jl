type twoxBC <: BCType
end
function call(obj::twoxBC, params::AdvectionEquationMod.ParamType, u, coords, nrm, t)
  u_bc = 2*coords[1]
  bndryflux = AdvectionEquationMod.RoeSolver(params, u, u_bc, nrm)
  return bndryflux
end

# input file used by all tests
global const test_lowlevel_inputfile = "input_vals_channel.jl"

"""
  Test that the mesh counts things correctly, and that coordinates and metrics
  are calculated correctly.  This makes it dependent on the element ordering
"""
function test_lowlevel_mesh(mesh, sbp, eqn, opts)
  facts("--- Testing Mesh --- ") do

    @fact mesh.numVert --> 4
    @fact mesh.numEdge --> 5
    @fact mesh.numEl --> 2
    @fact mesh.order --> 1
    @fact mesh.numDof --> 4
    @fact mesh.numNodes --> 4
    @fact mesh.numDofPerNode --> 1
    @fact mesh.numBoundaryFaces --> 4
    @fact mesh.numInterfaces --> 1
    @fact mesh.numNodesPerElement --> 3
    @fact mesh.numNodesPerType --> [1, 0 , 0]

    @fact mesh.bndry_funcs[1] --> AdvectionEquationMod.x5plusy5BC()
    @fact mesh.bndryfaces[1].element --> 1
    @fact mesh.bndryfaces[1].face --> 3
    @fact mesh.bndryfaces[2].element --> 1
    @fact mesh.bndryfaces[2].face --> 1
    @fact mesh.bndryfaces[3].element --> 2
    @fact mesh.bndryfaces[3].face --> 2
    @fact mesh.bndryfaces[4].element --> 2
    @fact mesh.bndryfaces[4].face --> 1

    @fact mesh.interfaces[1].elementL --> 1
    @fact mesh.interfaces[1].elementR --> 2
    @fact mesh.interfaces[1].faceL --> 2
    @fact mesh.interfaces[1].faceR --> 3

    jac_fac = 0.25
    fac = 2
    @fact mesh.coords[:, :, 2] --> roughly([4 4 0; 0 4 4.0])
    @fact mesh.coords[:, :, 1] --> roughly([0.0 4 0; 0 0 4])

    @fact mesh.dxidx[:, :, 1, 2] --> roughly(fac*[1 1; -1 0.0], atol=1e-14)
    @fact mesh.dxidx[:, :, 2, 2] --> roughly(fac*[1 1; -1 0.0], atol=1e-14)
    @fact mesh.dxidx[:, :, 3, 2] --> roughly(fac*[1 1; -1 0.0], atol=1e-14)

    @fact mesh.dxidx[:, :, 1, 1] --> roughly(fac*[1 0; 0 1.0], atol=1e-14)
    @fact mesh.dxidx[:, :, 2, 1] --> roughly(fac*[1 0; 0 1.0], atol=1e-14)
    @fact mesh.dxidx[:, :, 3, 1] --> roughly(fac*[1 0; 0 1.0], atol=1e-14)

    @fact mesh.jac --> roughly(jac_fac*ones(3,2))

  end  # end fact block

  return nothing
end

#test_lowlevel_mesh(mesh, sbp, eqn, opts)
add_func2!(AdvectionTests, test_lowlevel_mesh, test_lowlevel_inputfile, [TAG_SHORTTEST])

"""
  Test some basic functionality: dissassmble solution, mass matrix, etc.
"""
function test_lowlevel_core(mesh, sbp, eqn, opts)
  facts("--- Testing Functions Within AdvectionData_--- ") do
    Tsol = Float64
    u_vec = Tsol[1,2,3,4]
    u  = zeros(Tsol, 1, 3, 2)

    # checking disassembleSolution
    disassembleSolution(mesh, sbp, eqn, opts, u, u_vec)
    @fact u[1,1,2] --> roughly(2.0)
    @fact u[1,2,2] --> roughly(4.0)
    @fact u[1,3,2] --> roughly(3.0)
    @fact u[1,1,1] --> roughly(1.0)
    @fact u[1,2,1] --> roughly(2.0)
    @fact u[1,3,1] --> roughly(3.0)

    #checking assembleSolution
    fill!(u_vec, 0.0)
    assembleSolution(mesh, sbp, eqn, opts, u, u_vec)
    @fact u_vec --> roughly([1.0,4.0,6.0,4.0])

    # check mass matrix
    # just for testing, make jac != 1
    w_val = sbp.w[1]  # all of sbp.w entries are the same
    jac_val = 0.25
    M1 = w_val/jac_val
    M2 = 2*w_val/jac_val
    M3 = 2*w_val/jac_val
    M4 = w_val/jac_val
    M_test = [M1, M2, M3, M4]
    M_code = AdvectionEquationMod.calcMassMatrix(mesh, sbp, eqn)
    @fact M_code --> roughly(M_test, atol=1e-13)

    Minv_test = 1./M_test
    Minv_test = AdvectionEquationMod.calcMassMatrixInverse(mesh, sbp, eqn)
    @fact Minv_test --> roughly(Minv_test, atol=1e-13)

    arr = rand(1, 3, 2)
    arr_orig = copy(arr)
    AdvectionEquationMod.matVecA0inv(mesh, sbp, eqn, opts, arr)
    @fact arr --> roughly(arr_orig, atol=1e-14)

    AdvectionEquationMod.matVecA0(mesh, sbp, eqn, opts, arr)
    @fact arr --> roughly(arr_orig, atol=1e-14)
  end

  return nothing
end

#test_lowlevel_core(mesh, sbp, eqn, opts)
add_func2!(AdvectionTests, test_lowlevel_core, test_lowlevel_inputfile, [TAG_SHORTTEST])

"""
  Test things from common_funcs.jl
"""
function test_lowlevel_common(mesh, sbp, eqn, opts)
  facts("--- Testing common functions ---") do

    x = 1.
    y = 2.
    coords = [x, y]
    eqn.params.alpha_x = 0.0
    eqn.params.alpha_y = 0.0
    t = 0.0
    val = AdvectionEquationMod.calc_x5plusy5(eqn.params, coords, t)
    @fact val --> roughly(x^5 + y^5, atol=1e-14)
   
    val = AdvectionEquationMod.calc_exp_xplusy(eqn.params, coords, t)
    @fact val --> roughly(exp(x + y), atol=1e-14)
  end  # end facts block

  return nothing
end

#test_lowlevel_common(mesh, sbp, eqn, opts)
add_func2!(AdvectionTests, test_lowlevel_common, test_lowlevel_inputfile, [TAG_SHORTTEST])

"""
  Test Roe solver used for boundary conditions and some hand calculated
  boundary flux values.
"""
function test_lowlevel_bc(mesh, sbp, eqn, opts)
  facts("--- Testing Boundary Function ---") do

    # testing choice of u or u_bc
    u = 5.0
    u_bc = 2.5
    eqn.params.alpha_x = 1.5
    eqn.params.alpha_y = 0.5
    dxidx = [1. 0; 0 1]
    nrm = [1., 0]
    nrm2 = zeros(nrm)
    calcBCNormal(eqn.params, dxidx, nrm, nrm2)


    val = AdvectionEquationMod.RoeSolver(eqn.params, u, u_bc, nrm2)

    @fact val --> roughly(u*eqn.params.alpha_x, atol=1e-14)

    nrm = [-1.0, 0]
    calcBCNormal(eqn.params, dxidx, nrm, nrm2)
    val = AdvectionEquationMod.RoeSolver(eqn.params, u, u_bc, nrm2)
    @fact val --> roughly(-u_bc*eqn.params.alpha_x, atol=1e-14)


    nrm = [0, 1.0]
    calcBCNormal(eqn.params, dxidx, nrm, nrm2)
    val = AdvectionEquationMod.RoeSolver(eqn.params, u, u_bc, nrm2)
    @fact val --> roughly(u*eqn.params.alpha_y, atol=1e-14)

    nrm = [0, -1.0]
    calcBCNormal(eqn.params, dxidx, nrm, nrm2)
    val = AdvectionEquationMod.RoeSolver(eqn.params, u, u_bc, nrm2)
    @fact val --> roughly(-u_bc*eqn.params.alpha_y, atol=1e-14)

    # now test rotation using dxidx
    
    function get_rotation_matrix(theta)
      return [cos(theta) -sin(theta); sin(theta) cos(theta)]
    end

    theta = 30*pi/180  # angle between x axis and xi axis
    flow_direction = 25*pi/180
    alpha_mag = 2.0
    eqn.params.alpha_x = alpha_mag*cos(flow_direction)
    eqn.params.alpha_y = alpha_mag*sin(flow_direction)

    dxidx = get_rotation_matrix( theta)

    # nrm is in the xi-eta space
    # the xi axis is 30 degrees clockwise of the x axis, so the alpha vector
    # which is already at 25 degrees from the x axis needs to go 30 more
    # to be in the xi-eta coordinate system
    # this is still an outflow for nrm = [1, 0]
    nrm = [1., 0]  # flow is in the xi direction
    calcBCNormal(eqn.params, dxidx, nrm, nrm2)
    angle_diff = theta + flow_direction
    alpha_eff = alpha_mag*cos(angle_diff)  # effective alpha in the wall normal
                                           # direction
    val_exp = alpha_eff*u
    val = AdvectionEquationMod.RoeSolver(eqn.params, u, u_bc, nrm2)

    @fact val --> roughly(val_exp, atol=1e-14)

    # now check eta direction
    nrm = [0, 1.]
    calcBCNormal(eqn.params, dxidx, nrm, nrm2)
    alpha_eff = alpha_mag*sin(angle_diff)
    val_exp = u*alpha_eff
    val = AdvectionEquationMod.RoeSolver(eqn.params, u, u_bc, nrm2)

    @fact val --> roughly(val_exp, atol=1e-14)

    # now rotate the coordinate system so much that this becomes an inflow
    theta = 120*pi/180
    nrm = [1., 0]
    angle_diff = theta + flow_direction
    alpha_eff = alpha_mag*cos(angle_diff)
    val_exp = alpha_eff*u_bc

    dxidx = get_rotation_matrix( theta)
    calcBCNormal(eqn.params, dxidx, nrm, nrm2)
    val = AdvectionEquationMod.RoeSolver(eqn.params, u, u_bc, nrm2)
    @fact val --> roughly(val_exp, atol=1e-14)

    # check eta direction
    nrm = [0, 1]
    calcBCNormal(eqn.params, dxidx, nrm, nrm2)
    angle_diff = theta + flow_direction
    alpha_eff = alpha_mag*sin(angle_diff)
    val_exp = alpha_eff*u
    val = AdvectionEquationMod.RoeSolver(eqn.params, u, u_bc, nrm2)
    @fact val --> roughly(val_exp, atol=1e-14)
  end  # end facts block

  return nothing
end

#test_lowlevel_bc(mesh, sbp, eqn, opts)
add_func2!(AdvectionTests, test_lowlevel_bc, test_lowlevel_inputfile, [TAG_BC, TAG_SHORTTEST])

"""
  Test computing volume integrals over entire mesh
"""
function test_lowlevel_volumeintegrals()
  facts("--- Testing Volume Integrals ---")  do

    # use the 8 element mesh
    ARGS[1] = "input_vals_8el.jl"
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])

    println("precompute_boundary_flux = ", opts["precompute_boundary_flux"])
    Tmsh = eltype(mesh.dxidx)

    fill!(eqn.q, 0.0)
    eqn.params.alpha_x = 1.0
    eqn.params.alpha_y = 1.0
#    eqn.params.alpha_x = ones(Float64 ,1, mesh.numNodesPerElement, mesh.numEl)
    alpha_y = zero(eqn.params.alpha_x)

    fill!(eqn.res, 0.0)
    AdvectionEquationMod.evalVolumeIntegrals(mesh, sbp, eqn, opts)
    assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
    @fact eqn.res --> roughly(zeros(1, mesh.numNodesPerElement, mesh.numEl), atol=1e-12)


    # check that Qx.'*q = 0 when q = ones(3)
    fill!(eqn.q, 1.0)
    fill!(eqn.res, 0.0)
    dxidx1 = mesh.dxidx[:, :, :, 1]
    jac = mesh.jac[:, 1]
    dxidx_true = zeros(2, 2, 3)
    for i=1:3
      dxidx_true[:, :, i] = dxidx1[:, :, i]*jac[i]
    end

    Qx = sbp.Q[:, :, 1]*dxidx1[1, 1, 1] + sbp.Q[:, :, 2]*dxidx1[2, 1, 1]
    q = [eqn.q[1, 1, 1], eqn.q[1, 2, 1], eqn.q[1, 3, 1]]  # extract q values
    # sum reduces the vector to the value of the integral
    val_test = sum(Qx.'*q)  
    AdvectionEquationMod.evalVolumeIntegrals(mesh, sbp, eqn, opts)
    # extract residual values
    vec_code =  [eqn.res[1, 1, 1], eqn.res[1, 2, 1], eqn.res[1, 3, 1]]
    val_code = sum(vec_code)
    xi_flux = dxidx1[1,1,1]*1*q
    xi_component = sbp.Q[:, :, 1].'*xi_flux
    eta_flux = dxidx1[2, 1]*1*q
    eta_component = sbp.Q[:, :, 2].'*eta_flux
    val2_test = sum(xi_component + eta_component)
    @fact val_test --> roughly(val2_test, atol=1e-14)
    @fact val_test --> roughly(val_code, atol=1e-14)


    # test that the integral of qx dOmega^e when q = 2*x works
    # dq/dx = 2, thus the integral is = 1

    println("----- Checking q=2*x case -----")
    bc_func = twoxBC()
    mesh.bndry_funcs[1] = bc_func
    
    x1 = zeros(Tmsh, 3)
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        x = mesh.coords[1, j, i]
        eqn.q[1, j, i] = 2*x
        if i == 1
          x1[j] = x
        end

      end
    end

    fill!(eqn.res, 0.0)
    eqn.params.alpha_x = 1.0
    eqn.params.alpha_y = 0.0
    # check the boundry contribution
    AdvectionEquationMod.evalBoundaryIntegrals(mesh, sbp, eqn, opts)

    @fact sum(eqn.res[:, :, 1]) --> roughly(-2.0, atol=1e-14)
    @fact sum(eqn.res[:, :, 3]) --> roughly(-2.0, atol=1e-14)

    @fact sum(eqn.res[:, :, 6]) --> roughly(-2.0, atol=1e-14)
    @fact sum(eqn.res[:, :, 8]) --> roughly(-2.0, atol=1e-14)

    fill!(eqn.res, 0.0)
    AdvectionEquationMod.evalVolumeIntegrals(mesh, sbp, eqn, opts)
    for i=1:mesh.numEl
      Qx_i = sbp.Q[:, :, 1]*mesh.dxidx[1, 1, 1, i] + sbp.Q[:, :, 2]*mesh.dxidx[2, 1, 1, i]
      q_i = reshape(eqn.q[1, :, i], 3)
      val_test = Qx_i.'*q_i
      val_code = reshape(eqn.res[:, :, i], 3)
      @fact sum(val_code) --> roughly(0.0, atol=1e-14)  # proven by hand calc
      @fact val_code --> roughly(val_test, atol=1e-14)
    end

    println("----- Checking q=2*x^2 + 5 case -----")
    
    x1 = zeros(Tmsh, 3)
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        x = mesh.coords[1, j, i]
        eqn.q[1, j, i] = 2*x*x + 5
      end
    end

    fill!(eqn.res, 0.0)
    eqn.params.alpha_x = 1.0
    eqn.params.alpha_y = 0.0
    AdvectionEquationMod.evalVolumeIntegrals(mesh, sbp, eqn, opts)
    for i=1:mesh.numEl
      Qx_i = sbp.Q[:, :, 1]*mesh.dxidx[1, 1, 1, i] + sbp.Q[:, :, 2]*mesh.dxidx[2, 1, 1, i]
      q_i = reshape(eqn.q[1, :, i], 3)
      val_test = Qx_i.'*q_i
      val_code = reshape(eqn.res[:, :, i], 3)
      @fact val_code --> roughly(val_test, atol=1e-14)
    end


    ARGS[1] = "input_vals_8el_large.jl"
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])
    println("----- Checking q=2*x^2 + 5 case on large grid -----")
    
    x1 = zeros(Tmsh, 3)
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        x = mesh.coords[1, j, i]
        eqn.q[1, j, i] = 2*x*x + 5
      end
    end

    fill!(eqn.res, 0.0)
    eqn.params.alpha_x = 1.0
    eqn.params.alpha_y = 0.0
    AdvectionEquationMod.evalVolumeIntegrals(mesh, sbp, eqn, opts)
    for i=1:mesh.numEl
      Qx_i = sbp.Q[:, :, 1]*mesh.dxidx[1, 1, 1, i] + sbp.Q[:, :, 2]*mesh.dxidx[2, 1, 1, i]
      q_i = reshape(eqn.q[1, :, i], 3)
      val_test = Qx_i.'*q_i
      val_code = reshape(eqn.res[:, :, i], 3)
      @fact val_code --> roughly(val_test, atol=1e-14)
    end

    # back to the original mesh
    ARGS[1] = "input_vals_8el.jl"
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])

    println("----- Checking sinwave case -----")
    x1 = zeros(Tmsh, 3)
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        x = mesh.coords[1, j, i]
        alpha_x = eqn.params.alpha_x
        alpha_y = eqn.params.alpha_y
        eqn.q[1, j, i] = AdvectionEquationMod.calc_sinwave(eqn.params, mesh.coords[:, j, i], 0.25)
      end
    end

    mesh.bndry_funcs[1] = AdvectionEquationMod.sinwave_BC()
    fill!(eqn.res, 0.0)
    eqn.params.alpha_x = 1.0
    eqn.params.alpha_y = 0.0
    AdvectionEquationMod.evalVolumeIntegrals(mesh, sbp, eqn, opts)
    for i=1:mesh.numEl
      Qx_i = sbp.Q[:, :, 1]*mesh.dxidx[1, 1, 1, i] + sbp.Q[:, :, 2]*mesh.dxidx[2, 1, 1, i]
      q_i = reshape(eqn.q[1, :, i], 3)
      val_test = Qx_i.'*q_i
      val_code = reshape(eqn.res[:, :, i], 3)
      @fact val_code --> roughly(val_test, atol=1e-14)
    end

    assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
#TODO: uncomment when SBP boundaryintegrate is fixed
#= 
    ARGS[1] = "input_vals_channel_verylarge.jl"
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])

    println("----- Checking sinwave case on very large mesh -----")
    x1 = zeros(Tmsh, 3)
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        x = mesh.coords[1, j, i]
        alpha_x = eqn.params.alpha_x
        alpha_y, = eqn.params.alpha_y

        eqn.q[1, j, i] = AdvectionEquationMod.calc_sinwave(mesh.coords[:, j, i], eqn.params, 0.25)
      end
    end
    println("finished calculating sinwave")
    mesh.bndry_funcs[1] = AdvectionEquationMod.sinwave_BC()
    fill!(eqn.res, 0.0)
    eqn.params.alpha_x = 1.0
    eqn.params.alpha_y = 0.0
    AdvectionEquationMod.evalVolumeIntegrals(mesh, sbp, eqn, opts)
    println("called eval SCResidual")
    for i=1:mesh.numEl
      Qx_i = sbp.Q[:, :, 1]*mesh.dxidx[1, 1, 1, i] + sbp.Q[:, :, 2]*mesh.dxidx[2, 1, 1, i]
      q_i = reshape(eqn.q[1, :, i], 3)
      val_test = Qx_i.'*q_i
      val_code = reshape(eqn.res[:, :, i], 3)
      @fact val_code --> roughly(val_test, atol=1e-14)
    end
  =#
  end # end facts block

  return nothing
end

#test_lowlevel_volumeintegrals()
add_func1!(AdvectionTests, test_lowlevel_volumeintegrals, [TAG_VOLUMEINTEGRALS, TAG_SHORTTEST])
  #=
  context("--- Testing evalInteriorFlux ---") do
    
    fill!(eqn.res, 0.0)
    fill!(eqn.params.alpha_x, 1.0)
    fill!(eqn.params.alpha_y, 1.0)
    
    nbrnodeindex = Array(mesh.numNodesPerFace:-1:1)

  end
  =#




#=
 context("--- Testing dataPrep ---") do
 
  end



  context("--- Testing evalBoundaryIntegrals ---") do
  end

  context("--- Testing evalEuler --- ")  do

end # end facts block
=#
