include("../../src/solver/advection/startup_advection.jl")  # initialization and construction
fill!(eqn.res_vec, 0.0)
using ArrayViews

type twoxBC <: BCType
end
function call(obj::twoxBC, u, alpha_x, alpha_y, coords, dxidx, nrm, t)
  u_bc = 2*coords[1]
  bndryflux = AdvectionEquationMod.RoeSolver(u, u_bc, alpha_x, alpha_y, nrm, dxidx)
  return bndryflux
end




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

  arr = rand(1, 3, 2)
  arr_orig = copy(arr)
  AdvectionEquationMod.matVecA0inv(mesh, sbp, eqn, opts, arr)
  @fact arr => roughly(arr_orig, atol=1e-14)

  AdvectionEquationMod.matVecA0(mesh, sbp, eqn, opts, arr)
  @fact arr => roughly(arr_orig, atol=1e-14)
end


 context("--- Testing common functions ---") do

   x = 1.
   y = 2.
   coords = [x, y]
   val = AdvectionEquationMod.calc_x5plusy5(coords)
   @fact val => roughly(x^5 + y^5, atol=1e-14)
   
   val = AdvectionEquationMod.calc_exp_xplusy(coords)
   @fact val => roughly(exp(x + y), atol=1e-14)


 end



  context("--- Testing Boundary Function ---") do

    # testing choice of u or u_bc
    u = 5.0
    u_bc = 2.5
    alpha_x = 1.5
    alpha_y = 0.5
    dxidx = [1. 0; 0 1]
    nrm = [1., 0]

    val = AdvectionEquationMod.RoeSolver(u, u_bc, alpha_x, alpha_y, nrm, dxidx)

    @fact val => roughly(u*alpha_x, atol=1e-14)

    nrm = [-1.0, 0]
    val = AdvectionEquationMod.RoeSolver(u, u_bc, alpha_x, alpha_y, nrm, dxidx)
    @fact val => roughly(-u_bc*alpha_x, atol=1e-14)


    nrm = [0, 1.0]
    val = AdvectionEquationMod.RoeSolver(u, u_bc, alpha_x, alpha_y, nrm, dxidx)
    @fact val => roughly(u*alpha_y, atol=1e-14)

    nrm = [0, -1.0]
    val = AdvectionEquationMod.RoeSolver(u, u_bc, alpha_x, alpha_y, nrm, dxidx)
    @fact val => roughly(-u_bc*alpha_y, atol=1e-14)

    # now test rotation using dxidx
    
    function get_rotation_matrix(theta)
      return [cos(theta) -sin(theta); sin(theta) cos(theta)]
    end

    theta = 30*pi/180  # angle between x axis and xi axis
    flow_direction = 25*pi/180
    alpha_mag = 2.0
    alpha_x = alpha_mag*cos(flow_direction)
    alpha_y = alpha_mag*sin(flow_direction)

    dxidx = get_rotation_matrix( theta)

    # nrm is in the xi-eta space
    # the xi axis is 30 degrees clockwise of the x axis, so the alpha vector
    # which is already at 25 degrees from the x axis needs to go 30 more
    # to be in the xi-eta coordinate system
    # this is still an outflow for nrm = [1, 0]
    nrm = [1., 0]  # flow is in the xi direction
    angle_diff = theta + flow_direction
    alpha_eff = alpha_mag*cos(angle_diff)  # effective alpha in the wall normal
                                           # direction
    val_exp = alpha_eff*u
    val = AdvectionEquationMod.RoeSolver(u, u_bc, alpha_x, alpha_y, nrm, dxidx)
    println("val = ", val)
    println("val_exp = ", val_exp)

    @fact val => roughly(val_exp, atol=1e-14)

    # now check eta direction
    nrm = [0, 1.]
    alpha_eff = alpha_mag*sin(angle_diff)
    val_exp = u*alpha_eff
    val = AdvectionEquationMod.RoeSolver(u, u_bc, alpha_x, alpha_y, nrm, dxidx)

    @fact val => roughly(val_exp, atol=1e-14)

    # now rotate the coordinate system so much that this becomes an inflow
    theta = 120*pi/180
    nrm = [1., 0]
    angle_diff = theta + flow_direction
    alpha_eff = alpha_mag*cos(angle_diff)
    val_exp = alpha_eff*u_bc

    dxidx = get_rotation_matrix( theta)
    val = AdvectionEquationMod.RoeSolver(u, u_bc, alpha_x, alpha_y, nrm, dxidx)
    println("val = ", val)
    println("val_exp = ", val_exp)
    @fact val => roughly(val_exp, atol=1e-14)

    # check eta direction
    nrm = [0, 1]
    angle_diff = theta + flow_direction
    alpha_eff = alpha_mag*sin(angle_diff)
    val_exp = alpha_eff*u
    val = AdvectionEquationMod.RoeSolver(u, u_bc, alpha_x, alpha_y, nrm, dxidx)
    println("val = ", val)
    println("val_exp = ", val_exp)
    @fact val => roughly(val_exp, atol=1e-14)


 end

  context("--- Testing Volume Integrals ---")  do

    # use the 8 element mesh
    ARGS[1] = "input_vals_8el.jl"
    include(STARTUP_PATH)

    fill!(eqn.q, 0.0)
    alpha_x = ones(Float64 ,1, mesh.numNodesPerElement, mesh.numEl)
    alpha_y = zeros(alpha_x)
    println("typeof(alpha_x) = ", typeof(alpha_x))

    fill!(eqn.res, 0.0)
    AdvectionEquationMod.evalSCResidual(mesh, sbp, eqn, alpha_x, alpha_y)
    println("eqn.res = ", eqn.res)
    eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
    println("eqn.res_vec = \n", eqn.res_vec)
    @fact eqn.res => roughly(zeros(1, mesh.numNodesPerElement, mesh.numEl), atol=1e-14)


    # check that Qx.'*q = 0 when q = ones(3)
    fill!(eqn.q, 1.0)
    fill!(eqn.res, 0.0)
    dxidx1 = mesh.dxidx[:, :, :, 1]
    println("dxidx for element 1 = \n", dxidx1)
    println("|J| for element 1 = \n", mesh.jac[:, 1])
    println("coords for element 1 = \n", mesh.coords[:, :, 1])
    jac = mesh.jac[:, 1]
    dxidx_true = zeros(2, 2, 3)
    for i=1:3
      dxidx_true[:, :, i] = dxidx1[:, :, i]*jac[i]
    end

    Qx = sbp.Q[:, :, 1]*dxidx1[1, 1, 1] + sbp.Q[:, :, 2]*dxidx1[2, 1, 1]
    q = [eqn.q[1, 1, 1], eqn.q[1, 2, 1], eqn.q[1, 3, 1]]  # extract q values
    # sum reduces the vector to the value of the integral
    val_test = sum(Qx.'*q)  
    println("val_test = ", val_test)
    AdvectionEquationMod.evalSCResidual(mesh, sbp, eqn, alpha_x, alpha_y)
    # extract residual values
    vec_code =  [eqn.res[1, 1, 1], eqn.res[1, 2, 1], eqn.res[1, 3, 1]]
    val_code = sum(vec_code)
    println("val_code = ", val_code)
    xi_flux = dxidx1[1,1,1]*1*q
    xi_component = sbp.Q[:, :, 1].'*xi_flux
    eta_flux = dxidx1[2, 1]*1*q
    eta_component = sbp.Q[:, :, 2].'*eta_flux
    val2_test = sum(xi_component + eta_component)
    println("test xi_flux = ", xi_flux)
    println("test xi_component = ", xi_component)
    println("test eta_flux = ", eta_flux)
    println("test eta_component = ", eta_component)
    @fact val_test => roughly(val2_test, atol=1e-14)
    @fact val_test => roughly(val_code, atol=1e-14)


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

    println("eqn.q[:, :, 1] = ", eqn.q[:, :, 1])
    fill!(eqn.res, 0.0)
    fill!(eqn.alpha_x, 1.0)
    fill!(eqn.alpha_y, 0.0)
    # check the boundry contribution
    AdvectionEquationMod.evalBndry(mesh, sbp, eqn, alpha_x, alpha_y)

    println("eqn.res[:, :, 1] = ", eqn.res[:, :, 1])
    @fact sum(eqn.res[:, :, 1]) => roughly(-2.0, atol=1e-14)
    @fact sum(eqn.res[:, :, 3]) => roughly(-2.0, atol=1e-14)

    println("eqn.res[:, :, 6] = ", eqn.res[:, :, 6])
    @fact sum(eqn.res[:, :, 6]) => roughly(-2.0, atol=1e-14)
    @fact sum(eqn.res[:, :, 8]) => roughly(-2.0, atol=1e-14)

    fill!(eqn.res, 0.0)
    AdvectionEquationMod.evalSCResidual(mesh, sbp, eqn, alpha_x, alpha_y)
    for i=1:mesh.numEl
      println("element $i")
      Qx_i = sbp.Q[:, :, 1]*mesh.dxidx[1, 1, 1, i] + sbp.Q[:, :, 2]*mesh.dxidx[2, 1, 1, i]
      q_i = reshape(eqn.q[1, :, i], 3)
      val_test = Qx_i.'*q_i
      val_code = reshape(eqn.res[:, :, i], 3)
#      println("val_code = ", val_code)
      @fact sum(val_code) => roughly(0.0, atol=1e-14)  # proven by hand calc
      println("val_test = ", val_test)
      println("val_code = ", val_code)
      @fact val_code => roughly(val_test, atol=1e-14)
    end

    println("mesh.dofs = \n", mesh.dofs)
    eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
    val_test =[-1., 2/3, 1, -2, 1/3, -1, 1, 2/3, 1/3]
    println("eqn.res_vec = \n", eqn.res_vec)
    @fact eqn.res_vec => roughly(val_test, atol=1e-14)

    println("----- Checking q=2*x^2 + 5 case -----")
    
    x1 = zeros(Tmsh, 3)
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        x = mesh.coords[1, j, i]
        eqn.q[1, j, i] = 2*x*x + 5
      end
    end

    println("eqn.q[:, :, 1] = ", eqn.q[:, :, 1])
    fill!(eqn.res, 0.0)
    fill!(eqn.alpha_x, 1.0)
    fill!(eqn.alpha_y, 0.0)
    AdvectionEquationMod.evalSCResidual(mesh, sbp, eqn, alpha_x, alpha_y)
    for i=1:mesh.numEl
      println("element $i")
      Qx_i = sbp.Q[:, :, 1]*mesh.dxidx[1, 1, 1, i] + sbp.Q[:, :, 2]*mesh.dxidx[2, 1, 1, i]
      q_i = reshape(eqn.q[1, :, i], 3)
      val_test = Qx_i.'*q_i
      val_code = reshape(eqn.res[:, :, i], 3)
#      println("val_code = ", val_code)
      println("val_test = ", val_test)
      println("val_code = ", val_code)
      @fact val_code => roughly(val_test, atol=1e-14)
    end


    ARGS[1] = "input_vals_8el_large.jl"
    include(STARTUP_PATH)
    println("----- Checking q=2*x^2 + 5 case on large grid -----")
    
    x1 = zeros(Tmsh, 3)
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        x = mesh.coords[1, j, i]
        eqn.q[1, j, i] = 2*x*x + 5
      end
    end

    println("eqn.q[:, :, 1] = ", eqn.q[:, :, 1])
    fill!(eqn.res, 0.0)
    fill!(eqn.alpha_x, 1.0)
    fill!(eqn.alpha_y, 0.0)
    AdvectionEquationMod.evalSCResidual(mesh, sbp, eqn, alpha_x, alpha_y)
    for i=1:mesh.numEl
      println("element $i")
      Qx_i = sbp.Q[:, :, 1]*mesh.dxidx[1, 1, 1, i] + sbp.Q[:, :, 2]*mesh.dxidx[2, 1, 1, i]
      q_i = reshape(eqn.q[1, :, i], 3)
      val_test = Qx_i.'*q_i
      val_code = reshape(eqn.res[:, :, i], 3)
#      println("val_code = ", val_code)
      println("val_test = ", val_test)
      println("val_code = ", val_code)
      @fact val_code => roughly(val_test, atol=1e-14)
    end

    # back to the original mesh
    ARGS[1] = "input_vals_8el.jl"
    include(STARTUP_PATH)

    println("----- Checking sinwave case -----")
    x1 = zeros(Tmsh, 3)
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        x = mesh.coords[1, j, i]
        eqn.q[1, j, i] = AdvectionEquationMod.calc_sinwave(mesh.coords[:, j, i], 0.25)
      end
    end
    println("eqn.q = \n", eqn.q)

    mesh.bndry_funcs[1] = AdvectionEquationMod.sinwave_BC()
    println("eqn.q[:, :, 1] = ", eqn.q[:, :, 1])
    fill!(eqn.res, 0.0)
    fill!(eqn.alpha_x, 1.0)
    fill!(eqn.alpha_y, 0.0)
    AdvectionEquationMod.evalSCResidual(mesh, sbp, eqn, alpha_x, alpha_y)
    for i=1:mesh.numEl
      println("element $i")
      println("dxidx = ", mesh.dxidx[:, :, 1, i])
      Qx_i = sbp.Q[:, :, 1]*mesh.dxidx[1, 1, 1, i] + sbp.Q[:, :, 2]*mesh.dxidx[2, 1, 1, i]
      q_i = reshape(eqn.q[1, :, i], 3)
      val_test = Qx_i.'*q_i
      val_code = reshape(eqn.res[:, :, i], 3)
#      println("val_code = ", val_code)
      println("  val_test = ", val_test)
      println("  val_code = ", val_code)
      @fact val_code => roughly(val_test, atol=1e-14)
    end

    println("eqn.res = \n", eqn.res)
    eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
    println("eqn.res_vec = \n", eqn.res_vec)

    ARGS[1] = "input_vals_channel_verylarge.jl"
    include(STARTUP_PATH)

    println("----- Checking sinwave case on very large mesh -----")
    x1 = zeros(Tmsh, 3)
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        x = mesh.coords[1, j, i]
        eqn.q[1, j, i] = AdvectionEquationMod.calc_sinwave(mesh.coords[:, j, i], 0.25)
      end
    end
    println("eqn.q = \n", eqn.q)

    mesh.bndry_funcs[1] = AdvectionEquationMod.sinwave_BC()
    println("eqn.q[:, :, 1] = ", eqn.q[:, :, 1])
    fill!(eqn.res, 0.0)
    fill!(eqn.alpha_x, 1.0)
    fill!(eqn.alpha_y, 0.0)
    AdvectionEquationMod.evalSCResidual(mesh, sbp, eqn, alpha_x, alpha_y)
    for i=1:mesh.numEl
      println("element $i")
      Qx_i = sbp.Q[:, :, 1]*mesh.dxidx[1, 1, 1, i] + sbp.Q[:, :, 2]*mesh.dxidx[2, 1, 1, i]
      q_i = reshape(eqn.q[1, :, i], 3)
      val_test = Qx_i.'*q_i
      val_code = reshape(eqn.res[:, :, i], 3)
#      println("val_code = ", val_code)
      println("val_test = ", val_test)
      println("val_code = ", val_code)
      @fact val_code => roughly(val_test, atol=1e-14)
    end




  end




#=
 context("--- Testing dataPrep ---") do
 
  end



  context("--- Testing evalBoundaryIntegrals ---") do
  end

  context("--- Testing evalEuler --- ")  do

end # end facts block
=#
