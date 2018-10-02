
function test_bcs()

  fname = "input_vals_channel_bc.jl"
  mesh, sbp, eqn, opts = run_solver(fname)
  params = eqn.params
  gamma = params.gamma
  gamma_1 = params.gamma_1

  #----------------------------------------------------------------------------
  # Subsonic Inflow
  # construct a state such that the BC function will compute the boundary state
  # equal to the input state

  Ma = 0.5
  # these must match the values in the BC function
  pt = 102010.0/params.p_free
  Tt = 288.6/params.T_free

  println("pt = ", pt)
  println("Tt = ", Tt)

#  Ub = 1.0  #TODO: make sure this is subsonic
#  a = Ub/Ma

  operand = 1 + 0.5*gamma_1*Ma*Ma
  p = pt*(operand^(-gamma/gamma_1))
  T = Tt*(operand^(-1))

  rho = p/(params.R_ND*T)
  a = sqrt(gamma*p/rho)
  Ub = Ma*a
  E = p/params.gamma_1 + 0.5*rho*Ub*Ub

  q = zeros(mesh.numDofPerNode)
  q[1] = rho
  q[2] = 0
  q[3] = rho*Ub
  q[4] = E

  F = zeros(mesh.numDofPerNode)
  F2 = zeros(mesh.numDofPerNode)
  aux_vars = Float64[EulerEquationMod.calcPressure(params, q)]
  coords = mesh.coords_interface[:, 1, 1]
#  nrm_xy = mesh.nrm_face[:, 1, 1]
  nrm_xy = -Float64[0.0, 2.0]

  EulerEquationMod.calcEulerFlux(params, q, aux_vars, nrm_xy, F)
  obj = EulerEquationMod.BCDict["subsonicInflowBC"](mesh, eqn)
  obj(params, q, aux_vars, coords, nrm_xy, F2)

  @test isapprox( norm((F - F2)/norm(q)), 0.0) atol=1e-12


  #----------------------------------------------------------------------------
  # Subsonic Outflow

  p = 101300.0/params.p_free # must match the value in the BC function
  q[1] = 2
  q[2] = 0.2
  q[3] = 0.2
  q[4] = p/gamma_1 + 0.5*( q[2]*q[2] + q[3]*q[3])/q[1]

  nrm_xy = Float64[0.0, 2.0]

  EulerEquationMod.calcEulerFlux(params, q, aux_vars, nrm_xy, F)
  obj = EulerEquationMod.BCDict["subsonicOutflowBC"](mesh, eqn)
  obj(params, q, aux_vars, coords, nrm_xy, F2)

  @test isapprox( norm((F - F2)/norm(q)), 0.0) atol=1e-12




  return nothing
end

add_func1!(EulerTests, test_bcs, [TAG_SHORTTEST])

