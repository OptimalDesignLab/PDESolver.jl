# tests for shock capturing

import EulerEquationMod: AbstractShockSensor, AbstractShockCapturing

function test_shocksensorPP()

  @testset "Shock sensor PP" begin
    opts = read_input_file("input_vals_jac2d.jl")
    opts["order"] = 2
    delete!(opts, "calc_jac_explicit")
    opts["force_solution_complex"] = true
    mesh, sbp, eqn, opts = solvePDE(opts)

    Tsol = eltype(eqn.q); Tres = eltype(eqn.res)
    q = eqn.q[:, :, 1]
    jac = ones(Float64, mesh.numNodesPerElement)
    res = zeros(eltype(eqn.res), mesh.numDofPerNode, mesh.numNodesPerElement)

    sensor = eqn.params.sensor_pp
    capture = eqn.params.projection_shock_capturing
    # initial condition is constant, check the sensor reports no shock
    Se, ee = EulerEquationMod.getShockSensor(eqn.params, sbp, sensor, q, jac)

    @test abs(Se) < 1e-12
    @test ee == 0

    fill!(res, 0)
    EulerEquationMod.applyShockCapturing(eqn.params, sbp, sensor, capture, q, jac, res)
    @test maximum(abs.(res)) < 1e-13

    # test when a shock is present
    q[1, 3] += 5
    Se, ee = EulerEquationMod.getShockSensor(eqn.params, sbp, sensor, q, jac)

    @test abs(Se) > 1e-12
    @test ee > 0.99

    fill!(res, 0)
    w = copy(q)
    for i=1:mesh.numNodesPerElement
      w_i = sview(w, :, i)
      q_i = sview(q, :, i)
      EulerEquationMod.convertToIR(eqn.params, q_i, w_i)
    end
    EulerEquationMod.applyShockCapturing(eqn.params, sbp, sensor, capture, q, jac, res)

    @test sum(res .* w) < 0  # the term is negative definite

    # case 3: ee = 1
    test_shocksensor_diff(eqn.params, sbp, sensor, q, jac)
    test_shockcapturing_diff(eqn.params, sbp, sensor, capture, q, jac)

    # case 2: ee on sin wave
    q[1, 3] = 1.0105
    for i=1:mesh.numNodesPerElement
      for j=2:mesh.numDofPerNode
        q[j, i] += 0.1*(i + j)
      end
    end

    test_shocksensor_diff(eqn.params, sbp, sensor, q, jac)
    test_shockcapturing_diff(eqn.params, sbp, sensor, capture, q, jac)

  end  # end testset


  return nothing
end


"""
  Tests derivative of the shock sensor at a given state
"""
function test_shocksensor_diff(params, sbp, sensor::AbstractShockSensor, _q, jac)

  numDofPerNode, numNodesPerElement = size(_q)
  q = zeros(Complex128, numDofPerNode, numNodesPerElement)
  copy!(q, _q)
  Se_jac = zeros(q); Se_jac2 = zeros(q)
  ee_jac = zeros(q); ee_jac2 = zeros(q)

  h = 1e-20
  pert = Complex128(0, h)
  for i=1:numNodesPerElement
    for j=1:numDofPerNode
      q[j, i] += pert
      Se, ee = EulerEquationMod.getShockSensor(params, sbp, sensor, q, jac)
      Se_jac[j, i] = imag(Se)/h
      ee_jac[j, i] = imag(ee)/h
      q[j, i] -= pert
    end
  end

  EulerEquationMod.getShockSensor_diff(params, sbp, sensor, q, jac,
                                       Se_jac2, ee_jac2)

  @test maximum(abs.(Se_jac - Se_jac2)) < 1e-12
  @test maximum(abs.(ee_jac - ee_jac2)) < 1e-12

  # test vector mode
  q_dot = rand_realpart(size(q))
  q .+= pert*q_dot
  Se, ee = EulerEquationMod.getShockSensor(params, sbp, sensor, q, jac)
  Se_dot = imag(Se)/h
  ee_dot = imag(ee)/h
  q .-= pert*q_dot

  # run again to make sure intermediate arrays are zeroed out
  EulerEquationMod.getShockSensor_diff(params, sbp, sensor, q, jac,
                                       Se_jac2, ee_jac2)

  Se_dot2 = sum(Se_jac2 .* q_dot)
  ee_dot2 = sum(ee_jac2 .* q_dot)

  @test abs(Se_dot - Se_dot2) < 1e-12
  @test abs(ee_dot - ee_dot2) < 1e-12

  
  return nothing
end

#TODO: add abstract types
"""
  Tests Jacobian of shock capturing scheme
"""
function test_shockcapturing_diff(params, sbp, sensor::AbstractShockSensor,
                                  capture::AbstractShockCapturing,
                                   _q, jac)

  numDofPerNode, numNodesPerElement = size(_q)
  q = zeros(Complex128, size(_q))
  copy!(q, _q)
  q_dot = rand_realpart(size(q))
  #q_dot = zeros(numDofPerNode, numNodesPerElement)
  #q_dot[:, 5] = 2
  res = zeros(Complex128, size(q))
  
  h = 1e-20
  pert = Complex128(0, h)

  # complex step
  q .+= pert*q_dot
  EulerEquationMod.applyShockCapturing(params, sbp, sensor, capture, q, jac, res)
  res_dot = imag(res)./h
  q .-= pert*q_dot


  # AD
  res_jac = zeros(numDofPerNode, numDofPerNode, numNodesPerElement, numNodesPerElement)
  EulerEquationMod.applyShockCapturing_diff(params, sbp, sensor, capture, q,
                                            jac, res_jac)

  res_dot2 = zeros(Complex128, numDofPerNode, numNodesPerElement)
  for i=1:numNodesPerElement
    for j=1:numDofPerNode
      res_dot2 .+= res_jac[:, j, :, i] * q_dot[j, i]
    end
  end

  @test maximum(abs.(res_dot - res_dot2)) < 1e-12

  return nothing
end

add_func1!(EulerTests, test_shocksensorPP, [TAG_SHORTTEST, TAG_TMP])
