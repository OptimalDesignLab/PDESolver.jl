# functions for a homotopy that goes from one shock sensor to another

"""
  Data for doing homotopy from one shock capturing scheme to another.
  The homotopy function is:

  R_euler (1-lambda)*R_hard + lambda*R_easy

  where R_euler is the physics residual, excluding shock capturing, R_hard
  is the shock sensor that is difficult to converge, and R_easy is the shock
  sensor that is easy to converge.  It is recommended that the initial
  condition satisfy R_euler + R_easy.
"""
mutable struct SCHomotopy
  sensor_easy::AbstractShockSensor
  opts_orig::Dict{String, Bool}

  function SCHomotopy(mesh, sbp, eqn, opts)

    sensor_easy = createShockSensor(mesh, sbp, opts, opts["homotopy_shock_sensor"])
    opts_orig = Dict{String, Bool}

    return new(sensor_easy, opts_orig)
  end
end


const TermNames = ["addVolumeIntegrals", "addBoundaryIntegrals", "addStabilization", "addFaceIntegrals"]

"""
  Disable all terms of the residual evaluation except shock capturing

  **Inputs/Outputs**

   * obj: `HomotopyData`
   * opts

"""
function disableTerms(obj::HomotopyData, opts)

  for key in TermNames
    obj.opts_orig[key] = opts[key]
    opts[key] = false
  end

  return nothing
end


"""
  Reset the options dictionary to its original state

  **Inputs**

   * obj: `HomotopyData`
   * opts: the options dictionary
"""
function enableTerms(obj::HomotopyData, opts)

  for (key, val) in obj.opts_orig
    opts[key] = val
  end
end


"""
  Evaluates the homotopy function R_euler + (1-lambda)*R_hard + lambda*R_easy

  This function assumes parallel communication has already started
"""
function evalHomotopyFunction_sc(obj::HomotopyData, mesh, sbp, eqn, opts, t)
# this probably doesn't make sense if physics_func != evalResidual


  lambda = obj.lambda 

  # evaluate R_euler + (1-lambda)R_hard and store it in res_homotopy
  setShockSensorAlpha(eqn, 1-lambda)
  obj.physics_func(mesh, sbp, eqn, opts, t)
  setShockSensorAlpha(eqn, 1.0)  # just in case
  res_homotopy = copy(eqn.res); fill!(eqn.res, 0)

  # evaluate lambda * R_easy and store it in eqn.res
  sensor_orig = getShockSensor(eqn)
  sensor_easy = obj.sensor_easy
  setShockSensor(eqn, obj.sensor_easy)
  setShockSensorAlpha(eqn, lambda)
  disableTerms(obj, opts)
  obj.physics_func(mesh, sbp, eqn, opts, t)

  # reset things
  setShockSensorAlpha(eqn, 1)
  enableTerms(obj, opts)
  setShockSensor(eqn, sensor_orig)


  # combine res_homotopy and eqn.res
  @simd for i=1:length(eqn.res)
    eqn.res[i] += res_homotopy[i]
  end

  return nothing
end


"""
  Evaluates the derivative of the homotopy function with respect to lambda
"""
function calcdHdLambda_sc(obj::HomotopyData, mesh, sbp, eqn, opts, t, res_vec)


  fill!(eqn.res, 0)
  res_homotopy = copy(eqn.res)

  disableTerms(obj, opts)

  # evaluate -1*R_hard
  setShockSensorAlpha(eqn, -1)
  obj.physics_func(mesh, sbp, eqn, opts, t)
  copy!(res_homotopy, eqn.res)
  setShockSensorAlpha(eqn, 1)
  sensor_orig = getShockSensor(eqn)

  # evaluate 1*R_easy
  setShockSensor(eqn, obj.sensor_easy)
  setShockSensorAlpha(eqn, 1)
  obj.physics_func(mesh, sbp, eqn, opts, t)

  @simd for i=1:length(eqn.res)
    eqn.res[i] += res_homotopy[i]
  end

  array3DTo1D(mesh, sbp, eqn, opts, eqn.res, res_vec)

  # reset things
  setShockSensor(eqn, sensor_orig)
  enableTerms(obj, opts)

  return nothing
end


"""
  Computes Jacobian of homotopy function with respect to q
"""
function evalJacobian_homotopy_sc(lo::HomotopyMatLO, mesh, sbp, eqn, opts,
                                  ctx_residual, t)

  hdata = lo.hdata
  lambda = hdata.lambda

  # evaluate jacobian of R_euler + (1-lambda)*R_hard
  setShockSensorAlpha(eqn, 1-lambda)
  calcLinearOperator(lo.lo_inner, mesh, sbp, eqn, opts, ctx_residual, t)
  setShockSensorAlpha(eqn, 1)

  # evaluate Jacobian of lambda*R_easy
  disableTerms(hdata, opts)
  sensor_orig = getShockSensor(eqn)
  setShockSensor(eqn, hdata.sensor_easy)
  setShockSensorAlpha(eqn, lambda)

  A = getBaseLO(lo).A
  physicsJac(mesh, sbp, eqn, opts, A, ctx_residual, t, zero_jac=false)

  # reset things
  setShockSensorAlpha(eqn, 1)
  setShockSensor(eqn, sensor_orig)
  enableTerms(hdata, opts)

  return nothing
end
