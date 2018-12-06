

"""
  calcDragTimeAverage:
    Calculates the time average of drag. Reads this data from 'drag.dat' in the current directory.

  Input:
    mesh, sbp, eqn, opts: standard
    delta_t: Time step size
    itermax_fromnlsolver: The last iteration's number from the NL Solver

  Outputs: (return values)
    Cd:     Coefficient of drag
            Cd = <D>/(0.5*M^2)
    dCddM:  Derivative of the coefficient of drag. This derivative is in a partial sense.
            This is calculated by taking the partial derivative of Cd wrt M:
            dCd/dM = (-2<D>)/(0.5*M^3)

"""
function calcDragTimeAverage(mesh, sbp, eqn, opts, delta_t, itermax_fromnlsolver)

  # println(BSTDOUT, "------------------------------- in calcDragTimeAverage")

  dt = delta_t

  Ma = eqn.params.Ma
  data = readdlm("drag.dat")

  itermax_fromdata = size(data, 1)
  if itermax_fromdata > (itermax_fromnlsolver + 1)
    println(BSTDOUT, " itermax_fromdata: ", itermax_fromdata)
    println(BSTDOUT, " itermax_fromnlsolver + 1: ", itermax_fromnlsolver + 1)
    println(BSTDOUT, "--------------------------------")
    flush(BSTDOUT)
    error("HELLOOOO You forgot to delete drag.dat, or there's some problem with finaliter")
  end

  itermax = itermax_fromdata

  iter = round.(Int64, data[1:itermax, 1])
  drag = data[1:itermax, 2]

  # iter = iter - 1     # because iter starts at 2      ---- Now commented out bc of IC inclusion

  drag_timeavg = 0.0
  maxtime = dt*itermax - dt        # needs to have the minus dt here, because the IC doesn't count as its own time step

  # trapezoid rule
  for i = 1:itermax
    quad_weight = calcQuadWeight(i, delta_t, itermax)
    drag_timeavg += quad_weight * drag[i]
  end

  drag_timeavg = drag_timeavg * 1.0/maxtime

  # Cd calculations
  Cd = drag_timeavg/(0.5*Ma^2)
  println(" Cd = <D>/(0.5*M^2) = ", Cd)

  dCddM = (-4.0*drag_timeavg)/(Ma^3)
  # comes from dCd/dM = (-2<D>)/(0.5*M^3)
  println(" dCddM = (-4<D>)/(M^3) = ", dCddM)

  return Cd, dCddM

end     # end function calcDragTimeAverage

"""
  calcFinalIter:
    Simple function to decide which is actually the final iteration's number,
      because there are two ways that is determined.

  Input:
    t_steps:  max number of time steps, set by tmax
    itermax:  max number of time steps, set by user option itermax

  Output: (return value)
    finaliter: the actual number of the final iteration

"""
function calcFinalIter(t_steps, itermax)

  finaliter_setby_tmax = (t_steps + 1)
  finaliter_setby_itermax = (itermax + 1)
  if finaliter_setby_tmax <= finaliter_setby_itermax
    finaliter = finaliter_setby_tmax
  else
    finaliter = finaliter_setby_itermax
  end

  return finaliter

end     # end function calcFinalIter

"""
  calcQuadWeight:
    Determines the quadrature weight based upon the trapezoid rule. Handles the
    special case if # time steps is less than 3.
  
  Inputs:
    i:        iter number to calculate quadrature weight at
    delta_t:  time step size
    itermax:  max number of time steps. This version is agnostic to where that number comes from

  Output: (return value)
    quad_weight: The quadrature weight for this time step. This should be applied as a 
                 factor if, for example, you are forming a sum over all time steps
                 and adding into a field each time step.
"""
function calcQuadWeight(i, delta_t, finaliter)

  if (i == 1 || i == finaliter)
    quad_weight = delta_t/2.0             # first & last time step, trapezoid rule quadrature weight
  else
    quad_weight = delta_t                 # all other timesteps
  end

  if finaliter < 2        # if 1 or 2 timesteps, shift to regular rectangular rule. 
                          # this check is against 2, not 3, because the IC is not counted in this sequence of i's
    quad_weight = delta_t/2.0
  end

  return quad_weight
end     # end function calcQuadWeight

