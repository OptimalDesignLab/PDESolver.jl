
# function loadDragTimeAverage(mesh, sbp, eqn, opts, useArray::Bool, filename="drag.dat")

  # if useArray == false
  # if useArray == true


"""
  calcDragTimeAverage:
    Calculates the time average of drag. Reads this data from 'drag.dat' in the current directory.

  Input:
    mesh, sbp, eqn, opts: standard
    delta_t: Time step size
    itermax_fromnlsolver: The last iteration's number from the NL Solver
    useArray: Bool. If true, calculate the drag data from an array of values. 
              This array should be one element per time step's drag value.
              If false, calculate the drag data from drag.dat.

  Outputs: (return values)
    Cd:     Coefficient of drag
            Cd = <D>/(0.5*M^2)
    dCddM:  Derivative of the coefficient of drag. This derivative is in a partial sense.
            This is calculated by taking the partial derivative of Cd wrt M:
            dCd/dM = (-2<D>)/(0.5*M^3)

"""
function calcDragTimeAverage(mesh, sbp, eqn, opts, delta_t, itermax_fromnlsolver; useArray=false, drag_array=[0.0])

  # println(BSTDOUT, "------------------------------- in calcDragTimeAverage")
  myrank = mesh.myrank

  dt = delta_t

  Ma = eqn.params.Ma
  # data = readdlm("drag.dat")

  #=
  if opts["is_restart"]                 # if it's a restart, need to use log_cleaner, which saves the 
                                        # cleaned output in filename_cleaned.dat
    drag_filename = "drag_cleaned.dat"
  else
    drag_filename = "drag.dat"
  end
  =#

  if useArray == false
    println(BSTDOUT, "calculating drag from file")
    data = readdlm("drag.dat")
    itermax_fromdata = size(data, 1)
    itermax = itermax_fromdata
    drag = data[1:itermax, 2]
  else
    println(BSTDOUT, "calculating drag from array")
    if length(drag_array) == 1
      error(" drag_array not specified")
    end
    drag = drag_array
    itermax_fromdata = size(drag, 1)
    itermax = itermax_fromdata
  end


  println(BSTDOUT, " ----- Reading drag.dat -----")
  println(BSTDOUT, "  itermax_fromdata: ", itermax_fromdata)
  println(BSTDOUT, "  itermax_fromnlsolver: ", itermax_fromnlsolver)
  println(BSTDOUT, " ----------------------------")

  if itermax_fromdata > (itermax_fromnlsolver + 1)
    println(BSTDOUT, " itermax_fromdata: ", itermax_fromdata)
    println(BSTDOUT, " itermax_fromnlsolver + 1: ", itermax_fromnlsolver + 1)
    println(BSTDOUT, "--------------------------------")
    flush(BSTDOUT)
    println(BSTDOUT, "HELLOOOO Problem with length of drag data compared to what it should be.")
    println(BSTDOUT, "         Could be that duplicate entries exist due to restarting,")
    println(BSTDOUT, "         you forgot to delete drag.dat, or there's some problem with finaliter")
    println(BSTDOUT, "         Remember that we're appending to drag.dat now.")
    flush(BSTDOUT)
  end

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
  @mpi_master println(" in drag_timeaverage: Cd = <D>/(0.5*M^2) = ", Cd)

  dCddM = (-4.0*drag_timeavg)/(Ma^3)
  # comes from dCd/dM = (-2<D>)/(0.5*M^3)
  @mpi_master println(" in drag_timeaverage: dCddM = (-4<D>)/(M^3) = ", dCddM)

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

