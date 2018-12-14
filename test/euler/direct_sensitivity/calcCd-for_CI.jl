# Calculates the drag coefficient for a completed PDESolver run.
# This is done by reading in drag values from drag.dat and using the
#   trapezoid rule to time-average them.
#
# Input:
#   dir:    directory (relative) of the completed run
#
# Note: must have 'drag.dat', 'Ma.dat', and 'delta_t.dat' present
#
function calcCd(dir)

  file = string(dir, "/", "drag.dat")
  data = readdlm(file)
  println(" reading drag from ", file)

  file = string(dir, "/", "Ma.dat")
  Ma_arr = readdlm(file)
  Ma = Ma_arr[1,1]
  file = string(dir, "/", "delta_t.dat")
  dt_arr = readdlm(file)
  dt = dt_arr[1,1]

  maxiter = size(data,1)
  # println(" maxiter: ", maxiter)
  # maxiter = 20

  iter = round.(Int64, data[1:maxiter, 1])

  if iter[end] < (maxiter - 2)
    error("drag data appears to be too long; you might have forgotten to delete drag.dat")
  end

  drag = data[1:maxiter, 2]

  drag_timeavg = 0.0
  maxtime = dt*maxiter - dt       # needs to have the minus dt here, because the IC doesn't count as its own time step

  # println(" maxiter: ", maxiter)
  # println(" maxtime: ", maxtime)
  # println(" ")

  for i = 1:maxiter
    if (i == 1 || i == maxiter)
      quad_weight = dt/2.0
    else
      quad_weight = dt
    end
    if maxiter < 3        # if 1 or 2 timesteps, shift to regular rectangular rule
      quad_weight = dt/2.0
      println(" small maxiter; quad_weight = dt/2")
    end

    # println(" i: $i   quad_weight: $quad_weight   drag[i]: ", drag[i])
    drag_timeavg += quad_weight * drag[i]

  end

  drag_timeavg = drag_timeavg * 1.0/maxtime

  # println(" --- for dir: ", dir, " ---")
  # println(" drag_timeavg: ", drag_timeavg)

  Cd = drag_timeavg/(0.5*Ma^2)
  # println(" Cd = <D>/(0.5*M^2) = ", Cd)

  # D = drag_timeavg
  # println(" D = <D> = ", D)

  #--------- don't do this for FD.
  # dCddM = (-2.0*drag_timeavg)/(0.5*Ma^3)
  # println(" dCddM = (-2<D>)/(0.5*M^3) = ", dCddM)

  # println(" ")

  return Ma, Cd

end

export calcCd
