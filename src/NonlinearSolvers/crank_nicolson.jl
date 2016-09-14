# crank_nicolson.jl
# Crank-Nicolson implicit solver for PDEs

export crank_nicolson

@doc """
crank_nicolson

  This function performs Crank-Nicolson implicit time solution, using a function 
  of the form du/dt = f(u, t)
  

"""->
function crank_nicolson(f::Function, h::AbstractFloat, t_max::AbstractFloat,
                        q_vec::AbstractVector, res_vec::AbstractVector, pre_func,
                        post_func, ctx, opts)

  myrank = MPI.Comm_rank(MPI.COMM_WORLD)
  fstdout = BufferedIO(STDOUT)
  if myrank == 0
    println(fstdout, "\nEntered Crank-Nicolson")
    println(fstdout, "res_tol = ", res_tol)
  end

  for i = 2:(t_steps + 1)

    # assign q_i-1 and q_i
    # form R vector, then enter newton
    # Newton somehow needs Jac and r' inv
    # use newton to drive norm R below newton_tol

  end

  #returns t?

end


