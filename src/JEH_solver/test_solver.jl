module TestSolver

include("solver.jl")

using SummationByParts
using .Solver

export runtest

function runtest(;degree::Int=1, Nx::Int=2, Ny::Int=2, cfl::Float64=0.5,
                 output::Bool=false, energyhist::Bool=false, matrix::Bool=false)
  solver = Solver.EulerSolver{Float64}(degree, Nx, Ny, edgestab=true)
  u = ones((4, solver.mesh.numnodes))
  u[2,:] = 0.5
  u[3,:] = 0.5

#  println("coords = ", solver.mesh.x)

  (numrow, numcol) = size(solver.mesh.x)
  println("size(mesh.x) = ", size(solver.mesh.x))
  for i=1:numcol
    println(i, "| ", solver.mesh.x[:,i])
  end

  dt = cfl # approximation!

  # set the initial condition to be the exact solution
  for i = 1:solver.mesh.numnodes
    tmp = zeros(4)
    Solver.calcIsentropicVortex!(solver.mesh.x[:,i], tmp)
    u[:,i] = tmp
  end

  resnorm = Solver.timemarchRK!(solver, u, dt)  
  iter = 0
  while resnorm > 1e-5
    resnorm = Solver.timemarchRK!(solver, u, dt)
    if mod(iter, 1) == 0
      println("iter ",iter,": res. norm = ", resnorm)
    end
    iter += 1
  end

  L2err, maxerr = Solver.calcerror(solver, u)
  println("L2err = ",L2err)

  if output
    # output solution
    Solver.printsolution(solver, u, "euler.dat")
  end

end

end
