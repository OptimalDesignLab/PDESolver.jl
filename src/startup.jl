# this script invokes the solver for any physics

using PDESolver

mesh, sbp, eqn, opts = run_solver(ARGS[1])
