# script for invoking the solver using the input file specified on the
# command line

using PDESolver
using EulerEquationMod

mesh, sbp, eqn, opts = run_solver(ARGS[1])
