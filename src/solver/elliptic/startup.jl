# script for invoking the solver using the input file specified on the
# command line

using PDESolver
using EllipticEquationMod
using EllipticEquationMod.run_elliptic

mesh, sbp, eqn, opts = run_elliptic(ARGS[1])
