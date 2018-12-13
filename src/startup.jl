# this script invokes the solver for any physics

using PDESolver
using EulerEquationMod
using AdvectionEquationMod
using SimpleODEMod
using EllipticEquationMod

mesh, sbp, eqn, opts = run_solver(ARGS[1])
