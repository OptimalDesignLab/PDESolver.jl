using Documenter, PDESolver, ODLCommonTools

# need to use the physics modules here too
using AdvectionEquationMod
using EulerEquationMod
using SimpleODEMod
include("../test/TestSystem.jl")
using TestSystem

# some abbreviations
advec = "solver/advection"
euler = "solver/euler"
ode = "solver/simpleODE"

makedocs(
  format = :html,
  sitename = "PDESolver.jl",
  pages = Any["PDESolver Introduction" => "index.md"
              "PDESolver Concepts" => Any[
                    "Intro" => "concepts/intro.md"
                    "PUMI" => "concepts/pumi.md"
                    "SBP" => "concepts/sbp.md"
                   ]
              "Building PDESolver" => Any[
                    "build.md",
                    "Build Options" => "deps_readme.md"
                   ]
              "DOC To be broken up or organized" => Any[
                    "Code Interfaces" => "interfaces.md"
                    "Code Parallelization" => "parallel.md"
                   ]
              "Examples" => Any[
                    "Isentropic Vortex" => "examples/isentropic.md"
                    "Unsteady Vortex" => "examples/unsteady.md"
                   ]
              "Frontend" => Any[
                    "Introduction" => "pdesolver.md"
                    "PDESolver User Interface" => "pdesolver_user.md"
                    "PDESolver Physics Interface" => "pdesolver_physics.md"
                    "PDESolver Structure" => "pdesolver_structure.md"
                   ]
              "Invocation" => Any[
                    "Calling PDESolver" => "invocation/calling.md"
                    "Interactive Session (experimental)" => "invocation/interactive.md"
                   ]
              "Solver" => Any[
                  "solver/Readme.md"
                  "solver/misc.md"
                  "solver/SolverCommon.md"
                  "Advection" => Any[
                      "Introduction" => "$advec/advection.md"
                      "Datatypes" => "$advec/types.md"
                      "Volume Integrals" => "$advec/volume.md"
                      "Face Integrals" => "$advec/flux.md"
                      "Boundary Integrals" => "$advec/bc.md"
                      "Initial Condition" => "$advec/ic.md"
                      "Source Term" => "$advec/source.md"
                      "Common Functions" => "$advec/common.md"
                      "Adjoint" => "$advec/adjoint.md"
                      "Boundary Functional" => "$advec/boundary_functional.md"
                     ]
                  "Euler" => Any[
                      "Introduction" => "$euler/euler.md"
                      "Datatypes" => "$euler/types.md"
                      "Volume Integrals" => "$euler/volume.md"
                      "Volume Integrals Jacobian" => "$euler/volume_diff.md"
                      "Face Integrals" => "$euler/flux.md"
                      "Face Integrals Jacobian" => "$euler/flux_diff.md"
                      "Face Element Integrals" => "$euler/faceElementIntegrals.md"
                      "Boundary Integrals" => "$euler/bc.md"
                      "Boundary Integrals Jacobian" => "$euler/bc_diff.md"
                      "Initial Conditions" => "$euler/ic.md"
                      "Source Term" => "$euler/source.md"
                      "Common Functions" => "$euler/common.md"
                      "Conversion" => "$euler/conversion.md"
                      "Numerical Flux Functions" => "$euler/flux_functions.md"
                      "Numerical Flux Functions Jacobian" => "$euler/flux_functions_diff.md"
                      "Stabilization" => "$euler/stabilization.md"
                      "Adjoint" => "$euler/adjoint.md"
                      "Boundary Functional" => "$euler/boundary_functional.md"
                      "Misc" => "$euler/misc.md"
                      "Homotopy" => "$euler/homotopy.md"
                      "Homotopy Jacobian" => "$euler/homotopy_diff.md"
                      "Eigensystem" => "$euler/eigensystem.md"
                      "Startup" => "$euler/startup.md"
                     ]
                  "Simple ODE" => Any[
                      "Main" => "$ode/simpleODE.md"
                     ]

                 ]  # end Solver
              "Input" => Any[
                  "Introduction" => "input/input.md"
                  "Important Keys" => "input/keys.md"
                 ]
              "LinearSolvers" => Any[
                "Introduction" => "linearsolvers/linearsolvers.md"
                "Preconditioners" => "linearsolvers/pc.md"
                "Linear Operators" => "linearsolvers/lo.md"
                "Linear Solvers" => "linearsolvers/ls.md"
                ]  # end LinearSolvers
              "NonlinearSolvers" => Any[
                  "Introduction" => "NonlinearSolvers/nonlinearsolvers.md"
                  "Steady" => "NonlinearSolvers/steady.md"
                  "Unsteady" => Any[
                      "Intro" => "NonlinearSolvers/unsteady/intro.md"
                      "Runge-Kutta" => "NonlinearSolvers/unsteady/rk4.md"
                      "LSERK" => "NonlinearSolvers/unsteady/lserk.md"
                      "Crank-Nicolson" => "NonlinearSolvers/unsteady/cn.md"
                      "Crank-Nicolson: Unsteady Adjoint" => "NonlinearSolvers/unsteady/cn_uadj.md"
                    ]
                  "Newton's Method" => "NonlinearSolvers/newton.md"
                  "Jacobian Calculation" => "NonlinearSolvers/jacobian.md"
                  "Jacobian Freezing" => "NonlinearSolvers/jac_recalc.md"
                  "Residual Evalution" => "NonlinearSolvers/residual_evaluation.md"
                  "Matrix Interface" => "NonlinearSolvers/matrix.md"
                  "Newton Inner" => "NonlinearSolvers/newton_inner.md"
                 ]
              "Utils" => Any[
                  "Main" => "Utils/Utils.md"
                  "Parallel Constructs" => "Utils/parallel.md"
                  "Projections" => "Utils/projections.md"
                  "Logging" => "Utils/logging.md"
                  "Input/Output" => "Utils/io.md"
                  "Checkpointing" => "Utils/checkpoint.md"
                  "Misccellaneous" => "Utils/misc.md"
                 ]
                "Testing" => Any[
                    "Introduction" => "test/Testing.md"
                    "Local Testing" => "test/Readme.md"
                    "CI Testing" => "test/Travis.md"
                    "Test API" => "test/TestSystem.md"
                  ]
             ] # end Home
)  # end mkdocs
