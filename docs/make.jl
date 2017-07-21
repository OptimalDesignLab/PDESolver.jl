using Documenter, PDESolver

# some abbreviations
advec = "solver/advection"
euler = "solver/euler"
ode = "solver/simpleODE"

makedocs(
  format = :html,
  sitename = "PDESolver.jl",
  pages = Any["PDESolver Introduction" => "index.md"
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
                    "PDESolver PhysicsInterface" => "pdesolver_physics.md"
                   ]
              "Invocation" => Any[ 
                    "Calling PDESolver" => "invocation/calling.md"
                    "Interactive Session (experimental)" => "invocation/interactive.md"
                   ]
              "Solver" => Any[
                  "solver/Readme.md"
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
                      "Face Integrals" => "$euler/flux.md"
                      "Boundary Integrals" => "$euler/bc.md"
                      "Initial Conditions" => "$euler/ic.md"
                      "Source Term" => "$euler/source.md"
                      "Common Functions" => "$euler/common.md"
                      "Conversion" => "$euler/conversion.md"
                      "Numerical Flux Functions" => "$euler/flux_functions.md"
                      "Stabilization" => "$euler/stabilization.md"
                      "Adjoint" => "$euler/adjoint.md"
                      "Boundary Functions" => "$euler/boundary_functional.md"
                      "Misc" => "$euler/misc.md"
                     ]
                  "Simple ODE" => Any[ 
                      "Main" => "$ode/simpleODE.md"
                     ]

                 ]  # end Solver
              "Input" => Any[
                  "Introduction" => "input/input.md"
                 ]
              "NonlinearSolvers" => Any[
                  "Introduction" => "NonlinearSolvers/nonlinearsolvers.md"
                  "Steady" => "NonlinearSolvers/steady.md"
                  "Unsteady" => "NonlinearSolvers/unsteady.md"
                 ]
              "Utils" => Any[ 
                  "Main" => "Utils/Utils.md"
                  "Parallel Constructs" => "Utils/parallel.md"
                  "Projections" => "Utils/projections.md"
                  "Logging" => "Utils/logging.md"
                  "Input/Output" => "Utils/io.md"
                 ]
             ] # end Home
)  # end mkdocs

