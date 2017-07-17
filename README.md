# PDESolver
PDESolver is a multi-physics solver primarily focused on Computational Fluid
Dynamics.  It has been designed from ground up to support optimization,
robust simulation, and parallel scalability.

## Installation
The obtain this package, do `Pkg.clone(url)`, then `Pkg.build("PDESolver")`.  This will install all dependences, including those not listed in Metadata.

The dependencies not listed in Metadata are:

ODLCommonTools

SummationByParts

PumiInterface

Petsc

PumiInterface and Petsc are the only packages with non-trivial installation requirements, although for most cases the defaults will work just fine.

See those packages for details.


## Running the Code
Before running the code, you must `source` the shell script
`PumiInterface/src/use_julialib.sh`.  This enables Julia to find and use Pumi.

The code takes an input file that defines all options for the solver, including
the which physics to solve, which mesh to use, initial conditions,
boundary conditions, and discretization.  The file `src/input/input_vals.txt`
describes the input file format and all valid options.

### Simple Mode
Once an input file is prepared, the solver can be invoked with

```julia
mpirun -np x julia /path/to/PDESolver/src/startup.jl "input_file_name"
```

This will run the solver in parallel with `x` processes, solving the physics
specified in the input file.

### Advanced Mode
The code also provides a scripting interface.  This enables users to supply
initial conditions and boundary conditions, as well as do advanced
post-processing if desired.

The template script is:

```julia
using PDESolver  # load the interface to the code
using AdvectionEquationMod  # load the desired physics module

# register ICs and BCs here

input_fname = ARGS[1]  # use the first command line argument as the input file
                       # name

# solve the equation defined in the input file
mesh, sbp, eqn, opts = run_solver(input_fname)

# do post-processing using mesh, sbp, eqn, opts here

```

`mesh, sbp, eqn, opts` are the `AbstractMesh` object, the SummationByParts operator, the `AbstractSolutionData`, and the options dictionary, respectively.

The `AbstractMesh` object contains all the data about the mesh, including
coordinates and mapping jacobian information.  The `sbp` object is used for
numerically approximating derivatives. The `AbstractSolutionData` object contains
all the data about the solution and the residual of the equation.  The input
file gets parsed into the options dictionary, which is used to by the rest of
the code.

The PDESolver module exports several functions for use by the scripting
interface.  The first two `registerIC` and `registerBC` allow users to supply
initial conditions and boundary conditions.  See the documentation of the
functions for details.  The functions `printICNames` and `printBCNames` print
the names of all currently registers initial conditions and boundary conditions.

For example a user can, either in a Julia script or interactive (REPL) session:
```julia
using PDESolver
using AdvectionEquationMod

printBCNames(AdvectionEquationMod)
```

And this will print the names of all the boundary conditions currently known
to the Advection Equation module.

## Developer Guide

The abstraction hierarchy of the code is described in two documents.
`doc/interfaces.md` describes the required implementation and rationale for
the `AbstractMesh` and `AbstractSolutionData` types, as well as how different
parts of the code interact (ie. how the NonlinearSolvers interact with the
physics modules etc.).  `src/Readme.md` describes the user-facing interface
for the physics modules as well as what functions each physics module must
implement in order to be usable.

The code is organized as follows:

The code for evaluating the residual for an equation is in the src/solver directory.  
There are subdirectories for each physics.  Currently src/solver/euler is the most well developed physics.

The NonlinarSolvers module (`src/NonlinearSolvers`) implements both
time stepping methods for unsteady equations and non-linear root finding methods for steady problems.

The Utils module (`src/Utils`) contains auxiliary functions used by all
 physics modules.

The Input module (`src/input`) parses input files, provides default values,
and checks for unrecognized keys.  `src/input/input_vals.txt` lists all
possible user-supplied keys. `src/input/input_vals_internal.txt` list keys
used internally by the solver.  Users should never specify these keys in an
input file.

The `src/mesh_files` contains meshes used by the tests.  Please use the
smallest meshes possible for tests.

The directory `src/simple_mesh` contains the files to create simple structured meshes.

The file `src/estimateMem.jl` reads the file counts.txt, which is written by the mesh initialization, and estimates the steady state memory usage of the solver.
While this is not an upper bounds, it has proven to be quite accurate.  It does not include memory used by any factorizations the linear solver might do.


## Visualization
Paraview is used for visualization.  High order elements (p>2) are automatically
interpolated onto either a linear or quadratic mesh for visualization.
Empirically, this has been shown to produce better images than
sub-triangulating each element and doing an exact projection onto the finer
mesh.

# Version History
v0.1: the old master version (DG supported, I think)
v0.2: the work branch after the test system/frontend rewrite
v0.3: curvilinear entropy stable works (requires PumiInterface v0.3, SBP tag ticon_broken)
v0.4: curvilinear works (both entropy stable and Roe scheme)

[![Build Status](https://travis-ci.org/OptimalDesignLab/PDESolver.jl.svg)](https://travis-ci.org/OptimalDesignLab/PDESolver.jl)
[![Documentation Status](https://readthedocs.org/projects/pdesolverjl/badge/?version=latest)](https://readthedocs.org/projects/pdesolverjl/?badge=latest)
[![codecov](https://codecov.io/gh/OptimalDesignLab/PDESolver.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/OptimalDesignLab/PDESolver.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](http://www.optimaldesignlab.com/PDESolver.jl)

[//]: # (for some reason the webite is hosted at optimaldesignlab.com instead of github)

