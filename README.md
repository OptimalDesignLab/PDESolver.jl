# PDESolver

## Installation
The obtain this package, do `Pkg.clone(url)`, then `Pkg.build("PDESolver")`.  This will install all dependences, including those not listed in Metadata.

The dependencies not listed in Metadata are:

PDESolverCommon

SummationByParts

PumiInterface

PumiInterface is the only package with non-trivial installation requirements.  
See that package for details

## Developer Guide
The code for evaluating the residual for an equation is in the src/solver directory.  
There are subdirectories for each physics.  Currently src/solver/euler is the most well developed physics.

### Structure of each solver
The following describes the structure of the Euler solver.  Other solvers use a similar structure.

EulerEquationMod.jl contains the defintion for the EulerData object, which stores all data needed to solve the equation related to the physics.
euler.jl contains most of the functions needed to evaluate the equation.
bc.jl contains the functions that calculate the boundary flux for the supported boundary conditions.
bc_solvers.jl contains the boundary flux solvers such as the Roe solver.
ic.jl contains the functions to populate eqn.SL0 with the initial condition
common_funcs.jl contains functions that evaluate different conditions used by initial and boundary conditions.
stabilization.jl contains the functions to add stabilization.
euler_macros.jl contains macros to access the values stored in eqn.aux_vars

The nl_solvers directory contains the files to solve the non-linear equations.
It defines a module, nl_solvers, that export all the currently available methods, including Newtons method and RK4

The input directory contains the files to read and process the input arguments dictionary.

The directory simple_mesh contains the files to create simple structured meshes.


[![Build Status](https://travis-ci.org/OptimalDesignLab/PDESolver.jl.svg)](https://travis-ci.org/OptimalDesignLab/PDESolver.jl)
[![Coverage Status](https://coveralls.io/repos/OptimalDesignLab/PDESolver.jl/badge.png)](https://coveralls.io/r/OptimalDesignLab/PDESolver.jl)
[![Documentation Status](https://readthedocs.org/projects/pdesolverjl/badge/?version=latest)](https://readthedocs.org/projects/pdesolverjl/?badge=latest)
