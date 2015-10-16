# PDESolver
 
## Installation
The obtain this package, do `Pkg.clone(url)`, then `Pkg.build("PDESolver")`.  This will install all dependences, including those not listed in Metadata.

The dependencies not listed in Metadata are:

ODLCommonTools

SummationByParts

PumiInterface

Petsc

PumiInterface and Petsc are the only packages with non-trivial installation requirements, although for most cases the defaults will work just fine.

See those packages for details.

## Known issues:
Need to build PumiInterface every time:

cd ~/.julia/v0.4/PumiInterface/src/
source use_julialib.sh
./build_shared.scorec.sh6
cd PDESolver/tests
julia runtests.jl

## Running the Code
To run the code, execute `julia /path/to/startup.jl "input_file_name"`, where startup.jl is located in the subdirctory for the physics you are solving.  Each physics has its own startup.jl


## Developer Guide
The code for evaluating the residual for an equation is in the src/solver directory.  
There are subdirectories for each physics.  Currently src/solver/euler is the most well developed physics.

The nl_solvers directory contains the files to solve non-linear equations.
It defines a module, nl_solvers, that export all the currently available methods, including Newtons method and RK4

The input directory contains the files to read and process the input arguments dictionary.  
The input_vals.txt file details the recognized arguments.

The directory simple_mesh contains the files to create simple structured meshes.

The file src/estimateMem.jl reads the file counts.txt, which is written by the mesh initialization, and estimates the steady state memory usage of the solver.
While this is not an upper bounds, it has proven to be quite accurate.  It does not include memory usage by the linear solver.

### Mesh
Any implementation of the abstract type AbstractMesh must provide the data required by the solver about the mesh.
For performance reasons, all data used by the solver should be stored in the mesh object, typically in arrays.  During a resdual evaluation, the underlying mesh software should not be accessed.
More details about this are forthcoming.


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


#### Data Structures
EulerData.q stores the conservative variables in a 3 dimensional array, numDofPerNode x numNodesPerElement x numEl.
Although this contains some duplicate information compared to a column vector, it is much more convienient for doing computations.
EulerData.res stores the non-linear residual, and is the same shape ad EulerData.q.
The workflow of the solver is to take EulerData.q, do some operation on it, and store the residual in EulerData.res.
The evalEuler function is the driver function for evaluating the residual for a given EulerData object and mesh.
The input to the function is the q array, and the output is the res array.
Functions are provided for to take a column vector, usually EulerData.SL0, and disassmble it into EulerData.q, and conversely, to take EulerData.res and assemble it into a column vector, usually eqn.SL.

#### Precalculation of Quantities
Some quantities are preallocated and stored in the EulerData object.
The equation flux is stored for all equations.  
The 'aux_vars' field is used to store any additional variables.
The macros in euler_macros.jl provide the means to access the variables stored there, so the user does not have to do the bookkeeping.
If using Newtons method to solve the non-linear equation, the cost of storing the Jacobian dominates the total memory cost, to it is usually worth it to precalculate and store variables for speed.


#### Volume Integrals
The volume integrals are computed using the stored flux using the function evalVolumeIntegrals.


#### Boundary Conditions
Boundary integrals are computed using evalBoundaryIntegrals.  
The boundary flux must be calculated first using the function getBCFluxes.

Boundary conditions are implimented with functors.  
For each type of boundary condition, a dummy type is defined and the call function for it.
The dummy type should be suffixed with 'BC', and be a subtype of 'BCTypes'.  
All call functions should take the same arguments, even if they are not needed.
The function should calculate the boundary flux at a node, where flux going into the element is positive.
An instance of the dummy type must be added  to the dictionary BCDict.  
The function getBCFunctors is called during initialization to get the dummy types from the dictionary and store
them in the mesh object for use by calcBoundaryFlux, which calculates and stores the boundary flux.

#### Stabilization 
The addStabilization function add the stabilization term to the residual.  It calls functions in the stabilization.jl file.




## Visualization
Paraview is used for visualization.  High order elements (p>2) are automatically
subtriangulated into a set of linear elements and nodal solutions are copied 
from the high order mesh to the linear mesh.  The field faceNums_old_1 on the 
linear mesh shows the the element number of the high order element the linear 
element was created from.

[![Build Status](https://travis-ci.org/OptimalDesignLab/PDESolver.jl.svg)](https://travis-ci.org/OptimalDesignLab/PDESolver.jl)
[![Coverage Status](https://coveralls.io/repos/OptimalDesignLab/PDESolver.jl/badge.png)](https://coveralls.io/r/OptimalDesignLab/PDESolver.jl)
[![Documentation Status](https://readthedocs.org/projects/pdesolverjl/badge/?version=latest)](https://readthedocs.org/projects/pdesolverjl/?badge=latest)
