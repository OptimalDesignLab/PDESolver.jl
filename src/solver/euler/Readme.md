# EulerEquationMod

This module implements the Euler physics.  The internal structure of the
code is described below.

## Mesh

Any implementation of the abstract type AbstractMesh must provide the data required by the solver about the mesh.
For performance reasons, all data used by the solver should be stored in the mesh object, typically in arrays.  During a resdual evaluation, the underlying mesh software should not be accessed.
More details about this are forthcoming.


## Structure of each solver

The following describes the structure of the Euler solver.  Other solvers use a similar structure.

EulerEquationMod.jl contains the defintion for the EulerData object, which stores all data needed to solve the equation related to the physics.
euler.jl contains most of the functions needed to evaluate the equation.
bc.jl contains the functions that calculate the boundary flux for the supported boundary conditions.
bc_solvers.jl contains the boundary flux solvers such as the Roe solver.
ic.jl contains the functions to populate eqn.q_vec with the initial condition
common_funcs.jl contains functions that evaluate different conditions used by initial and boundary conditions.
stabilization.jl contains the functions to add stabilization.
euler_macros.jl contains macros to access the values stored in eqn.aux_vars


### Data Structures

EulerData.q stores the conservative variables in a 3 dimensional array, numDofPerNode x numNodesPerElement x numEl.
Although this contains some duplicate information compared to a column vector, it is much more convienient for doing computations.
EulerData.res stores the non-linear residual, and is the same shape ad EulerData.q.
The workflow of the solver is to take EulerData.q, do some operation on it, and store the residual in EulerData.res.
The evalResidual function is the driver function for evaluating the residual for a given EulerData object and mesh.
The input to the function is the q array, and the output is the res array.


### Precalculation of Quantities

Some quantities are preallocated and stored in the EulerData object.
The equation flux is stored for all equations.  
The 'aux_vars' field is used to store any additional variables.
The macros in euler_macros.jl provide the means to access the variables stored there, so the user does not have to do the bookkeeping.
If using Newtons method to solve the non-linear equation, the cost of storing the Jacobian dominates the total memory cost so the size of the data arrays
is not a limiting factor, however it may be faster to compute quantities on
the fly (and thus reduce the size of the working data set) than precalculate
and access the data again later.


### Volume Integrals

The volume integrals are computed using the stored flux using the function evalVolumeIntegrals.


### Boundary Conditions

Boundary integrals are computed using evalBoundaryIntegrals.  
The boundary flux must be calculated first using the function getBCFluxes.

Boundary conditions are implemented with functors.  
For each type of boundary condition, a dummy type is defined and the call function for it.
The dummy type should be suffixed with 'BC', and be a subtype of 'BCTypes'.  
All call functions should take the same arguments, even if they are not needed.
The function should calculate the boundary flux at a node, where flux going into the element is positive.
An instance of the dummy type must be added  to the dictionary BCDict.  
The function getBCFunctors is called during initialization to get the dummy types from the dictionary and store
them in the mesh object for use by calcBoundaryFlux, which calculates and stores the boundary flux.

### Stabilization

The addStabilization function add the stabilization term to the residual.  It calls functions in the stabilization.jl file.
