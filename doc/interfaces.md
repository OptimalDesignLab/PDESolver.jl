# Interfaces in PDESolver
PDESolver depends on the the three main objects, the AbstractSolutionData object,  AbstractMesh object, and the SBP object implementing certain interfaces.
This document describes what the interfaces are, and gives some hints for how to implement them.

Before doing so, a short description of what general Julia interfaces look like is in order.  
The paradigm of Julia code is that of "objects with associated functions", where a new Type is defined, and then functions that take the Type as an argument are defined.
The functions define the interface to the Type.
The Type holds data (ie. state), and the functions perform operations on that state (ie. behavior).
Perhaps counter-intuitively, it is generally not recommended for users of a type to  access the fields directly.
Instead, any needed operations on the data that the Type holds should be provided through functions.
The benefit of this convention is that it imposes no requirements on how the Type stores its data or implements its behavior.
This is important because a user of the Type should not be concerned with these things.
The user needs to know what behavior the Type has, but not how it is implemented.
This distinction becomes even more important when there are multiple implementations certain functionality.
The user should be able to seamlessly transition between different implementations.
This requires all implementations have the same interface.

The question of how to enforce interfaces, and how strongly to do so, is still an open question in Julia.
Some relevant Github issues:
* [5](https://github.com/JuliaLang/julia/issues/5)
* [4935](https://github.com/JuliaLang/julia/issues/4935)
* [6975](https://github.com/JuliaLang/julia/issues/6975)


One of the strongest arguments against the "functions as interfaces" idea is that for many applications, just storing data in an array is best.
Creating interface functions for the Type to implement the array interface would be a lot of extra code with no benefit.
For this reason, it makes sense to directly access the fields of some Types, to avoid trivial get/set methods.
We do this extensively in PDESolver, because arrays are the natural choice for storing the kind of data used in PDESolver.

## AbstractSolutionData
ODLCommonTools defines:

`abstract AbstractSolutionData{Tsol}`.

The purpose of an `AbstractSolutionData` is to hold all the data related to the solution of an equation.
This includes the solution at every node and any auxiliary quantities.
The storage for any quantity that is calculated over the entire mesh should be allocated as part of this object, in order to avoid repeatedly reallocated the array for every residual evaluation.
In general, there should never be a need to allocate a vector longer than the number of degrees of freedom at a node (or a matrix similarly sized matrix) during a residual evaluation.
Structuring code such that it conforms with this requirement has significant performance benefits because it reduces memory allocation/deallocation.

The static parameter `Tsol` is the datatype of the solution variables.

### Required Fields
The required fields of an `AbstractSolutionData are`:
```
  q::AbstractArray{Tsol, 3}
  q_vec::AbstractArray{Tsol, 1}
  res::AbstractArray{Tres, 3}
  res_vec::AbstractArray{Tres, 1}
  M::AbstractArray{Float64, 1}
  Minv::AbstractArray{Float64, 1}
  disassembleSolution::Function
  assembleSolution::Function
  multiplyA0inv::Function
  majorIterationCallback::Function
```

The purpose of these fields are:

`q`: to hold the solution variables in an element-based array. 
     This array should be `numDofPerNode` x `numNodesPerElement` x `numEl`.
     The residual evaluation *only* uses `q`, never `q_vec`

`q_vec`: to hold the solution variables as a vector, used for any linear algebra operations and time stepping.
This array should have a length equal to the total number of degrees of freedom in the mesh.
Even though this vector is not used by the residual evaluation, it is needed for many other operations, so it is allocated here so the memory can be reused.
There are functions to facilitate the scattering of values from `q_vec` to `q`.
Note that for Continuous Galerkin type discretization (as opposed to Discontinuous Galerkin discretizations), there is not a corresponding "gather" operation (ie. `q` -> `q_vec`).

`res`: similar to `q`, except that the residual evaluation function populates it with the residual values.  
       As with `q`, the residual evaluation function only interacts with this array, never with `res_vec`.

`res_vec`: similar to `q_vec`.  Unlike `q_vec` there are functions to perform an additive reduction (basically a "gather") of `res` to `res_vec`.  For continuous Galerkin discretizations, the corresponding "scatter" (ie. `res_vec` -> res`) may not exist.

`M`:  The mass matrix of the entire mesh.  Because SBP operators have diagonal mass matrices, this is a vector.  Length numDofPerNode x numNodes (where numNodes is the number of nodes in the entire mesh).

`Minv`:  The inverse of the mass matrix.

`disassembleSolution`:  Function that takes the a vector such as `q_vec` and scatters it to an array such as `q`.
                        This function must have the signature:
                        `disassembleSolution(mesh::AbstractMesh, sbp, eqn::AbstractSolutionData, opts, q_arr:AbstractArray{T, 3}, q_vec::AbstractArray{T, 1}`
                        Because this variable is a field of a type, it will be dynamically dispatched.
                        Although this is slower than compile-time dispatch, the cost is insignificant compared to the cost of evaluating the residual, so the added flexibility of having this function as a field is worth the cost.

`assembleSolution`:  Function that takes an array such as `res` and performs an additive reduction to a vector such as `res_vec`.
                     This function must have the signature:

                     `assembleSolution(mesh::AbstractMesh, sbp, eqn::AbstractSolutionData, opts, res_arr::AbstractArray{T, 3}, res_vec::AbstractArray{T, 1}, zero_resvec=true)`
                     The argument `zero_resvec` determines whether `res_vec` is zeroed before the reduction is performed.
                     Because it is an additive reduction, elements of the vector are only added to, never overwritten, so forgetting to zero out the vector could cause strange results.
                     Thus the default is true.

`multiplyA0inv`:  Multiplies the solution values at each node in an array such as `res` by the inverse of the coefficient matrix of the time term of the equation.
                  This function is used by time marching methods.
                  For some equations, this matrix is the identity matrix, so it can be a no-op, while for others might not be.
                  The function must have the signature:

                  `multiplyA0inv(mesh::AbstractMesh, sbp, eqn::AbstractSolutionData, opts, res_arr::AbstractArray{Tsol, 3})`


`majorIterationCallback`:  function called before every step of Newton's method or stage of an explicit time marching scheme.
This function is used to do output and logging.
The function must have the signature:

`function majorIterationCallback(itr, mesh::AbstractMesh, sbp::SBPOperator, eqn::AbstractEulerData, opts)`

##AbstractMesh
ODLCommonTools defines:

`abstract AbstractMesh{Tmsh}`.

The purpose of an `AbstractMesh` is to hold all the mesh related data that the solver will need.
It also serves to establish an interface between the solver and whatever mesh software is used.
By storing all data in the fields of the `AbstractMesh` object, the details of how the mesh software stores and allows retrieval of data are not needed by the solver.
This should make it easy to accommodate different mesh software without making any changes to the solver.

The static parameter `Tmsh` is used to enable differentiation with respect to the mesh variable in the future.

###Required Fields
```
  # counts
  numVert::Integer
  numEl::Integer
  numNodes::Integer
  numDof::Integer
  numDofPerNode::Integer
  numNodesPerElement::Integer
  order::Integer

  # mesh data
  coords::AbstractArray{Tmsh, 3}
  dxidx::AbstractArray{Tmsh, 4}
  jac::AbstractArray{Tmsh, 2}

  # boundary condition data
  numBC::Integer
  numBoundaryEdges::Integer
  bndryfaces::AbstractArray{Boundary, 1}
  bndry_offsets::AbstractArray{Integer, 1}
  bndry_funcs::AbstractArray{BCType, 1}
  
  # interior edge data
  numInterfaces::Integer
  interfaces::AbstractArray{Interface, 1}

  # degree of freedom number data
  dofs::AbstractArray{Integer, 2}
  sparsity_bnds::AbstractArray{Integer, 2}
  sparsity_nodebnds::AbstractArray{Integer, 2}

  # mesh coloring data
  numColors::Integer
  color_masks::AbstractArray{ AbstractArray{Number, 1}, 1}, 
  neighbor_nums::AbstractArray{Integer, 2}
  pertNeighborEls::AbstractArray{Integer, 2}
```
####Counts

`numVert`:  number of vertices in the mesh

`numEl`:  number of elements in the mesh

`numNodes`: number of nodes in the mesh

`numDof`:  number of degrees of freedom in the mesh (= `numNodes` * `numDofPerNode`)

`numDofPerNode`:  number of degrees of freedom on each node.

`numNodesPerElement`:  number of nodes on each element.

`order`:  order of the discretization (ie. first order, second order...), where an order `p` discretization should have a convergence rate of `p+1`.


####Mesh Data
`coords`: `n` x `numNodesPerElement` x `numEl` array, where `n` is the dimensionality of   the equation being solved (2D or 3D typically).  `coords[:, nodenum, elnum] = [x, y, z]` coordinates of node `nodenum` of element `elnum`.

`dxidx`:  `n` x `n` x `numNodesPerElement` x `numEl`, where `n` is defined above.
It stores the mapping jacobian scaled by `( 1/det(jac) dxi/dx )` where `xi` are the parametric coordinates, `x` are the physical (x,y,z) coordinates, and `jac` is the determinant of the mapping jacobian `dxi/ dx`.

`jac`  : `numNodesPerElement` x `numEl` array, holding the determinant of the mapping jacobian `dxi/dx` at each node of each element.


####Boundary Condition Data
The mesh object stores data related to applying boundary conditions.
Boundary conditions are imposed weakly, so there is no need to remove degrees of freedom from the mesh when Dirichlet boundary conditions are applied.
In order to accommodate any combination of boundary conditions, an array of functors are stored as part of the mesh object, along with lists of which mesh edges (or faces in 3D) should have which boundary condition applied to them


`numBC`: number of different types of boundary conditions used.

`numBoundaryEdges`: number of mesh edges that have boundary conditions applied to them.

`bndryfaces`:  array of Boundary objects (which contain the element number and the local index of the edge), of length `numBoundaryEdges`.

`bndry_offsets`:  array of length numBC+1, where `bndry_offsets[i]` is the index  in `bndryfaces` where the edges that have boundary condition `i` applied to them start.
The final entry in `bndry_offsets` should be `numBoundaryEdges + 1`.
Thus `bndryfaces[ bndry_offsets[i]:(bndry_offsets[i+1] - 1) ]` contains all the boundary edges that have boundary condition `i` applied to them.

`bndry_funcs`:  array of boundary functors, length `numBC`.  All boundary functors are subtypes of `BCType`.  Because `BCType` is an abstract type, the elements of this array should not be used directly, but passed as an argument to another function, to avoid type instability.

####Interior Edge Data
Data about interior mesh edges (or faces in 3D) is stored to enable use of edge stabilization or Discontinuous Galerkin type discretizations.
Only data for edges (faces) that are shared by two elements are stored (ie. boundary edges are not considered).

`numInterfaces`:  number of interior edges

`interfaces`:  array of Interface types (which contain the element numbers for the two elements sharing the edge, and the local index of the edge from the perspective of the two elements, and an indication of the relative edge orientation).
The two element are referred to as `elementL` and `elementR`, but the choice of which element is `elementL` and which is `elementR` is arbitrary.
The length of the array is numInterfaces.
Unlike `bndryfaces`, the entries in the array do not have to be in any particular order.

####Degree of Freedom Numbering Data
`dofs`:  `numDofPerNode` x `numNodesPerElement` x `numEl` array.
Holds the degree of freedom number of each degree of freedom.

`sparsity_bnds`:  2 x `numDof` array.
`sparsity_bnds[:, i]` holds the maximum, minimum degree of freedom numbers associated with degree of freedom `i`.
In this context, degrees of freedom `i` and `j` are associated if entry `(i,j)` of the jacobian is non-zero.
In actuality, `sparsity_bnds` need only define upper and lower bounds for degree of freedom associations (ie. they need not be tight bounds).
This array is used to to define the sparsity pattern of the jacobian matrix.

`sparsity_nodebnds`:  2 x numNodes array.
`sparsity_bnds[:, i]` holds the maximum, minimum node associated with node `i`, similar the information stored in `sparsity_bnds` for degrees of freedom.


####Mesh Coloring Data
The NonlinearSolvers module uses algorithmic differentiation to compute the Jacobian.
Doing so efficiently requires perturbing multiple degrees of freedom simultaneously, but perturbing associated degrees of freedom at the same time leads to incorrect results.
Mesh coloring assigns each element of the mesh to a group (color) such that every degree of freedom on each element is not associated with any other degree of freedom on any other element of the same color.
An important aspect of satisfying this condition is the use of the element-based arrays (all arrays that store data for a quantity over the entire mesh are `ncomp` x `numNodesPerElement` x `numEl`).
In such an array, any node that is part of 2 or more elements has one entry for each element.
When performing algorithmic differentiation, this enables perturbing a degree of freedom on one element without perturbing it on the other elements that share the degree of freedom.

For example, consider a node that is shared by two elements.
Let us say it is node 2 of element 1 and node 3 of element 2.
This means `AbstractSolutionData.q[:, 2, 1]` stores the solution variables for this node on the first element, and `AbstractSolutionData.q[:, 3, 2]` stores the solution variables for the second element.
Because these are different entries in the array `AbstractSolutionData.q`, they can be perturbed independently.
Because `AbstractSolutionData.res` has the same format, the perturbations to `AbstractSolutionData.q[:, 2, 1] are mapped to `AbstractSolutionData.res[:, 2, 1]` for a typical continuous Galerkin type discretization.
This is a direct result of having an element-based discretization.

There are some discretizations, however, that are not strictly element-based.
Edge stabilization, for example, causes all the degrees of freedom of one element to be associated with any elements it shares an edge with.
To deal with this, we use the idea of a distance-n coloring.
A distance-n coloring is a coloring where there are n elements in between two element of the same color. 
For element-based discretizations with element-based arrays, every element in the mesh can be the same color.
This is a distance-0 coloring.
For an edge stabilization discretization, a distance-1 coloring is required, where every element is a different color than any neighbors it shares and edge with.
(As a side node, the algorithms that perform a distance-1 coloring are rather complicated, so in practice we use a distance-2 coloring instead).

In order to do algorithmic differentiation, the `AbstractMesh` object must store the information that determines which elements are perturbed for which colors, and, for the edge stabilization case, how to relate a perturbation in the output `AbstractSolutionData.res` to the degree of freedom in `AbstractSolutionData.q` in O(1) time.
Each degree of freedom on an element is perturbed independently of the other degrees of freedom on the element, so the total number of residual evaluations is the number of colors times the number of degrees of freedom on an element.

The fields required are:

`numColors`:  The number of colors in the mesh.

`color_masks`:  array of length `numColors`.  Each entry in the array is itself an array of length `numEl`.  Each entry of the inner array is either a 1 or a 0, indicating if the current element is perturbed or not for the current color.
For example, in `color_mask_i = color_masks[i]; mask_elj = color_mask_i[j]`, the variable `mask_elj` is either a 1 or a zero, determining whether or not element `j` is perturbed as part of color `i`.


`neighbor_nums`:  `numEl` x `numColors` array.  `neighbor_nums[i,j]` is the element number of of the element whose perturbation is affected element `i` when color `j` is being perturbed, or zero if element `i` is not affected by any perturbation.  

###Other functions
The mesh module must also define and export the functions


```
saveSolutionToMesh(mesh::MeshImplementationType, vec::AbstractVector)
writeVisFiles(mesh::MeshImplementationType, fname::ASCIIString)`

where the first function takes a vector of length `numDof` and saves it to the mesh, and the second writes Paraview files for the mesh, including the solution field.


##Physics Module
For every new physics to be solved, a new module should be created.
The purpose of this module is to evaluate the equation:

`dq/dt = f(q)`

For steady problems, `dq/dt = 0` and the module evaluates the residual.
For unsteady problems, the form `dq/dt = f(q)` is suitable for explicit time marching.

###AbstractSolutionData Implementation

Each physics module should define and export a subtype of `AbstractSolutionData{Tsol}`.
The implementation of `AbstractSolutionData{Tsol}` must inherit the `Tsol` static parameter, and may have additional static parameters as well.
It may also be helpful to define additional abstract types within the physics module to provide different levels of abstractions.
For example, the Euler physics module defines:

```
abstract AbstractEulerData{Tsol} <: AbstractSolutionData{Tsol}
abstract EulerData {Tsol, Tdim, Tres, var_type} <: AbstractEulerData{Tsol}
type EulerData_{Tsol, Tres, Tdim, Tmsh, var_type} <: EulerData{Tsol, Tdim, Tres, var_type}
```

The first line is effectively just a name change and may not be necessary.
The second line adds the static parameters `Tdim`, `Tres`, and `var_type` while inheriting the `Tsol` type from `AbstractEulerData`.
`Tdim` is the dimensionality of the equation, `Tres` is the datatype of the residual variables, and `var_type` is a symbol indicating whether the equation is being solved with conservative or entropy variables.
The third line defines a concrete type that implements all the features required of an `AbstractSolutionData`, and adds a static parameter `Tmsh`, the datatype of the mesh variables.  
The additional static parameter is necessary because one field of `EulerData_` has type `Tmsh`.
Note that there could be multiple implementations of `AbstractSolutionData` for the Euler equations, perhaps with different fields to store certain data or not.
All these implementations will need to have the static parameters `Tsol`, `Tdim`, `Tres`, and `var_type`, so `EulerData` is defined as an abstract type, so all implementations can inherit from it.
All high level function involved in evaluating the residual will take in an argument of type `EulerData`.
Only when low level functions need to dispatch based on which implementation is used would it take in an `EulerData_` or another implementation.

###Residual Evaluation
In addition to the `AbstractSolutionData` implementation, each physics module must define and export a function with the signature

`function_name(mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerData, opts)`

where `EulerData` should be replaced with an appropriate type name for the physics module implementation of `AbstractSolutionData`.
This function should use the `eqn.q` array to populate `eqn.res`.
This must be done in such a way that algorithmic differentiation will work as described above.
The physics module should *never* use `q_vec` or `res_vec`.

###Initialization
The physics module should also define and export an initialization function with the signature:

`init(mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerData, opts, pmesh=mesh)`

that performs all initialization required by the module.
This function is called before the first residual evaluation, just after the `mesh`, `sbp`, and `eqn` objects are constructed.
It must also do any initialization for the objects that could not be done during their construction.
Currently, this only includes populating the `bndry_funcs` field of the `mesh` object.
This cannot be done when the mesh is constructed because the mesh object should be independent of the physics.
When using iterative solvers, there might be a second mesh object used for residual evaluations relating to the preconditioner.
This mesh is called `pmesh` and should be initialized in the same way as the regular `mesh` object.

###Initial Condition
The physics module must export a dictionary called ICDict that maps strings to the functions that apply the initial condition to the entire mesh.
These functions must have the signature:

`fname(mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerData, opts, u0::AbstractVector)`

where `u0` is a vector of length `numDof` that is populated with the initial condition.


###Callbacks
The physics module must have the `majorIterationCallback` function described in the `AbstractSolutionData` section.
This function need not be exported

A minor iteration callback function is being considered.



 

##NonlinearSolvers
The NonlinearSolvers module does not use any abstract types, but it does place some requirements on the other components of PDESolver.
It has two main components: a 4th order Runge-Kutta (RK4) explicit time marching method and a Newton's method.
The RK4 method can do either pseudo-time stepping for steady problems or real time stepping for unsteady problems.
Newton's method can solve steady nonlinear problems.

The RK4 function has the signature:

`rk4(f::Function, h::FloatingPoint, t_max::FloatingPoint, mesh::AbstractMesh, sbp, eqn::AbstractSolutionData, opts; kwargs...)`

where `f` is the residual evaluation function described above.
For RK4, `Tsol=Tres=Tmsh=Float64` because RK4 does not do any kind of algorithmic differentiation.

The Newton's method function has the signature:

`newton(func::Function, mesh::AbstractMesh, sbp, eqn::AbstractSolutionData, opts, pmesh=mesh; kwargs...)`

where `pmesh` is the preconditioning mesh object, used for constructing the preconditioner for an iterative solver.
The values of the static parameters for the `mesh` and `eqn` objects depend on the method used for calculating the jacobian `dR/dq` where `R` is the residual.
If finite differences are used, then `Tsol=Tres=Tmsh=Float64` is required.
If the complex step method is used, then `Tsol=Tres=Complex128` and `Tmsh=Float64` is required.
The startup.jl script that drives the simulation run should determine these values automatically based on the options dictionary `opts`.


##Input Options
Many of the components of PDESolver have different options that control how they work and what they do.
In order to  provide a unified method of specifying these options, an dictionary of type `Dict{Any, Any}` is read in from a disk file.
This dictionary (called `opts` in function signatures), is passed to all high and mid level function so they can use values in the dictionary to determine their control flow.
Low level functions need to be extremely efficient, so they cannot have conditional logic, therefore they are not passed the dictionary.
Note that retrieving values from a dictionary is very slow compared to accessing the fields of a type, so all values that are accessed repeatedly should be stored as the field of a type.

##Initialization of a Simulation 
Now that all the individual pieces have been described, here is how a simulation is launched.

First, the options dictionary is read in.  
Default values are supplied for any key that is not specified, if a reasonable default value exists.

Second, the `sbp` operator is constructed.

Third, the `mesh` object is constructed, using the options dictionary and the `sbp` operator.  Some of the options in the dictionary are used to determine how the mesh gets constructed.  For example, the options dictionary specifies what kind of mesh coloring to do.

Fourth, the `eqn` object is constructed, using the `mesh`, `sbp`, and `opts` objects.  

Fifth, the physics module is `init` function is called, which initializes the physics module and finishes any initialization that `mesh` and `eqn` objects require.

At this point, the procedure becomes a bit more complicated because there are optional steps.
Only the required steps are listed below.

Sixth, the initial condition is applied to `eqn.q_vec`.

Seventh, a nonlinear solver is called.  Which solver is called and what parameters is uses are determined by the options dictionary.

Eighth, post-processing is done, if required by the options dictionary.


