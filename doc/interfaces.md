# Interfaces in PDESolver
PDESolver depends on the the three main objects, the AbstractSolutionData object,  AbstractMesh object, and the SBP object implementing certain interfaces.
This document describes what the interfaces are, and gives some hints for how to implement them.

Before doing so, a short description of what general Julia interfaces look like is in order.  
The paradigm of Julia code is that of "objects with associated functions", where a new Type is defined, and then functions that take the Type as an argument are defined.
The functions define the interface to the Type.
The Type holds data (ie. state), and the functions perform operations on that state (ie. behavior).
Perhaps counter-intuitively, it is generally not recommended for users of a type to  acess the fields directly.
Instead, any needed operations on the data that the Type holds should be provided through functions.
The benefit of this convention is that it imposes no requirements on how the Type stores its data or the implementations its behavior.
This is important because a user of the Type should not be concerned with these things.
The user needs to know what behavior the Type has, but not how it is implemented.
This becomes even more important when there are multiple implementations certain functionality.
The user should be able to seemlessly transition between different implementations.
This requires all implementations have the same interface.

The question of how to enforce interfaces, and how strongly to do so, is still an open question in Julia.
Some relevent Github issues:
5
4935
6975

One of the strongest arguments against the "functions as interfaces" idea is that for many applications, just storing things in an array is best.
Creating interface functions for the Type to implement the array interface would be a lot of extra code with no benefit.
For this reason, it makes sense to directly access the fields of some Types, to avoid trivial get/set methods.
We do this extensively in PDESolver, because arrays are the natural choice for storing data.

## AbstractSolutionData
ODLCommonTools defines:

`abstract AbstractSolutionData{Tsol}`.

The purpose of an `AbstractSolutionData` is to hold all the data related to the solution of an equation.
This includes the solution at every node and any auxiliary quantities.
The storage for any quantity that is calculated over the entire mesh should be allocated as part of this object, in order to avoid repeatedly reallocated the array for every residual evaluation.
In general, there should never be a need to allocate a vector longer than the number of degrees of freedom at a node (or a matrix similarly sized matrix) during a residual evaluation.
If it seems like this is necessary, reconsider the structure of the code you are writing.
There are significiant performance benefits to this rule because it reduces memory allocation.

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
```

The purpose of these fields are:
`q`: to hold the solution variables in an element based array. 
     This array should be numDofPerNode x numNodesPerElement x numEl.
     The residual evaluation *only* uses `q`, never `q_vec`

`q_vec`: to hold the solution variables as a vector, used for any linear algebra operations and timestepping.
         This array should have a length equal to the total number of degrees of freedom in the mesh.
         There are functions to facilitate the scattering of values from `q_vec` to `q`.
         Note that for Continuous Galerkin type discretization (as opposed to Discontinuous Galerkin discretizations), there is not a corresponding "gather" operation (ie. `q` -> `q_vec`).

`res`: similar to `q`, except that the residual evaluation function populates it with the residual values.  
       As with `q`, the residual evaluation function only interacts with this array, never with `res_vec`.

`res_vec`: similar to `q_vec`.  Unlike `q_vec` there are functions to perform an additive reduction (basically a "gather") of `res` to `res_vec`.  For continuous Galerkin discretizations, the corresponding "scatter" (ie. `res_vec` -> res`) may not exist.

`M`:  The mass matrix of the entire mesh.  Because SBP operators have diagonal mass matrices, this is a vector.  Length numDofPerNode x numNodes (where numNodes is the number of nodes in the entire mesh).

`Minv`:  The inverse of the mass matrix.

`disassembleSolution`:  Function that takes the a vector such as `q_vec` and scatters it to an array such as `q`.
                        This function must have the signature:
                        `disassembleSolution(mesh::AbstractMesh, sbp, eqn::AbstractSolutionData, opts, q_arr:AbstractAray{T, 3}, q_vec::AbstractArray{T, 1}`
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



##AbstractMesh
ODLCommonTools defines:

`abstract AbstractMesh{Tmsh}`.

The purpoose of an `AbstractMesh` is to hold all the mesh related data that the solver will need.
It also serves to establish an interface between the solver and whatever mesh software is used.
By storing all data in the fields of the `AbstractMesh` object, the details of how the mesh software stores and allows retrieval of data are not needed by the solver.
This should make it easy to accomodate different mesh softwares without making any changes to the solver.

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
  neighbor_nums::AbstractArray{Integer, 2}
  color_masks::AbstractArray{ AbstractArray{Number, 1}, 1}, 
  pertNeighborEls::AbstractArray{Integer, 2}
```
####Counts
`numVert`:  number of vertices in the mesh
`numEl`:  number of elements in the mesh
`numNodes`: number of nodes in the mesh
`numDof`:  number of degrees of freedom in the mesh (= numNodes * numDofPerNode)
`numDofPerNode`:  number of degrees of freedom on each node.
`numNodesPerElement`:  number of nodes on each element.
`order`:  order of the discretization (ie. first order, second order...), where an order `p` discretization should have a convergence rate of `p+1`.


####Mesh Data
`coords`: n x numNodesPerElement x numEl array, where n is the dimensionality of   the equation being solved (2D or 3D typically).  `coords[:, nodenum, elnum] = [x, y, z]` coordinates.

`dxidx`:  n x n x numNodesPerElement x numEl, where n is defined above.  It stores the mapping jacobian scaled by `( 1/det(jac) dxi/dx )` where xi are the parametric coordinates, x are the physical (x,y,z) coordinates, and jac is the determinant of the mapping jacobian `dxi/ dx`.

`jac`  : numNodesPerElement x numEl array, holding the determinant of the mapping jacobian `dxi/dx` at each node of each element.


####Boundary Condition Data
The mesh object stores data related to applying boundary conditions.  Boundary conditions are imposed weakly, so there is no need to remove degrees of freedom from the mesh when Dirchlet boundary conditions are applied.
In order to accomodate any combination of boundary conditions, an array of functors are as part of the mesh object, along with lists of which mesh edges (or faces in 3D) should have which boundary condition applied to them


`numBC`: number of different types of boundary conditions used.
`numBoundaryEdges`: number of mesh edges that have boundary conditions applied to them.

`bndryfaces`:  array of Boundary objects (which contain the element number and the local index of the edge), of length `numBoundaryEdges`.
`bndry_offsets`:  array of length numBC+1, where `bndry_offsets[i]` is the index  in `bndryfaces` where the edges that have boundary condition `i` applied to them start.  The final entry should be `numBoundaryEdges + 1`.  Thus `bndryfaces[ bndry_offsets[i]:(bndry_offsets[i+1] - 1) ]` contains all the boundary edges that have boundary condition `i` applied to them.

`bndry_funcs`:  array of boundary functors, length numBC.  All boundary functors are subtypes of `BCType`.  Because `BCType` is an abstract type, the elements of this array should not be used directly, but passed as an argument to another function, to avoid type instability.

####Interior Edge Data
Data about interior mesh edges (or faces in 3D) is stored to enable use of edge stabilization or Discontinuous Galerkin type discretization.  Only data for edges (faces) that are shared by two elements is stored (ie. boundary edges are not considered).

`numInterfaces`:  number of interior edges
`interfaces`:  array of Interface types (which contain the element numbers for the two elements sharing the edge, and the local index of the edge from the perspective of the two elements, and an indication of the relative edge orientation).  The two element are referred to as `elementL` and `elementR`, but the choice of which element is `elementL` and which is `elementR` is arbitrary.  The length of the array is numInterfaces.  Unlike `bndryfaces`, the entries in the array do not have to be in any particular order.

####Degree of Freedom Numbering Data
`dofs`:  numDofPerNode x numNodesPerElement x numEl array.  Holds the degree of freedom number of each degree of freedom.

`sparsity_bnds`:  2 x numDof array.  `sparsity_bnds[:, i]` holds the maximum, minimum degree of freedom numbers assoicated with degree of freedom `i`.  In this context, degrees of freedom `i` and `j` are associated if entry `(i,j)` of the jacobian is non-zero.  In actuality, `sparsity_bnds` need only define upper and lower bounds for degree of freedom associations (ie. they need not be tight bounds).  This array is used to to define the sparsity pattern of the jacobian matrix.

`sparsity_nodebnds`:  2 x numNodes array.  `sparsity_bnds[:, i]` holds the maximum, minimum node associated with node `i`, similar the information stored in `sparsity_bnds` for degrees of freedom.


####Mesh Coloring Data
The NonlinearSolvers module uses algorithmic differentiation to compute the Jacobian.  Doing so efficiently requires perturbing multiple degrees of freedom simultaniously, but perturbing associated degrees of freedom at the same time leads to incorrect results.  Mesh coloring assigns each element of the mesh to a group (color) such that every degrees of freedom on each element is not associated with any other degree of freedom on any other element.
