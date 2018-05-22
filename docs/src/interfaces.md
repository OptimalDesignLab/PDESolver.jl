# Interfaces in PDESolver

```@meta
  CurrentModule = ODLCommonTools
```

PDESolver depends on the the three main objects, the `AbstractSolutionData`
object,  `AbstractMesh` object, and the SBP object implementing certain interfaces.
This document describes what the interfaces are, and gives some hints for how to
 implement them.

Before doing so, a short description of what general Julia interfaces look like
is in order.  
The paradigm of Julia code is that of "objects with associated functions", where
 a new Type is defined, and then functions that take the Type as an argument are
defined.
The functions define the interface to the Type.
The Type holds data (ie. state), and the functions perform operations on that
state (ie. behavior).
Perhaps counter-intuitively, it is generally not recommended for users of a type
to access the fields directly.
Instead, any needed operations on the data that the Type holds should be provided
through functions.
The benefit of this convention is that it imposes no requirements on how the Type
stores its data or implements its behavior.
This is important because a user of the Type should not be concerned with these
things.
The user needs to know what behavior the Type has, but not how it is implemented.
This distinction becomes even more important when there are multiple
implementations certain functionality.
The user should be able to seamlessly transition between different implementations.
This requires all implementations have the same interface.

The question of how to enforce interfaces, and how strongly to do so, is still an open question in Julia.
Some relevant Github issues:

* [5](https://github.com/JuliaLang/julia/issues/5)
* [4935](https://github.com/JuliaLang/julia/issues/4935)
* [6975](https://github.com/JuliaLang/julia/issues/6975)


One of the strongest arguments against the "functions as interfaces" idea is
that for many applications, just storing data in an array is best.
Creating interface functions for the Type to implement the array interface would
be a lot of extra code with no benefit.
For this reason, it makes sense to directly access the fields of some Types,
to avoid trivial get/set methods.
We do this extensively in PDESolver, because arrays are the natural choice for
storing the kind of data used in PDESolver.

## AbstractSolutionData
ODLCommonTools defines:

```@docs
AbstractSolutionData
```

The purpose of an `AbstractSolutionData` is to hold all the data related to the
solution of an equation.
This includes the solution at every node and any auxiliary quantities.
The storage for any quantity that is calculated over the entire mesh should be
allocated as part of this object, in order to avoid repeatedly reallocated the
array for every residual evaluation.
In general, there should never be a need to allocate a vector longer than the
number of degrees of freedom at a node (or a matrix similarly sized matrix)
during a residual evaluation.
Structuring code such that it conforms with this requirement has significant
performance benefits because it reduces memory allocation/deallocation.

The static parameter `Tsol` is the datatype of the solution variables and `Tres`
is the datatype of the residual (when computing the Jacobian with finite
differences or algorithmic differentiation, these will be the same).

### Required Fields
The required fields of an `AbstractSolutionData` are:
```
  q::AbstractArray{Tsol, 3}
  q_vec::AbstractArray{Tsol, 1}
  shared_data::AbstractArray{SharedFaceData{Tsol}, 1}
  res::AbstractArray{Tres, 3}
  res_vec::AbstractArray{Tres, 1}
  res_edge::AbstractArray{Tres, 4}
  M::AbstractArray{Float64, 1}
  Minv::AbstractArray{Float64, 1}
  multiplyA0inv::Function
  majorIterationCallback::Function
  params{Tsol..., Tdim}::AbstractParamType{Tdim}
```

### Optional Fields
```
  file_dict::ASCIIString, IO}
```

### Field Meanings

```@meta
  CurrentModule = Utils
```


The purpose of these fields are:

`q`: to hold the solution variables in an element-based array.
     This array should be `numDofPerNode` x `numNodesPerElement` x `numEl`.
     The residual evaluation *only* uses `q`, never `q_vec`

`q_vec`: to hold the solution variables as a vector, used for any linear algebra
operations and time stepping.
This array should have a length equal to the total number of degrees of freedom
in the mesh.
Even though this vector is not used by the residual evaluation, it is needed for
many other operations, so it is allocated here so the memory can be reused.
There are functions to facilitate the scattering of values from `q_vec` to `q`.
Note that for Continuous Galerkin type discretization (as opposed to
Discontinuous Galerkin discretizations), there is not a corresponding "gather"
operation (ie. `q` -> `q_vec`).  See [`Utils.array1DTo3D`](@ref) 
and [`Utils.array3DTo1D`](@ref).

`shared_data` is a vector of length `npeers`.  Each element contains the data
              needed send and receive the `q` variables to/from other
              the corresponding MPI rank listed in `mesh.peer_parts`.
              The precise contents of `SharedFaceData` is documented in the
              `Utils` module, however they do include the send and receive
              buffers.

`res`: similar to `q`, except that the residual evaluation function populates
       it with the residual values.  
       As with `q`, the residual evaluation function only interacts with this array,
       never with `res_vec`.

`res_vec`: similar to `q_vec`.  Unlike `q_vec` there are functions to perform an
           additive reduction (basically a "gather") of `res` to `res_vec`.  
           For continuous Galerkin discretizations, the corresponding "scatter"
           (ie. `res_vec` -> `res`) may not exist.

`res_edge`: a 4 dimensional array sometimes used for edge-based data structure.
            This feature does not work.  This array must exist and should be
            of size 0.
`M`:  The mass matrix of the entire mesh.  Because SBP operators have diagonal
      mass matrices, this is a vector.  Length numDofPerNode x numNodes (where
      numNodes is the number of nodes in the entire mesh).

`Minv`:  The inverse of the mass matrix.

`multiplyA0inv`:  Multiplies the solution values at each node in an array such as `res` by the inverse of the coefficient matrix of the time term of the equation.
                  This function is used by time marching methods.
                  For some equations, this matrix is the identity matrix, so it can be a no-op, while for others might not be.
                  The function must have the signature:

`multiplyA0inv(mesh::AbstractMesh, sbp, eqn::AbstractSolutionData, opts, res_arr::AbstractArray{Tsol, 3})`


`majorIterationCallback`:  function called before every step of Newton's method
or stage of an explicit time marching scheme.
This function is used to do output and logging.
The function must have the signature:

`function majorIterationCallback(itr, mesh::AbstractMesh, sbp::AbstractSBP, eqn::AbstractEulerData, opts)`

`params`:  user defined type that inherits from `AbstractParamType`:

```@meta
  CurrentModule = ODLCommonTools
```


```@docs
AbstractParamType
```

The purpose of this type is to store any variables that need to be quickly accessed or updated.
The only required fields are:

 * `t::Float64`: hold the current time value
 * `order`: order of accuracy of the discretization (same as `AbstractMesh.order`)
 * `x_design`: vector of design variables that don't fit into any other catagory
              (ie. shape variables or other aerodynamic variables like angle of
              attack)
 *  `time::Timings`: an object to record how long different parts of the code take,
  defined in the Utils module.
 * `f`: a file handle that prints to the file `log_myrank.dat`, where `myrank`
        is the MPI rank of the current process, but only in debug mode.  In
        non-debug mode, this file should not be opened.

`file_dict`: dictionary that maps from the file name to a file handle.  This
             field is not required, but if it is present, all copies of the
             equation object must share the same dictionary (to avoid problems
             with buffering).

## AbstractMesh

ODLCommonTools defines:

```@docs
AbstractMesh
```

The purpose of an `AbstractMesh` is to hold all the mesh related data that the
 solver will need.  It also serves to establish an interface between the solver
 and whatever mesh software is used.  By storing all data in the fields of the
 `AbstractMesh` object, the details of how the mesh software stores and allows
 retrieval of data are not needed by the solver.  This should make it easy to
 accommodate different mesh software without making any changes to the solver.

The static parameter `Tmsh` is used to enable differentiation with respect to
the mesh variable in the future.

### Required Fields
```
  # counts
  numVert::Integer
  numEl::Integer
  numNodes::Integer
  numDof::Integer
  numDofPerNode::Integer
  numNodesPerElement::Integer
  order::Integer
  numNodesPerFace::Int

  # parallel counts
  npeers::Int
  numGlobalEl::Int
  numSharedEl::Int
  peer_face_counts::Array{Int, 1}
  local_element_counts::Array{Int, 1}
  remote_element_counts::array{Int, 1}

  # MPI Info
  comm::MPI.Comm
  myrank::Int
  commsize::Int
  peer_parts::Array{Int, 1}

  # Discretization type
  isDG::Bool
  isInterpolated::bool

  # mesh data
  coords::AbstractArray{Tmsh, 3}
  dxidx::AbstractArray{Tmsh, 4}
  jac::AbstractArray{Tmsh, 2}

  # boundary data
  coords_bndry::Array{Tmsh, 3}
  dxidx_bndry::Array{Tmsh, 4}
  jac_bndry::Array{Tmsh, 2}
  nrm_bndry::Array{Tmsh, 3}

  # interior face data
  coords_interface::Array{Tmsh, 3}
  dxidx_face::Array{Tmsh, 4}
  jac_face::Array{Tmsh, 2}
  nrm_face::Array{Tmsh, 3}

  # parallel data
  coords_sharedface::Array{Array{Tmsh, 3}, 1}
  dxidx_sharedface::Array{Array{Tmsh, 4}, 1}
  jac_sharedface::Array{Array{Tmsh, 2}, 1}  
  nrm_sharedface::Array{Array{Tmsh, 3}, 1}

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
  dofs::AbstractArray{Integer, 3}
  dof_offset::Int
  sparsity_bnds::AbstractArray{Integer, 2}
  sparsity_nodebnds::AbstractArray{Integer, 2}

  # mesh coloring data
  numColors::Integer
  maxColors::Integer
  color_masks::AbstractArray{ AbstractArray{Number, 1}, 1}
  shared_element_colormasks::Array{Array{BitArray{1}, 1}, 1}
  pertNeighborEls::AbstractArray{Integer, 2}

  # parallel bookkeeping
  bndries_local::Array{Array{Boundary, 1}, 1}
  bndries_remote::Array{Array{Boundary, 1}, 1}
  shared_interfaces::Array{Array{Interface, 1}, 1}
  shared_element_offsets::Array{Int, 1}
  local_element_lists::Array{Array{Int, 1}, 1}

```

TODO: update with vert_coords etc

#### Counts

`numVert`:  number of vertices in the mesh

`numEl`:  number of elements in the mesh

`numNodes`: number of nodes in the mesh

`numDof`:  number of degrees of freedom in the mesh (= `numNodes` * `numDofPerNode`)

`numDofPerNode`:  number of degrees of freedom on each node.

`numNodesPerElement`:  number of nodes on each element.

`order`:  order of the discretization (ie. first order, second order...), where
          an order `p` discretization should have a convergence rate of `p+1`.

`numNodesPerFace`: number of nodes on an edge in 2D or face in 3D.  For interpolated
                   meshes it is the number of interpolation points on the face

#### Parallel Counts
`npeers`: number of processes that have elements that share a face with the current
         process

`numGlobalEl`: number of locally owned elements + number of non-local elements that
               share a face with a locally owned element

`numSharedEl`: number of non-local elements that share a face with a locally owned
               element

`local_element_counts`: array of length `npeers`, number of local elements
                        that share a face with each remote process

`remote_element_counts`: array of length `npeers`, number of remote elements
                         that share a face with with current process

#### MPI Info
`comm`: the MPI Communicator the mesh is defined on
`myrank`: rank of current MPI process (0-based)
`commsize`: number of MPI processes on this communicator
`peer_parts`: array of MPI proccess ranks for each process that has elements that
              share a face with the current process

#### Discretization Type
`isDG`: true if mesh is a DG type discretization
`isInterpolated`: true if mesh requires data to be interpolated to the faces
                  of an element

#### Mesh Data
`coords`: `n` x `numNodesPerElement` x `numEl` array, where `n` is the dimensionality of
           the equation being solved (2D or 3D typically).  
           `coords[:, nodenum, elnum] = [x, y, z]` coordinates of node `nodenum`
           of element `elnum`.

`dxidx`:  `n` x `n` x `numNodesPerElement` x `numEl`, where `n` is defined above.
It stores the mapping jacobian scaled by `( 1/det(jac) dxi/dx )` where `xi` are
the parametric coordinates, `x` are the physical (x,y,z) coordinates, and
`jac` is the determinant of the mapping jacobian `dxi/ dx`.

`jac`  : `numNodesPerElement` x `numEl` array, holding the determinant of the
         mapping jacobian `dxi/dx` at each node of each element.

#### Boundary Data
This data is used for interpolated mesh only.

`coords_bndry`: coordinates of nodes on the boundary of the mesh,
                2 x `numFaceNodes` x `numBoundaryEdges`

`dxidx_bndry`: 2 x 2 x `numFaceNodes` x `numBoundary edges array of `dxidx`
               interpolated to the boundary of the mesh.  This fields is
               deprecated, use `nrm_bndry` instead.

`jac_bndry`: `numFaceNodes` x `numBoundaryEdges` array of `jac` interpolated
              to the boundary.  This fields is deprecated, use `nrm_bndry`
              instead.

`nrm_bndry`: `dim` x `numNodesPerFace` x `numBoundaryFaces` array containing the
             scaled face normal vector at each face node of each boundary face.
             The scaling factor is believed to be `1/|J|`, where `|J|` is the
             determinant of the mapping jacobian, however this needs to be
             verified.  The normal vector is oriented outwards.

#### Interior Face Data

`coords_face`: a `mesh.dim` x `mesh.numNodesPerFace` x `mesh.numInterfaces`
               array containing he coordinates at each face node of each
               interface.

`dxidx_face`: 2 x 2 x `numFaceNodes` x `numInterfaces` array of `dxidx`
              interpolated to the face shared between two elements.

`jac_face`: `numNodesPerFace` x `numInterfaces` array of `jac` interpolated
             to the face shared between two element

`nrm_face`: `dim` x `numNodesPerFace` x `numInterfaces` containing the scaled
            face normal vector at each face node as each interface.  The scaling
            factor is believed to be `1/|J|` where `|J|` is the determinant of 
            the mapping jacobian, however this needs to be verified.  The normal
            vector is oriented outwards from the perspective of `elementL` of
            the corresponding entry of `interfaces`.

#### Parallel Data
This data is required for parallelizing interpolated DG meshes

`coords_sharedface`: array of arrays, one array for each peer process,
                     containing the coordinates of the nodes on the faces
                     shared between a local element on a non-local element.
                     Each array is 2 x `numFaceNodes` x number of faces shared
                     with this process.
`dxidx_sharedface`: similar to `coords_sharedface`, `dxidx` interpolated to
                    faces between elements in different processes, each
                    array is 2 x 2 x `numFaceNodes` x number of faces shared
                    with this process.  This field is deprecated, use 
                    nrm_sharedface` instead.

`jac_sharedface`: similar to `coords_sharedface`, `jac` interpolated to faces
                  between a local element and a non-local element. Each array
                  is `numFaceNodes` x number of faces shared with this process.
                  This field is deprecated, use `nrm_sharedface`.

`nrm_sharedface`; similar to `coords_sharedface`, the inner array contains the
                  scaled normal vector at each node of each shared face.  Each
                  inner array is `dim` x `numNodesPerFace` x number of faces 
                  shared with this process.  The normal vector is oriented
                  outward from the perspective of the *local* element.

#### Boundary Condition Data
The mesh object stores data related to applying boundary conditions.
Boundary conditions are imposed weakly, so there is no need to remove degrees of
freedom from the mesh when Dirichlet boundary conditions are applied.
In order to accommodate any combination of boundary conditions, an array of
functors are stored as part of the mesh object, along with lists of which mesh
edges (or faces in 3D) should have which boundary condition applied to them


`numBC`: number of different types of boundary conditions used.

`numBoundaryEdges`: number of mesh edges that have boundary conditions applied
                    to them.

`bndryfaces`:  array of Boundary objects (which contain the element number and
               the local index of the edge), of length `numBoundaryEdges`.

`bndry_offsets`:  array of length numBC+1, where `bndry_offsets[i]` is the index
                  in `bndryfaces` where the edges that have boundary condition
                  `i` applied to them start.
                  The final entry in `bndry_offsets` should be `numBoundaryEdges + 1`.
                  Thus `bndryfaces[ bndry_offsets[i]:(bndry_offsets[i+1] - 1) ]`
                  contains all the boundary edges that have boundary condition
                  `i` applied to them.

`bndry_funcs`:  array of boundary functors, length `numBC`.  All boundary
                functors are subtypes of `BCType`.  Because `BCType` is an
                abstract type, the elements of this array should not be used
                directly, but passed as an argument to another function, to
                 avoid type instability.

#### Interior Edge Data
Data about interior mesh edges (or faces in 3D) is stored to enable use of
edge stabilization or Discontinuous Galerkin type discretizations.
Only data for edges (faces) that are shared by two elements are stored
(ie. boundary edges are not considered).

`numInterfaces`:  number of interior edges

`interfaces`:  array of Interface types (which contain the element numbers for
               the two elements sharing the edge, and the local index of the
               edge from the perspective of the two elements, and an indication
               of the relative edge orientation).
               The two element are referred to as `elementL` and `elementR`,
               but the choice of which element is `elementL` and which is
               `elementR` is arbitrary.
               The length of the array is numInterfaces.
              Unlike `bndryfaces`, the entries in the array do not have to be in
              any particular order.

#### Degree of Freedom Numbering Data
`dofs`:  `numDofPerNode` x `numNodesPerElement` x `numEl` array.
Holds the local degree of freedom number of each degree of freedom.
Although the exact method used to assign dof numbers is not critical, all
degrees of freedom on a node must be numbered sequentially.

`dof_offset`: offset added to the local dof number to make it a global dof
              number.

`sparsity_bnds`:  2 x `numDof` array.
`sparsity_bnds[:, i]` holds the maximum, minimum degree of freedom numbers
associated with degree of freedom `i`.
In this context, degrees of freedom `i` and `j` are associated if entry `(i,j)`
of the jacobian is non-zero.
In actuality, `sparsity_bnds` need only define upper and lower bounds for degree
of freedom associations (ie. they need not be tight bounds).
This array is used to to define the sparsity pattern of the jacobian matrix.

`sparsity_nodebnds`:  2 x numNodes array.
`sparsity_bnds[:, i]` holds the maximum, minimum node associated with node `i`,
similar the information stored in `sparsity_bnds` for degrees of freedom.


#### Mesh Coloring Data
The NonlinearSolvers module uses algorithmic differentiation to compute the
Jacobian.
Doing so efficiently requires perturbing multiple degrees of freedom
simultaneously, but perturbing associated degrees of freedom at the same time
leads to incorrect results.
Mesh coloring assigns each element of the mesh to a group (color) such that
every degree of freedom on each element is not associated with any other degree
 of freedom on any other element of the same color.
An important aspect of satisfying this condition is the use of the
 element-based arrays (all arrays that store data for a quantity over the entire
 mesh are `ncomp` x `numNodesPerElement` x `numEl`).
In such an array, any node that is part of 2 or more elements has one entry for
each element.
When performing algorithmic differentiation, this enables perturbing a degree of
 freedom on one element without perturbing it on the other elements that share
 the degree of freedom.

For example, consider a node that is shared by two elements.
Let us say it is node 2 of element 1 and node 3 of element 2.
This means `AbstractSolutionData.q[:, 2, 1]` stores the solution variables for
this node on the first element, and `AbstractSolutionData.q[:, 3, 2]` stores the
solution variables for the second element.
Because these are different entries in the array `AbstractSolutionData.q`,
they can be perturbed independently.
Because `AbstractSolutionData.res` has the same format, the perturbations to
`AbstractSolutionData.q[:, 2, 1]` are mapped to
`AbstractSolutionData.res[:, 2, 1]` for a typical continuous Galerkin type
discretization.
This is a direct result of having an element-based discretization.

There are some discretizations, however, that are not strictly element-based.
Edge stabilization, for example, causes all the degrees of freedom of one
element to be associated with any elements it shares an edge with.
To deal with this, we use the idea of a distance-n coloring.
A distance-n coloring is a coloring where there are n elements in between two
element of the same color.
For element-based discretizations with element-based arrays, every element in
the mesh can be the same color.  This is a distance-0 coloring.
For an edge stabilization discretization, a distance-1 coloring is required,
where every element is a different color than any neighbors it shares and edge
with.
(As a side node, the algorithms that perform a distance-1 coloring are rather
complicated, so in practice we use a distance-2 coloring instead).

In order to do algorithmic differentiation, the `AbstractMesh` object must store
the information that determines which elements are perturbed for which colors,
and, for the edge stabilization case, how to relate a perturbation in the output
 `AbstractSolutionData.res` to the degree of freedom in `AbstractSolutionData.q`
in O(1) time.
Each degree of freedom on an element is perturbed independently of the other
degrees of freedom on the element, so the total number of residual evaluations
is the number of colors times the number of degrees of freedom on an element.

The fields required are:

`numColors`:  The number of colors in the mesh.
`maxColors`: the maximum number of colors on any process

`color_masks`:  array of length `numColors`.  Each entry in the array is itself
                an array of length `numEl`.  Each entry of the inner array is
                either a 1 or a 0, indicating if the current element is
                perturbed or not for the current color.
For example, in `color_mask_i = color_masks[i]; mask_elj = color_mask_i[j]`,
the variable `mask_elj` is either a 1 or a zero, determining whether or not
element `j` is perturbed as part of color `i`.

`shared_element_colormasks`: array of BitArrays controlling when to perturb
                             non-local elements.  There are `npeers` arrays,
                             each of length number of non-local elements shared
                             with this process

`pertNeighborEls`:  `numEl` x `numColors` array.  `neighbor_nums[i,j]` is the
element number of of the element whose perturbation is affected element `i`
when color `j` is being perturbed, or zero if element `i` is not affected by any
 perturbation.  

#### Parallel Bookkeeping
`bndries_local`: array of arrays of `Boundary`s describing faces shared
                 with non-local elements from the local side (ie. the
                 local element number and face).  The number of arrays is
                 `npeers`.

`bndries_remote`: similar to `bndries_local`, except describing the faces
                  from the non-local side.  Note that the element numbers are
                  from the *remote* process

`shared_interfaces`: array of arrays of `Interface`s describing the faces
                     shared between local and non-local elements.  The local
                     element is *always* `elementL`.  The remote element is
                     assigned a local number greater than `numEl`.

`shared_element_offsets`: array of length npeers+1 that contains the first
                          element number assigned to elements on the shared
                          interface.  The last entry is equal to the highest
                          element number on the last peer process + 1.  The
                          elements numbers assigned to a given peer must form
                          a contiguous range.

`local_element_lists`: array of arrays containing the element numbers of the
                       elements that share a face with each peer process
### Other Functions
The mesh module must also defines and exports the functions


```
saveSolutionToMesh(mesh::MeshImplementationType, vec::AbstractVector)
writeVisFiles(mesh::MeshImplementationType, fname::ASCIIString)
```

where the first function takes a vector of length `numDof` and saves it to the
mesh, and the second writes Paraview files for the mesh, including the solution
field.

### Subtypes

```@docs
AbstractCGMesh
AbstractDGMesh
```

## Physics Module
For every new physics to be solved, a new module should be created.
The purpose of this module is to evaluate the equation:

`M dq/dt = f(q)`

where `M` is the mass matrix.
For steady problems, `dq/dt = 0` and the module evaluates the residual.
For unsteady problems, the form `M dq/dt = f(q)` is suitable for explicit time
marching.

Every physics module must define a boundary condition called `defaultBC`.  If
the user does not specify a boundary condition on any geometric edges, the
mesh constructor will add a new boundary condition and assign all mesh edges
classified on the unspecified geometric edges are assigned to it.  This boundary
condition can be a no-op if that is correct for the physics, or it can be
the face integral that should be done for every element (recall that boundary
faces do not appear in `mesh.interfaces`).

### Interface to NonlinearSolvers
The `evalResidual` function and the fields `eqn.q` and `eqn.res` are the
interface between the NonlinearSolvers and the physics modules.  The Nonlinear
solvers populate `eqn.q`, and use `evalResidual` to populate `eqn.res`, from
which the next value if `eqn.q` is calculated.  The algorthmic differentiation
mechanism described above is uses several residual evaluations to compute the
Jacobian if needed for a given method.  Some methods, such as RK4, are better
expressed in terms of `eqn.q_vec` and `eqn.res_vec` rather than `eqn.q` and
`eqn.res`.  `eqn.array3DTo1D` and `eqn.disassmbleSolution` exist to
transform back and forth between the vector and 3D array forms.  In order to
compute the Jacobian efficiently, it is necessary to work with the 3D arrays.
For this reason, `evalResidual` must work only with `eqn.q` and `eqn.res`
and let the caller decide whether or not to transform into the vector form.

Newton's method supports both finite differences and complex step for
calculating the Jacobian, and the static parameters need to be set accordingly.
If finite differences are used, then `Tsol=Tres=Tmsh=Float64` is required.  If
complex step is used, then `Tsol=Tres=Complex128` and `Tmsh = Float64` is needed.

### Interface to users
The interface to users is described in `src/Readme.md`

### Interface to Summation-by-Parts
The physics modules use the interface provided by the Summation-by-Parts package
to approximate derivatives numerically.  The reason for passing around the
`sbp` object is that the SBP interface requires it to be passed in along with
the data arrays it operates on.


## Functional Programming
An important aspect of the use of the `mesh`, `sbp`, `eqn`, and `opts` to define
interfaces is that the physics modules and nonlinear solvers are written in a
purely functional programming style, which is to say that the behavior of every
function is determined entirely by the arguments to the function, and the only
effects of the function are to modify an argument or return a value.

This property is important for writing generic, reusable code.
For example, when using iterative solvers, a preconditioner is usually required,
and constructing the preconditioner requires doing residual evaluations.
In some cases, the preconditioner will use a different mesh or a different mesh
coloring.
Because the physics modules are written in a functional style, they can be used
to evaluate the residual on a different mesh simply by passing the residual
evaluation function a different mesh object.

A critical aspect of function programming is that there is *no global state*.
The state of the solver is determined entirely by the state of the objects that
are passed around.


## Variable Naming Conventions
In an attempt to make code more uniform and readable, certain variable names are
reserved for certain uses.

* `mesh`:  object that implements `AbstractMesh`
* `pmesh`:  mesh used for preconditioning
* `sbp`:  Summation-by-Parts object
* `eqn`:  object that implements `AbstractSolutionData`
* `opts`: options dictionary
* `params`: parameter object (used for values that might be in `opts` but need to be accessed quickly)
* `x`: the first real coordinate direction
* `y`: the second real coordinate direction
* `xi`: the first parametric coordinate direction
* `eta`: the second parametric coordinate direction
* `h`:  mesh spacing
* `hx`: mesh spacing in x direction
* `hy`: mesh spacing in y direction
* `p`: pressure at a node
* `a`; speed of sound at a node
* `s`: entropy at a node
* `gamma`: specific heat ratio
* `gamma_1`: `gamma` - 1
* `R`: specific gas constant in ideal gas law (units J/(Kg * K) in SI)
* `delta_t`: time step
* `t`: current time
* `nrm`: a normal vector of some kind
* `A0`: the coefficient matrix of the time term at a node
* `Axi`: the flux jacobian at a node in the `xi` direction
* `Aeta`: the flux jacobian at a node in the `eta` direction
* `Ax`: the flux jacobian in the `x` direction at a node
* `Ay`: the flux jacobian in the `y` direction at a node
