var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "PDESolver Introduction",
    "title": "PDESolver Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#PDESolver-Documentation-1",
    "page": "PDESolver Introduction",
    "title": "PDESolver Documentation",
    "category": "section",
    "text": "Welcome to the PDESolver documentation.  These documents provide an overview of PDESolver, a Julia based solver for partial differential equations. This page will describe the form of the equations and the Summation-By-Parts operators used to discretize them.  The Table of Contents links to the components of PDESolver that calculate each term."
},

{
    "location": "index.html#Form-of-the-Equation-1",
    "page": "PDESolver Introduction",
    "title": "Form of the Equation",
    "category": "section",
    "text": "PDESolver discretizes equation in the form:fracpartial qpartial t = mathcalR(q t)where q is the variable being solved for, mathcalR(u t) is called the residual and t is time. The residual contains all the spatial derivatives.  For example, the residual for the 1D advection equation is $ a\\frac{\\partial q}{\\partial x} where a is the advection velocity. The code is capable of solving both unsteady problems, in the form shown above, and steady problem of the formmathcalR(q) = 0The code is structured such that the physics modules are responsible for  evaluating the residual and the Nonlinear Solvers are responsible for the time term in the unsteady case, or solving the nonlinear rootfinding problem for steady cases."
},

{
    "location": "index.html#index_sbp-1",
    "page": "PDESolver Introduction",
    "title": "Summation-by-Parts operators",
    "category": "section",
    "text": "Summation-by-Parts (SBP) operators are used to discretize the spatial derivatives in the residual. In particular, we use the multi-dimensional Summation-by-Parts operators first defined in   Hicken, J.E., Del Rey Fernandez, D.C, and Zingg, D.W., \"Multi-dimensional \n  Summation-by-Parts Operators: General Theory and Application to \n  Simplex Elements\", SIAM Journal on Scientific Computing, Vol, 38, No. 4, 2016,\n   pp. A1935-A1958See the introductor PDF (to be posted shortly) for an introduction to the operators."
},

{
    "location": "index.html#Table-of-Contents-1",
    "page": "PDESolver Introduction",
    "title": "Table of Contents",
    "category": "section",
    "text": "Pages = [\"invocation/calling.md\",\n         \"solver/Readme.md\",\n         \"solver/advection/advection.md\",\n         \"solver/euler/euler.md\",\n         \"solver/simpleODE/simpleODE.md\",\n         \"NonlinearSolvers/nonlinearsolvers.md\",\n         \"input/input.md\",\n         \"Utils/Utils.md\"\n        ]\nDepth=1"
},

{
    "location": "build.html#",
    "page": "Building PDESolver",
    "title": "Building PDESolver",
    "category": "page",
    "text": ""
},

{
    "location": "build.html#Building-PDESolver-1",
    "page": "Building PDESolver",
    "title": "Building PDESolver",
    "category": "section",
    "text": "This page describes how to build the PDESolver package. Julia has two kinds of dependencies, system packages and Julia packages. System packages can either be installed using a Linux package manager (such as apt-get on Debian-based system) or built from source. Julia packages are install using the Julia package manager."
},

{
    "location": "build.html#System-Dependencies-1",
    "page": "Building PDESolver",
    "title": "System Dependencies",
    "category": "section",
    "text": "Before installing PDESolver, you should already have the following installed:Julia v0.4\nC/C++/Fortran compilers\nAn MPI implementation (including MPI compiler wrappers)\nCMake v3.0.0 or later\nBLAS and LAPACKThe build process for Julia itself can be found here.If you are on a Debian based system, the following command will install the remaining dependencies:  sudo apt-get install build-essential gfortran cmake libblas-dev liblapack-dev mpich"
},

{
    "location": "build.html#PDESolver-Installation-1",
    "page": "Building PDESolver",
    "title": "PDESolver Installation",
    "category": "section",
    "text": "PDESolver is not listed in Julia's METADATA, so you will have to clone the repository and then build the code.After installing the system System Dependencies, run the following Julia commands to install the Julia dependencies and and PDESolver itself:  Pkg.clone(\"https://github.com/OptimalDesignLab/PDESolver.jl.git\")\n  Pkg.resolve()  # install all packages in PDESolvers REQUIRE file\n  Pkg.build(\"PDESolver\")  # install all packages not in REQUIRE file and\n                          # build PDESolverThis will install PDESolver and all Julia dependencies into the directory specified by the JULIA_PKGDIR environmental variable. If this variable is not specified, it will install to ~/.julia/v0.4/PDESolver.If there are errors building any Julia dependencies, see the Installation page for methods of installing particular versions of the dependencies.After installation it is recommended to run the test suite. To do so, run the following commands in the terminal:  cd /path/to/pdesolver # (for example ~/.julia/v0.4/PDESolver)\n  cd ./test\n  ./runtests_fast.shIf the tests complete without error, then the package is properly installed.TODO: link to examples page"
},

{
    "location": "deps_readme.html#",
    "page": "Build Options",
    "title": "Build Options",
    "category": "page",
    "text": ""
},

{
    "location": "deps_readme.html#Installation-1",
    "page": "Build Options",
    "title": "Installation",
    "category": "section",
    "text": "PDESolver depends on several packages, some registered Julia packages, others  not."
},

{
    "location": "deps_readme.html#Regular-Packages-1",
    "page": "Build Options",
    "title": "Regular Packages",
    "category": "section",
    "text": "Like any Julia package, most of the dependencies that are registered packages are listed in the REQUIRE file along with version numbers known to work. The Julia package manager is used to resolve any dependency conflicts with  other installed packages."
},

{
    "location": "deps_readme.html#Non-Standard-Packages-1",
    "page": "Build Options",
    "title": "Non-Standard Packages",
    "category": "section",
    "text": "The dependencies that are not registered packages are installed by manually  cloning their repositories and building them.  Commit hashes are used to  identify versions known to work.  The commits are checked out into a branch called pdesolver_version.If a package is already installed, by default it will left alone (ie. the  specified commit will not be checked out).  This behavior can be overrided  by declaring certain environmental variables in the shell before launching  Julia (see below)"
},

{
    "location": "deps_readme.html#Manual-Installation-Process-1",
    "page": "Build Options",
    "title": "Manual Installation Process",
    "category": "section",
    "text": "The packages listed in the REQUIRE file can be installed manually, by  cloning the repository and checkout out a particular commit known to work. This will force the checking out of he commit and re-building of the package  even if it is already installed."
},

{
    "location": "deps_readme.html#Environmental-Variables-1",
    "page": "Build Options",
    "title": "Environmental Variables",
    "category": "section",
    "text": "PDESOLVER_INSTALL_DEPS_MANUAL: install the packages in the REQUIRE file manually, forcefully. PDESOLVER_FORCE_DEP_INSTALL_ALL: forces the checkout and re-installation   of the non-standard packages, even if already installed PDESOLVER_FORCE_DEP_INSTALL_pkg_name: force the checkout and re-installation of the package named pkg_name.For all these environmental variables, the value is not important, only the  existance of the variable."
},

{
    "location": "deps_readme.html#Offline-Installation-1",
    "page": "Build Options",
    "title": "Offline Installation",
    "category": "section",
    "text": "If the internet is not accessible from the machine where the code is to be  installed, it is possible to download the packages to a specified  directory on another machine and then copy them to the machine where they will be installed.  To do this, set the envinronmental variable  PDESOLVER_BUNDLE_DEPS to tell the build script to perform bundling instead of  installation, and set the path to the directory where the deps should be stored in the variable PDESOLVER_PKGDIR.After running the build script, copy the contents of the PDESOLVER_PKGDIR to the package directory on the target machine, set PDESOLVER_UNBUNDLE_DEPS, and run the PDESolver/build.jl."
},

{
    "location": "deps_readme.html#Logging-1",
    "page": "Build Options",
    "title": "Logging",
    "category": "section",
    "text": "The file deps/install.log is created and written to with the progress of the checkout process, including error information if a build fails."
},

{
    "location": "interfaces.html#",
    "page": "Code Interfaces",
    "title": "Code Interfaces",
    "category": "page",
    "text": ""
},

{
    "location": "interfaces.html#Interfaces-in-PDESolver-1",
    "page": "Code Interfaces",
    "title": "Interfaces in PDESolver",
    "category": "section",
    "text": "  CurrentModule = ODLCommonToolsPDESolver depends on the the three main objects, the AbstractSolutionData object,  AbstractMesh object, and the SBP object implementing certain interfaces. This document describes what the interfaces are, and gives some hints for how to  implement them.Before doing so, a short description of what general Julia interfaces look like is in order.   The paradigm of Julia code is that of \"objects with associated functions\", where  a new Type is defined, and then functions that take the Type as an argument are defined. The functions define the interface to the Type. The Type holds data (ie. state), and the functions perform operations on that state (ie. behavior). Perhaps counter-intuitively, it is generally not recommended for users of a type to access the fields directly. Instead, any needed operations on the data that the Type holds should be provided through functions. The benefit of this convention is that it imposes no requirements on how the Type stores its data or implements its behavior. This is important because a user of the Type should not be concerned with these things. The user needs to know what behavior the Type has, but not how it is implemented. This distinction becomes even more important when there are multiple implementations certain functionality. The user should be able to seamlessly transition between different implementations. This requires all implementations have the same interface.The question of how to enforce interfaces, and how strongly to do so, is still an open question in Julia. Some relevant Github issues:5\n4935\n6975One of the strongest arguments against the \"functions as interfaces\" idea is that for many applications, just storing data in an array is best. Creating interface functions for the Type to implement the array interface would be a lot of extra code with no benefit. For this reason, it makes sense to directly access the fields of some Types, to avoid trivial get/set methods. We do this extensively in PDESolver, because arrays are the natural choice for storing the kind of data used in PDESolver."
},

{
    "location": "interfaces.html#ODLCommonTools.AbstractSolutionData",
    "page": "Code Interfaces",
    "title": "ODLCommonTools.AbstractSolutionData",
    "category": "Type",
    "text": "This abstract type is the supertype for all the objects that store the    solution data. Every physics module should implement its own subtype.\n\nStatic parameters:\n\nTsol: datatype of solution variables\nTres: datatype of the mesh variables\n\nSee the AbstractSolutionData for the description of everything this   type must implement.\n\n\n\n"
},

{
    "location": "interfaces.html#AbstractSolutionData-1",
    "page": "Code Interfaces",
    "title": "AbstractSolutionData",
    "category": "section",
    "text": "ODLCommonTools defines:AbstractSolutionDataThe purpose of an AbstractSolutionData is to hold all the data related to the solution of an equation. This includes the solution at every node and any auxiliary quantities. The storage for any quantity that is calculated over the entire mesh should be allocated as part of this object, in order to avoid repeatedly reallocated the array for every residual evaluation. In general, there should never be a need to allocate a vector longer than the number of degrees of freedom at a node (or a matrix similarly sized matrix) during a residual evaluation. Structuring code such that it conforms with this requirement has significant performance benefits because it reduces memory allocation/deallocation.The static parameter Tsol is the datatype of the solution variables and Tres is the datatype of the residual (when computing the Jacobian with finite differences or algorithmic differentiation, these will be the same)."
},

{
    "location": "interfaces.html#ODLCommonTools.AbstractParamType",
    "page": "Code Interfaces",
    "title": "ODLCommonTools.AbstractParamType",
    "category": "Type",
    "text": "This abstract type is the supertype for all Param objects, which hold values    needed for the computation in a place that is fast to access.\n\nThe Param type is also useful for dispatching to low level functions which     the AbstractSolutionData might not be passed (depending on the organization     of the physics module.\n\n\n\n"
},

{
    "location": "interfaces.html#Required-Fields-1",
    "page": "Code Interfaces",
    "title": "Required Fields",
    "category": "section",
    "text": "The required fields of an AbstractSolutionData are:  q::AbstractArray{Tsol, 3}\n  q_vec::AbstractArray{Tsol, 1}\n  shared_data::AbstractArray{SharedFaceData{Tsol}, 1}\n  res::AbstractArray{Tres, 3}\n  res_vec::AbstractArray{Tres, 1}\n  M::AbstractArray{Float64, 1}\n  Minv::AbstractArray{Float64, 1}\n  disassembleSolution::Function\n  assembleSolution::Function\n  multiplyA0inv::Function\n  majorIterationCallback::Function\n  params{Tsol..., Tdim}::AbstractParamType{Tdim}The purpose of these fields are:q: to hold the solution variables in an element-based array.      This array should be numDofPerNode x numNodesPerElement x numEl.      The residual evaluation only uses q, never q_vecq_vec: to hold the solution variables as a vector, used for any linear algebra operations and time stepping. This array should have a length equal to the total number of degrees of freedom in the mesh. Even though this vector is not used by the residual evaluation, it is needed for many other operations, so it is allocated here so the memory can be reused. There are functions to facilitate the scattering of values from q_vec to q. Note that for Continuous Galerkin type discretization (as opposed to Discontinuous Galerkin discretizations), there is not a corresponding \"gather\" operation (ie. q -> q_vec).shared_data is a vector of length npeers.  Each element contains the data               needed send and receive the q variables to/from other               the corresponding MPI rank listed in mesh.peer_parts.               The precise contents of SharedFaceData is documented in the               Utils module, however they do include the send and receive               buffers.res: similar to q, except that the residual evaluation function populates        it with the residual values.          As with q, the residual evaluation function only interacts with this array,        never with res_vec.res_vec: similar to q_vec.  Unlike q_vec there are functions to perform an            additive reduction (basically a \"gather\") of res to res_vec.              For continuous Galerkin discretizations, the corresponding \"scatter\"            (ie. res_vec -> res`) may not exist.M:  The mass matrix of the entire mesh.  Because SBP operators have diagonal       mass matrices, this is a vector.  Length numDofPerNode x numNodes (where       numNodes is the number of nodes in the entire mesh).Minv:  The inverse of the mass matrix.disassembleSolution:  Function that takes the a vector such as q_vec and                         scatters it to an array such as q.                         This function must have the signature:                         disassembleSolution(mesh::AbstractMesh, sbp, eqn::AbstractSolutionData, opts, q_arr:AbstractArray{T, 3}, q_vec::AbstractArray{T, 1})                         Because this variable is a field of a type, it will be dynamically dispatched.                         Although this is slower than compile-time dispatch, the cost is insignificant compared to the cost of evaluating the residual, so the added flexibility of having this function as a field is worth the cost.assembleSolution:  Function that takes an array such as res and performs an additive reduction to a vector such as res_vec.                      This function must have the signature:                      assembleSolution(mesh::AbstractMesh, sbp, eqn::AbstractSolutionData, opts, res_arr::AbstractArray{T, 3}, res_vec::AbstractArray{T, 1}, zero_resvec=true)                      The argument zero_resvec determines whether res_vec is zeroed before the reduction is performed.                      Because it is an additive reduction, elements of the vector are only added to, never overwritten, so forgetting to zero out the vector could cause strange results.                      Thus the default is true.multiplyA0inv:  Multiplies the solution values at each node in an array such as res by the inverse of the coefficient matrix of the time term of the equation.                   This function is used by time marching methods.                   For some equations, this matrix is the identity matrix, so it can be a no-op, while for others might not be.                   The function must have the signature:multiplyA0inv(mesh::AbstractMesh, sbp, eqn::AbstractSolutionData, opts, res_arr::AbstractArray{Tsol, 3})majorIterationCallback:  function called before every step of Newton's method or stage of an explicit time marching scheme. This function is used to do output and logging. The function must have the signature:function majorIterationCallback(itr, mesh::AbstractMesh, sbp::AbstractSBP, eqn::AbstractEulerData, opts)params:  user defined type that inherits from AbstractParamType:AbstractParamTypeThe purpose of this type is to store any variables that need to be quickly accessed or updated. The only required fields are: * t::Float64: hold the current time value * order: order of accuracy of the discretization (same as AbstractMesh.order) *  time::Timings: an object to record how long different parts of the code take,   defined in the Utils module."
},

{
    "location": "interfaces.html#ODLCommonTools.AbstractMesh",
    "page": "Code Interfaces",
    "title": "ODLCommonTools.AbstractMesh",
    "category": "Type",
    "text": "This abstract type is the supertype for all mesh objects.  Every interface to   a mesh software should define its own implementation.\n\nStatic parameters:\n\nTmsh: datatype of the mesh data (coordinates, mapping to/from parametric\n      space, mapping jacobian).\n\nSee the AbstractMesh for the description of everything this   type must implement.\n\n\n\n"
},

{
    "location": "interfaces.html#AbstractMesh-1",
    "page": "Code Interfaces",
    "title": "AbstractMesh",
    "category": "section",
    "text": "ODLCommonTools defines:AbstractMeshThe purpose of an AbstractMesh is to hold all the mesh related data that the  solver will need.  It also serves to establish an interface between the solver  and whatever mesh software is used.  By storing all data in the fields of the  AbstractMesh object, the details of how the mesh software stores and allows  retrieval of data are not needed by the solver.  This should make it easy to  accommodate different mesh software without making any changes to the solver.The static parameter Tmsh is used to enable differentiation with respect to the mesh variable in the future."
},

{
    "location": "interfaces.html#Required-Fields-2",
    "page": "Code Interfaces",
    "title": "Required Fields",
    "category": "section",
    "text": "  # counts\n  numVert::Integer\n  numEl::Integer\n  numNodes::Integer\n  numDof::Integer\n  numDofPerNode::Integer\n  numNodesPerElement::Integer\n  order::Integer\n  numNodesPerFace::Int\n\n  # parallel counts\n  npeers::Int\n  numGlobalEl::Int\n  numSharedEl::Int\n  peer_face_counts::Array{Int, 1}\n  local_element_counts::Array{Int, 1}\n  remote_element_counts::array{Int, 1}\n\n  # MPI Info\n  comm::MPI.Comm\n  myrank::Int\n  commsize::Int\n  peer_parts::Array{Int, 1}\n\n  # Discretization type\n  isDG::Bool\n  isInterpolated::bool\n\n  # mesh data\n  coords::AbstractArray{Tmsh, 3}\n  dxidx::AbstractArray{Tmsh, 4}\n  jac::AbstractArray{Tmsh, 2}\n\n  # interpolated data\n  coords_bndry::Array{Tmsh, 3}\n  dxidx_bndry::Array{Tmsh, 4}\n  jac_bndry::Array{T1, 2}\n  dxidx_face::Array{Tmsh, 4}\n  jac_face::Array{Tmsh, 2}\n\n  # parallel data\n  coords_sharedface::Array{Array{Tmsh, 3}, 1}\n  dxidx_sharedface::Array{Array{Tmsh, 4}, 1}\n  jac_sharedface::Array{Array{Tmsh, 2}, 1}  \n\n  # boundary condition data\n  numBC::Integer\n  numBoundaryEdges::Integer\n  bndryfaces::AbstractArray{Boundary, 1}\n  bndry_offsets::AbstractArray{Integer, 1}\n  bndry_funcs::AbstractArray{BCType, 1}\n\n  # interior edge data\n  numInterfaces::Integer\n  interfaces::AbstractArray{Interface, 1}\n\n  # degree of freedom number data\n  dofs::AbstractArray{Integer, 2}\n  dof_offset::Int\n  sparsity_bnds::AbstractArray{Integer, 2}\n  sparsity_nodebnds::AbstractArray{Integer, 2}\n\n  # mesh coloring data\n  numColors::Integer\n  maxColors::Integer\n  color_masks::AbstractArray{ AbstractArray{Number, 1}, 1}\n  shared_element_colormasks::Array{Array{BitArray{1}, 1}, 1}\n  pertNeighborEls::AbstractArray{Integer, 2}\n\n  # parallel bookkeeping\n  bndries_local::Array{Array{Boundary, 1}, 1}\n  bndries_remote::Array{Array{Boundary, 1}, 1}\n  shared_interfaces::Array{Array{Interface, 1}, 1}\n  shared_element_offsets::Array{Int, 1}\n  local_element_lists::Array{Array{Int, 1}, 1}\n"
},

{
    "location": "interfaces.html#Counts-1",
    "page": "Code Interfaces",
    "title": "Counts",
    "category": "section",
    "text": "numVert:  number of vertices in the meshnumEl:  number of elements in the meshnumNodes: number of nodes in the meshnumDof:  number of degrees of freedom in the mesh (= numNodes * numDofPerNode)numDofPerNode:  number of degrees of freedom on each node.numNodesPerElement:  number of nodes on each element.order:  order of the discretization (ie. first order, second order...), where           an order p discretization should have a convergence rate of p+1.numNodesPerFace: number of nodes on an edge in 2D or face in 3D.  For interpolated                    meshes it is the number of interpolation points on the face"
},

{
    "location": "interfaces.html#Parallel-Counts-1",
    "page": "Code Interfaces",
    "title": "Parallel Counts",
    "category": "section",
    "text": "npeers: number of processes that have elements that share a face with the current          processnumGlobalEl: number of locally owned elements + number of non-local elements that                share a face with a locally owned elementnumSharedEl: number of non-local elements that share a face with a locally owned                elementlocal_element_counts: array of length npeers, number of local elements                         that share a face with each remote processremote_element_counts: array of length npeers, number of remote elements                          that share a face with with current process"
},

{
    "location": "interfaces.html#MPI-Info-1",
    "page": "Code Interfaces",
    "title": "MPI Info",
    "category": "section",
    "text": "comm: the MPI Communicator the mesh is defined on myrank: rank of current MPI process (0-based) commsize: number of MPI processes on this communicator peer_parts: array of MPI proccess ranks for each process that has elements that               share a face with the current process"
},

{
    "location": "interfaces.html#Discretization-Type-1",
    "page": "Code Interfaces",
    "title": "Discretization Type",
    "category": "section",
    "text": "isDG: true if mesh is a DG type discretization isInterpolated: true if mesh requires data to be interpolated to the faces                   of an element"
},

{
    "location": "interfaces.html#Mesh-Data-1",
    "page": "Code Interfaces",
    "title": "Mesh Data",
    "category": "section",
    "text": "coords: n x numNodesPerElement x numEl array, where n is the dimensionality of            the equation being solved (2D or 3D typically).              coords[:, nodenum, elnum] = [x, y, z] coordinates of node nodenum            of element elnum.dxidx:  n x n x numNodesPerElement x numEl, where n is defined above. It stores the mapping jacobian scaled by ( 1/det(jac) dxi/dx ) where xi are the parametric coordinates, x are the physical (x,y,z) coordinates, and jac is the determinant of the mapping jacobian dxi/ dx.jac  : numNodesPerElement x numEl array, holding the determinant of the          mapping jacobian dxi/dx at each node of each element."
},

{
    "location": "interfaces.html#Interpolated-Data-1",
    "page": "Code Interfaces",
    "title": "Interpolated Data",
    "category": "section",
    "text": "This data is used for interpolated mesh only.coords_bndry: coordinates of nodes on the boundary of the mesh,                 2 x numFaceNodes x numBoundaryEdgesdxidx_bndry: 2 x 2 x numFaceNodes x numBoundary edges array of `dxidx                interpolated to the boundary of the meshjac_bndry: numFaceNodes x numBoundaryEdges array of jac interpolated               to the boundarydxidx_face: 2 x 2 x numFaceNodes x numInterfaces array of dxidx               interpolated to the face shared between two elementsjac_face: numNodesPerFace x numInterfaces array of jac interpolated              to the face shared between two element"
},

{
    "location": "interfaces.html#Parallel-Data-1",
    "page": "Code Interfaces",
    "title": "Parallel Data",
    "category": "section",
    "text": "This data is required for parallelizing interpolated DG meshescoords_sharedface: array of arrays, one array for each peer process,                      containing the coordinates of the nodes on the faces                      shared between a local element on a non-local element.                      Each array is 2 x numFaceNodes x number of faces shared                      with this process. dxidx_sharedface: similar to coords_sharedface, dxidx interpolated to                     faces between elements in different processes, each                     array is 2 x 2 x numFaceNodes x number of faces shared                     with this process. jac_sharedface: similar to coords_sharedface, jac interpolated to faces                   between a local element and a non-local element. Each array                   is numFaceNodes x number of faces shared with this process."
},

{
    "location": "interfaces.html#Boundary-Condition-Data-1",
    "page": "Code Interfaces",
    "title": "Boundary Condition Data",
    "category": "section",
    "text": "The mesh object stores data related to applying boundary conditions. Boundary conditions are imposed weakly, so there is no need to remove degrees of freedom from the mesh when Dirichlet boundary conditions are applied. In order to accommodate any combination of boundary conditions, an array of functors are stored as part of the mesh object, along with lists of which mesh edges (or faces in 3D) should have which boundary condition applied to themnumBC: number of different types of boundary conditions used.numBoundaryEdges: number of mesh edges that have boundary conditions applied                     to them.bndryfaces:  array of Boundary objects (which contain the element number and                the local index of the edge), of length numBoundaryEdges.bndry_offsets:  array of length numBC+1, where bndry_offsets[i] is the index                   in bndryfaces where the edges that have boundary condition                   i applied to them start.                   The final entry in bndry_offsets should be numBoundaryEdges + 1.                   Thus bndryfaces[ bndry_offsets[i]:(bndry_offsets[i+1] - 1) ]                   contains all the boundary edges that have boundary condition                   i applied to them.bndry_funcs:  array of boundary functors, length numBC.  All boundary                 functors are subtypes of BCType.  Because BCType is an                 abstract type, the elements of this array should not be used                 directly, but passed as an argument to another function, to                  avoid type instability."
},

{
    "location": "interfaces.html#Interior-Edge-Data-1",
    "page": "Code Interfaces",
    "title": "Interior Edge Data",
    "category": "section",
    "text": "Data about interior mesh edges (or faces in 3D) is stored to enable use of edge stabilization or Discontinuous Galerkin type discretizations. Only data for edges (faces) that are shared by two elements are stored (ie. boundary edges are not considered).numInterfaces:  number of interior edgesinterfaces:  array of Interface types (which contain the element numbers for                the two elements sharing the edge, and the local index of the                edge from the perspective of the two elements, and an indication                of the relative edge orientation).                The two element are referred to as elementL and elementR,                but the choice of which element is elementL and which is                elementR is arbitrary.                The length of the array is numInterfaces.               Unlike bndryfaces, the entries in the array do not have to be in               any particular order."
},

{
    "location": "interfaces.html#Degree-of-Freedom-Numbering-Data-1",
    "page": "Code Interfaces",
    "title": "Degree of Freedom Numbering Data",
    "category": "section",
    "text": "dofs:  numDofPerNode x numNodesPerElement x numEl array. Holds the local degree of freedom number of each degree of freedom. Although the exact method used to assign dof numbers is not critical, all degrees of freedom on a node must be numbered sequentially.dof_offset: offset added to the local dof number to make it a global dof               number.sparsity_bnds:  2 x numDof array. sparsity_bnds[:, i] holds the maximum, minimum degree of freedom numbers associated with degree of freedom i. In this context, degrees of freedom i and j are associated if entry (i,j) of the jacobian is non-zero. In actuality, sparsity_bnds need only define upper and lower bounds for degree of freedom associations (ie. they need not be tight bounds). This array is used to to define the sparsity pattern of the jacobian matrix.sparsity_nodebnds:  2 x numNodes array. sparsity_bnds[:, i] holds the maximum, minimum node associated with node i, similar the information stored in sparsity_bnds for degrees of freedom."
},

{
    "location": "interfaces.html#Mesh-Coloring-Data-1",
    "page": "Code Interfaces",
    "title": "Mesh Coloring Data",
    "category": "section",
    "text": "The NonlinearSolvers module uses algorithmic differentiation to compute the Jacobian. Doing so efficiently requires perturbing multiple degrees of freedom simultaneously, but perturbing associated degrees of freedom at the same time leads to incorrect results. Mesh coloring assigns each element of the mesh to a group (color) such that every degree of freedom on each element is not associated with any other degree  of freedom on any other element of the same color. An important aspect of satisfying this condition is the use of the  element-based arrays (all arrays that store data for a quantity over the entire  mesh are ncomp x numNodesPerElement x numEl). In such an array, any node that is part of 2 or more elements has one entry for each element. When performing algorithmic differentiation, this enables perturbing a degree of  freedom on one element without perturbing it on the other elements that share  the degree of freedom.For example, consider a node that is shared by two elements. Let us say it is node 2 of element 1 and node 3 of element 2. This means AbstractSolutionData.q[:, 2, 1] stores the solution variables for this node on the first element, and AbstractSolutionData.q[:, 3, 2] stores the solution variables for the second element. Because these are different entries in the array AbstractSolutionData.q, they can be perturbed independently. Because AbstractSolutionData.res has the same format, the perturbations to AbstractSolutionData.q[:, 2, 1] are mapped to AbstractSolutionData.res[:, 2, 1] for a typical continuous Galerkin type discretization. This is a direct result of having an element-based discretization.There are some discretizations, however, that are not strictly element-based. Edge stabilization, for example, causes all the degrees of freedom of one element to be associated with any elements it shares an edge with. To deal with this, we use the idea of a distance-n coloring. A distance-n coloring is a coloring where there are n elements in between two element of the same color. For element-based discretizations with element-based arrays, every element in the mesh can be the same color.  This is a distance-0 coloring. For an edge stabilization discretization, a distance-1 coloring is required, where every element is a different color than any neighbors it shares and edge with. (As a side node, the algorithms that perform a distance-1 coloring are rather complicated, so in practice we use a distance-2 coloring instead).In order to do algorithmic differentiation, the AbstractMesh object must store the information that determines which elements are perturbed for which colors, and, for the edge stabilization case, how to relate a perturbation in the output  AbstractSolutionData.res to the degree of freedom in AbstractSolutionData.q in O(1) time. Each degree of freedom on an element is perturbed independently of the other degrees of freedom on the element, so the total number of residual evaluations is the number of colors times the number of degrees of freedom on an element.The fields required are:numColors:  The number of colors in the mesh. maxColors: the maximum number of colors on any processcolor_masks:  array of length numColors.  Each entry in the array is itself                 an array of length numEl.  Each entry of the inner array is                 either a 1 or a 0, indicating if the current element is                 perturbed or not for the current color. For example, in color_mask_i = color_masks[i]; mask_elj = color_mask_i[j], the variable mask_elj is either a 1 or a zero, determining whether or not element j is perturbed as part of color i.shared_element_colormasks: array of BitArrays controlling when to perturb                              non-local elements.  There are npeers arrays,                              each of length number of non-local elements shared                              with this processpertNeighborEls:  numEl x numColors array.  neighbor_nums[i,j] is the element number of of the element whose perturbation is affected element i when color j is being perturbed, or zero if element i is not affected by any  perturbation.  "
},

{
    "location": "interfaces.html#Parallel-Bookkeeping-1",
    "page": "Code Interfaces",
    "title": "Parallel Bookkeeping",
    "category": "section",
    "text": "bndries_local: array of arrays of Boundarys describing faces shared                  with non-local elements from the local side (ie. the                  local element number and face).  The number of arrays is                  npeers.bndries_remote: similar to bndries_local, except describing the faces                   from the non-local side.  Note that the element numbers are                   from the remote processshared_interfaces: array of arrays of Interfaces describing the faces                      shared between local and non-local elements.  The local                      element is always elementL.  The remote element is                      assigned a local number greater than numEl.shared_element_offsets: array of length npeers+1 that contains the first                           element number assigned to elements on the shared                           interface.  The last entry is equal to the highest                           element number on the last peer process + 1.  The                           elements numbers assigned to a given peer must form                           a contiguous range.local_element_lists: array of arrays containing the element numbers of the                        elements that share a face with each peer process"
},

{
    "location": "interfaces.html#Other-Functions-1",
    "page": "Code Interfaces",
    "title": "Other Functions",
    "category": "section",
    "text": "The mesh module must also defines and exports the functionssaveSolutionToMesh(mesh::MeshImplementationType, vec::AbstractVector)\nwriteVisFiles(mesh::MeshImplementationType, fname::ASCIIString)where the first function takes a vector of length numDof and saves it to the mesh, and the second writes Paraview files for the mesh, including the solution field."
},

{
    "location": "interfaces.html#ODLCommonTools.AbstractCGMesh",
    "page": "Code Interfaces",
    "title": "ODLCommonTools.AbstractCGMesh",
    "category": "Type",
    "text": "ODLCommonTools.AbstractCGMesh\n\nThe abstrac type is the supertype of all continuous Galerkin meshes\n\n\n\n"
},

{
    "location": "interfaces.html#ODLCommonTools.AbstractDGMesh",
    "page": "Code Interfaces",
    "title": "ODLCommonTools.AbstractDGMesh",
    "category": "Type",
    "text": "ODLCommonTools.AbstractDGGMesh\n\nThe abstrac type is the supertype of all discontinuous Galerkin meshes\n\n\n\n"
},

{
    "location": "interfaces.html#Subtypes-1",
    "page": "Code Interfaces",
    "title": "Subtypes",
    "category": "section",
    "text": "AbstractCGMesh\nAbstractDGMesh"
},

{
    "location": "interfaces.html#Physics-Module-1",
    "page": "Code Interfaces",
    "title": "Physics Module",
    "category": "section",
    "text": "For every new physics to be solved, a new module should be created. The purpose of this module is to evaluate the equation:M dq/dt = f(q)where M is the mass matrix. For steady problems, dq/dt = 0 and the module evaluates the residual. For unsteady problems, the form M dq/dt = f(q) is suitable for explicit time marching.Every physics module must define a boundary condition called defaultBC.  If the user does not specify a boundary condition on any geometric edges, the mesh constructor will add a new boundary condition and assign all mesh edges classified on the unspecified geometric edges are assigned to it.  This boundary condition can be a no-op if that is correct for the physics, or it can be the face integral that should be done for every element (recall that boundary faces do not appear in mesh.interfaces)."
},

{
    "location": "interfaces.html#Interface-to-NonlinearSolvers-1",
    "page": "Code Interfaces",
    "title": "Interface to NonlinearSolvers",
    "category": "section",
    "text": "The evalResidual function and the fields eqn.q and eqn.res are the interface between the NonlinearSolvers and the physics modules.  The Nonlinear solvers populate eqn.q, and use evalResidual to populate eqn.res, from which the next value if eqn.q is calculated.  The algorthmic differentiation mechanism described above is uses several residual evaluations to compute the Jacobian if needed for a given method.  Some methods, such as RK4, are better expressed in terms of eqn.q_vec and eqn.res_vec rather than eqn.q and eqn.res.  eqn.assembleSolution and eqn.disassmbleSolution exist to transform back and forth between the vector and 3D array forms.  In order to compute the Jacobian efficiently, it is necessary to work with the 3D arrays. For this reason, evalResidual must work only with eqn.q and eqn.res and let the caller decide whether or not to transform into the vector form.Newton's method supports both finite differences and complex step for calculating the Jacobian, and the static parameters need to be set accordingly. If finite differences are used, then Tsol=Tres=Tmsh=Float64 is required.  If complex step is used, then Tsol=Tres=Complex128 and Tmsh = Float64 is needed."
},

{
    "location": "interfaces.html#Interface-to-users-1",
    "page": "Code Interfaces",
    "title": "Interface to users",
    "category": "section",
    "text": "The interface to users is described in src/Readme.md"
},

{
    "location": "interfaces.html#Interface-to-Summation-by-Parts-1",
    "page": "Code Interfaces",
    "title": "Interface to Summation-by-Parts",
    "category": "section",
    "text": "The physics modules use the interface provided by the Summation-by-Parts package to approximate derivatives numerically.  The reason for passing around the sbp object is that the SBP interface requires it to be passed in along with the data arrays it operates on."
},

{
    "location": "interfaces.html#Functional-Programming-1",
    "page": "Code Interfaces",
    "title": "Functional Programming",
    "category": "section",
    "text": "An important aspect of the use of the mesh, sbp, eqn, and opts to define interfaces is that the physics modules and nonlinear solvers are written in a purely functional programming style, which is to say that the behavior of every function is determined entirely by the arguments to the function, and the only effects of the function are to modify an argument or return a value.This property is important for writing generic, reusable code. For example, when using iterative solvers, a preconditioner is usually required, and constructing the preconditioner requires doing residual evaluations. In some cases, the preconditioner will use a different mesh or a different mesh coloring. Because the physics modules are written in a functional style, they can be used to evaluate the residual on a different mesh simply by passing the residual evaluation function a different mesh object.A critical aspect of function programming is that there is no global state. The state of the solver is determined entirely by the state of the objects that are passed around.##Variable Naming Conventions In an attempt to make code more uniform and readable, certain variable names are reserved for certain uses.mesh\n:  object that implements \nAbstractMesh\npmesh\n:  mesh used for preconditioning\nsbp\n:  Summation-by-Parts object\neqn\n:  object that implements \nAbstractSolutionData\nopts\n: options dictionary\nparams\n: parameter object (used for values that might be in \nopts\n but need to be accessed quickly)\nx\n: the first real coordinate direction\ny\n: the second real coordinate direction\nxi\n: the first parametric coordinate direction\neta\n: the second parametric coordinate direction\nh\n:  mesh spacing\nhx\n: mesh spacing in x direction\nhy\n: mesh spacing in y direction\np\n: pressure at a node\na\n; speed of sound at a node\ns\n: entropy at a node\ngamma\n: specific heat ratio\ngamma_1\n: \ngamma\n - 1\nR\n: specific gas constant in ideal gas law (units J/(Kg * K) in SI)\ndelta_t\n: time step\nt\n: current time\nnrm\n: a normal vector of some kind\nA0\n: the coefficient matrix of the time term at a node\nAxi\n: the flux jacobian at a node in the \nxi\n direction\nAeta\n: the flux jacobian at a node in the \neta\n direction\nAx\n: the flux jacobian in the \nx\n direction at a node\nAy\n: the flux jacobian in the \ny\n direction at a node"
},

{
    "location": "parallel.html#",
    "page": "Code Parallelization",
    "title": "Code Parallelization",
    "category": "page",
    "text": ""
},

{
    "location": "parallel.html#Parallel-Overview-1",
    "page": "Code Parallelization",
    "title": "Parallel Overview",
    "category": "section",
    "text": "This document describes how PDEs are solved in parallel.In general, the mesh is partitioned and each part is assigned to a different MPI process. Each element is owned by exactly one process.  During  initialization, the mesh constructor on each process figures out which other processes have elements that share a face (edge in 2D) with local elements. It counts how many faces and elements are shared (a single element could have multiple faces on the parallel boundary), and assigns local number to both the elements and the degrees of freedom on the elements.  The non-local elements  are given numbers greater than numEl, the number of local elements.   The degrees of freedom are re-numbered such that newly assigned dof number plus the dof_offset for the current process equals the global dof number, which is defined by the local dof number assigned by the process that owns the element.  As a  result, dof numbers for elements that live on processes with lower ranks will be negative.As part of counting the number of shared faces and elements, 3 arrays are formed: bndries_local, bndries_remote, and shared_interfaces which  describe the shared faces from the local side, the remote side, and a  unified view of the interface, respectively.  This allows treating the shared faces like either a special kind of boundary condition or a proper  interface, similar to the interior interfaces.There are 2 modes of parallel operation, one for explicit time marching and the other for Newton's method."
},

{
    "location": "parallel.html#Explicit-Time-Marching-1",
    "page": "Code Parallelization",
    "title": "Explicit Time Marching",
    "category": "section",
    "text": "In this mode, each process each process sends the solution values at the  shared faces to the other processes.  Each process then evaluates the residual using the received values and updates the solution.The function exchangeFaceData is designed to perform the sending and  receiving of data.  Non-blocking communications are used, and the function does not wait for the communication to finish before returning.  The  MPI_Requests for the sends and receives are stored in the appropriate fields of the mesh.  It is the responsibility of each physics module call  exchangeFaceData and to wait for the communication to finish before using the data.  Because the receives could be completed in any order, it is  recommended to use MPI_Waitany to wait for the first receive to complete,  do as many computations as possible on the data, and then call MPI_Waitany again for the next receive."
},

{
    "location": "parallel.html#Newton's-Method-1",
    "page": "Code Parallelization",
    "title": "Newton's Method",
    "category": "section",
    "text": "For Newton's method, each process sends the solution values for all the  elements on the shared interface at the beginning of a Jacobian calculation.  Each process is then responsible for perturbing the solutions values of both  the local and non-local elements.  The benefit of this is that parallel  communication is required once per Jacobian calculation, rather than once  per residual evaluation as with the explicit time marching mode.The function exchangeElementData copies the data from the shared elements into the send buffer and sends it, and also posts the corresponding receives. It does not wait for the communications to finish before returning.   The function is called by Newton's method after a new solution is calculated, so the physics module does not have to do it, but the physics module does  have to wait for the receives to finish before using the data.  This is  necessary to allow overlap of the communication with computation.  As  with the explicit time marching mode, use of MPI_Waitany is recommended. "
},

{
    "location": "examples/isentropic.html#",
    "page": "Isentropic Vortex",
    "title": "Isentropic Vortex",
    "category": "page",
    "text": ""
},

{
    "location": "examples/isentropic.html#Example-1:-Steady-Isentropic-Vortex-1",
    "page": "Isentropic Vortex",
    "title": "Example 1: Steady Isentropic Vortex",
    "category": "section",
    "text": "Save the following as input_vals.jl.arg_dict = Dict{Any, Any}(\n# specify the physics and SBP operator\n\"physics\" => \"Euler\",  # specify physics to run\n\"dimensions\" => 2,  # this is a two dimensional mesh\n\n# specify temporal and spatial discretization\n\"run_type\" => 5,  # steady Newton (5 historically was complex-stepped Newton, as opposed to 4 being FD)\n\"jac_method\" => 2,  # complex-step Newton Jacobian calculation\n\"jac_type\" => 1,  # store the Jacobian as a Julia sparse matrix\n\"t_max\" => 10.0,  # make time\n\"operator_type\" => \"SBPOmega\",  # specify SBP operator\n\"order\" => 1,  # p = 1 operator\n\"use_DG\" => true,  # use discontinuous galerkin solver\n\"Flux_name\" => \"RoeFlux\",  # numerical flux function used in face integrals\n\"CFL\" => 0.10,  # CFL number\n\"itermax\" => 20,\n\n# specify the problem itself\n\"IC_name\" => \"ICIsentropicVortex\",  # initial condtiion\n\"numBC\" => 1,  # number of boundary conditions\n\"BC1\" => [0, 1, 2, 3],  # geometric edges to apply the BC to\n\"BC1_name\" => \"isentropicVortexBC\",   # name of boundary condition\n\n# specify mesh\n\"smb_name\" => \"SRCMESHES/squarevortex_small.smb\",\n\n# misc options\n\"write_vis\" => true,  # write paraview files\n\"do_postproc\" => true,  # calculate error at end of run\n\"exact_soln_func\" => \"ICIsentropicVortex\",  # use this function for the exact soltuion (to calculate error)\n)Run the case with julia ~/.julia/v0.4/PDESolver/src/solver/euler/startup.jl input_vals.jl."
},

{
    "location": "examples/unsteady.html#",
    "page": "Unsteady Vortex",
    "title": "Unsteady Vortex",
    "category": "page",
    "text": ""
},

{
    "location": "examples/unsteady.html#Example-2:-Unsteady-Advecting-Isentropic-Vortex-1",
    "page": "Unsteady Vortex",
    "title": "Example 2: Unsteady Advecting Isentropic Vortex",
    "category": "section",
    "text": ""
},

{
    "location": "pdesolver.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "pdesolver.html#Documentation-of-the-PDESolver-Module-1",
    "page": "Introduction",
    "title": "Documentation of the PDESolver Module",
    "category": "section",
    "text": "PDESolver is the primary module of the solver. Its purpose is to load all the components of PDESolver and provide an interface for running the code. It also provides an interface for registering new initial conditions and boundary conditions with the physics modules.  Pages = [ \"pdesolver_user.md\"\n            \"pdesolver_physics.md\"\n          ]\n  Depth = 1"
},

{
    "location": "pdesolver_user.html#",
    "page": "PDESolver User Interface",
    "title": "PDESolver User Interface",
    "category": "page",
    "text": ""
},

{
    "location": "pdesolver_user.html#PDESolver-User-Interface-Interface-1",
    "page": "PDESolver User Interface",
    "title": "PDESolver User Interface Interface",
    "category": "section",
    "text": "CurrentModule = PDESolver"
},

{
    "location": "pdesolver_user.html#PDESolver.run_solver",
    "page": "PDESolver User Interface",
    "title": "PDESolver.run_solver",
    "category": "Function",
    "text": "This function provides a way to invoke any physics solver based on the   specification of the physics in the input file.   This requires loading the input file twice, once to figure out the physics,   and a second time when the physics-specific startup function is called\n\nThe physics module must have already been registered using register_physics\n\nInputs:\n\ninput_file: an AbstractString specifying the path to the input file\n\nOutputs:\n\nmesh: the AbstractMesh object used during the solve\nsbp: the SBP operator used by the solver\neqn: the AbstractSolutionData object during the solve.  At exit,\n     eqn.q_vec should have the final solution in it\nopts: the options dictionary\n\n\n\n"
},

{
    "location": "pdesolver_user.html#Invoking-the-Code-1",
    "page": "PDESolver User Interface",
    "title": "Invoking the Code",
    "category": "section",
    "text": "This page describes the API for running the code, as well as utility scripts that provide simple operation of the code.The entry point for running the code is the function run_solverrun_solverEach physics module is required to register itself when using PhysicsModule is run. When using this API to invoke the code, users should be sure to have loaded the required physics module before calling run_solver.In order to facilitate running the code, the script src/startup.jl loads all the physics modules and passes the first command line argument as the input file name"
},

{
    "location": "pdesolver_user.html#PDESolver.registerIC",
    "page": "PDESolver User Interface",
    "title": "PDESolver.registerIC",
    "category": "Function",
    "text": "This function registers a new initial condition function with the specified   physics module.  The function must have the signature:\n\n```    ICfunc(mesh::AbstractMesh, sbp::AbstractSBP, eqn:AbstractSolutionData{Tsol},                 opts:Dict, q_vec::AbstractVector{Tsol})   ```\n\nwhere q_vec is a vector of length numDof that will be populated with the   initial condition.  If the function is used to provide the exact solution   for an unsteady problem (for calculating the error via the exact_soln_func   key in the options dictionary), then it should use eqn.params.t as the   current time value.\n\nThis function does not attempt to verify that the functor has the correct   signature, so the user should take care to ensure it is correct for the   physics module.\n\nInputs:\n\nmod: physics module to register the function with\nfname: name associated with the function, used as the value for the\n       for any key in the options dictionary that specifies an initial\n       condition function\nfunc: the function itself\n\nOutputs:\n\nnone\n\n\n\n"
},

{
    "location": "pdesolver_user.html#PDESolver.registerBC",
    "page": "PDESolver User Interface",
    "title": "PDESolver.registerBC",
    "category": "Function",
    "text": "This function registers a new boundary condition functor with the specified   physics module.  The exact signature of the functor varies by physics module.   The purpose of the functor is to calculate the flux across a boundary at   a particular node. It usually takes in the solution values at that node, as   well as the coordinates and mapping jacobian.  The params type is typically   used to dispatch to different methods for 2D or 3D, if needed.  Generally,   this functor will calculate the boundary state and then call a numerical   flux function of some kind to compute the flux.\n\nThis function does not attempt to verify that the functor has the correct   signature, so the user should take care to ensure it is correct for the   physics module.\n\nInputs:\n\nmod: module to register the the functor with\nfname: the name associated with this function, used as the value for any\n       key in the options dictionary that specifies a boundary condition,\n       for example `BC1_name`\nfunc: the functor itself\n\nOutputs:\n\nnone\n\n\n\n"
},

{
    "location": "pdesolver_user.html#Registering-ICs-and-BCs-1",
    "page": "PDESolver User Interface",
    "title": "Registering ICs and BCs",
    "category": "section",
    "text": "These functions allow users to register new initial conditions and boundary conditions without modifying the source coderegisterIC\nregisterBC"
},

{
    "location": "pdesolver_user.html#PDESolver.printICNames",
    "page": "PDESolver User Interface",
    "title": "PDESolver.printICNames",
    "category": "Function",
    "text": "Print all initial condition names for a given physics module\n\nInputs:\n\nmod: the module to print the ICs for\nf: the IO to print to, defaults to STDOUT\n\nOutputs:     none\n\n\n\n"
},

{
    "location": "pdesolver_user.html#PDESolver.printBCNames",
    "page": "PDESolver User Interface",
    "title": "PDESolver.printBCNames",
    "category": "Function",
    "text": "Like printICNames, but for boundary conditions, see that function for details\n\n\n\n"
},

{
    "location": "pdesolver_user.html#PDESolver.printPhysicsModules",
    "page": "PDESolver User Interface",
    "title": "PDESolver.printPhysicsModules",
    "category": "Function",
    "text": "Prints the name of all currently registered physics modules and the module   name itself\n\n\n\n"
},

{
    "location": "pdesolver_user.html#Interactive-Functions-1",
    "page": "PDESolver User Interface",
    "title": "Interactive Functions",
    "category": "section",
    "text": "These functions are useful for inspecting the state of the front-end.printICNames\nprintBCNames\nprintPhysicsModules"
},

{
    "location": "pdesolver_physics.html#",
    "page": "PDESolver PhysicsInterface",
    "title": "PDESolver PhysicsInterface",
    "category": "page",
    "text": ""
},

{
    "location": "pdesolver_physics.html#Documentation-of-the-PDESolver-Physics-Module-Interface-1",
    "page": "PDESolver PhysicsInterface",
    "title": "Documentation of the PDESolver Physics Module Interface",
    "category": "section",
    "text": "CurrentModule = PDESolverThe PDESolver Module provides some structure for physics modules to plug into. See Interfaces in PDESolver for details. Each physics module should extend evalResidual with a new method to which takes the AbstractSolutionData defined by the physics module as an argument. This structure allows evalResidual to be visible to the NonlinearSolver module while being defined by the physics modules."
},

{
    "location": "pdesolver_physics.html#PDESolver.evalResidual",
    "page": "PDESolver PhysicsInterface",
    "title": "PDESolver.evalResidual",
    "category": "Function",
    "text": "This function evalutes dq/dt = R(q).  For steady problems it evalutes R(q)   at some state q.  The state is stored in eqn.q, and eqn.res is populated with   R(q).  Note that these are 3 dimensional arrays.  The physics modules only   interact with the 3 dimensional arrays, never the vectors eqn.q_vec and   eqn.res_vec.  Each physics module must implement this function for its   own subtype of AbstractSolutionData (ie. with a more specific type for   the eqn argument and equallty specific types for the other arguments).   This is important because evalResidual is common to all physics modules,   so a user or some other part of the code can call evalResidual(mesh, sbp   eqn, opts), and Julia's multiple dispatch will figure out the right method   to call based on the type of the eqn argument.\n\nThe evaluation of the residual R(q) should depend only on the data stored in   mesh, sbp, eqn, and opts, and any data that depend on q should be recalculated   every time the function is called.  This function is used as a building block   by other parts of the solver, particularly the NonlinearSolvers.  See   interfaces.md for details\n\nInputs:     mesh: an AbstractMesh describing the mesh on which to solve the physics     sbp: an SBP operator     eqn: a subtype of AbstractSolution data, used to store all of the data used          by the physics module     opts: the options dictionary     t: the current time value, defaults to 0.0\n\n\n\nAdvectionEquationMod.evalResidual\n\nThis function evaluates the Advection equation.\n\nInputs\n\n \nmesh\n : Abstract mesh object\n \nsbp\n  : Summation-by-parts operator\n \neqn\n  : Advection equation object\n \nopts\n : Options dictionary\n \nt\n    :\n\nEffectively updates eqn.res – not eqn.res_vec. To make them consistent, use assembleSolution on eqn.res and eqn.res_vec\n\nOutputs\n\n None\n\n\n\nEulerEquationMod.evalResidual\n\nThis function drives the evaluation of the EulerEquations.   It is agnostic to the dimension of the equation. and the types the arguments   are paramaterized on.\n\nThe function calls only high level functions, all of which take the same   four arguments.  Mid level function also take the same arguments.\n\nThe input/output variables are eqn.q and eqn.res, respectively.   eqn.q_vec and eqn.res_vec exist for reusable storage outside the residual   evaluation.  They should never be used inside the residual evaluation.\n\nThe function disassembleSolution takes q_vec and puts it into eqn.q   The function assembleSolution takes eqn.res and puts it into res_vec\n\nArguments:     * mesh  : a mesh object     * sbp   : SBP operator object     * eqn   : an EulerData object     * opts  : options dictionary\n\nThe optional time argument is used for unsteady equations\n\n\n\nSimpleODEMod.evalResidual\n\nThis function evaluates the simple ODE equation.\n\n** Inputs **\n\n \nmesh\n : Abstract mesh object\n \nsbp\n  : Summation-by-parts operator\n \neqn\n  : Simple ODE equation object\n \nopts\n : Options dictionary\n \nt\n    :\n\nOutputs\n\n None\n\n\n\n"
},

{
    "location": "pdesolver_physics.html#Evaluating-the-Physics-1",
    "page": "PDESolver PhysicsInterface",
    "title": "Evaluating the Physics",
    "category": "section",
    "text": "evalResidual"
},

{
    "location": "pdesolver_physics.html#PDESolver.register_physics",
    "page": "PDESolver PhysicsInterface",
    "title": "PDESolver.register_physics",
    "category": "Function",
    "text": "This function registered a new physics module with the global list of all   known physics modules.  Every physics module should do this as part of   module initialization.  The name, the module, and the startup function must   be unique (ie. they must not already exist in the list).  This function   throws and exception if they are not.\n\nInputs:\n\nmodname:  an ASCIIString name for this entry in the list.  It is used\n          to retrieve the module and startup function in the \n          retrieve_physics function. Typically the name is capitalized.\nmod:  the Module itself\nstartup_func: the function for running the physics.  It must have signature\n              startup_func(fname::ASCIIString), where fname is the name of\n              an input file\n\nOutputs:\n\nnone\n\n\n\n"
},

{
    "location": "pdesolver_physics.html#PDESolver.retrieve_physics",
    "page": "PDESolver PhysicsInterface",
    "title": "PDESolver.retrieve_physics",
    "category": "Function",
    "text": "Retrieves the physics module and function registered using register_physics\n\nInput:\n\nmodname: an ASCIIString containing the name of the module supplied to\n         `register_physics`\n\nOutputs:\n\nmod: the physics Module\nfunc: the function to evaluate the physics\n\n\n\n"
},

{
    "location": "pdesolver_physics.html#Registration-Functions-1",
    "page": "PDESolver PhysicsInterface",
    "title": "Registration Functions",
    "category": "section",
    "text": "These function provide the basic API for registering and retrieving physics modules with the PDESolver module.register_physics\nretrieve_physicsSee also registerIC and registerBC."
},

{
    "location": "invocation/calling.html#",
    "page": "Calling PDESolver",
    "title": "Calling PDESolver",
    "category": "page",
    "text": ""
},

{
    "location": "invocation/calling.html#Running-the-Code-1",
    "page": "Calling PDESolver",
    "title": "Running the Code",
    "category": "section",
    "text": "Before running the code, you must source the shell script PumiInterface/src/use_julialib.sh.  This enables Julia to find and use Pumi.The code takes an input file that defines all options for the solver, including the which physics to solve, which mesh to use, initial conditions, boundary conditions, and discretization.  The file src/input/input_vals.txt describes the input file format and all valid options."
},

{
    "location": "invocation/calling.html#Simple-Mode-1",
    "page": "Calling PDESolver",
    "title": "Simple Mode",
    "category": "section",
    "text": "Once an input file is prepared, the solver can be invoked withmpirun -np x julia /path/to/PDESolver/src/startup.jl \"input_file_name\"This will run the solver in parallel with x processes, solving the physics specified in the input file."
},

{
    "location": "invocation/calling.html#Advanced-Mode-1",
    "page": "Calling PDESolver",
    "title": "Advanced Mode",
    "category": "section",
    "text": "The code also provides a scripting interface.  This enables users to supply initial conditions and boundary conditions, as well as do advanced post-processing if desired.The template script is:using PDESolver  # load the interface to the code\nusing AdvectionEquationMod  # load the desired physics module\n\n# register ICs and BCs here\n\ninput_fname = ARGS[1]  # use the first command line argument as the input file\n                       # name\n\n# solve the equation defined in the input file\nmesh, sbp, eqn, opts = run_solver(input_fname)\n\n# do post-processing using mesh, sbp, eqn, opts here\nmesh, sbp, eqn, opts are the AbstractMesh object, the SummationByParts operator, the AbstractSolutionData, and the options dictionary, respectively.The AbstractMesh object contains all the data about the mesh, including coordinates and mapping jacobian information.  The sbp object is used for numerically approximating derivatives. The AbstractSolutionData object contains all the data about the solution and the residual of the equation.  The input file gets parsed into the options dictionary, which is used to by the rest of the code.The PDESolver module exports several functions for use by the scripting interface.  The first two registerIC and registerBC allow users to supply initial conditions and boundary conditions.  See the documentation of the functions for details.  The functions printICNames and printBCNames print the names of all currently registers initial conditions and boundary conditions.For example a user can, either in a Julia script or interactive (REPL) session:using PDESolver\nusing AdvectionEquationMod\n\nprintBCNames(AdvectionEquationMod)And this will print the names of all the boundary conditions currently known to the Advection Equation module."
},

{
    "location": "invocation/interactive.html#",
    "page": "Interactive Session (experimental)",
    "title": "Interactive Session (experimental)",
    "category": "page",
    "text": ""
},

{
    "location": "invocation/interactive.html#Running-PDESolver-in-an-interactive-session-using-Julia's-REPL-1",
    "page": "Interactive Session (experimental)",
    "title": "Running PDESolver in an interactive session using Julia's REPL",
    "category": "section",
    "text": "Julia's REPL is a powerful tool for debugging.  Being able to run commands interactively can reveal important information   that would otherwise be available only through cumbersome print statements and  error logs.This page describes the first steps involved for running PDESolver interactively. Calling PDESolver from Julia's REPL in the fashion described here is an experimental    method, and is not tested consistently. Additionally, not all parts of PDESolver are usable after the steps below. Instead, this document is intended to provide a springing-off point for the developer    who wishes to run PDESolver interactively, and can adapt the commands below to their needs."
},

{
    "location": "invocation/interactive.html#Script-1",
    "page": "Interactive Session (experimental)",
    "title": "Script",
    "category": "section",
    "text": "using PDESolver\nusing AdvectionEquationMod\nusing ArrayViews\nusing ODLCommonTools\nusing SummationByParts\nusing PdePumiInterface\nusing NonlinearSolvers\nusing ForwardDiff\nusing Utils\nimport ODLCommonTools.sview\nusing MPI\nusing Input\n\n# after here, hopefully this is the absolute minimum commands needed\nMPI.Init()\n\n# example for input file; substitute your own\nopts = read_input(\"sine_adv_input_CN_adjoint.jl\")\n\nTdim = opts[\"dimensions\"]\n# note: funny character in opts, last entry when loaded in REPL\nopts[\"Tdim\"] = Tdim\ndofpernode = 1\n\nsbp, mesh, pmesh, Tsol, Tres, Tmsh, mesh_time = createMeshAndOperator(opts, dofpernode)\n\nvar_type = opts[\"variable_type\"]\n\n# initialize the eqn object with the desired physics\n# eqn = EulerData_{Tsol, Tres, Tdim, Tmsh, var_type}(mesh, sbp, opts)\neqn = AdvectionData_{Tsol, Tres, Tdim, Tmsh}(mesh, sbp, opts)\n\ninit(mesh, sbp, eqn, opts)\n\n# sample commands for populating some fields of eqn\nrand!(eqn.q)\nrand!(eqn.res)\n# q_vec and res_vec are populated after this too"
},

{
    "location": "invocation/interactive.html#REPL-debugging-commands-1",
    "page": "Interactive Session (experimental)",
    "title": "REPL debugging commands",
    "category": "section",
    "text": "Here are some helpful Julia commands for debugging/investigation within the REPL at this stage. Refer to the Julia documentation for details, or press '?' within the REPL    and type the command you wish to see help for.size()\nlength()\ntypeof()\nfieldnames()\nisimmutable()\nwhos()"
},

{
    "location": "solver/Readme.html#",
    "page": "Overview of Physics Modules",
    "title": "Overview of Physics Modules",
    "category": "page",
    "text": ""
},

{
    "location": "solver/Readme.html#Overview-of-Physics-Modules-1",
    "page": "Overview of Physics Modules",
    "title": "Overview of Physics Modules",
    "category": "section",
    "text": ""
},

{
    "location": "solver/Readme.html#AbstractSolutionData-and-Physics-Module-Implementation-1",
    "page": "Overview of Physics Modules",
    "title": "AbstractSolutionData and Physics Module Implementation",
    "category": "section",
    "text": "This document describes some best practices for implementing a physics module. These practices are not required, but have proven to be useful for producing organized, readable, and reusable code."
},

{
    "location": "solver/Readme.html#Levels-of-Functions-1",
    "page": "Overview of Physics Modules",
    "title": "Levels of Functions",
    "category": "section",
    "text": "It is useful to divide functions into 3 catagories, high, mid, and low level functions.  The purpose of high level functions is to decide which method of performing an operation should be used and call other functions to do it. For example, the Euler physics modules has evalVolumeIntegrals and evalBoundaryIntegrals as high level functions.  There are several different ways of calculating both the volume and boundary integrals.  The options dictionary is used to decide what mid level function to call.  Each mid level function implements a different way of doing the calculation.The purpose of mid level functions is to loop over the mesh and call a low level function for each node.  For example, the function getEulerFlux loops over the nodes of the mesh and calls a function to calculate the Euler flux at each node.  Mid level function names usually start with get to indicate that their purpose is to calculate some quantity but don't do the calculation themselves.Low level functions calculate a quantity at a node.  For example, calcEulerFlux calculates the Euler flux at a single node.  Low level function names usually start with calc to indicate that they perform a specific calculation. Often, different discretizations use the same structure of loops, but do a slightly different calculation at each node.  Low level functions are called inside the innermost loop of the code, so it would be too expensive to have if statements to select which low level function to call, so various tricks involving Julia's multiple dispatch system are used to get the compiler to decide which low level function to call.  These will be described later in this document.It is often useful to dispatch to low level functions based on Tdim and var_type.  For this reason the Euler equation implementation of AbstractParamType istype ParamType{Tdim, var_type, Tsol, Tres, Tmsh} <: AbstractParamType{Tdim}The other static parameters are necessary because ParamType has fields of those datatypes."
},

{
    "location": "solver/Readme.html#AbstractSolutionData-implementation-1",
    "page": "Overview of Physics Modules",
    "title": "AbstractSolutionData implementation",
    "category": "section",
    "text": "Each physics module should define and export a subtype of AbstractSolutionData{Tsol, Tres}. The implementation of AbstractSolutionData{Tsol, Tres} must inherit the Tsol and Tres static parameters, and may have additional static parameters as well. It may also be helpful to define additional abstract types within the physics module to provide different levels of abstractions. For example, the Euler physics module defines:abstract AbstractEulerData{Tsol, Tres} <: AbstractSolutionData{Tsol, Tres}\nabstract EulerData {Tsol, Tdim, Tres, var_type} <: AbstractEulerData{Tsol, Tres}\ntype EulerData_{Tsol, Tres, Tdim, Tmsh, var_type} <: EulerData{Tsol, Tdim, Tres, var_type}The first line is effectively just a name change and may not be necessary. The second line adds the static parameters Tdim, and var_type while inheriting the Tsol and Tres types from AbstractEulerData. Tdim is the dimensionality of the equation, Tres is the datatype of the residual variables, and var_type is a symbol indicating whether the equation is being solved with conservative or entropy variables. The third line defines a concrete type that implements all the features required of an AbstractSolutionData, and adds a static parameter Tmsh, the datatype of the mesh variables.   The additional static parameter is necessary because one field of EulerData_ has type Tmsh. Note that there could be multiple implementations of AbstractSolutionData for the Euler equations, perhaps with different fields to store certain data or not. All these implementations will need to have the static parameters Tsol, Tdim, Tres, and var_type, so EulerData is defined as an abstract type,  allowing all implementations to inherit from it. All high level functions involved in evaluating the residual will take in an argument of type EulerData. Only when low level functions need to dispatch based on which implementation is  used would it take in an EulerData_ or another implementation."
},

{
    "location": "solver/Readme.html#Variable-Conversion-1",
    "page": "Overview of Physics Modules",
    "title": "Variable Conversion",
    "category": "section",
    "text": "Some equations can be written in different variables, and need to convert between them.  To do this, it is function convertFromNaturalToWorkingVars{Tsol}(params::ParamType{2, :var_type},                qc::AbstractArray{Tsol,1}, qe::AbstractArray{Tsol,1})that converts from the \"natural\" variables in which to write an equation to some other set of variables at a node.  For the Euler equations, the \"natural\" variables would be the conservative variables, and one example of \"other\" variables would be the entropy variables.It is also sometimes useful to define the opposite conversion, ie. from the working variables to the natural variables."
},

{
    "location": "solver/Readme.html#Input-Options-1",
    "page": "Overview of Physics Modules",
    "title": "Input Options",
    "category": "section",
    "text": "Many of the components of PDESolver have different options that control how they work and what they do. In order to  provide a unified method of specifying these options, an dictionary  of type Dict{ASCIIString, Any} is read in from a disk file. This dictionary (called opts in function signatures), is passed to all high and mid level functions so they can use values in the dictionary to determine their  control flow. Low level functions need to be extremely efficient, so they cannot have conditional logic, therefore they are not passed the dictionary. Note that retrieving values from a dictionary is very slow compared to accessing the fields of a type, so all values that are accessed repeatedly should be stored  as the field of a type."
},

{
    "location": "solver/Readme.html#Functors-1",
    "page": "Overview of Physics Modules",
    "title": "Functors",
    "category": "section",
    "text": "Functors are a trick used to get Julia's dispatch system to make decisions at compile time rather than runtime.  This is particularly useful for boundary conditions, where the list of mesh faces that have boundary conditions applied is determined at runtime, but having conditional statements that execute for every node on the mesh boundary would be slow.  Instead a construct is used as follows:type myBC <: BCType  # create a singleton type\nend\n\nfunction call(obj::myBC, q::AbstractVector, bndryflux::AbstractVector)\n  # calculate boundary flux here\nendThis defines a datatype and adds a method to the call function for that type. The call function is what makes a datatype callable like a function.  This method is called as follows:functor = myBC()  # construct and object of type myBC\nq = rand(4)\nbndryflux = zeros(4)\nfunctor(q, bndryflux)  # the Julia compiler turns this into call(functor, q, bndryflux)  The way this is used for boundary conditions is through a two level construct where an outer function passes a functor to an inner function.  Julia's JIT will generate a method of the inner function that is specialized to the functor (this is why it is important that the functor is a datatype).  For example:function getBCFluxes(mesh, sbp, eqn, opts)\n\n  for i=1:mesh.numBC  # loop over different boundary conditions\n    functor_i = mesh.bndry_functor[i]  # get the functor for this boundary condition\n    start_index = mesh.bndry_offsets[i]\n    end_index = mesh.bndry_offsets[i+1] - 1\n    # get data for boundary faces start_index:end_index\n\n    calcBoundaryFlux(functor_i, data for boundary faces start_index:end_index)\n  end\nend  # end function\n\n  function calcBoundaryFlux(functor_i::BCType, data for boundary faces start_index:end_index)\n    for i=1:length(start_index:end_index)\n      for j=1:num_nodes_on_face\n        # get data for this boundary face node\n        functor_i(data for this boundary face node)\n      end\n    end\n\n  end  # end functionThe benefit of this arrangement is that mesh.numBC different version of calcBoundaryFlux get compiled, one for each functor, and each version knows about the call method that was defined for the functor it is passed.  This two level scheme allows the compiler to make all the decisions about what function to call (ie. the call method of the functor), avoiding any conditional logic at runtimeThis idea is also applicable to the flux functions used by DG methods."
},

{
    "location": "solver/Readme.html#Initialization-of-a-Simulation-1",
    "page": "Overview of Physics Modules",
    "title": "Initialization of a Simulation",
    "category": "section",
    "text": "This section lists an outline of how a simulation gets launched After step 4, the procedure becomes a bit more complicated because there are optional steps. Only the required steps are listed below.The options dictionary is read in.  Default values are supplied for any key that is not specified, if a reasonable default value exists.\nSecond, the \nsbp\n operator is constructed.\nThe \nmesh\n object is constructed, using the options dictionary and the \nsbp\n operator.  Some of the options in the dictionary are used to determine how the mesh gets constructed.  For example, the options dictionary specifies what kind of mesh coloring to do.\nThe \neqn\n object is constructed, using the \nmesh\n, \nsbp\n, and \nopts\n objects\nThe physics module \ninit\n function is called, which initializes the physics module and finishes any initialization that \nmesh\n and \neqn\n objects require.\nThe initial condition is applied to \neqn.q_vec\n.\nA nonlinear solver is called.  Which solver is called and what parameters it uses are determined by the options dictionary.\nPost-processing is done, if required by the options dictionary."
},

{
    "location": "solver/misc.html#",
    "page": "Assorted Function and Types",
    "title": "Assorted Function and Types",
    "category": "page",
    "text": ""
},

{
    "location": "solver/misc.html#Assorted-Function-and-Types-1",
    "page": "Assorted Function and Types",
    "title": "Assorted Function and Types",
    "category": "section",
    "text": "  CurrentModule = ODLCommonToolsThis page contains several types and functions that are used throughout the physics modules.  Many of these are defined in ODLCommonTools"
},

{
    "location": "solver/misc.html#ODLCommonTools.BCType",
    "page": "Assorted Function and Types",
    "title": "ODLCommonTools.BCType",
    "category": "Type",
    "text": "Abstract supertype of all boundary condition functors\n\n\n\n"
},

{
    "location": "solver/misc.html#ODLCommonTools.BCType_revm",
    "page": "Assorted Function and Types",
    "title": "ODLCommonTools.BCType_revm",
    "category": "Type",
    "text": "Abstract supertype of all boundary condition functors that compute the   reverse mode with respect to the metrics\n\n\n\n"
},

{
    "location": "solver/misc.html#ODLCommonTools.SRCType",
    "page": "Assorted Function and Types",
    "title": "ODLCommonTools.SRCType",
    "category": "Type",
    "text": "Abstract supertype of all source term functors\n\n\n\n"
},

{
    "location": "solver/misc.html#ODLCommonTools.FluxType",
    "page": "Assorted Function and Types",
    "title": "ODLCommonTools.FluxType",
    "category": "Type",
    "text": "Abstract supertype of all numerical flux functions used by standard DG face   integrals\n\n\n\n"
},

{
    "location": "solver/misc.html#ODLCommonTools.FluxType_revm",
    "page": "Assorted Function and Types",
    "title": "ODLCommonTools.FluxType_revm",
    "category": "Type",
    "text": "Abstract supertype of all numerical flux functions used by standard DG   face integral that compute the reverse mode with respect to the metrics\n\n\n\n"
},

{
    "location": "solver/misc.html#Abstract-Functor-Types-1",
    "page": "Assorted Function and Types",
    "title": "Abstract Functor Types",
    "category": "section",
    "text": "Abstract types are provided for commonly used Functors:BCType\nBCType_revm\nSRCType\nFluxType\nFluxType_revm"
},

{
    "location": "solver/misc.html#ODLCommonTools.Boundary",
    "page": "Assorted Function and Types",
    "title": "ODLCommonTools.Boundary",
    "category": "Type",
    "text": "ODLCommonTools.Boundary\n\nUsed to identify boundary faces in a finite-element grid.\n\nFields\n\nelement\n : index of the element to which the boundary face belongs\nface\n : the face index of the boundary (local index to the element)\n\nExample\n\nTo mark face 2 of element 7 to be a boundary face, use Boundary(7,2)\n\n\n\n"
},

{
    "location": "solver/misc.html#ODLCommonTools.Interface",
    "page": "Assorted Function and Types",
    "title": "ODLCommonTools.Interface",
    "category": "Type",
    "text": "ODLCommonTools.Interface\n\nUsed to identify interfaces between elements in a finite-element grid.\n\nFields\n\nelementL\n : index of the so-called left element in the pair\nelementR\n : index of the so-called right element in the pair\nfaceL\n : the face index of the interface with respect to the left element\nfaceR\n : the face index of the interface with respect to the right element\norient\n : orientation of the 'right' element relative to the 'left'\n\nExample\n\nConsider an interface between elements 2 and 5.  Suppose the interface is on face 1 of element 2 and face 3 of element 5.  Furthermore, suppose element 5 has orientation 1 relative to element 1 (defintion of orientation TBD).  This can be indicated as Interface(2,5,1,3,1)\n\n\n\n"
},

{
    "location": "solver/misc.html#ODLCommonTools.getElementL",
    "page": "Assorted Function and Types",
    "title": "ODLCommonTools.getElementL",
    "category": "Function",
    "text": "This function returns either the element field of a Boundary or the   elementL field of an interface.\n\n\n\n"
},

{
    "location": "solver/misc.html#ODLCommonTools.getFaceL",
    "page": "Assorted Function and Types",
    "title": "ODLCommonTools.getFaceL",
    "category": "Function",
    "text": "This function returns either the face field of a Boundary or the   faceL field of an Interface\n\n\n\n"
},

{
    "location": "solver/misc.html#Base.show-Tuple{IO,ODLCommonTools.Boundary}",
    "page": "Assorted Function and Types",
    "title": "Base.show",
    "category": "Method",
    "text": "Show method for Boundary objects\n\n\n\n"
},

{
    "location": "solver/misc.html#Base.show-Tuple{IO,ODLCommonTools.Interface}",
    "page": "Assorted Function and Types",
    "title": "Base.show",
    "category": "Method",
    "text": "Show method for Interface objects\n\n\n\n"
},

{
    "location": "solver/misc.html#Boundaries-and-Interfaces-1",
    "page": "Assorted Function and Types",
    "title": "Boundaries and Interfaces",
    "category": "section",
    "text": "All physics modules need to apply boundary conditions and all DG schemes and some CG schemes need to do face integrals. The types and functions described here assist in identifying them:Boundary\nInterface\ngetElementL\ngetFaceL\nshow(::IO, ::Boundary)\nshow(::IO, ::Interface)"
},

{
    "location": "solver/advection/advection.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "solver/advection/advection.html#Advection-Physics-Documentation-1",
    "page": "Introduction",
    "title": "Advection Physics Documentation",
    "category": "section",
    "text": "  CurrentModule = AdvectionEquationModThis module evaluates the residual for the constant-coefficient advection equationfracpartial qpartial t = - boldsymbola cdot fracpartial qpartial boldsymbolx + Swhere boldsymbola is the vector of advection velocities in the x, y, and z directions, boldsymbolx is the vector and x, y, and z directions, and S is the source term.This weak form of the equation is discretized as described in the Introduction.  Pages = [ \"types.md\"\n            \"volume.md\"\n            \"flux.md\"\n            \"bc.md\"\n            \"ic.md\"\n            \"source.md\"\n            \"common.md\"\n            \"adjoint.md\"\n            \"boundary_functional.md\"\n          ]\n  Depth = 1"
},

{
    "location": "solver/advection/types.html#",
    "page": "Datatypes",
    "title": "Datatypes",
    "category": "page",
    "text": ""
},

{
    "location": "solver/advection/types.html#AdvectionEquationMod.ParamType",
    "page": "Datatypes",
    "title": "AdvectionEquationMod.ParamType",
    "category": "Type",
    "text": "Subtype of AbstractParamType.\n\nStatic Parameters:     Tsol     Tres     Tdim\n\nThis is a container passed to all low level function, useful for storing   miscellaneous parameters or constants\n\n\n\n"
},

{
    "location": "solver/advection/types.html#AdvectionEquationMod.ParamType2",
    "page": "Datatypes",
    "title": "AdvectionEquationMod.ParamType2",
    "category": "Constant",
    "text": "Convenient alias for all 2D ParamTypes\n\n\n\n"
},

{
    "location": "solver/advection/types.html#AdvectionEquationMod.ParamType3",
    "page": "Datatypes",
    "title": "AdvectionEquationMod.ParamType3",
    "category": "Constant",
    "text": "Convenient alias for all 3D ParamTypes\n\n\n\n"
},

{
    "location": "solver/advection/types.html#AdvectionEquationMod.AbstractAdvectionData",
    "page": "Datatypes",
    "title": "AdvectionEquationMod.AbstractAdvectionData",
    "category": "Type",
    "text": "Direct subtype of AbstractSolutionData, inheriting Tsol and   Tres as static parameter\n\n\n\n"
},

{
    "location": "solver/advection/types.html#AdvectionEquationMod.AdvectionData",
    "page": "Datatypes",
    "title": "AdvectionEquationMod.AdvectionData",
    "category": "Type",
    "text": "Subtype of AbstractAdvectionData, inheriting its static parameters   and adding Tdim.\n\n\n\n"
},

{
    "location": "solver/advection/types.html#AdvectionEquationMod.AdvectionData_",
    "page": "Datatypes",
    "title": "AdvectionEquationMod.AdvectionData_",
    "category": "Type",
    "text": "AdvectionEquationMod.AdvectionData_\n\nThis type is an implementation of the abstract AdvectionData.  It is   parameterized by Tsol, the datatype of the solution variables and Tmsh,   the datatype of the mesh variables.   Tres is the 'maximum' type of Tsol and Tmsh.   Tdim is the dimensionality of the equation being solve (2 or 3).\n\nThis type is (ultimately) a subtype of AbstractSolutionData    and contains all the required fields.\n\n\n\n"
},

{
    "location": "solver/advection/types.html#AdvectionEquationMod.PhysicsName",
    "page": "Datatypes",
    "title": "AdvectionEquationMod.PhysicsName",
    "category": "Constant",
    "text": "This physics is named Advection\n\n\n\n"
},

{
    "location": "solver/advection/types.html#Advection-Types-1",
    "page": "Datatypes",
    "title": "Advection Types",
    "category": "section",
    "text": "  CurrentModule = AdvectionEquationModThis page provides the documentation for the DataTypes defined in the Advection moduleParamType\nParamType2\nParamType3\nAbstractAdvectionData\nAdvectionData\nAdvectionData_\nPhysicsName"
},

{
    "location": "solver/advection/volume.html#",
    "page": "Volume Integrals",
    "title": "Volume Integrals",
    "category": "page",
    "text": ""
},

{
    "location": "solver/advection/volume.html#AdvectionEquationMod.evalVolumeIntegrals",
    "page": "Volume Integrals",
    "title": "AdvectionEquationMod.evalVolumeIntegrals",
    "category": "Function",
    "text": "AdvectionEquationMod.evalVolumeIntegrals\n\nEvaluates the volume integrals of the weak form.  eqn.res is updated   with the result.  Both the precompute_volume_flux and non    precompute_volume_flux versions are contained within this function\n\nInputs:\n\nmesh : mesh type    sbp  : Summation-by-parts operator    eqn  : Advection equation object    opts : options dictionary\n\nOutputs\n\nNone\n\n\n\n"
},

{
    "location": "solver/advection/volume.html#AdvectionEquationMod.calcAdvectionFlux",
    "page": "Volume Integrals",
    "title": "AdvectionEquationMod.calcAdvectionFlux",
    "category": "Function",
    "text": "Populates eqn.flux_parametric.  Repeatedly calls the other method of this   function.\n\nInputs:\n\nmesh\nsbp\neqn\nopts\n\n\n\nCalculates the advection flux in the parametric directions at a node.\n\nInputs:\n\nparams: a ParamType object\nq: the solution value at the node\nalphas_xy: the advection velocities in the x-y directions, vector of length\n           Tdim\ndxidx: scaled mapping jacobian at the node, Tdim x Tdim matrix\n\nInputs/Outputs:\n\nflux: vector of length Tdim to populate with the flux in the parametric\n      directions\n\n\n\n"
},

{
    "location": "solver/advection/volume.html#Advection-Volume-Integrals-1",
    "page": "Volume Integrals",
    "title": "Advection Volume Integrals",
    "category": "section",
    "text": "  CurrentModule = AdvectionEquationModThis page describes the different functions involved in computing the volume integralsevalVolumeIntegrals\ncalcAdvectionFlux"
},

{
    "location": "solver/advection/flux.html#",
    "page": "Face Integrals",
    "title": "Face Integrals",
    "category": "page",
    "text": ""
},

{
    "location": "solver/advection/flux.html#AdvectionEquationMod.evalFaceIntegrals",
    "page": "Face Integrals",
    "title": "AdvectionEquationMod.evalFaceIntegrals",
    "category": "Function",
    "text": "AdvectionEquationMod.evalFaceIntegrals\n\nThis function evaluates the interior face integrals for DG methods, using   the flux function from eqn.flux_func.  The solution variables are interpolated   to the faces, the flux computed, and then interpolated back to the   solution points.\n\nThis function also logs some quantities to disk (TODO: move this to   Utils/logging)\n\nInputs:     mesh:  an AbstractDGMesh     sbp     eqn     opts\n\n\n\n"
},

{
    "location": "solver/advection/flux.html#AdvectionEquationMod.calcFaceFlux",
    "page": "Face Integrals",
    "title": "AdvectionEquationMod.calcFaceFlux",
    "category": "Function",
    "text": "This function calculates the DG flux between a specified set of faces,   using the solution data at the faces stored in eqn.q_face.   Note that the flux is negated because the face integrals have a    negative sign in the weak form.\n\nInputs:\n\nmesh\nsbp\neqn\nfunctor: the functor that calculates the flux at a node\ninterfaces: an array of type Interface that specifies which interfaces\n            to calculate the flux for\n\nInputs/Outputs:\n\nface_flux: array to store the flux in, numDofPerNode x nnodesPerFace\n           x length(interfaces)\n\nThe functor must have the signature:   func( uL, qR, alpha_x, alpha_y, dxidx, nrm, params)\n\nwhere uL and uR are the solution values for a node on the left and right   elements, alpha_x and alpha_y are the x and y advection velocities,   dxidx is the scaled mapping jacobian for elementL, and nrm is the face   normal in reference space.  params is eqn.params\n\n\n\n"
},

{
    "location": "solver/advection/flux.html#AdvectionEquationMod.calcFaceIntegrals_nopre",
    "page": "Face Integrals",
    "title": "AdvectionEquationMod.calcFaceIntegrals_nopre",
    "category": "Function",
    "text": "Compute the face integrals without using eqn.q_face or eqn.flux_face.   The integral is computed directly and res is updated\n\nInputs:\n\nmesh\nsbp\neqn\nopts\nflux_func: the flux functor that computes the face flux at a node\n\n\n\n"
},

{
    "location": "solver/advection/flux.html#AdvectionEquationMod.calcSharedFaceIntegrals",
    "page": "Face Integrals",
    "title": "AdvectionEquationMod.calcSharedFaceIntegrals",
    "category": "Function",
    "text": "Thin wrapper around calcSharedFaceIntegrals_inner.  This function is passed   to finishDataExchange, and internally calls calcSharedFaceIntegrals_inner.   See finishExchangeData for details on the interface and    calcSharedFaceIntegrals_inner for the integral that is computed.\n\n\n\n"
},

{
    "location": "solver/advection/flux.html#AdvectionEquationMod.calcSharedFaceIntegrals_inner",
    "page": "Face Integrals",
    "title": "AdvectionEquationMod.calcSharedFaceIntegrals_inner",
    "category": "Function",
    "text": "AdvectionEquationMod.calcSharedFaceIntegrals\n\nThis function calculates the shared face integrals for the faces shared   with a single peer process.  This function is for   opts[\"parallel_type\"] == \"face\" and regular face integrals (ie. not the   entropy-stable face integrals) only.\n\nInputs:\n\nmesh\nsbp\neqn\nopts:\ndata: a SharedFaceData specifying which faces to compute\nfunctor: the FluxType to use for the face flux\n\n\n\n"
},

{
    "location": "solver/advection/flux.html#AdvectionEquationMod.calcSharedFaceIntegrals_inner_nopre",
    "page": "Face Integrals",
    "title": "AdvectionEquationMod.calcSharedFaceIntegrals_inner_nopre",
    "category": "Function",
    "text": "Like calcSharedFaceIntegrals_inner_nopre, but it computes the integral one   face at a time rather than computing all the flux, storing it in   eqn.flux_sharedface and then doing the integral\n\nSee calcSharedFaceIntegrals_inner for a description of the arguments\n\n\n\n"
},

{
    "location": "solver/advection/flux.html#AdvectionEquationMod.calcSharedFaceIntegrals_element",
    "page": "Face Integrals",
    "title": "AdvectionEquationMod.calcSharedFaceIntegrals_element",
    "category": "Function",
    "text": "Thin wrapper around calcSharedFaceIntegrals_inner.  This function is passed   to finishDataExchange, and internally calls calcSharedFaceIntegrals_inner.   See finishDataExchange for details on the interface and    calcSharedFaceIntegrals_inner for the integral that is computed.\n\n\n\n"
},

{
    "location": "solver/advection/flux.html#AdvectionEquationMod.calcSharedFaceIntegrals_element_inner",
    "page": "Face Integrals",
    "title": "AdvectionEquationMod.calcSharedFaceIntegrals_element_inner",
    "category": "Function",
    "text": "Like calcSharedFaceIntegrals_inner, but for the case when   opts[\"parallel_data\"] == element.  This effectively means it has to   interpolate the solution from the elements to the faces and then do the   integral\n\n\n\n"
},

{
    "location": "solver/advection/flux.html#AdvectionEquationMod.calcSharedFaceIntegrals_element_inner_nopre",
    "page": "Face Integrals",
    "title": "AdvectionEquationMod.calcSharedFaceIntegrals_element_inner_nopre",
    "category": "Function",
    "text": "Like calcSharedFaceIntegrals_element_inner, but computes the   integral one   face at a time instead of computing the entire flux and then integrating.\n\n\n\n"
},

{
    "location": "solver/advection/flux.html#Advection-Face-Integrals-1",
    "page": "Face Integrals",
    "title": "Advection Face Integrals",
    "category": "section",
    "text": "  CurrentModule = AdvectionEquationModThis page describes the functions that compute the face integrals used by DG schemes.evalFaceIntegrals\ncalcFaceFlux\ncalcFaceIntegrals_nopre\ncalcSharedFaceIntegrals\ncalcSharedFaceIntegrals_inner\ncalcSharedFaceIntegrals_inner_nopre\ncalcSharedFaceIntegrals_element\ncalcSharedFaceIntegrals_element_inner\ncalcSharedFaceIntegrals_element_inner_nopre"
},

{
    "location": "solver/advection/bc.html#",
    "page": "Boundary Integrals",
    "title": "Boundary Integrals",
    "category": "page",
    "text": ""
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.evalBoundaryIntegrals",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.evalBoundaryIntegrals",
    "category": "Function",
    "text": "Evaluate boundary integrals for advection equation, updating eqn.res with the result.\n\nInputs\n\n \nmesh\n : Abstract mesh type\n \nsbp\n  : Summation-by-parts operator\n \neqn\n  : Advection equation object\n \nopts\n : options dictionary\n\nOutputs\n\n None\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.calcBoundaryFlux",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.calcBoundaryFlux",
    "category": "Function",
    "text": "This function calculates the boundary flux for the portion of the boundary   with a particular boundary condition.  The eqn.q are converted to    conservative variables if needed.  For the DG version, eqn.q_bndry must   already be populated with the q variables interpolated to the boundary\n\nInputs:\n\nmesh : AbstractMesh   sbp : AbstractSBP   eqn : AdvectionEquation   functor : a callable object that calculates the boundary flux at a node   idx_range: the Range describing which Boundaries have the current BC   bndry_facenums:  An array with elements of type Boundary that tell which                    element faces have the boundary condition   Outputs:\n\nbndryflux : the array to store the boundary flux, corresponds to                bndry_facenums\n\nnote that bndry_facenums and bndryflux must be only the portion of the    their parent arrays that correspond to the Boundaries that have the    current boundary condition applied.\n\nThe functor must have the signature:   functor( q, aux_vars, x, dxidx, nrm, bndryflux_i, eqn.params)   where q are the conservative variables.   where all arguments (except params and nrm) are vectors of values at a node.\n\nparams is the ParamType associated with the the EulerEquation object   nrm = mesh.sbpface.normal[:, current_node]\n\nThis is a mid level function.\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.calcBoundaryFlux_nopre",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.calcBoundaryFlux_nopre",
    "category": "Function",
    "text": "This function computes the boundary integrals (and should probably be renamed)   without using eqn.q_bndry of eqn.bndryflux.  eqn.res is updated with   the results.\n\nSee calcBoundaryFlux for the meaning of the arguments\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.BCDict",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.BCDict",
    "category": "Constant",
    "text": "AdvectionEquationMod.BCDict\n\nIt stores all the possible boundary condition dictionary options. Whenever a  new boundary condition is created, it should get added to BCDict.\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.getBCFunctors",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.getBCFunctors",
    "category": "Function",
    "text": "AdvectionEquationMod.getBCFunctors\n\nThis function uses the opts dictionary to populate mesh.bndry_funcs with the the functors\n\nThis is a high level function.\n\nInputs\n\n \nmesh\n : Abstract mesh type\n \nsbp\n  : Summation-by-parts operator\n \neqn\n  : Advection equation object\n \nopts\n : Input dictionary options\n\nOutputs\n\n None\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#Advection-Boundary-Integrals-1",
    "page": "Boundary Integrals",
    "title": "Advection Boundary Integrals",
    "category": "section",
    "text": "  CurrentModule = AdvectionEquationModThis page describes the functions that compute the boundary integrals.  evalBoundaryIntegrals\n  calcBoundaryFlux\n  calcBoundaryFlux_nopre\n  BCDict\n  getBCFunctors"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.defaultBC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.defaultBC",
    "category": "Type",
    "text": "Default BC to calculate the boundary face integral (no numerical flux   functions)\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.exp2xplus2yBC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.exp2xplus2yBC",
    "category": "Type",
    "text": "AdvectionEquationMod.exp2xplus2yBC\n\nUses the Roe solver to calculate the boundary flux using calc_exp2xplus2y to get the boundary state.\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.exp3xplusyBC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.exp3xplusyBC",
    "category": "Type",
    "text": "AdvectionEquationMod.exp3xplusyBC\n\nUses the Roe solver to calculate the boundary flux using calc_exp3xplusy to get the boundary state.\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.exp5xplus4yplus2BC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.exp5xplus4yplus2BC",
    "category": "Type",
    "text": "AdvectionEquationMod.exp5xplus4yplus2BC\n\nUses the Roe solver to calculate the boundary flux using calc_exp5xplus4yplus2  to get the boundary state.\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.exp5xplusyBC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.exp5xplusyBC",
    "category": "Type",
    "text": "AdvectionEquationMod.exp5xplusyBC\n\nUses the Roe solver to calculate the boundary flux using calc_exp5xplusy to get the boundary state.\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.exp_xplusyBC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.exp_xplusyBC",
    "category": "Type",
    "text": "AdvectionEquationMod.exp_xplusyBC\n\nCalculates q at the boundary which is equal to exp(x+y). It is a nodal  level function.\n\nInputs\n\n \nu\n : Advection variable (eqn.q)\n \nalpha_x\n & \nalpha_y\n : velocities in the X & Y directions\n \ncoords\n : Nodal coordinates\n \ndxidx\n  : Mapping Jacobian\n \nnrm\n    : SBP face-normal vectors\n \nbndryflux\n : Flux at the boundary\n\nOutputs\n\n None\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.exp_xyBC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.exp_xyBC",
    "category": "Type",
    "text": "AdvectionEquationMod.exp_xyBC\n\nUses the Roe solver to calculate the boundary flux using calc_exp_xy to get the boundary state\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.mms1BC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.mms1BC",
    "category": "Type",
    "text": "AdvectionEquationMod.mms1BC\n\nUses the Roe solver to calculate the boundary flux using calc_mms1 to get   the boundary state.\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.p1BC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.p1BC",
    "category": "Type",
    "text": "AdvectionEquationMod.p1BC\n\nUses the Roe solver to calculate the boundary flux using calc_p1 to   get the boundary state\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.p2BC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.p2BC",
    "category": "Type",
    "text": "AdvectionEquationMod.p2BC\n\nUses the Roe solver to calculate the boundary flux using calc_p2 to   get the boundary state\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.p3BC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.p3BC",
    "category": "Type",
    "text": "AdvectionEquationMod.p3BC\n\nUses the Roe solver to calculate the boundary flux using calc_p3 to   get the boundary state\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.p4BC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.p4BC",
    "category": "Type",
    "text": "AdvectionEquationMod.p4BC\n\nUses the Roe solver to calculate the boundary flux using calc_p4 to   get the boundary state.\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.p5BC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.p5BC",
    "category": "Type",
    "text": "AdvectionEquationMod.p5BC\n\nUses the Roe solver to calculate the boundary flux using calc_p5 to   get the boundary state.\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.sinwave_BC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.sinwave_BC",
    "category": "Type",
    "text": "AdvectionEquationMod.sinwave_BC\n\nUses the Roe solver to calculate the boundary flux using calc_sinewave to   get the boundary state\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.sinwavey_BC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.sinwavey_BC",
    "category": "Type",
    "text": "AdvectionEquationMod.sinwavey_BC\n\nUses the Roe solver to calculate the boundary flux using calc_sinewavey to   get the boundary state\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.sinwavey_pertBC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.sinwavey_pertBC",
    "category": "Type",
    "text": "AdvectionEquationMod.sinwavey_pertBC\n\nUses the Roe solver to calculate the boundary flux using calc_sinewave_pert to   get the boundary state\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.unsteadymmsBC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.unsteadymmsBC",
    "category": "Type",
    "text": "BC for unsteadymms\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.unsteadypolyBC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.unsteadypolyBC",
    "category": "Type",
    "text": "BC for unsteadypoly\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.x4BC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.x4BC",
    "category": "Type",
    "text": "AdvectionEquationMod.x4BC\n\nUses the Roe solver to calculate the boundary flux using calc_x4 to   get the boundary state.\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.x5plusy5BC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.x5plusy5BC",
    "category": "Type",
    "text": "AdvectionEquationMod.x5plusy5BC\n\nCalculates q at the boundary which is equal to x^5 + y^5. It is a nodal  level function.\n\nInputs\n\n \nu\n : Advection variable (eqn.q)\n \nparams\n: the equation ParamType\n \ncoords\n : Nodal coordinates\n \nnrm_scaled\n    : scaled face normal vector in x-y space\n \nt\n:  current time value\n\nOutputs *  bndryflux : Flux at the boundary\n\n None\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.xplusyBC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.xplusyBC",
    "category": "Type",
    "text": "AdvectionEquationMod.xplusyBC\n\nUses Roe solver to calculate the boundary flux using calc_xplusy to get the  boundary state\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#Boundary-Conditions-1",
    "page": "Boundary Integrals",
    "title": "Boundary Conditions",
    "category": "section",
    "text": "This section describes all boundary conditions currently available  Modules = [AdvectionEquationMod]\n  Pages = [\"advection/bc.jl\"]\n  Order = [:type]"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.RoeSolver-Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2},Tsol,Any,Any}",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.RoeSolver",
    "category": "Method",
    "text": "AdvectionEquationMod.RoeSolver\n\nRoe solver for the advection equations. It determines the boundary flux on  each boundary. It is called at the nodal level\n\nInputs\n\n \nu\n    : Solution of advection equation at a particular node\n \nu_bc\n : Prescribed solution value at the boundary\n \nparams\n: the equation ParamType object\n \nnrm\n  : scaled face normal vector in x-y space\n\nOutputs\n\n \nbndryflux\n : Boundary flux at the particular node\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.flux1-Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2},Any,Any,Any,Any}",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.flux1",
    "category": "Method",
    "text": "flux1\n\nCalculates the boundary flux for the advection equation. It works at the nodal level.\n\nInputs\n\n \nu_sbp_\n: The entry from u_sbp for this node\n \ndxidx\n : The jacobian for this node\n \nnrm\n   : nrm is the normal vector\n \nnet_flux\n:\n \nparams\n: the equation ParamType\n\nOutputs\n\n None\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#Numerical-Flux-Functions-1",
    "page": "Boundary Integrals",
    "title": "Numerical Flux Functions",
    "category": "section",
    "text": "This section lists the numerical flux functions used to impose the boundary conditions weakly.  Modules = [AdvectionEquationMod]\n  Pages = [\"advection/bc_solvers.jl\"]"
},

{
    "location": "solver/advection/ic.html#",
    "page": "Initial Condition",
    "title": "Initial Condition",
    "category": "page",
    "text": ""
},

{
    "location": "solver/advection/ic.html#AdvectionEquationMod.ICDict",
    "page": "Initial Condition",
    "title": "AdvectionEquationMod.ICDict",
    "category": "Constant",
    "text": "Dictionary that maps IC names to functions.  Every new IC should be added   to the list\n\n\n\n"
},

{
    "location": "solver/advection/ic.html#AdvectionEquationMod.ICFile-Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP{Tsbp},AdvectionEquationMod.AdvectionData{Tsol,Tres,Tdim},Any,AbstractArray{Tsol,1}}",
    "page": "Initial Condition",
    "title": "AdvectionEquationMod.ICFile",
    "category": "Method",
    "text": "AdvectionEquationMod.ICFile\n\nThis function reads a vector from a file on disk and set the solution to it. The vector must contain the same number of entries as there are degrees of  freedom in the mesh. \n\nThis function is useful for things like restarting from a checkpoint. In this case, the file should be the output of writedlm(eqn.q).  The degree  of freedom number must be the same for both simulation for this to work (the  file contains no degree of freedom number information).\n\nInputs\n\n \nmesh\n \nsbp\n \neqn\n \nopts\n\nInputs/Outputs\n\n \nu0\n: vector to populate with the solution\n\nAliasing restrictions: none.\n\n\n\n"
},

{
    "location": "solver/advection/ic.html#AdvectionEquationMod.ICexp_xplusy-Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP{Tsbp},AdvectionEquationMod.AdvectionData{Tsol,Tres,2},Any,AbstractArray{Tsol,N}}",
    "page": "Initial Condition",
    "title": "AdvectionEquationMod.ICexp_xplusy",
    "category": "Method",
    "text": "AdvectionEquationMod.ICexp_xplusy\n\nComputes the initial conditions for the state variable     u = exp(x + y + z), z is omitted in 2d\n\nInputs\n\n \nmesh\n : AbstractMesh type\n \nsbp\n  : Summation-by-parts operator\n \neqn\n  : Advection equation object\n \nopts\n : Options dictionary\n \nu0\n   : Array that stores inital state variables\n\nOutputs\n\n None\n\n\n\n"
},

{
    "location": "solver/advection/ic.html#AdvectionEquationMod.ICx5plusy5-Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP{Tsbp},AdvectionEquationMod.AdvectionData{Tsol,Tres,2},Any,AbstractArray{Tsol,N}}",
    "page": "Initial Condition",
    "title": "AdvectionEquationMod.ICx5plusy5",
    "category": "Method",
    "text": "AdvectionEquationMod.ICx5plusy5\n\nComputes the initial conditions for the state variable     u = x^5 + y^5 + z^5; z is omitted it 2d\n\nInputs\n\n \nmesh\n : AbstractMesh type\n \nsbp\n  : Summation-by-parts operator\n \neqn\n  : Advection equation object\n \nopts\n : Options dictionary\n \nu0\n   : Array that stores inital state variables\n\nOutputs\n\n None\n\n\n\n"
},

{
    "location": "solver/advection/ic.html#Advection-Initial-Conditions-1",
    "page": "Initial Condition",
    "title": "Advection Initial Conditions",
    "category": "section",
    "text": "This pages describes the functions that apply initial conditions  CurrentModule = AdvectionEquationMod  Modules = [AdvectionEquationMod]\n  Pages = [\"advection/ic.jl\"]"
},

{
    "location": "solver/advection/source.html#",
    "page": "Source Term",
    "title": "Source Term",
    "category": "page",
    "text": ""
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCDict",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCDict",
    "category": "Constant",
    "text": "AdvectionEquationMod.SRCDict\n\nIt stores all the possible boundary condition dictionary options. Whenever a    new boundary condition is created, it should get added to BCDict.\n\nAll functors must have the signature:\n\nsrc_func(params, coords::ParamType, t)\n\nwhere coords is the vector of length 2 containing the x and y coordinates   of the node, params.alpha_x and params.alpha_y are the advection velocities in the x an y   directions, and t is the current time\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRC0",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRC0",
    "category": "Type",
    "text": "AdvectionEquationMod.SRC0\n\nThis is the zero source term.  This is the default of source term   is specified\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRC1",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRC1",
    "category": "Type",
    "text": "AdvectionEquationMod.SRC1\n\nThis source term returns 1 everywhere.\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRC2",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRC2",
    "category": "Type",
    "text": "AdvectionEquationMod.SRC2\n\nThis source term returns 1 everywhere.\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCexp2xplus2y",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCexp2xplus2y",
    "category": "Type",
    "text": "AdvectionEquationMod.SRCexp2xplus2y\n\nCalculates the source term for q = exp(2x + 2y)\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCexp3xplusy",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCexp3xplusy",
    "category": "Type",
    "text": "AdvectionEquationMod.SRCexp3xplusy\n\nCalculates the source term for q = exp(3*x + y)\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCexp5xplus4yplus2",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCexp5xplus4yplus2",
    "category": "Type",
    "text": "AdvectionEquationMod.SRCexp5xplus4yplus2\n\nCalculates the source term for q = exp(5x + 4y +2)\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCexp5xplusy",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCexp5xplusy",
    "category": "Type",
    "text": "AdvectionEquationMod.SRCexp5xplusy\n\nCalculates the source term for q = exp(5*x + y)\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCexp_xplusy",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCexp_xplusy",
    "category": "Type",
    "text": "AdvectionEquationMod.SRCexp_xplusy\n\nThis is a source term that returns a source term for e^(x+y)\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCexp_xy",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCexp_xy",
    "category": "Type",
    "text": "AdvectionEquationMod.SRCexp_xy\n\nCalculates the source term for q = exp(x*y)\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCmms1",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCmms1",
    "category": "Type",
    "text": "AdvectionEquationMod.SRC1\n\nThis source term that returns: the derivative of mms1\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCp1",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCp1",
    "category": "Type",
    "text": "AdvectionEquationMod.SRCp1\n\nThis source term that returns: the source term for a manufactured solution   using a 1st order polynomial\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCp2",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCp2",
    "category": "Type",
    "text": "AdvectionEquationMod.SRCp3\n\nThis source term that returns: the source term for a manufactured solution   using a 2nd order polynomial\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCp3",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCp3",
    "category": "Type",
    "text": "AdvectionEquationMod.SRCp3\n\nThis source term that returns: the source term for a manufactured solution   using a 3rd order polynomial\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCp4",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCp4",
    "category": "Type",
    "text": "AdvectionEquationMod.SRCp4\n\nThis source term that returns: the source term for a manufactured solution   using a 4th order polynomial\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCp5",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCp5",
    "category": "Type",
    "text": "AdvectionEquationMod.SRCp5\n\nThis source term that returns: the source term for a manufactured solution   using a 5th order polynomial\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCunsteadymms",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCunsteadymms",
    "category": "Type",
    "text": "Source term for unsteady mms\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCunsteadypoly",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCunsteadypoly",
    "category": "Type",
    "text": "Source term for unsteady poly\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCx",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCx",
    "category": "Type",
    "text": "AdvectionEquationMod.SRCx\n\nThis source term that returns: f = x\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCx4",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCx4",
    "category": "Type",
    "text": "AdvectionEquationMod.SRCx4\n\nThis source term that returns: the source term for a manufactured solution   using a 4th order polynomial\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCx5plusy5",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCx5plusy5",
    "category": "Type",
    "text": "AdvectionEquationMod.SRCx5plusy5\n\nThis is a source term that returns a source term for e^(x+y)\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCxplusy",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCxplusy",
    "category": "Type",
    "text": "AdvectionEquationMod.SRCxplusy\n\ncalculates the source term for q = x + y\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.getSRCFunctors-Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP{T<:Number},AdvectionEquationMod.AdvectionData{Tsol,Tres,Tdim},Any}",
    "page": "Source Term",
    "title": "AdvectionEquationMod.getSRCFunctors",
    "category": "Method",
    "text": "This function gets the functor specified by opts[\"SRCname\"] and stores   it to the equation object.  Currently one 1 source functor is allowed.\n\n\n\n"
},

{
    "location": "solver/advection/source.html#Advection-Source-Term-1",
    "page": "Source Term",
    "title": "Advection Source Term",
    "category": "section",
    "text": "This pages describes the functions that apply source terms  CurrentModule = AdvectionEquationMod  Modules = [AdvectionEquationMod]\n  Pages = [\"advection/source.jl\"]"
},

{
    "location": "solver/advection/common.html#",
    "page": "Common Functions",
    "title": "Common Functions",
    "category": "page",
    "text": ""
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_exp2xplus2y-Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_exp2xplus2y",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_exp2xplus2y\n\nCalculates and return the expression u = exp(2x + 2y)\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_exp3xplusy-Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_exp3xplusy",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_exp3xplusy\n\nCalculates and return the expression u = exp(3*x + y)\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_exp5xplus4yplus2-Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_exp5xplus4yplus2",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_exp5xplus4yplus2\n\nCalculates and returns the expression u = exp(5x + 4y +2)\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_exp5xplusy-Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_exp5xplusy",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_exp5xplusy\n\nCalculates and return the expression u = exp(5*x + y)\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_exp_xplusy-Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_exp_xplusy",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_exp_xplusy\n\nCalculates and returns e^(x + y)\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_exp_xy-Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_exp_xy",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_exp_xy\n\nCalculates and returns u = exp(x*y)\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_mms1-Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_mms1",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_mms1\n\nCalculates and returns the value of the solution for doing Method of    Manufactured solutions.  This is for debugging only, and could change   at any time.\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_mms1dx-Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_mms1dx",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_mms1dx\n\nCalculates and returns the x derivative of calc_mms1\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p1-Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p1",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_p1\n\nCalculates and returns a 1st order polynomial of x and y\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p1dx-Tuple{Union{AdvectionEquationMod.ParamType{Tsol,Tres,2},AdvectionEquationMod.ParamType{Tsol,Tres,3}},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p1dx",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_p1dx\n\nCalculates and returns the x derivative of calc_p1\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p1dy-Tuple{Union{AdvectionEquationMod.ParamType{Tsol,Tres,2},AdvectionEquationMod.ParamType{Tsol,Tres,3}},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p1dy",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_p1dy\n\nCalculates and returns the y derivative of calc_p1\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p1dz-Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,3},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p1dz",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_p1dz\n\nCalculates and returns the z derivative of calc_p1\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p2-Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p2",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_p2\n\nCalculates and returns a 2nd order polynomial in x and y\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p2dx-Tuple{Union{AdvectionEquationMod.ParamType{Tsol,Tres,2},AdvectionEquationMod.ParamType{Tsol,Tres,3}},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p2dx",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_p2dx\n\nCalculates and returns a the x derivative of calc_p2\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p2dy-Tuple{Union{AdvectionEquationMod.ParamType{Tsol,Tres,2},AdvectionEquationMod.ParamType{Tsol,Tres,3}},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p2dy",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_p2dy\n\nCalculates and returns the y derivative of calc_p2\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p2dz-Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,3},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p2dz",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_p2dz\n\nCalculates and returns the z derivative of calc_p2\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p3-Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p3",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_p3\n\nCalculates and returns a 3rd order polynomial in x and y (and z in 3d)\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p3dx-Tuple{Union{AdvectionEquationMod.ParamType{Tsol,Tres,2},AdvectionEquationMod.ParamType{Tsol,Tres,3}},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p3dx",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_p3dx\n\nCalculates and returns the x derivataive of calc_p3\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p3dy-Tuple{Union{AdvectionEquationMod.ParamType{Tsol,Tres,2},AdvectionEquationMod.ParamType{Tsol,Tres,3}},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p3dy",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_p3dy\n\nCalculates and returns the y derivative of calc_p3\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p3dz-Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,3},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p3dz",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_p3dz\n\nCalculates and returns the z derivative of calc_p3\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p4-Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p4",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_p4\n\nCalculates and returns a 4th order polynomial in x and y\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p4dx-Tuple{Union{AdvectionEquationMod.ParamType{Tsol,Tres,2},AdvectionEquationMod.ParamType{Tsol,Tres,3}},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p4dx",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_p4x\n\nCalculates and returns the x derivative of calc_p4\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p4dy-Tuple{Union{AdvectionEquationMod.ParamType{Tsol,Tres,2},AdvectionEquationMod.ParamType{Tsol,Tres,3}},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p4dy",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_p4dy\n\nCalculates and returns the y derivative of calc_p4\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p4dz-Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,3},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p4dz",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_p4dz\n\nCalculates and returns the z derivative of calc_p4\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p5-Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p5",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_p5\n\nCalculates and returns a 5th order polynomial in x and y (and z in 3d)\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p5dx-Tuple{Union{AdvectionEquationMod.ParamType{Tsol,Tres,2},AdvectionEquationMod.ParamType{Tsol,Tres,3}},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p5dx",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_p5dx\n\nCalculates and returns the x derivative of calc_p5\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p5dy-Tuple{Union{AdvectionEquationMod.ParamType{Tsol,Tres,2},AdvectionEquationMod.ParamType{Tsol,Tres,3}},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p5dy",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_p5y\n\nCalculates and returns the y derivative of calc_p5\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p5dz-Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,3},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p5dz",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_p5z\n\nCalculates and returns the z derivative of calc_p5\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_sinwave-Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_sinwave",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_sinwave\n\nCalculates and returns sin(-x + t)\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_sinwavey-Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_sinwavey",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_sinwavey\n\nCalculates and returns sin(y)^2 + 5*sin(y) + 3/sin(y)\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_sinwavey_pert-Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_sinwavey_pert",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_sinwavey_pert\n\nCalculates and returns 1000sin(x)calc_sinwavey\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_unsteadymms-Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_unsteadymms",
    "category": "Method",
    "text": "u = exp(x + y + z + t) in 3d (z = 0 in 2d)\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_x4-Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_x4",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_x4\n\nCalculates and returns a 4th order polynomial in x\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_x4der-Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_x4der",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_x5plusy5\n\nCalculates and returns the x derivative of calc_x4\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_x5plusy5-Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_x5plusy5",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_x5plusy5\n\nCalculates and returns x^5 + y^5\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_xplusy-Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2},AbstractArray{Tmsh,N},Any}",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_xplusy",
    "category": "Method",
    "text": "AdvectionEquationMod.calc_x5plusy5\n\nCalculates and returns u = x+y\n\n\n\n"
},

{
    "location": "solver/advection/common.html#Advection-Common-Functions-1",
    "page": "Common Functions",
    "title": "Advection Common Functions",
    "category": "section",
    "text": "This page describes functions that are used to calculated quantities that are needed for IC, BCs, and source terms. They are particularly useful for defining manufactures solutions.  CurrentModule = AdvectionEquationMod  Modules = [AdvectionEquationMod]\n  Pages = [\"advection/common_funcs.jl\"]"
},

{
    "location": "solver/advection/adjoint.html#",
    "page": "Adjoint",
    "title": "Adjoint",
    "category": "page",
    "text": ""
},

{
    "location": "solver/advection/adjoint.html#Advection-Adjoint-1",
    "page": "Adjoint",
    "title": "Advection Adjoint",
    "category": "section",
    "text": ""
},

{
    "location": "solver/advection/boundary_functional.html#",
    "page": "Boundary Functional",
    "title": "Boundary Functional",
    "category": "page",
    "text": ""
},

{
    "location": "solver/advection/boundary_functional.html#Advection-Boundary-Functional-1",
    "page": "Boundary Functional",
    "title": "Advection Boundary Functional",
    "category": "section",
    "text": ""
},

{
    "location": "solver/euler/euler.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "solver/euler/euler.html#Euler-Physics-Documentation-1",
    "page": "Introduction",
    "title": "Euler Physics Documentation",
    "category": "section",
    "text": "Describe the equation being solved here  Pages = [ \"advection.md\"\n            \"types.md\"\n            \"volume.md\"\n            \"flux.md\"\n            \"bc.md\"\n            \"ic.md\"\n            \"source.md\"\n            \"common.md\"\n            \"conversion.md\"\n            \"flux_functions.md\"\n            \"stabilization.md\"\n            \"adjoint.md\"\n            \"boundary_functional.md\"\n            \"misc.md\"\n          ]\n  Depth = 1"
},

{
    "location": "solver/euler/types.html#",
    "page": "Datatypes",
    "title": "Datatypes",
    "category": "page",
    "text": ""
},

{
    "location": "solver/euler/types.html#Euler-Datatype-Documentation-1",
    "page": "Datatypes",
    "title": "Euler Datatype Documentation",
    "category": "section",
    "text": ""
},

{
    "location": "solver/euler/volume.html#",
    "page": "Volume Integrals",
    "title": "Volume Integrals",
    "category": "page",
    "text": ""
},

{
    "location": "solver/euler/volume.html#Volume-Integrals-1",
    "page": "Volume Integrals",
    "title": "Volume Integrals",
    "category": "section",
    "text": ""
},

{
    "location": "solver/euler/flux.html#",
    "page": "Face Integrals",
    "title": "Face Integrals",
    "category": "page",
    "text": ""
},

{
    "location": "solver/euler/flux.html#Face-Integral-1",
    "page": "Face Integrals",
    "title": "Face Integral",
    "category": "section",
    "text": ""
},

{
    "location": "solver/euler/bc.html#",
    "page": "Boundary Integrals",
    "title": "Boundary Integrals",
    "category": "page",
    "text": ""
},

{
    "location": "solver/euler/bc.html#Boundary-Integrals-1",
    "page": "Boundary Integrals",
    "title": "Boundary Integrals",
    "category": "section",
    "text": ""
},

{
    "location": "solver/euler/ic.html#",
    "page": "Initial Conditions",
    "title": "Initial Conditions",
    "category": "page",
    "text": ""
},

{
    "location": "solver/euler/ic.html#Initial-Conditions-1",
    "page": "Initial Conditions",
    "title": "Initial Conditions",
    "category": "section",
    "text": ""
},

{
    "location": "solver/euler/source.html#",
    "page": "Source Term",
    "title": "Source Term",
    "category": "page",
    "text": ""
},

{
    "location": "solver/euler/source.html#Source-Term-1",
    "page": "Source Term",
    "title": "Source Term",
    "category": "section",
    "text": ""
},

{
    "location": "solver/euler/common.html#",
    "page": "Common Functions",
    "title": "Common Functions",
    "category": "page",
    "text": ""
},

{
    "location": "solver/euler/common.html#Common-Functions-1",
    "page": "Common Functions",
    "title": "Common Functions",
    "category": "section",
    "text": ""
},

{
    "location": "solver/euler/conversion.html#",
    "page": "Conversion",
    "title": "Conversion",
    "category": "page",
    "text": ""
},

{
    "location": "solver/euler/conversion.html#Conversion-Between-Different-Variables-1",
    "page": "Conversion",
    "title": "Conversion Between Different Variables",
    "category": "section",
    "text": ""
},

{
    "location": "solver/euler/flux_functions.html#",
    "page": "Numerical Flux Functions",
    "title": "Numerical Flux Functions",
    "category": "page",
    "text": ""
},

{
    "location": "solver/euler/flux_functions.html#Numerical-Flux-Functions-1",
    "page": "Numerical Flux Functions",
    "title": "Numerical Flux Functions",
    "category": "section",
    "text": "bc_solvers.jl should be renamed to this"
},

{
    "location": "solver/euler/stabilization.html#",
    "page": "Stabilization",
    "title": "Stabilization",
    "category": "page",
    "text": ""
},

{
    "location": "solver/euler/stabilization.html#Stabilization-Documentation-1",
    "page": "Stabilization",
    "title": "Stabilization Documentation",
    "category": "section",
    "text": ""
},

{
    "location": "solver/euler/adjoint.html#",
    "page": "Adjoint",
    "title": "Adjoint",
    "category": "page",
    "text": ""
},

{
    "location": "solver/euler/adjoint.html#Adjoint-1",
    "page": "Adjoint",
    "title": "Adjoint",
    "category": "section",
    "text": ""
},

{
    "location": "solver/euler/boundary_functional.html#",
    "page": "Boundary Functions",
    "title": "Boundary Functions",
    "category": "page",
    "text": ""
},

{
    "location": "solver/euler/boundary_functional.html#Boundary-Functional-1",
    "page": "Boundary Functions",
    "title": "Boundary Functional",
    "category": "section",
    "text": ""
},

{
    "location": "solver/euler/misc.html#",
    "page": "Misc",
    "title": "Misc",
    "category": "page",
    "text": ""
},

{
    "location": "solver/euler/misc.html#Miscellaneous-Function-1",
    "page": "Misc",
    "title": "Miscellaneous Function",
    "category": "section",
    "text": "A bunch of the things in euler_funcs.jl"
},

{
    "location": "solver/simpleODE/simpleODE.html#",
    "page": "Main",
    "title": "Main",
    "category": "page",
    "text": ""
},

{
    "location": "solver/simpleODE/simpleODE.html#Simple-ODE-Documentation-1",
    "page": "Main",
    "title": "Simple ODE Documentation",
    "category": "section",
    "text": ""
},

{
    "location": "input/input.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "input/input.html#Input-1",
    "page": "Introduction",
    "title": "Input",
    "category": "section",
    "text": "A PDESolver execution is controlled by an options dictonary.  The user must supply this dictonary in the form of a Julia source file that declares a Dict{Any, Any} called arg_dict. The file path (relative to the users pwd must be passed to the solver as the the first optional argument.  For example, to launch a PDESolver run, executejulia /path/to/startup.jl \"path/to/dictonary\"Descriptions of the valid keys and values can be found in the file input_vals.txt. PDESolver supplies default values for all keys for which there is a sane default.  Values that might be unused are initialized to -1 or \"none\"."
},

{
    "location": "input/input.html#Conventions-1",
    "page": "Introduction",
    "title": "Conventions",
    "category": "section",
    "text": "Physics modules generally use the majorIterationCallback function to log important quantities to files.  Such logging should be controlled by two keys, \"write_outname\" where outname is the name of the quantity, which has a boolean value, and \"write_outname_fname\" that has a string value containing the name of the file to write (including extension).  Examples out things that can be logged are entropy and kinetic energy.  Both these keys should have default values, and users should generally not need to modify the second one."
},

{
    "location": "input/input.html#Input.make_input-Tuple{Dict{K,V},AbstractString}",
    "page": "Introduction",
    "title": "Input.make_input",
    "category": "Method",
    "text": "Make an input file from an options dictionary.\n\nInputs:     dict: the options dictionary     fname: the file name (without extension)\n\n\n\n"
},

{
    "location": "input/input.html#Input.read_input-Tuple{AbstractString}",
    "page": "Introduction",
    "title": "Input.read_input",
    "category": "Method",
    "text": "PDESolver.read_input\n\nThis function reads a file which must  be a julia source file that declares   a dictionary of option keywords and values for the options named arg_dict.   See the documention on input variables for valid keywords.\n\nread_input() returns the dictionary after doing some sanity checks and   supplying default values for any unspecified keys.\n\nAfter supplying default values, it prints the dictonary to arg_dict_output.jl,   which is a valid julia source file and can be read in to re-run a simulation.\n\nThis function checks whether the keys in arg_dict are recognized keywords   and prints a warning to STDERR if an unrecognized key is found.  The list of   known keys is read from the julia source file known_keys.jl\n\nInputs:     * fname : name of file to read\n\nOutputs:     arg_dict: a Dict{Any, Any} containing the option keywords and values\n\n\n\n"
},

{
    "location": "input/input.html#Input.checkForIllegalOptions_post-Tuple{Any}",
    "page": "Introduction",
    "title": "Input.checkForIllegalOptions_post",
    "category": "Method",
    "text": "Check the user supplied options for errors after supplying default options.\n\n\n\n"
},

{
    "location": "input/input.html#Input.checkForIllegalOptions_pre-Tuple{Any}",
    "page": "Introduction",
    "title": "Input.checkForIllegalOptions_pre",
    "category": "Method",
    "text": "Check the user supplied options for errors before supplying default values\n\n\n\n"
},

{
    "location": "input/input.html#Input.checkKeys-Tuple{Any,Any}",
    "page": "Introduction",
    "title": "Input.checkKeys",
    "category": "Method",
    "text": "PDESolver.checkKeys\n\nThis function verifies all the keys in the first argument are also keys   of the second argument and prints a warning to STDERR if they are not.\n\nInputs     arg_dict: first dictonary     known_keys: second dictonary\n\nOutputs:     cnt: number of unrecognized keys\n\n\n\n"
},

{
    "location": "input/input.html#Key-Validation-1",
    "page": "Introduction",
    "title": "Key Validation",
    "category": "section",
    "text": "After supplying default values, PDESolver checks that all keys in the dictonary are recognized keys.  It does this by comparing against the list of keys documented in the input_vals.txt file.  A warning of is printed to STDERR if an unrecognized key is found.The mechanics of the key validation are as follows.  The extract_keys.jl script reads the input_vals.txt file (and the input_vals_internal.txt), looking for any string that is surrounded by double quotes and starts in the first character of a line.  All keys are written to a dictonary in the file known_keys.jl.  This file is included by PDESolver.The script extract_keys.jl is run as part of PDESolver installation.  If new keys are added it must be run again to silence warnings during key validation.  Modules = [Input]\n  Pages = [\"input/extract_keys.jl\",\n           \"input/Input.jl\",\n           \"input/make_input.jl\",\n           \"input/read_input.jl\"]"
},

{
    "location": "NonlinearSolvers/nonlinearsolvers.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "NonlinearSolvers/nonlinearsolvers.html#NonlinearSolvers-Documentation-1",
    "page": "Introduction",
    "title": "NonlinearSolvers Documentation",
    "category": "section",
    "text": ""
},

{
    "location": "NonlinearSolvers/steady.html#",
    "page": "Steady",
    "title": "Steady",
    "category": "page",
    "text": ""
},

{
    "location": "NonlinearSolvers/steady.html#Steady-NonlinearSolver-Documentation-1",
    "page": "Steady",
    "title": "Steady NonlinearSolver Documentation",
    "category": "section",
    "text": ""
},

{
    "location": "NonlinearSolvers/unsteady/intro.html#",
    "page": "Intro",
    "title": "Intro",
    "category": "page",
    "text": ""
},

{
    "location": "NonlinearSolvers/unsteady/intro.html#Unsteady-NonlinearSolver-Basics-1",
    "page": "Intro",
    "title": "Unsteady NonlinearSolver Basics",
    "category": "section",
    "text": "PDESolver contains multiple time-stepping methods. They are selectable with the \"run_type\" input option."
},

{
    "location": "NonlinearSolvers/unsteady/rk4.html#",
    "page": "Runge-Kutta",
    "title": "Runge-Kutta",
    "category": "page",
    "text": ""
},

{
    "location": "NonlinearSolvers/unsteady/rk4.html#Fourth-Order-Runge-Kutta-(RK4)-1",
    "page": "Runge-Kutta",
    "title": "Fourth Order Runge-Kutta (RK4)",
    "category": "section",
    "text": "RK4 is a very widely used explicit time stepping method. See the Wikipedia page for the mathematical details."
},

{
    "location": "NonlinearSolvers/unsteady/lserk.html#",
    "page": "LSERK",
    "title": "LSERK",
    "category": "page",
    "text": ""
},

{
    "location": "NonlinearSolvers/unsteady/lserk.html#Low-Storage-Explicit-Runge-Kutta-(LSERK)-1",
    "page": "LSERK",
    "title": "Low Storage Explicit Runge-Kutta (LSERK)",
    "category": "section",
    "text": ""
},

{
    "location": "NonlinearSolvers/unsteady/cn.html#",
    "page": "Crank-Nicolson",
    "title": "Crank-Nicolson",
    "category": "page",
    "text": ""
},

{
    "location": "NonlinearSolvers/unsteady/cn.html#Crank-Nicolson-(CN)-1",
    "page": "Crank-Nicolson",
    "title": "Crank-Nicolson (CN)",
    "category": "section",
    "text": "Crank-Nicolson is an implicit time marching scheme. While it adds a measure of complexity, it is unconditionally stable.  Note that this does not mean that accuracy is unaffected by time step size,  degree of the operator (p = 1 or p = 2, etc.), or element size."
},

{
    "location": "NonlinearSolvers/unsteady/cn.html#Selecting-CN-and-Problem-Setup-1",
    "page": "Crank-Nicolson",
    "title": "Selecting CN & Problem Setup",
    "category": "section",
    "text": "CN is selected with \"run_type\" => 20 in your problem's input options file. As this is an unsteady problem, you will also need to specify your    time step size (\"delta_t\" option) and your final solution time (\"t_max\" option)"
},

{
    "location": "NonlinearSolvers/unsteady/cn.html#Forward-in-time-1",
    "page": "Crank-Nicolson",
    "title": "Forward-in-time",
    "category": "section",
    "text": "\\begin{aligned} \\underbrace{\\mathbf{I} - \\frac{\\Delta t}{2} \\underbrace{\\frac{\\partial R\\left(u^{i+1}\\right)}{\\partial u^{i+1}}}_{\\text{jac from \\texttt{physicsJac}}} \n}_{\\text{\\texttt{cnJac}}} \n= \\\n\\underbrace{- \\left( \n  \\underbrace{u^{i+1}}_{\\text{\\texttt{eqn_nextstep.q_vec}}} - \\      \\frac{\\Delta t}{2} \\underbrace{R\\left(u^{i+1}\\right)}_{\\text{\\texttt{eqn_nextstep.res_vec}}} - \\    \\underbrace{u^i}_{\\text{\\texttt{eqn.q_vec}}} - \\      \\frac{\\Delta t}{2} \\underbrace{R\\left(u^i\\right)}_{\\text{\\texttt{eqn.res_vec}}} \n\\right)}_{\\text{\\texttt{cnRhs}}} \\end{aligned}This equation is solved with PDESolver's Newton's method."
},

{
    "location": "NonlinearSolvers/unsteady/cn.html#Backward-in-time-(unsteady-adjoint)-1",
    "page": "Crank-Nicolson",
    "title": "Backward-in-time (unsteady adjoint)",
    "category": "section",
    "text": ""
},

{
    "location": "Utils/Utils.html#",
    "page": "Main",
    "title": "Main",
    "category": "page",
    "text": ""
},

{
    "location": "Utils/Utils.html#Utilties-1",
    "page": "Main",
    "title": "Utilties",
    "category": "section",
    "text": "This module contains functions and types that are useful for the solver but independent of the equation being solved. Additional utility functions are located in the ODLCommontools. The functions defined in the Utils module are useful in the context of PDESolver and depend on the functions and datatypes defined in the other parts of the solver. The functions defined in ODLCommonTools are more general in nature and usable independent of PDESolver.  Pages = [\"io.md\"\n           \"logging.md\"\n           \"projections.md\"\n           \"parallel.md\"\n          ]\n  Depth = 1"
},

{
    "location": "Utils/parallel.html#",
    "page": "Parallel Constructs",
    "title": "Parallel Constructs",
    "category": "page",
    "text": ""
},

{
    "location": "Utils/parallel.html#Parallel-Constructs-Documentations-1",
    "page": "Parallel Constructs",
    "title": "Parallel Constructs Documentations",
    "category": "section",
    "text": "  CurrentModule = UtilsThese function define the primative operations used by the physics modules to exchange data in parallel. When using these functions, the should not have to make any MPI calls directly, they should all be encapsulated within the provided functions.TODO: crossref to physics module documentationThe Types and Basic API section describes the SharedFaceData datatype and the basic functions that operate on it. The Parallel Data Exchange section describes the functions used by the physics modules that start and finish parallel parallel communication."
},

{
    "location": "Utils/parallel.html#Utils.SharedFaceData",
    "page": "Parallel Constructs",
    "title": "Utils.SharedFaceData",
    "category": "Type",
    "text": "This type holds all the data necessary to perform MPI communication with   a given peer process that shared mesh edges (2D) or faces (3D) with the   current process.\n\nFields:\n\npeernum: the MPI rank of the peer process\npeeridx: the index of this peer in mesh.peer_parts\nmyrank: MPI rank of the current process\ncomm: MPI communicator used to define the above\n\nq_send: the send buffer, a 3D array of n x m x d.  While these dimensions\n        are arbitrary, there are two commonly used case.  If\n        opts[\"parallel_type\"] == face, then m is mesh.numNodesPerFace and\n        d is the number of faces shared with peernum.\n        If opts[\"parallel_type\"] == element, then \n        m = mesh.numNodesPerElement and d is the number of elements that\n        share faces with peernum.\nq_recv: the receive buffer.  Similar to q_send, except the size needs to\n        to be the number of entities on the *remote* process.\n\nsend_waited: has someone called MPI.Wait() on send_req yet?  Some MPI\n             implementations complain if Wait() is called on a Request\n             more than once, so use this field to avoid doing so.\nsend_req: the MPI.Request object for the Send/Isend/whatever other type of\n          Send\nsend_status: the MPI.Status object returned by calling Wait() on send_req\n\nrecv_waited: like send_waited, but for the receive\nrecv_req: like send_req, but for the receive\nrecv_status: like send_status, but for the receive\n\nbndries_local: Vector of Boundaries describing the faces from the local\n               side of the interface\nbndries_remote: Vector of Boundaries describing the facaes from the remote\n                side (see the documentation for PdePumiInterface before\n                using this field)\ninterfaces: Vector of Interfaces describing the faces from both sides (see\n            the documentation for PdePumiInterfaces, particularly the\n            mesh.shared_interfaces field, before using this field\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Utils.SharedFaceData-Tuple{ODLCommonTools.AbstractMesh{Tmsh},Int64,Array{T,3},Array{T,3}}",
    "page": "Parallel Constructs",
    "title": "Utils.SharedFaceData",
    "category": "Method",
    "text": "Outer constructor for SharedFaceData.\n\nInputs:\n\nmesh: a mesh object\npeeridx: the index of a peer in mesh.peer_parts\nq_send: the send buffer\nq_recv: the receive buffer\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Utils.getSharedFaceData-Tuple{Type{Tsol},ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP{T<:Number},Any}",
    "page": "Parallel Constructs",
    "title": "Utils.getSharedFaceData",
    "category": "Method",
    "text": "This function returns a vector of SharedFaceData objects, one for each   peer processes the current process shares mesh edge (2D) or face (3D) with.   This function is intended to be used by the AbstractSolutionData constructors,   although it can be used to create additional vectors of SharedFaceData   objects.\n\nif opts[\"parallel_data\"] == \"face\", then the send and receive buffers   are numDofPerNode x numNodesPerFace x number of shared faces.\n\nif opts[\"parallel_data\"] == \"element\", the send and receive buffers are     numDofPerNode x numNodesPerElement x number of elements that share the     faces.  Note that the number of elements that share the faces can be     different for the send and receive buffers.\n\nInputs:\n\nTsol: element type of the arrays\nmesh: an AbstractMesh object\nsbp: an SBP operator\nopts: the options dictonary\n\nOutputs:\n\ndata_vec: Vector{SharedFaceData}.  data_vec[i] corresponds to \n          mesh.peer_parts[i]\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Base.copy!-Tuple{Utils.SharedFaceData{T},Utils.SharedFaceData{T}}",
    "page": "Parallel Constructs",
    "title": "Base.copy!",
    "category": "Method",
    "text": "In-place copy for SharedFaceData.  This copies the buffers, but does not   retain the state of the Request and Status fields.  Instead they are   initialized the same as the constructor.\n\nThis function may only be called after receiving is complete,   otherwise an exception is thrown.\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Base.copy-Tuple{Utils.SharedFaceData{T}}",
    "page": "Parallel Constructs",
    "title": "Base.copy",
    "category": "Method",
    "text": "Copy function for SharedFaceData.  Note that this does not retain the   send_req/send_status (and similarly for the recceive) state   of the original object.  Instead, they are initialized the same as the   constructor.\n\nThis function may only be called after receiving is complete,   otherwise an exception is thrown.\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Utils.assertReceivesConsistent-Tuple{Array{Utils.SharedFaceData{T},1}}",
    "page": "Parallel Constructs",
    "title": "Utils.assertReceivesConsistent",
    "category": "Method",
    "text": "Like assertSendsConsistent, but for the receives\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Utils.assertReceivesWaited-Tuple{Array{Utils.SharedFaceData{T},1}}",
    "page": "Parallel Constructs",
    "title": "Utils.assertReceivesWaited",
    "category": "Method",
    "text": "This function verifies all the receives have been waited on for the    supplied SharedFaceData objects\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Utils.assertSendsConsistent-Tuple{Array{Utils.SharedFaceData{T},1}}",
    "page": "Parallel Constructs",
    "title": "Utils.assertSendsConsistent",
    "category": "Method",
    "text": "Verify either all or none of the sends have been waited on.  Throw an   exception otherwise.\n\nInputs:\n\nshared_data: Vector of SharedFaceData objects\n\nOutput:\n\nval: number of receives that have been waited on\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Utils.waitAllReceives-Tuple{Array{Utils.SharedFaceData{T},1}}",
    "page": "Parallel Constructs",
    "title": "Utils.waitAllReceives",
    "category": "Method",
    "text": "This function is like MPI.Waitall, operating on the recvs of a vector of    SharedFaceData objects\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Utils.waitAllSends-Tuple{Array{Utils.SharedFaceData{T},1}}",
    "page": "Parallel Constructs",
    "title": "Utils.waitAllSends",
    "category": "Method",
    "text": "This function is like MPI.Waitall, operating on the sends of a vector of    SharedFaceData objects\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Utils.waitAnyReceive-Tuple{Array{Utils.SharedFaceData{T},1}}",
    "page": "Parallel Constructs",
    "title": "Utils.waitAnyReceive",
    "category": "Method",
    "text": "Like MPI.WaitAny, but operates on the receives of  a vector of SharedFaceData.   Only the index of the Request that was waited on is returned,    the Status and recv_waited fields of hte SharedFaceData are updated internally\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Types-and-Basic-API-1",
    "page": "Parallel Constructs",
    "title": "Types and Basic API",
    "category": "section",
    "text": "  Modules = [Utils]\n  Pages = [\"Utils/parallel_types.jl\"]"
},

{
    "location": "Utils/parallel.html#Utils.startSolutionExchange",
    "page": "Parallel Constructs",
    "title": "Utils.startSolutionExchange",
    "category": "Function",
    "text": "This function is a thin wrapper around exchangeData().  It is used for the   common case of sending and receiving the solution variables to other processes.   It uses eqn.shared_data to do the parallel communication.   eqn.shared_data must be passed into the corresponding finishDataExchange   call.\n\nInputs:     mesh: an AbstractMesh     sbp: an SBP operator     eqn: an AbstractSolutionData     opts: options dictionary\n\nKeyword arguments:     tag: MPI tag to use for communication, defaults to TAG_DEFAULT     wait: wait for sends and receives to finish before exiting\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Utils.exchangeData",
    "page": "Parallel Constructs",
    "title": "Utils.exchangeData",
    "category": "Function",
    "text": "This function posts the MPI sends and receives for a vector of SharedFaceData.  It works for both opts[\"parallel_data\"] == \"face\" or \"element\".  The only   difference between these two cases is the populate_buffer() function.\n\nThe previous receives using these SharedFaceData objects should have   completed by the time this function is called.  An exception is throw   if this is not the case.\n\nThe previous sends are likely to have completed by the time this function   is called, but they are waited on if not.  This function might not perform   well if the previous sends have not completed.   #TODO: fix this using WaitAny\n\nInputs:     mesh: an AbstractMesh     sbp: an SBPOperator     eqn: an AbstractSolutionData     opts: the options dictionary     populate_buffer: function with the signature:                      populate_buffer(mesh, sbp, eqn, opts, data::SharedFaceData)                      that populates data.q_send   Inputs/Outputs:     shared_data: vector of SharedFaceData objects representing the parallel                  communication to be done\n\nKeyword Arguments:     tag: MPI tag to use for this communication, defaults to TAG_DEFAULT          This tag is typically used by the communication of the solution          variables to other processes.  Other users of this function should          provide their own tag\n\nwait: wait for the sends and receives to finish before returning.  This\n      is a debugging option only.  It will kill parallel performance.\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Utils.finishExchangeData",
    "page": "Parallel Constructs",
    "title": "Utils.finishExchangeData",
    "category": "Function",
    "text": "This is the counterpart of exchangeData.  This function finishes the   receives started in exchangeData.\n\nThis function (efficiently) waits for a receive to finish and calls   a function to do calculations for on that data. If opts[\"parallel_data\"]   == \"face\", it also permutes the data in the receive buffers to agree   with the ordering of elementL.  For opts[\"parallel_data\"] == \"element\",   users should call SummationByParts.interiorFaceInterpolate to interpolate   the data to the face while ensuring proper permutation.\n\nInputs:     mesh: an AbstractMesh     sbp: an SBPOperator     eqn: an AbstractSolutionData     opts: the options dictonary     calc_func: function that does calculations for a set of shared faces                described by a single SharedFaceData.  It must have the signature                calc_func(mesh, sbp, eqn, opts, data::SharedFaceData)\n\nInputs/Outputs:     shared_data: vector of SharedFaceData, one for each peer process that                  needs to be communicated with.  By the time calc_func is                  called, the SharedFaceData passed to it has its q_recv field                  populated.  See note above about data permutation.\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Utils.TAG_DEFAULT",
    "page": "Parallel Constructs",
    "title": "Utils.TAG_DEFAULT",
    "category": "Constant",
    "text": "Default MPI tag used for sending and receiving solution variables.\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Parallel-Data-Exchange-1",
    "page": "Parallel Constructs",
    "title": "Parallel Data Exchange",
    "category": "section",
    "text": "The functions in this section are used to start sending data in parallel and finish receiving it. All functions operate on a Vector of SharedFaceData that define what data to send to which peer processes.  See Parallel Overview for a high-level overview of how the code is parallelized.Sending the data to the other processes is straight-forward.  Receiving it (efficiently) is not. In particular, [finishExchangeData] waits to receive data from one peer process, calls a user supplied callback function to do calculations involving the received data, and then waits for the next receive to finish. This is significantly more efficient than waiting for all receives to finish and then doing computations on all the data.This section describes the API the physics modules use to do parallel  communication.  The Internals section describes the helper functions used in the implementation.startSolutionExchange\nexchangeData\nfinishExchangeData\nTAG_DEFAULT"
},

{
    "location": "Utils/parallel.html#Utils.verifyReceiveCommunication",
    "page": "Parallel Constructs",
    "title": "Utils.verifyReceiveCommunication",
    "category": "Function",
    "text": "Utils.verifyCommunication\n\nThis function checks the data provided by the Status object to verify a    communication completed successfully.  The sender's rank and the number of   elements is checked agains the expected sender and the buffer size\n\nInputs:     data: a SharedFaceData\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Utils.getSendDataFace",
    "page": "Parallel Constructs",
    "title": "Utils.getSendDataFace",
    "category": "Function",
    "text": "This function populates the send buffer from eqn.q for    opts[\"parallle_data\"]  == \"face\"\n\nInputs:     mesh: a mesh     sbp: an SBP operator     eqn: an AbstractSolutionData     opts: options dictonary\n\nInputs/Outputs:     data: a SharedFaceData.  data.q_send will be overwritten\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Utils.getSendDataElement",
    "page": "Parallel Constructs",
    "title": "Utils.getSendDataElement",
    "category": "Function",
    "text": "This function populates the send buffer from eqn.q for    opts[\"parallle_data\"]  == \"element\"\n\nInputs:\n\nmesh: a mesh\nsbp: an SBP operator\neqn: an AbstractSolutionData\nopts: options dictonary\n\nInputs/Outputs:\n\ndata: a SharedFaceData.  data.q_send will be overwritten\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Utils.@mpi_master",
    "page": "Parallel Constructs",
    "title": "Utils.@mpi_master",
    "category": "Macro",
    "text": "Utils.mpi_master\n\nThis macro introduces an if statement that causes the expression to be    executed only if the variable myrank is equal to zero.  myrank must exist   in the scope of the caller\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Utils.@time_all",
    "page": "Parallel Constructs",
    "title": "Utils.@time_all",
    "category": "Macro",
    "text": "Utils.time_all\n\nThis macro returns the value produced by the expression as well as    the execution time, the GC time, and the amount of memory allocated\n\n\n\n"
},

{
    "location": "Utils/parallel.html#utils_parallel_internals-1",
    "page": "Parallel Constructs",
    "title": "Internals",
    "category": "section",
    "text": "These helper functions are used by the functions in Parallel Data Exchange.verifyReceiveCommunication\ngetSendDataFace\ngetSendDataElement\n@mpi_master\n@time_all"
},

{
    "location": "Utils/projections.html#",
    "page": "Projections",
    "title": "Projections",
    "category": "page",
    "text": ""
},

{
    "location": "Utils/projections.html#Projections-documentation-1",
    "page": "Projections",
    "title": "Projections documentation",
    "category": "section",
    "text": "  CurrentModule = UtilsThe functions here project the vector of conservative variables back and  forth between x-y-z and n-t-b (normal-tangential-binormal) cordinates."
},

{
    "location": "Utils/projections.html#Utils.calcLength-Tuple{ODLCommonTools.AbstractParamType{2},AbstractArray{T,1}}",
    "page": "Projections",
    "title": "Utils.calcLength",
    "category": "Method",
    "text": "Calculates the length of a vector, using a manually unrolled loop.   Methods are available for 2D and 3D.\n\nInputs:\n\nparams: an AbstractParamType{Tdim}, used to dispatch to the right method\nnrm: the vector to calculate the length of.  Must have length 2 in 2D\n     and 3 in 3D\n\nOutputs:\n\nlength of nrm\n\n\n\n"
},

{
    "location": "Utils/projections.html#Utils.getProjectionMatrix-Tuple{ODLCommonTools.AbstractParamType{2},AbstractArray{T,1},AbstractArray{T,2}}",
    "page": "Projections",
    "title": "Utils.getProjectionMatrix",
    "category": "Method",
    "text": "This function populates the matrix P that project from x-y-z to n-t-b.   Methods are available for 2D and 3D, determined from the AbstractParamType   object.  This is a somewhat specialized routine in that it only works for   vectors like the state vector of the Euler or Navier-Stokes equations,   where the first and last equations are coordinate system invarient and   only the middle 2 (in 2D) or 3 (in 3D) equations need to be rotated.   Specialized multiplication routines are provied as projectToXY and   projectToNT.\n\nInputs:\n\nparams:  An AbstractParamType{Tdim}\nnrm: the normal direction in x-y coordinate system.  Must be a unit vector\n\nInputs/Outputs\n\nP:  matrix to be populated with projection matrix\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "Utils/projections.html#Utils.projectToNT-Tuple{ODLCommonTools.AbstractParamType{2},AbstractArray{T,2},AbstractArray{T,1},AbstractArray{T,1}}",
    "page": "Projections",
    "title": "Utils.projectToNT",
    "category": "Method",
    "text": "This function projects a vector x from x-y coordinates to normal-tangential   (n-t) coordinates, where A is a projection matrix from x-y to n-t, obtained   from getProjectionMatrix.  This function is a specialized matrix-vector   multiplication routine and only works for matrices with the particular   sparsity pattern created by getProjectionMatrix.   Methods are available for 2D and 3D\n\nInputs:\n\nparams: an AbstractParamType{{Tdim} used to dispatch to the 2D or 3D\n        method\nP: the projection matrix\nx: vector to be projected\n\nInputs/Outputs:\n\nb: the result of the projection\n\nAliasing restrictions: x and b cannot alias\n\n\n\n"
},

{
    "location": "Utils/projections.html#Utils.projectToXY-Tuple{ODLCommonTools.AbstractParamType{2},AbstractArray{T,2},AbstractArray{T,1},AbstractArray{T,1}}",
    "page": "Projections",
    "title": "Utils.projectToXY",
    "category": "Method",
    "text": "This function is similar to projectToNT, except it project from n-t   to x-y.  Note that P is still the projection matrix from getProjectionMatrix   that projects from x-y to n-t.\n\n\n\n"
},

{
    "location": "Utils/projections.html#Utils.getBinormalVector-Tuple{Number,Number,Number,Number,Number,Number}",
    "page": "Projections",
    "title": "Utils.getBinormalVector",
    "category": "Method",
    "text": "This function computes a vector normal to the 2 supplied vectors and   returns the components.\n\nInputs:\n\nn1, n2, n3: the components of the first vector\nt1, t2, t3: the components of the second vector\n\nOutputs:\n\nb1, b2, b3: the components of the binormal vector\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "Utils/projections.html#Utils.getOrthogonalVector-Tuple{ODLCommonTools.AbstractParamType{2},AbstractArray{T,1}}",
    "page": "Projections",
    "title": "Utils.getOrthogonalVector",
    "category": "Method",
    "text": "This function generates a unit vector orthogonal to the input vector nrm and    returns its components.  Methods are availble for 2D and 3D\n\nInputs:\n\nparams: an AbstractParamType, used to dispatch to the 2D or 3D method\nnrm: the input vector\n\nOutputs:\n\nt1, t2, (and t3 in 3D): the components of the unit vector orthogonal to\n                        nrm\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "Utils/projections.html#Detailed-Documentation-1",
    "page": "Projections",
    "title": "Detailed Documentation",
    "category": "section",
    "text": "  Modules = [Utils]\n  Pages = [\"Utils/projections.jl\"]"
},

{
    "location": "Utils/logging.html#",
    "page": "Logging",
    "title": "Logging",
    "category": "page",
    "text": ""
},

{
    "location": "Utils/logging.html#Logging-Documentation-1",
    "page": "Logging",
    "title": "Logging Documentation",
    "category": "section",
    "text": "  CurrentModule = UtilsThe functions define here assist with debugging by logging some results to disk."
},

{
    "location": "Utils/logging.html#Utils.sharedFaceLogging-Tuple{Any,Any,ODLCommonTools.AbstractSolutionData{Tsol,Tres},Any,Utils.SharedFaceData{T},Any,Any}",
    "page": "Logging",
    "title": "Utils.sharedFaceLogging",
    "category": "Method",
    "text": "Utils.sharedFaceLogging\n\nThis function writes files for the function that evalutes the shared face   integrals.  This function should only be called in debug mode.\n\nIf opts[\"writeqface\"] is true,  writes the q values at the face to a file    q_sharedface_i_myrank.dat, where i is the peer process index (not MPI rank)   and myrank is the MPI rank of the current process. \n\nIf opts[\"write_fluxface\"] is true, it writes eqn.flux_sharedface to a file   flux_sharedface_i_myrank.dat\n\nInputs:     mesh     sbp     eqn: and AbstractSolutionData     opts: options dictonary     data: SharedFaceData object     qL_arr: an array holding solution values at face nodes on the              left side of the interface,  numDofPerNode x numfacenodes x              numsharedfaces on this partition boundary     qR_arr: solution values at face nodes on the right side of the interface.             same shape as qL_arr\n\nalso, the eqn.flux_shared is used to write the face flux.  It is the same   shape as qL_arr\n\nAliasing restrictions: qL_arr, qR_arr, and eqn.flux_sharedface must not alias.\n\n\n\n"
},

{
    "location": "Utils/logging.html#Detailed-Documentation-1",
    "page": "Logging",
    "title": "Detailed Documentation",
    "category": "section",
    "text": "  Modules = [Utils]\n  Pages = [\"Utils/logging.jl\"]"
},

{
    "location": "Utils/io.html#",
    "page": "Input/Output",
    "title": "Input/Output",
    "category": "page",
    "text": ""
},

{
    "location": "Utils/io.html#Input/Output-1",
    "page": "Input/Output",
    "title": "Input/Output",
    "category": "section",
    "text": "  CurrentModule = UtilsThe functions and types defined here facilitate writing to STDOUT and STDERR as well as files. In particular, performance tests have shown that buffering is important when running in parallel. To this end, the BufferedIO type is introduced that stores output in an in-memory buffer before writing to an underlying stream. This can create some difficulty when debuging, because output written to a buffered stream is not immediately written to the underlying stream. In cases where precise control of output is needed, users should call flush to make sure all output is written to the underlying stream"
},

{
    "location": "Utils/io.html#Utils.BSTDERR",
    "page": "Input/Output",
    "title": "Utils.BSTDERR",
    "category": "Constant",
    "text": "Buffered version of STDERR.  This should always be used instead of STDERR\n\n\n\n"
},

{
    "location": "Utils/io.html#Utils.BSTDOUT",
    "page": "Input/Output",
    "title": "Utils.BSTDOUT",
    "category": "Constant",
    "text": "Buffered version of STDOUT.  This should always be used instead of STDOUT\n\n\n\n"
},

{
    "location": "Utils/io.html#Utils.BufferedIO",
    "page": "Input/Output",
    "title": "Utils.BufferedIO",
    "category": "Type",
    "text": "Utils.BufferedIO\n\nThis type provides a means to buffer IO in memory before writing it to a file.   Data written to the object is stored in the IOBuffer until flush() is called,    when the buffer is (efficiently) dumped into the file\n\n\n\n"
},

{
    "location": "Utils/io.html#Utils.BufferedIO",
    "page": "Input/Output",
    "title": "Utils.BufferedIO",
    "category": "Type",
    "text": "Utils.BufferedIO\n\nConstructor for BufferedIO type.  Takes an IOStream and creates an IOBuffer.   If no IOStream is given, a dummy stream is created.  If a dummy stream is   used, flush must never be called on the returned BufferedIO object\n\nInputs:\n\nf: an IOStream object, defaults to a dummy stream\n\nOutputs:\n\na BufferedIO object\n\n\n\n"
},

{
    "location": "Utils/io.html#Utils.BufferedIO",
    "page": "Input/Output",
    "title": "Utils.BufferedIO",
    "category": "Type",
    "text": "Alternative constructor for BufferedIO, emulating the open() function.   This function creates the underlying file using open() and then creates   a BufferedIO around it.\n\nInputs:\n\nfname: AbstractString, name of file to open\nmode: file open mode, see documentation of open(), defaults to append\n\nOutputs:\n\na BufferedIO object\n\n\n\n"
},

{
    "location": "Utils/io.html#Base.close-Tuple{Utils.BufferedIO{T<:IO}}",
    "page": "Input/Output",
    "title": "Base.close",
    "category": "Method",
    "text": "Base function close extended for BufferedIO\n\n\n\n"
},

{
    "location": "Utils/io.html#Base.flush-Tuple{Utils.BufferedIO{T<:IO}}",
    "page": "Input/Output",
    "title": "Base.flush",
    "category": "Method",
    "text": "Base function flush extended for BufferedIO\n\n\n\n"
},

{
    "location": "Utils/io.html#Detailed-Documentation-1",
    "page": "Input/Output",
    "title": "Detailed Documentation",
    "category": "section",
    "text": "  Modules = [Utils]\n  Pages = [\"Utils/io.jl\"]"
},

]}
