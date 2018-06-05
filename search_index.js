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
    "text": "Summation-by-Parts (SBP) operators are used to discretize the spatial derivatives in the residual. In particular, we use the multi-dimensional Summation-by-Parts operators first defined in   Hicken, J.E., Del Rey Fernandez, D.C, and Zingg, D.W., \"Multi-dimensional \n  Summation-by-Parts Operators: General Theory and Application to \n  Simplex Elements\", SIAM Journal on Scientific Computing, Vol, 38, No. 4, 2016,\n   pp. A1935-A1958See the introductory PDF (to be posted shortly) for an introduction to the operators."
},

{
    "location": "index.html#Table-of-Contents-1",
    "page": "PDESolver Introduction",
    "title": "Table of Contents",
    "category": "section",
    "text": "Pages = [\"invocation/calling.md\",\n         \"solver/Readme.md\",\n         \"solver/advection/advection.md\",\n         \"solver/euler/euler.md\",\n         \"solver/simpleODE/simpleODE.md\",\n         \"NonlinearSolvers/nonlinearsolvers.md\",\n         \"input/input.md\",\n         \"Utils/Utils.md\"\n        ]\nDepth=1"
},

{
    "location": "concepts/intro.html#",
    "page": "Intro",
    "title": "Intro",
    "category": "page",
    "text": ""
},

{
    "location": "concepts/intro.html#Intro-1",
    "page": "Intro",
    "title": "Intro",
    "category": "section",
    "text": ""
},

{
    "location": "concepts/pumi.html#",
    "page": "PUMI",
    "title": "PUMI",
    "category": "page",
    "text": ""
},

{
    "location": "concepts/pumi.html#Pumi-1",
    "page": "PUMI",
    "title": "Pumi",
    "category": "section",
    "text": ""
},

{
    "location": "concepts/sbp.html#",
    "page": "SBP",
    "title": "SBP",
    "category": "page",
    "text": ""
},

{
    "location": "concepts/sbp.html#SBP-1",
    "page": "SBP",
    "title": "SBP",
    "category": "section",
    "text": ""
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
    "text": "PDESolver is not listed in Julia\'s METADATA, so you will have to clone the repository and then build the code.After installing the system System Dependencies, run the following Julia commands to install the Julia dependencies and and PDESolver itself:  Pkg.clone(\"https://github.com/OptimalDesignLab/PDESolver.jl.git\")\n  Pkg.resolve()  # install all packages in PDESolvers REQUIRE file\n  Pkg.build(\"PDESolver\")  # install all packages not in REQUIRE file and\n                          # build PDESolverThis will install PDESolver and all Julia dependencies into the directory specified by the JULIA_PKGDIR environmental variable. If this variable is not specified, it will install to ~/.julia/v0.4/PDESolver.If there are errors building any Julia dependencies, see the Installation page for methods of installing particular versions of the dependencies.After installation it is recommended to run the test suite. To do so, run the following commands in the terminal:  # need to configure environment before using PumiInterface\n  cd /path/to/PumiInterface/src # for example ~/.julia/v0.4/PumiInterface/src)\n  source ./use_julialib.sh\n  \n  # run the PDESolver tests\n  cd /path/to/pdesolver # (for example ~/.julia/v0.4/PDESolver)\n  cd ./test\n  ./runtests_fast.shIf the tests complete without error, then the package is properly installed.TODO: link to examples page"
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
    "text": "PDESOLVER_INSTALL_DEPS_MANUAL: install the packages in the REQUIRE filemanually.  If the package is already installed, does nothing (does not check out the specified version).PDESOLVER_FORCE_DEP_INSTALL_ALL: forces the checkout and re-installation of the non-standard packages, even if already installedPDESOLVER_FORCE_DEP_INSTALL_pkg_name: force the checkout and re-installationof the package named pkg_name.For all these environmental variables, the value is not important, only the  existance of the variable."
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
    "text": "The file deps/install.log is created and written to with the progress of the checkout process, including error information if a build fails. Because Pkg.build swallows errors, failures inside Pkg.build are not logged."
},

{
    "location": "deps_readme.html#Known-Issues-1",
    "page": "Build Options",
    "title": "Known Issues",
    "category": "section",
    "text": ""
},

{
    "location": "deps_readme.html#Stack-Smashing-1",
    "page": "Build Options",
    "title": "Stack Smashing",
    "category": "section",
    "text": "On Intel Skylake, there is a problem with reusing the symbolic factorization when using the sparse direct linar solver.  The error message looks like:*** stack smashing detected ***: julia terminatedAn excerpt of the stack trace is:umfdl_local_search at /mnt/scratch/common/julia/0.4/bin/../lib/julia/libumfpack.so (unknown line)\numfdl_kernel at /mnt/scratch/common/julia/0.4/bin/../lib/julia/libumfpack.so (unknown line)\numfpack_dl_numeric at /mnt/scratch/common/julia/0.4/bin/../lib/julia/libumfpack.so (unknown line)\numfpack_numeric! at sparse/umfpack.jl:168\n_linearSolve at /home/creanj/pkgtest/v0.4/PDESolver/src/linearsolvers/ls_standard.jl:292\nlinearSolve at /home/creanj/pkgtest/v0.4/PDESolver/src/Utils/parallel.jl:312\nnewtonInner at util.jl:179To avoid this problem, see the SKYLAKE_STACKSMASH the src/linearsolver/ls_standard.jl file."
},

{
    "location": "deps_readme.html#Building-Julia-Offline-1",
    "page": "Build Options",
    "title": "Building Julia Offline",
    "category": "section",
    "text": "Building Julia on a system that cannot access the internet is a little bit different than build on a regular computer.  The approach isOn a computer with internet accessclone the Julia repository\ncheckout the version of Julia you want to build\nrun make -C deps getall to download the source tarballs for all the things Julia builds during its build process\ncreate a tarball of the julia directory\ncopy the julia tarball to the offline systemOn the system without internet accessextract the julia tarball\nadd the line prefix=/path/to/installation/directory to the file Make.user (create the file if it does not exist\nrun make (or parallel make with the -j switch\nmake the directories julia/doc/_build/html/en\ncreate the file julia/doc/_build/html/en/index.html\nrun make installThe purpose of creating the index.html file is that it will trick the Julia build system into not building the documentation.  This is necessary because build documentation requires downloading several Julia packages from the internet.After make install has finished, the julia directory can be deleted, which can save a lot of disk space on the filesystems used on clusers/supercomputers."
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
    "category": "type",
    "text": "This abstract type is the supertype for all the objects that store the    solution data. Every physics module should implement its own subtype.\n\nStatic parameters\n\nTsol: datatype of solution variables\nTres: datatype of the mesh variables\n\nSee the AbstractSolutionData for the description of everything this   type must implement.\n\n\n\n"
},

{
    "location": "interfaces.html#AbstractSolutionData-1",
    "page": "Code Interfaces",
    "title": "AbstractSolutionData",
    "category": "section",
    "text": "ODLCommonTools defines:AbstractSolutionDataThe purpose of an AbstractSolutionData is to hold all the data related to the solution of an equation. This includes the solution at every node and any auxiliary quantities. The storage for any quantity that is calculated over the entire mesh should be allocated as part of this object, in order to avoid repeatedly reallocated the array for every residual evaluation. In general, there should never be a need to allocate a vector longer than the number of degrees of freedom at a node (or a matrix similarly sized matrix) during a residual evaluation. Structuring code such that it conforms with this requirement has significant performance benefits because it reduces memory allocation/deallocation.The static parameter Tsol is the datatype of the solution variables and Tres is the datatype of the residual (when computing the Jacobian with finite differences or algorithmic differentiation, these will be the same)."
},

{
    "location": "interfaces.html#Required-Fields-1",
    "page": "Code Interfaces",
    "title": "Required Fields",
    "category": "section",
    "text": "The required fields of an AbstractSolutionData are:  q::AbstractArray{Tsol, 3}\n  q_vec::AbstractArray{Tsol, 1}\n  shared_data::AbstractArray{SharedFaceData{Tsol}, 1}\n  res::AbstractArray{Tres, 3}\n  res_vec::AbstractArray{Tres, 1}\n  res_edge::AbstractArray{Tres, 4}\n  M::AbstractArray{Float64, 1}\n  Minv::AbstractArray{Float64, 1}\n  multiplyA0inv::Function\n  majorIterationCallback::Function\n  params{Tsol..., Tdim}::AbstractParamType{Tdim}"
},

{
    "location": "interfaces.html#Optional-Fields-1",
    "page": "Code Interfaces",
    "title": "Optional Fields",
    "category": "section",
    "text": "  file_dict::String, IO}"
},

{
    "location": "interfaces.html#ODLCommonTools.AbstractParamType",
    "page": "Code Interfaces",
    "title": "ODLCommonTools.AbstractParamType",
    "category": "type",
    "text": "This abstract type is the supertype for all Param objects, which hold values    needed for the computation in a place that is fast to access.\n\nThe Param type is also useful for dispatching to low level functions which     the AbstractSolutionData might not be passed (depending on the organization     of the physics module.\n\nStatic Parameters:\n\nTdim: the dimensionality of the equation being solved (2d or 3d usually)\n\n\n\n"
},

{
    "location": "interfaces.html#Field-Meanings-1",
    "page": "Code Interfaces",
    "title": "Field Meanings",
    "category": "section",
    "text": "  CurrentModule = UtilsThe purpose of these fields are:q: to hold the solution variables in an element-based array.      This array should be numDofPerNode x numNodesPerElement x numEl.      The residual evaluation only uses q, never q_vecq_vec: to hold the solution variables as a vector, used for any linear algebra operations and time stepping. This array should have a length equal to the total number of degrees of freedom in the mesh. Even though this vector is not used by the residual evaluation, it is needed for many other operations, so it is allocated here so the memory can be reused. There are functions to facilitate the scattering of values from q_vec to q. Note that for Continuous Galerkin type discretization (as opposed to Discontinuous Galerkin discretizations), there is not a corresponding \"gather\" operation (ie. q -> q_vec).  See Utils.array1DTo3D  and Utils.array3DTo1D.shared_data is a vector of length npeers.  Each element contains the data               needed send and receive the q variables to/from other               the corresponding MPI rank listed in mesh.peer_parts.               The precise contents of SharedFaceData is documented in the               Utils module, however they do include the send and receive               buffers.res: similar to q, except that the residual evaluation function populates        it with the residual values.          As with q, the residual evaluation function only interacts with this array,        never with res_vec.res_vec: similar to q_vec.  Unlike q_vec there are functions to perform an            additive reduction (basically a \"gather\") of res to res_vec.              For continuous Galerkin discretizations, the corresponding \"scatter\"            (ie. res_vec -> res) may not exist.res_edge: a 4 dimensional array sometimes used for edge-based data structure.             This feature does not work.  This array must exist and should be             of size 0. M:  The mass matrix of the entire mesh.  Because SBP operators have diagonal       mass matrices, this is a vector.  Length numDofPerNode x numNodes (where       numNodes is the number of nodes in the entire mesh).Minv:  The inverse of the mass matrix.multiplyA0inv:  Multiplies the solution values at each node in an array such as res by the inverse of the coefficient matrix of the time term of the equation.                   This function is used by time marching methods.                   For some equations, this matrix is the identity matrix, so it can be a no-op, while for others might not be.                   The function must have the signature:multiplyA0inv(mesh::AbstractMesh, sbp, eqn::AbstractSolutionData, opts, res_arr::AbstractArray{Tsol, 3})majorIterationCallback:  function called before every step of Newton\'s method or stage of an explicit time marching scheme. This function is used to do output and logging. The function must have the signature:function majorIterationCallback(itr, mesh::AbstractMesh, sbp::AbstractSBP, eqn::AbstractEulerData, opts)params:  user defined type that inherits from AbstractParamType:  CurrentModule = ODLCommonToolsAbstractParamTypeThe purpose of this type is to store any variables that need to be quickly accessed or updated. The only required fields are:t::Float64: hold the current time value\norder: order of accuracy of the discretization (same as AbstractMesh.order)\nx_design: vector of design variables that don\'t fit into any other catagory            (ie. shape variables or other aerodynamic variables like angle of            attack)\ntime::Timings: an object to record how long different parts of the code take,defined in the Utils module.f: a file handle that prints to the file log_myrank.dat, where myrank      is the MPI rank of the current process, but only in debug mode.  In      non-debug mode, this file should not be opened.file_dict: dictionary that maps from the file name to a file handle.  This              field is not required, but if it is present, all copies of the              equation object must share the same dictionary (to avoid problems              with buffering)."
},

{
    "location": "interfaces.html#ODLCommonTools.AbstractMesh",
    "page": "Code Interfaces",
    "title": "ODLCommonTools.AbstractMesh",
    "category": "type",
    "text": "This abstract type is the supertype for all mesh objects.  Every interface to   a mesh software should define its own implementation.\n\nStatic parameters\n\nTmsh: datatype of the mesh data (coordinates, mapping to/from parametric      space, mapping jacobian).\n\nSee the AbstractMesh for the description of everything this   type must implement.\n\n\n\n"
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
    "text": "  # counts\n  numVert::Integer\n  numEl::Integer\n  numNodes::Integer\n  numDof::Integer\n  numDofPerNode::Integer\n  numNodesPerElement::Integer\n  order::Integer\n  numNodesPerFace::Int\n\n  # parallel counts\n  npeers::Int\n  numGlobalEl::Int\n  numSharedEl::Int\n  peer_face_counts::Array{Int, 1}\n  local_element_counts::Array{Int, 1}\n  remote_element_counts::array{Int, 1}\n\n  # MPI Info\n  comm::MPI.Comm\n  myrank::Int\n  commsize::Int\n  peer_parts::Array{Int, 1}\n\n  # Discretization type\n  isDG::Bool\n  isInterpolated::bool\n\n  # mesh data\n  coords::AbstractArray{Tmsh, 3}\n  dxidx::AbstractArray{Tmsh, 4}\n  jac::AbstractArray{Tmsh, 2}\n\n  # boundary data\n  coords_bndry::Array{Tmsh, 3}\n  dxidx_bndry::Array{Tmsh, 4}\n  jac_bndry::Array{Tmsh, 2}\n  nrm_bndry::Array{Tmsh, 3}\n\n  # interior face data\n  coords_interface::Array{Tmsh, 3}\n  dxidx_face::Array{Tmsh, 4}\n  jac_face::Array{Tmsh, 2}\n  nrm_face::Array{Tmsh, 3}\n\n  # parallel data\n  coords_sharedface::Array{Array{Tmsh, 3}, 1}\n  dxidx_sharedface::Array{Array{Tmsh, 4}, 1}\n  jac_sharedface::Array{Array{Tmsh, 2}, 1}  \n  nrm_sharedface::Array{Array{Tmsh, 3}, 1}\n\n  # boundary condition data\n  numBC::Integer\n  numBoundaryEdges::Integer\n  bndryfaces::AbstractArray{Boundary, 1}\n  bndry_offsets::AbstractArray{Integer, 1}\n  bndry_funcs::AbstractArray{BCType, 1}\n\n  # interior edge data\n  numInterfaces::Integer\n  interfaces::AbstractArray{Interface, 1}\n\n  # degree of freedom number data\n  dofs::AbstractArray{Integer, 3}\n  dof_offset::Int\n  sparsity_bnds::AbstractArray{Integer, 2}\n  sparsity_nodebnds::AbstractArray{Integer, 2}\n\n  # mesh coloring data\n  numColors::Integer\n  maxColors::Integer\n  color_masks::AbstractArray{ AbstractArray{Number, 1}, 1}\n  shared_element_colormasks::Array{Array{BitArray{1}, 1}, 1}\n  pertNeighborEls::AbstractArray{Integer, 2}\n\n  # parallel bookkeeping\n  bndries_local::Array{Array{Boundary, 1}, 1}\n  bndries_remote::Array{Array{Boundary, 1}, 1}\n  shared_interfaces::Array{Array{Interface, 1}, 1}\n  shared_element_offsets::Array{Int, 1}\n  local_element_lists::Array{Array{Int, 1}, 1}\nTODO: update with vert_coords etc"
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
    "location": "interfaces.html#Boundary-Data-1",
    "page": "Code Interfaces",
    "title": "Boundary Data",
    "category": "section",
    "text": "This data is used for interpolated mesh only.coords_bndry: coordinates of nodes on the boundary of the mesh,                 2 x numFaceNodes x numBoundaryEdgesdxidx_bndry: 2 x 2 x numFaceNodes x numBoundary edges array ofdxidxinterpolated to the boundary of the mesh.  This fields is                deprecated, usenrm_bndry` instead.jac_bndry: numFaceNodes x numBoundaryEdges array of jac interpolated               to the boundary.  This fields is deprecated, use nrm_bndry               instead.nrm_bndry: dim x numNodesPerFace x numBoundaryFaces array containing the              scaled face normal vector at each face node of each boundary face.              The scaling factor is believed to be 1/|J|, where |J| is the              determinant of the mapping jacobian, however this needs to be              verified.  The normal vector is oriented outwards."
},

{
    "location": "interfaces.html#Interior-Face-Data-1",
    "page": "Code Interfaces",
    "title": "Interior Face Data",
    "category": "section",
    "text": "coords_face: a mesh.dim x mesh.numNodesPerFace x mesh.numInterfaces                array containing he coordinates at each face node of each                interface.dxidx_face: 2 x 2 x numFaceNodes x numInterfaces array of dxidx               interpolated to the face shared between two elements.jac_face: numNodesPerFace x numInterfaces array of jac interpolated              to the face shared between two elementnrm_face: dim x numNodesPerFace x numInterfaces containing the scaled             face normal vector at each face node as each interface.  The scaling             factor is believed to be 1/|J| where |J| is the determinant of              the mapping jacobian, however this needs to be verified.  The normal             vector is oriented outwards from the perspective of elementL of             the corresponding entry of interfaces."
},

{
    "location": "interfaces.html#Parallel-Data-1",
    "page": "Code Interfaces",
    "title": "Parallel Data",
    "category": "section",
    "text": "This data is required for parallelizing interpolated DG meshescoords_sharedface: array of arrays, one array for each peer process,                      containing the coordinates of the nodes on the faces                      shared between a local element on a non-local element.                      Each array is 2 x numFaceNodes x number of faces shared                      with this process. dxidx_sharedface: similar to coords_sharedface, dxidx interpolated to                     faces between elements in different processes, each                     array is 2 x 2 x numFaceNodes x number of faces shared                     with this process.  This field is deprecated, use                      nrm_sharedface` instead.jac_sharedface: similar to coords_sharedface, jac interpolated to faces                   between a local element and a non-local element. Each array                   is numFaceNodes x number of faces shared with this process.                   This field is deprecated, use nrm_sharedface.nrm_sharedface; similar to coords_sharedface, the inner array contains the                   scaled normal vector at each node of each shared face.  Each                   inner array is dim x numNodesPerFace x number of faces                    shared with this process.  The normal vector is oriented                   outward from the perspective of the local element."
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
    "text": "The mesh module must also defines and exports the functionssaveSolutionToMesh(mesh::MeshImplementationType, vec::AbstractVector)\nwriteVisFiles(mesh::MeshImplementationType, fname::String)where the first function takes a vector of length numDof and saves it to the mesh, and the second writes Paraview files for the mesh, including the solution field."
},

{
    "location": "interfaces.html#ODLCommonTools.AbstractCGMesh",
    "page": "Code Interfaces",
    "title": "ODLCommonTools.AbstractCGMesh",
    "category": "type",
    "text": "ODLCommonTools.AbstractCGMesh\n\nThe abstrac type is the supertype of all continuous Galerkin meshes\n\n\n\n"
},

{
    "location": "interfaces.html#ODLCommonTools.AbstractDGMesh",
    "page": "Code Interfaces",
    "title": "ODLCommonTools.AbstractDGMesh",
    "category": "type",
    "text": "ODLCommonTools.AbstractDGGMesh\n\nThe abstract type is the supertype of all discontinuous Galerkin meshes\n\n\n\n"
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
    "text": "The evalResidual function and the fields eqn.q and eqn.res are the interface between the NonlinearSolvers and the physics modules.  The Nonlinear solvers populate eqn.q, and use evalResidual to populate eqn.res, from which the next value if eqn.q is calculated.  The algorthmic differentiation mechanism described above is uses several residual evaluations to compute the Jacobian if needed for a given method.  Some methods, such as RK4, are better expressed in terms of eqn.q_vec and eqn.res_vec rather than eqn.q and eqn.res.  eqn.array3DTo1D and eqn.disassmbleSolution exist to transform back and forth between the vector and 3D array forms.  In order to compute the Jacobian efficiently, it is necessary to work with the 3D arrays. For this reason, evalResidual must work only with eqn.q and eqn.res and let the caller decide whether or not to transform into the vector form.Newton\'s method supports both finite differences and complex step for calculating the Jacobian, and the static parameters need to be set accordingly. If finite differences are used, then Tsol=Tres=Tmsh=Float64 is required.  If complex step is used, then Tsol=Tres=Complex128 and Tmsh = Float64 is needed."
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
    "text": "An important aspect of the use of the mesh, sbp, eqn, and opts to define interfaces is that the physics modules and nonlinear solvers are written in a purely functional programming style, which is to say that the behavior of every function is determined entirely by the arguments to the function, and the only effects of the function are to modify an argument or return a value.This property is important for writing generic, reusable code. For example, when using iterative solvers, a preconditioner is usually required, and constructing the preconditioner requires doing residual evaluations. In some cases, the preconditioner will use a different mesh or a different mesh coloring. Because the physics modules are written in a functional style, they can be used to evaluate the residual on a different mesh simply by passing the residual evaluation function a different mesh object.A critical aspect of function programming is that there is no global state. The state of the solver is determined entirely by the state of the objects that are passed around."
},

{
    "location": "interfaces.html#Variable-Naming-Conventions-1",
    "page": "Code Interfaces",
    "title": "Variable Naming Conventions",
    "category": "section",
    "text": "In an attempt to make code more uniform and readable, certain variable names are reserved for certain uses.mesh:  object that implements AbstractMesh\npmesh:  mesh used for preconditioning\nsbp:  Summation-by-Parts object\neqn:  object that implements AbstractSolutionData\nopts: options dictionary\nparams: parameter object (used for values that might be in opts but need to be accessed quickly)\nx: the first real coordinate direction\ny: the second real coordinate direction\nxi: the first parametric coordinate direction\neta: the second parametric coordinate direction\nh:  mesh spacing\nhx: mesh spacing in x direction\nhy: mesh spacing in y direction\np: pressure at a node\na; speed of sound at a node\ns: entropy at a node\ngamma: specific heat ratio\ngamma_1: gamma - 1\nR: specific gas constant in ideal gas law (units J/(Kg * K) in SI)\ndelta_t: time step\nt: current time\nnrm: a normal vector of some kind\nA0: the coefficient matrix of the time term at a node\nAxi: the flux jacobian at a node in the xi direction\nAeta: the flux jacobian at a node in the eta direction\nAx: the flux jacobian in the x direction at a node\nAy: the flux jacobian in the y direction at a node"
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
    "text": "This document describes how PDEs are solved in parallel.In general, the mesh is partitioned and each part is assigned to a different MPI process. Each element is owned by exactly one process.  During  initialization, the mesh constructor on each process figures out which other processes have elements that share a face (edge in 2D) with local elements. It counts how many faces and elements are shared (a single element could have multiple faces on the parallel boundary), and assigns local number to both the elements and the degrees of freedom on the elements.  The non-local elements  are given numbers greater than numEl, the number of local elements.   The degrees of freedom are re-numbered such that newly assigned dof number plus the dof_offset for the current process equals the global dof number, which is defined by the local dof number assigned by the process that owns the element.  As a  result, dof numbers for elements that live on processes with lower ranks will be negative.As part of counting the number of shared faces and elements, 3 arrays are formed: bndries_local, bndries_remote, and shared_interfaces which  describe the shared faces from the local side, the remote side, and a  unified view of the interface, respectively.  This allows treating the shared faces like either a special kind of boundary condition or a proper  interface, similar to the interior interfaces.There are 2 modes of parallel operation, one for explicit time marching and the other for Newton\'s method."
},

{
    "location": "parallel.html#Explicit-Time-Marching-1",
    "page": "Code Parallelization",
    "title": "Explicit Time Marching",
    "category": "section",
    "text": "In this mode, each process each process sends the solution values at the  shared faces to the other processes.  Each process then evaluates the residual using the received values and updates the solution.The function exchangeFaceData is designed to perform the sending and  receiving of data.  Non-blocking communications are used, and the function does not wait for the communication to finish before returning.  The  MPI_Requests for the sends and receives are stored in the appropriate fields of the mesh.  It is the responsibility of each physics module call  exchangeFaceData and to wait for the communication to finish before using the data.  Because the receives could be completed in any order, it is  recommended to use MPI_Waitany to wait for the first receive to complete,  do as many computations as possible on the data, and then call MPI_Waitany again for the next receive."
},

{
    "location": "parallel.html#Newton\'s-Method-1",
    "page": "Code Parallelization",
    "title": "Newton\'s Method",
    "category": "section",
    "text": "For Newton\'s method, each process sends the solution values for all the  elements on the shared interface at the beginning of a Jacobian calculation.  Each process is then responsible for perturbing the solutions values of both  the local and non-local elements.  The benefit of this is that parallel  communication is required once per Jacobian calculation, rather than once  per residual evaluation as with the explicit time marching mode.The function exchangeElementData copies the data from the shared elements into the send buffer and sends it, and also posts the corresponding receives. It does not wait for the communications to finish before returning.   The function is called by Newton\'s method after a new solution is calculated, so the physics module does not have to do it, but the physics module does  have to wait for the receives to finish before using the data.  This is  necessary to allow overlap of the communication with computation.  As  with the explicit time marching mode, use of MPI_Waitany is recommended. "
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
    "location": "pdesolver_user.html#PDESolver-User-Interface-1",
    "page": "PDESolver User Interface",
    "title": "PDESolver User Interface",
    "category": "section",
    "text": "CurrentModule = PDESolver"
},

{
    "location": "pdesolver_user.html#PDESolver.run_solver",
    "page": "PDESolver User Interface",
    "title": "PDESolver.run_solver",
    "category": "function",
    "text": "This function provides a way to invoke any physics solver based on the   specification of the physics in the input file.\n\nThe physics module must have already been registered using register_physics\n\nThis function is equivalent to solvePDE, it is maintained here   for legacy purposes only.\n\nInputs\n\ninput_file: an AbstractString specifying the path to the input file\n\nOutputs\n\nmesh: the AbstractMesh object used during the solve\nsbp: the SBP operator used by the solver\neqn: the AbstractSolutionData object during the solve.  At exit,      eqn.q_vec should have the final solution in it\nopts: the options dictionary\n\n\n\n"
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
    "category": "function",
    "text": "This function registers a new initial condition function with the specified   physics module.  The function must have the signature:\n\nICfunc(mesh::AbstractMesh, sbp::AbstractSBP, eqn:AbstractSolutionData{Tsol},                 opts:Dict, q_vec::AbstractVector{Tsol})\n\nwhere q_vec is a vector of length numDof that will be populated with the   initial condition.  If the function is used to provide the exact solution   for an unsteady problem (for calculating the error via the exact_soln_func   key in the options dictionary), then it should use eqn.params.t as the   current time value.\n\nThis function does not attempt to verify that the functor has the correct   signature, so the user should take care to ensure it is correct for the   physics module.\n\nThe physics module must have an associative container called ICDict that   accepts Strings as keys and functions as values.\n\nIt is an error to register two functions under the same name or to register   one function under two different names.\n\nInputs:\n\nmod: physics module to register the function with\nfname: name associated with the function, used as the value for the\n       for any key in the options dictionary that specifies an initial\n       condition function\nfunc: the function itself\n\nOutputs:\n\nnone\n\n\n\n"
},

{
    "location": "pdesolver_user.html#PDESolver.registerBC",
    "page": "PDESolver User Interface",
    "title": "PDESolver.registerBC",
    "category": "function",
    "text": "This function registers a new boundary condition functor with the specified   physics module.  The exact signature of the functor varies by physics module.   The purpose of the functor is to calculate the flux across a boundary at   a particular node. It usually takes in the solution values at that node, as   well as the coordinates and mapping jacobian.  The params type is typically   used to dispatch to different methods for 2D or 3D, if needed.  Generally,   this functor will calculate the boundary state and then call a numerical   flux function of some kind to compute the flux.\n\nThis function does not attempt to verify that the functor has the correct   signature, so the user should take care to ensure it is correct for the   physics module.\n\nIt is an error to register two functors with the same name or the same   functor under two different names.  An exception is made for \"reanalysisBC\",   because it is special.\n\nInputs\n\nmod: module to register the the functor with\nfname: the name associated with this function, used as the value for any        key in the options dictionary that specifies a boundary condition,        for example BC1_name\nfunc: the functor itself\n\nOutputs\n\nnone\n\nNotes For Implementers\n\nThis function works by utilizing mod.BCDict, which must be an associative   collection mapping BC names to functors.\n\n\n\n"
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
    "category": "function",
    "text": "Print all initial condition names for a given physics module\n\nInputs:\n\nmod: the module to print the ICs for\nf: the IO to print to, defaults to STDOUT\n\nOutputs:     none\n\n\n\n"
},

{
    "location": "pdesolver_user.html#PDESolver.printBCNames",
    "page": "PDESolver User Interface",
    "title": "PDESolver.printBCNames",
    "category": "function",
    "text": "Like printICNames, but for boundary conditions, see that function for details\n\n\n\n"
},

{
    "location": "pdesolver_user.html#PDESolver.printPhysicsModules",
    "page": "PDESolver User Interface",
    "title": "PDESolver.printPhysicsModules",
    "category": "function",
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
    "page": "PDESolver Physics Interface",
    "title": "PDESolver Physics Interface",
    "category": "page",
    "text": ""
},

{
    "location": "pdesolver_physics.html#Documentation-of-the-PDESolver-Physics-Module-Interface-1",
    "page": "PDESolver Physics Interface",
    "title": "Documentation of the PDESolver Physics Module Interface",
    "category": "section",
    "text": "CurrentModule = PDESolverThe PDESolver Module provides some structure for physics modules to plug into. See Interfaces in PDESolver for details. Each physics module should create a new subtype of AbstractSolutionData. This physics module should then extend the required methods, specializing one of the arguments with the new (physics-specific) type. This structure allows the API for the physics module to be visible to other parts of the code (such as NonlinearSolvers), while being defined  within the physics module.One exception to this pattern of defining generic functions and then specializing one argument is the initialization functions. These functions run before the AbstractSolutionData object is created, so they can not dispatch to different methods. Instead, each physics module has to register itself and provide the required initialization functions. See the Registration Functions section."
},

{
    "location": "pdesolver_physics.html#PDESolver.createFunctional-Union{Tuple{I}, Tuple{ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,ODLCommonTools.AbstractSolutionData,Any,AbstractString,Array{I,1}}} where I<:Integer",
    "page": "PDESolver Physics Interface",
    "title": "PDESolver.createFunctional",
    "category": "method",
    "text": "Creates a functional object.\n\nEach physics modules should extend this function with a new method,   specializing the eqn argument.\n\nInputs\n\nmesh : Abstract PUMI mesh\nsbp  : Summation-by-parts operator\neqn  : AbstractSolutionData object\nopts : Options dictionary\nfunctional_name: the name of the functional\nfunctional_bcs: the boundary condition numbers the functional is                   computed on.\n\nOutputs\n\nfunctional: an AbstractFunctional.  Usually an             AbstractIntegralFunctional, although this is not             required.\n\n\n\n"
},

{
    "location": "pdesolver_physics.html#PDESolver.evalFunctional-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP,ODLCommonTools.AbstractSolutionData{Tsol,Tres} where Tres,Any,ODLCommonTools.AbstractFunctional}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "PDESolver Physics Interface",
    "title": "PDESolver.evalFunctional",
    "category": "method",
    "text": "High level function that evaluates the given functional   This function is agnostic to the type of the functional being   computed and calls a mid level functional-type specific function for the    actual evaluation.\n\nThe functional is evaluated at the state in eqn.q_vec.\n\nArguments\n\nmesh :  Abstract mesh object\nsbp  : Summation-By-Parts operator\neqn  : AbstractSolutionData object\nopts : Options dictionary\nfunctionalData : Object of type AbstractFunctional.                      This type determines the functional being computed and                     holds all the relevant data.\n\n\n\n"
},

{
    "location": "pdesolver_physics.html#PDESolver.evalFunctionalDeriv-Union{Tuple{ODLCommonTools.AbstractDGMesh{Tmsh},SummationByParts.AbstractSBP,ODLCommonTools.AbstractSolutionData{Tsol,Tres} where Tres,Any,ODLCommonTools.AbstractIntegralFunctional,AbstractArray{T,3} where T}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "PDESolver Physics Interface",
    "title": "PDESolver.evalFunctionalDeriv",
    "category": "method",
    "text": "Computes a 3D array of hte derivative of a functional wrt eqn.q.\n\nThe derivative is evaluated at the state in eqn.q_vec.\n\nInputs\n\nmesh\nsbp\neqn\nopts\nfunctionalData: AbstractIntegralFunctional to evaluate\n\nInputs/Outputs\n\nfunc_deriv_arr: array to hold derivative of function wrt eqn.q, same                 size as equation.q\n\nOptions Keys\n\nThis funciton is not compatible with precompute_q_bndry = false\n\n\n\n"
},

{
    "location": "pdesolver_physics.html#PDESolver.evalHomotopyJacobian-Tuple{ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,ODLCommonTools.AbstractSolutionData,Dict,PDESolver.AssembleElementData,Number}",
    "page": "PDESolver Physics Interface",
    "title": "PDESolver.evalHomotopyJacobian",
    "category": "method",
    "text": "Evaluates the Jacobian of the evalHomotopy multiplied by the   homotopy parameter lambda.  Conceptually similar to evalJacobian.\n\nInputs\n\nmesh: an AbstractMesh\nsbp: an AbstractSBP\neqn: an AbstractSolutionData (physics modules should specialize this      argument)\nopts: options dictionary\nassembler: an object used to assemble contributions into the Jacobian\nlambda: the homotopy parameter\n\n\n\n"
},

{
    "location": "pdesolver_physics.html#PDESolver.evalJacobian",
    "page": "PDESolver Physics Interface",
    "title": "PDESolver.evalJacobian",
    "category": "function",
    "text": "Similar function to evalResidual, but instead of computing the   residual, it computes and assembles the Jacobian of the residual with   respect to eqn.q into the system matrix.\n\nThe functions assembleElement and assembleInterface   in NonlinearSolvers should be used to assembles the Jacobians of   individual elements and interfaces into the system matrix.\n\nPhysics modules that can compute element Jacobians directly   should extend this function with a new method, otherwise the coloring   algorithm with be used with evalResidual to compute the   system matrix.  For physics modules that support multiple formuations,   this function should throw an exception if called with an unsupported   formulation.\n\nInputs\n\nmesh: an AbstractMesh\nsbp: an SBP operator\neqn: AbstractSolutionData (physics modules should specialize this      argument)\nopts: options dictionary\nassembler: object that must be passed to assembleElement and             assembleInterface\n\nKeyword Arguments\n\nstart_comm: whether or not to start parallel communication,             default false because the functions Nonlinear solvers             have already done parallel communication when this             function gets called.\n\nOptions Keys\n\nTODO: describe the key that controls whether this function gets called         or not.\n\n\n\n"
},

{
    "location": "pdesolver_physics.html#PDESolver.evalJacobianStrong",
    "page": "PDESolver Physics Interface",
    "title": "PDESolver.evalJacobianStrong",
    "category": "function",
    "text": "This function evaluates the Jacobian of the strong form of the spatial   residual.  Note that the residual is written\n\ndu/dt = -(Q * f) + SAT\n\nNote the negative sign.\n\nCurrently this function neglects the SAT terms (both interface and boundary   conditions)\n\nInputs\n\nmesh: an AbstractMesh\nsbp: an SBP operator\neqn: AbstractSolutionData (physics modules should specialize this      argument)\nopts: options dictionary\nassembler: object that must be passed to assembleElement and             assembleInterface\n\n\n\n"
},

{
    "location": "pdesolver_physics.html#PDESolver.evalResidual",
    "page": "PDESolver Physics Interface",
    "title": "PDESolver.evalResidual",
    "category": "function",
    "text": "This function evalutes dq/dt = R(q).  For steady problems it evalutes R(q)   at some state q.  The state is stored in eqn.q, and eqn.res is populated with   R(q).  Note that these are 3 dimensional arrays.  The physics modules only   interact with the 3 dimensional arrays, never the vectors eqn.q_vec and   eqn.res_vec.  Each physics module must implement this function for its   own subtype of AbstractSolutionData (ie. with a more specific type for   the eqn argument and equallty specific types for the other arguments).   This is important because evalResidual is common to all physics modules,   so a user or some other part of the code can call evalResidual(mesh, sbp   eqn, opts), and Julia\'s multiple dispatch will figure out the right method   to call based on the type of the eqn argument.\n\nThe evaluation of the residual R(q) should depend only on the data stored in   mesh, sbp, eqn, and opts, and any data that depend on q should be recalculated   every time the function is called.  This function is used as a building block   by other parts of the solver, particularly the NonlinearSolvers.  See   interfaces.md for details\n\nInputs\n\nmesh: an AbstractMesh describing the mesh on which to solve the physics\nsbp: an SBP operator\neqn: a subtype of AbstractSolution data, used to store all of the data used by the physics module\nopts: the options dictionary\nt: the current time value, defaults to 0.0\n\n#TODO: list required options keys\n\n\n\n"
},

{
    "location": "pdesolver_physics.html#PDESolver.solvePDE",
    "page": "PDESolver Physics Interface",
    "title": "PDESolver.solvePDE",
    "category": "function",
    "text": "This function takes fully initialzed objects and solves the partial   differential equation.\n\nEvery physics should extend this function with a new method, specialzing the   eqn argument type.\n\nObjects can be created with the createObjects function.\n\nThe functions in the SolverCommon module are helpful in writing   this function.\n\nInputs\n\nmesh: an AbstractMesh\nsbp: an SBP Operator\neqn: an AbstractSolutionData\nopts: the options dictionary\npmesh: mesh for preconditioning (optional)\n\nOutputs\n\nmesh\nsbp\neqn: on exit, eqn.q_vec should have the converged solution in it.\nopts\n\nOptions Keys\n\nRelfunc_name: also writes vtk files called \"solution_relfunc\"              if key not present, ignored               TODO: fix that\nIC_name\ncalc_error: also write vtk files called \"solution_error\"\ncalc_trunc_error\nperturb_ic\ncalc_dt\n\nFor options like calc_dt and Relfunc_name, it is very important that\nthe computed quantity be saved to the options dictionary for use later\nin the code (ie. and not passed directly to another function).  The\ncode won\'t restart correctly if this happens.\n\n\n\n"
},

{
    "location": "pdesolver_physics.html#PDESolver.updateMetricDependents-Tuple{ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,ODLCommonTools.AbstractSolutionData,Any}",
    "page": "PDESolver Physics Interface",
    "title": "PDESolver.updateMetricDependents",
    "category": "method",
    "text": "This function recalculates all quantities that depend on the mesh metrics.   This function handles changes in the metric values only (not changed in   mesh topology).\n\nThis function should be called after the mesh is warped (ie. using mesh   movement).\n\nEvery physics module should extend this function with a new method   specializing the eqn argument type.\n\nInputs\n\nmesh\nsbp\neqn\nopts\n\n\n\n"
},

{
    "location": "pdesolver_physics.html#Functions-to-be-Extended-1",
    "page": "PDESolver Physics Interface",
    "title": "Functions to be Extended",
    "category": "section",
    "text": "All the functions in this section should be extended by each physics module with a new method, specializing the type of the eqn argument.Modules = [PDESolver]\nPages = [\"src/interface.jl\"]"
},

{
    "location": "pdesolver_physics.html#PDESolver.createFunctional",
    "page": "PDESolver Physics Interface",
    "title": "PDESolver.createFunctional",
    "category": "function",
    "text": "Creates a functional object, using the data described in the options   dictionary\n\nArguments\n\nmesh : Abstract PUMI mesh\nsbp  : Summation-by-parts operator\neqn  : AbstractSolutionData\nopts : Options dictionary\nfunctional_number : which functional (of the functionals described                       in the options dictionary) to create\n\nOutputs\n\nsame as other method\n\n\n\n"
},

{
    "location": "pdesolver_physics.html#PDESolver.createObjects-Tuple{AbstractString}",
    "page": "PDESolver Physics Interface",
    "title": "PDESolver.createObjects",
    "category": "method",
    "text": "This function returns the fully initialized objects needed to solve an   equations.\n\nInputs\n\nfname: input file name\n\nOutputs\n\nmesh: an AbstractMesh of some kind (depending on input file)\nsbp: an SBP operator of some kind (depending on input file)\neqn: an AbstractSolutionData of some kind (depending on input file)\nopts: options dictionary\npmesh: mesh for computing preconditioning, may be the same object as mesh\n\n\n\n"
},

{
    "location": "pdesolver_physics.html#PDESolver.createObjects-Tuple{Dict}",
    "page": "PDESolver Physics Interface",
    "title": "PDESolver.createObjects",
    "category": "method",
    "text": "Method for creating objects from an options dictionary\n\nInputs\n\nopts: the options dictionary\n\nOutputs\n\nsee other method\n\n\n\n"
},

{
    "location": "pdesolver_physics.html#PDESolver.createObjects-Tuple{ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,Dict}",
    "page": "PDESolver Physics Interface",
    "title": "PDESolver.createObjects",
    "category": "method",
    "text": "This method constructs a new equation object for the given mesh, sbp   and opts.\n\nThis is useful for solving on a submesh.\n\nInputs\n\nmesh: an AbstractMesh\nsbp: an AbstractSBP\nopts: options dictionary\n\nOutputs\n\nsame as other method\n\n\n\n"
},

{
    "location": "pdesolver_physics.html#PDESolver.solvePDE-Tuple{Union{AbstractString, Dict}}",
    "page": "PDESolver Physics Interface",
    "title": "PDESolver.solvePDE",
    "category": "method",
    "text": "Additional methods for solvePDE, allows solving the PDE starting   with either a file name or options dictionary.\n\nSee the other method for details\n\nInputs\n\nopts: either a file name to load an options dictionary from or the       options dictionary itself\n\nOutputs\n\nsame as other method\n\n\n\n"
},

{
    "location": "pdesolver_physics.html#Additional-Functions-1",
    "page": "PDESolver Physics Interface",
    "title": "Additional Functions",
    "category": "section",
    "text": "These functions are available for all physics module, but do not need to be extended with a new method. They usually boil down to calling one of the functions in the previous section.Modules = [PDESolver]\nPages = [\"src/interface2.jl\"]"
},

{
    "location": "pdesolver_physics.html#PDESolver.register_physics",
    "page": "PDESolver Physics Interface",
    "title": "PDESolver.register_physics",
    "category": "function",
    "text": "This function registered a new physics module with the global list of all   known physics modules.  Every physics module should do this as part of   module initialization.  The name, the module, and the startup function must   be unique (ie. they must not already exist in the list).  This function   throws and exception if they are not.\n\nInputs\n\nmodname:  an String name for this entry in the list.  It is used            to retrieve the module and startup function in the             retrieve_physics function. Typically the name is capitalized.\nmod:  the Module itself\n_createObjects: function that creates the mesh, sbp, eqn, and opts                  objects. Must have a method with signature                  _createObjects(opt::Dict)\n             If this physics modules supports solving on a submesh,\n             this function should also have a method\n             `_createObjects(mesh::AbstractMesh, sbp::AbstractSBP, opts::Dict)`\n             to create a new equation object.  See [`createObjects`](@ref).\n_checkOptions: physics-specific function for supplying default options and               checking options.  Note that most of the input option               processing is done in read_input,               _checkOptions need only do the physics-specific part.               This function must have signature _checkOptions(opts::Dict).\n\nOutputs\n\nnone\n\n\n\n"
},

{
    "location": "pdesolver_physics.html#PDESolver.retrieve_physics",
    "page": "PDESolver Physics Interface",
    "title": "PDESolver.retrieve_physics",
    "category": "function",
    "text": "Retrieves the physics module and function registered using register_physics\n\nInput\n\nmodname: an String containing the name of the module supplied to         register_physics\n\nOutputs\n\nmod: the physics Module\n_createObjects: function that creates the solver objects\n_checkOptions: the options checking function\n\n\n\n"
},

{
    "location": "pdesolver_physics.html#Registration-Functions-1",
    "page": "PDESolver Physics Interface",
    "title": "Registration Functions",
    "category": "section",
    "text": "These function provide the basic API for registering and retrieving physics modules with the PDESolver module.register_physics\nretrieve_physicsSee also registerIC and registerBC."
},

{
    "location": "pdesolver_structure.html#",
    "page": "PDESolver Structure",
    "title": "PDESolver Structure",
    "category": "page",
    "text": ""
},

{
    "location": "pdesolver_structure.html#PDESolver.run_solver-Tuple{AbstractString}",
    "page": "PDESolver Structure",
    "title": "PDESolver.run_solver",
    "category": "method",
    "text": "This function provides a way to invoke any physics solver based on the   specification of the physics in the input file.\n\nThe physics module must have already been registered using register_physics\n\nThis function is equivalent to solvePDE, it is maintained here   for legacy purposes only.\n\nInputs\n\ninput_file: an AbstractString specifying the path to the input file\n\nOutputs\n\nmesh: the AbstractMesh object used during the solve\nsbp: the SBP operator used by the solver\neqn: the AbstractSolutionData object during the solve.  At exit,      eqn.q_vec should have the final solution in it\nopts: the options dictionary\n\n\n\n"
},

{
    "location": "pdesolver_structure.html#Starting-a-Simulation-1",
    "page": "PDESolver Structure",
    "title": "Starting a Simulation",
    "category": "section",
    "text": "  CurrentModule = PDESolverThis page describes functions located in the PDESolver module that tie together the physics modules and the Nonlinear solvers.    Modules = [PDESolver]\n  Pages = [\"src/startup_func.jl\", \"src/initialization.jl\"]To run a simulation, the following steps must be doneread the input dictionary\ncreate an AbstractMesh and AbstractSBP, and AbstractSolutionData\nLoad an initial condition\nCalculate various quantities\nInvoke a NonlinearSolver\nDo postprocessingSome of these steps are handled by the PDESolver module, and a few are handled by the physics module."
},

{
    "location": "pdesolver_structure.html#Input-Dictionary-1",
    "page": "PDESolver Structure",
    "title": "Input Dictionary",
    "category": "section",
    "text": "The first step is to read the input file.   Reading input files is split into two parts.  The first part is done by the Input module, which loads the file from disk and supplies default values. The second part is done by the physics-specific registered with register_physics. this  function that verifies the physics module supports the given options (especially checking for combinations of  CurrentModule = EulerEquationModoptions that might not be supported).  See, for example, checkOptions.  CurrentModule = PDESolver"
},

{
    "location": "pdesolver_structure.html#Creating-Objects-1",
    "page": "PDESolver Structure",
    "title": "Creating Objects",
    "category": "section",
    "text": "The next thing the physics module needs to do is create AbstractSBP and AbstractMesh, and AbstractSolutionData  objects. The function createMeshAndOperator should be used by all physics modules to create the first two.The AbstractSolutionData is created by the physics module itself. The details of how to do this are left up to the physics module, but the return values of createMeshAndOperator should be used for static parameter values.These creation of all three objects are performed by the _createObjects functions provided to register_physics."
},

{
    "location": "pdesolver_structure.html#Load-an-initial-condition-1",
    "page": "PDESolver Structure",
    "title": "Load an initial condition",
    "category": "section",
    "text": "The function solvePDE is extended by each physics module to is do the remaining operations.The details of how to load an initial condition are left up to the physics module, but the end result must be the initial condition is present in eqn.q_vec.Physics modules generally use a Dictionary to map IC names (which is how  ICs are referred to in the input file) to the function that applies the IC.  See registerIC for further description."
},

{
    "location": "pdesolver_structure.html#Various-calculations-1",
    "page": "PDESolver Structure",
    "title": "Various calculations",
    "category": "section",
    "text": "After loading the IC, the options dictionary may require the calculation of a few quantities.  See solvePDE for the list of options keys that must be supported."
},

{
    "location": "pdesolver_structure.html#Invoke-a-NonlinearSolver-1",
    "page": "PDESolver Structure",
    "title": "Invoke a NonlinearSolver",
    "category": "section",
    "text": "  CurrentModule = PDESolverThe next step is calling a Nonlinear Solver. The function call_nlsolver takes the objects already constructed and calls the appropriate nonlinear solver. Currently, there are no strong guarantees about where the solution is stored (eqn.q or eqn.q_vec) when this function returns (TODO: fix that)."
},

{
    "location": "pdesolver_structure.html#Do-Postprocessing-1",
    "page": "PDESolver Structure",
    "title": "Do Postprocessing",
    "category": "section",
    "text": "  CurrentModule = EulerEquationModThe options dictionary may require post-processing be done, for example calculating the solution error when the analytical solution is known. Each physics module usually defines a function to do this. See postproc for an example.  CurrentModule = PDESolver"
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
    "text": "Before running the code, you must source the shell script PumiInterface/src/use_julialib.sh.  This enables Julia to find and use Pumi.cd /path/to/PumiInterface/src # for example ~/.julia/v0.4/PumiInterface/src)\nsource ./use_julialib.sh  # note that you must be in the same directory as the script when sourcing itThe code takes an input file that defines all options for the solver, including the which physics to solve, which mesh to use, initial conditions, boundary conditions, and discretization.  The file src/input/input_vals.txt describes the input file format and all valid options."
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
    "location": "invocation/interactive.html#Running-PDESolver-in-an-interactive-session-using-Julia\'s-REPL-1",
    "page": "Interactive Session (experimental)",
    "title": "Running PDESolver in an interactive session using Julia\'s REPL",
    "category": "section",
    "text": "Julia\'s REPL is a powerful tool for debugging.  Being able to run commands interactively can reveal important information   that would otherwise be available only through cumbersome print statements and  error logs.This page describes the first steps involved for running PDESolver interactively. Calling PDESolver from Julia\'s REPL in the fashion described here is an experimental    method, and is not tested consistently. Additionally, not all parts of PDESolver are usable after the steps below. Instead, this document is intended to provide a springing-off point for the developer    who wishes to run PDESolver interactively, and can adapt the commands below to their needs."
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
    "text": "Here are some helpful Julia commands for debugging/investigation within the REPL at this stage. Refer to the Julia documentation for details, or press \'?\' within the REPL    and type the command you wish to see help for.size()\nlength()\ntypeof()\nfieldnames()\nisimmutable()\nwhos()"
},

{
    "location": "solver/Readme.html#",
    "page": "Overview of Physics Modules",
    "title": "Overview of Physics Modules",
    "category": "page",
    "text": ""
},

{
    "location": "solver/Readme.html#sec:physics_modules-1",
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
    "text": "It is useful to divide functions into 3 catagories, high, mid, and low level functions.  The purpose of high level functions is to decide which method of performing an operation should be used and call other functions to do it. For example, the Euler physics modules has evalVolumeIntegrals and evalBoundaryIntegrals as high level functions.  There are several different ways of calculating both the volume and boundary integrals.  The options dictionary is used to decide what mid level function to call.  Each mid level function implements a different way of doing the calculation.The purpose of mid level functions is to loop over the mesh and call a low level function for each node.  For example, the function getEulerFlux loops over the nodes of the mesh and calls a function to calculate the Euler flux at each node.  Mid level function names usually start with get to indicate that their purpose is to calculate some quantity but don\'t do the calculation themselves.Low level functions calculate a quantity at a node.  For example, calcEulerFlux calculates the Euler flux at a single node.  Low level function names usually start with calc to indicate that they perform a specific calculation. Often, different discretizations use the same structure of loops, but do a slightly different calculation at each node.  Low level functions are called inside the innermost loop of the code, so it would be too expensive to have if statements to select which low level function to call, so various tricks involving Julia\'s multiple dispatch system are used to get the compiler to decide which low level function to call.  These will be described later in this document.It is often useful to dispatch to low level functions based on Tdim and var_type.  For this reason the Euler equation implementation of AbstractParamType istype ParamType{Tdim, var_type, Tsol, Tres, Tmsh} <: AbstractParamType{Tdim}The other static parameters are necessary because ParamType has fields of those datatypes."
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
    "text": "Many of the components of PDESolver have different options that control how they work and what they do. In order to  provide a unified method of specifying these options, an dictionary  of type Dict{String, Any} is read in from a disk file. This dictionary (called opts in function signatures), is passed to all high and mid level functions so they can use values in the dictionary to determine their  control flow. Low level functions need to be extremely efficient, so they cannot have conditional logic, therefore they are not passed the dictionary. Note that retrieving values from a dictionary is very slow compared to accessing the fields of a type, so all values that are accessed repeatedly should be stored  as the field of a type."
},

{
    "location": "solver/Readme.html#Functors-1",
    "page": "Overview of Physics Modules",
    "title": "Functors",
    "category": "section",
    "text": "Functors are a trick used to get Julia\'s dispatch system to make decisions at compile time rather than runtime.  This is particularly useful for boundary conditions, where the list of mesh faces that have boundary conditions applied is determined at runtime, but having conditional statements that execute for every node on the mesh boundary would be slow.  Instead a construct is used as follows:mutable struct myBC <: BCType  # create a singleton type\nend\n\nfunction call(obj::myBC, q::AbstractVector, bndryflux::AbstractVector)\n  # calculate boundary flux here\nendThis defines a datatype and adds a method to the call function for that type. The call function is what makes a datatype callable like a function.  This method is called as follows:functor = myBC()  # construct and object of type myBC\nq = rand(4)\nbndryflux = zeros(4)\nfunctor(q, bndryflux)  # the Julia compiler turns this into call(functor, q, bndryflux)  The way this is used for boundary conditions is through a two level construct where an outer function passes a functor to an inner function.  Julia\'s JIT will generate a method of the inner function that is specialized to the functor (this is why it is important that the functor is a datatype).  For example:function getBCFluxes(mesh, sbp, eqn, opts)\n\n  for i=1:mesh.numBC  # loop over different boundary conditions\n    functor_i = mesh.bndry_functor[i]  # get the functor for this boundary condition\n    start_index = mesh.bndry_offsets[i]\n    end_index = mesh.bndry_offsets[i+1] - 1\n    # get data for boundary faces start_index:end_index\n\n    calcBoundaryFlux(functor_i, data for boundary faces start_index:end_index)\n  end\nend  # end function\n\n  function calcBoundaryFlux(functor_i::BCType, data for boundary faces start_index:end_index)\n    for i=1:length(start_index:end_index)\n      for j=1:num_nodes_on_face\n        # get data for this boundary face node\n        functor_i(data for this boundary face node)\n      end\n    end\n\n  end  # end functionThe benefit of this arrangement is that mesh.numBC different version of calcBoundaryFlux get compiled, one for each functor, and each version knows about the call method that was defined for the functor it is passed.  This two level scheme allows the compiler to make all the decisions about what function to call (ie. the call method of the functor), avoiding any conditional logic at runtimeThis idea is also applicable to the flux functions used by DG methods."
},

{
    "location": "solver/Readme.html#Creating-a-Manufactured-Solution-1",
    "page": "Overview of Physics Modules",
    "title": "Creating a Manufactured Solution",
    "category": "section",
    "text": "Using the Method of Manufactured Solutions is an effective way to verify the correctness of the code.  A guide to deriving the solution can be found here."
},

{
    "location": "solver/Readme.html#Solution-and-Derivative-Functions-1",
    "page": "Overview of Physics Modules",
    "title": "Solution and Derivative Functions",
    "category": "section",
    "text": "For simple equations such as linear advection, the general approach is to define a function that evalutes the manufactured solution and its derivatives at a given point, and from that create the required source term, initial condition, and boundary condition. Most physics modules have a file called common_funcs.jl where the solution and derivative evaluation functions go. Make sure the new function you create has the same signature as the other functions in the file. It is often useful to create different methods of the same function when creating related manufactured solutions for 2D and 3D. If the signature requires an AbstractParamType, its static parameter can be used to distinguish the methods. Creating two methods of the same function, rather than two different functions, will make it possible to define source terms, boundary conditions, and initial conditions that work for both 2D and 3D.For complicated equations such as Euler, it is tedious to construct the  source term, boundary condition, and initial condition from the solution and its derivatives. In this case, it is better to use a symbolic math program to generate expressions for the source term, boundary condition, and initial condition directly. Some symbolic math programs have the option to generate C or Fortran code, which can be easily converted to Julia code."
},

{
    "location": "solver/Readme.html#Source-Term-1",
    "page": "Overview of Physics Modules",
    "title": "Source Term",
    "category": "section",
    "text": "To create the source term functor, locate the file where the source terms are defined for the physics modules, usually called source.jl and create a new functor and associated call() method (see description of functors above). Make sure the functor object is a subtype of SRCType and the call() method has the same signature (except for the first argument) as the other call methods in the file. The name of the functor should be SRCfoo, where foo is the name of the  source term.For simple equations such as linear advection, the body of the call() function should construct the source term from the functions in common_funcs.jl. For more complicated equations, the code that evalutes the source term at a given point should be placed in the body of the call() function directly.Note that the purpose of this function is to evalute the value of the source term, not to do integration of any kind.Once the functor is created, it should be added to the list of source terms (usually a Dictionary located at the bottom of the file where the source terms are defined). Consult the physics module documentation for details."
},

{
    "location": "solver/Readme.html#Boundary-Condition-1",
    "page": "Overview of Physics Modules",
    "title": "Boundary Condition",
    "category": "section",
    "text": "Construction of a boundary term is similar to construction of a source term. Locate the file where the boundary conditions are defined, usually bc.jl, and add a new functor. Make sure the functor type is a subtype of BCType and the call() method has the same signature (except for the first argument) as the other call() methods in the file. The naming convention for BC functors is fooBC, where foo is the name of the  boundary condition. The body of the call() method should evalute the flux caused by the imposition of the boundary condition (because boundary conditions are imposed weakly). This is typically accomplished by calculating the boundary condition state and then calling a numerical flux function with both the current state and the boundary state.For simple equations, the boundary state should construct the boundary state by calling the functions in common_funcs.jl. For more complicated equations, the code to evalute the boundary state should be contained in the call() method body.Once the functor is created, it should be added to the list of boundary conditions, usually a dictionary located at the bottom of the file where the boundary conditions are defined. Consults the physical module documentation for details."
},

{
    "location": "solver/Readme.html#Initial-condition-1",
    "page": "Overview of Physics Modules",
    "title": "Initial condition",
    "category": "section",
    "text": "Initial conditions are a bit different than boundary conditions and source terms because they do not use functors (functors are unnecessary because ICs are evaluated infrequently). Locate the file where initial conditions are defined, typically ic.jl, and create a new function with the same signature as the the other functions in the file. This function should loop over all elements in the mesh, every node on the element, and use mesh.dofs to assign the solution to proper indices of the supplied vector. The naming convention for IC functions is ICfoo, where foo is the name of the initial condition.For simple equations, the solution should be calculated using the functions in common_funcs.jl, otherwise it should be calculated in the initial condition function.Initial condition functions are used to calculate errors during post-processing, so it is important for the initial condition function to evaluate the solution at the proper time for unsteady problems.After the initial condition function is created, it should be added to the list of initial conditions, usually a dictionary at the bottom of the file where the initial conditions are defined. See the physics module documentation for details."
},

{
    "location": "solver/Readme.html#Initialization-of-a-Simulation-1",
    "page": "Overview of Physics Modules",
    "title": "Initialization of a Simulation",
    "category": "section",
    "text": "This section lists an outline of how a simulation gets launched After step 4, the procedure becomes a bit more complicated because there are optional steps. Only the required steps are listed below.The options dictionary is read in.  Default values are supplied for any key that is not specified, if a reasonable default value exists.\nSecond, the sbp operator is constructed.\nThe mesh object is constructed, using the options dictionary and the sbp operator.  Some of the options in the dictionary are used to determine how the mesh gets constructed.  For example, the options dictionary specifies what kind of mesh coloring to do.\nThe eqn object is constructed, using the mesh, sbp, and opts objects\nThe physics module init function is called, which initializes the physics module and finishes any initialization that mesh and eqn objects require.\nThe initial condition is applied to eqn.q_vec.\nA nonlinear solver is called.  Which solver is called and what parameters it uses are determined by the options dictionary.\nPost-processing is done, if required by the options dictionary."
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
    "category": "type",
    "text": "Abstract supertype of all boundary condition functors\n\n\n\n"
},

{
    "location": "solver/misc.html#ODLCommonTools.BCType_revm",
    "page": "Assorted Function and Types",
    "title": "ODLCommonTools.BCType_revm",
    "category": "type",
    "text": "Abstract supertype of all boundary condition functors that compute the   reverse mode with respect to the metrics\n\n\n\n"
},

{
    "location": "solver/misc.html#ODLCommonTools.SRCType",
    "page": "Assorted Function and Types",
    "title": "ODLCommonTools.SRCType",
    "category": "type",
    "text": "Abstract supertype of all source term functors\n\n\n\n"
},

{
    "location": "solver/misc.html#ODLCommonTools.FluxType",
    "page": "Assorted Function and Types",
    "title": "ODLCommonTools.FluxType",
    "category": "type",
    "text": "Abstract supertype of all numerical flux functions used by standard DG face   integrals\n\n\n\n"
},

{
    "location": "solver/misc.html#ODLCommonTools.FluxType_revm",
    "page": "Assorted Function and Types",
    "title": "ODLCommonTools.FluxType_revm",
    "category": "type",
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
    "category": "type",
    "text": "ODLCommonTools.Boundary\n\nUsed to identify boundary faces in a finite-element grid.\n\nFields\n\nelement : index of the element to which the boundary face belongs\nface : the face index of the boundary (local index to the element)\n\nExample\n\nTo mark face 2 of element 7 to be a boundary face, use Boundary(7,2)\n\n\n\n"
},

{
    "location": "solver/misc.html#ODLCommonTools.Interface",
    "page": "Assorted Function and Types",
    "title": "ODLCommonTools.Interface",
    "category": "type",
    "text": "ODLCommonTools.Interface\n\nUsed to identify interfaces between elements in a finite-element grid.\n\nFields\n\nelementL : index of the so-called left element in the pair\nelementR : index of the so-called right element in the pair\nfaceL : the face index of the interface with respect to the left element\nfaceR : the face index of the interface with respect to the right element\norient : orientation of the \'right\' element relative to the \'left\'\n\nExample\n\nConsider an interface between elements 2 and 5.  Suppose the interface is on face 1 of element 2 and face 3 of element 5.  Furthermore, suppose element 5 has orientation 1 relative to element 1 (defintion of orientation TBD).  This can be indicated as Interface(2,5,1,3,1)\n\n\n\n"
},

{
    "location": "solver/misc.html#ODLCommonTools.getElementL",
    "page": "Assorted Function and Types",
    "title": "ODLCommonTools.getElementL",
    "category": "function",
    "text": "This function returns either the element field of a Boundary or the   elementL field of an interface.\n\n\n\n"
},

{
    "location": "solver/misc.html#ODLCommonTools.getFaceL",
    "page": "Assorted Function and Types",
    "title": "ODLCommonTools.getFaceL",
    "category": "function",
    "text": "This function returns either the face field of a Boundary or the   faceL field of an Interface\n\n\n\n"
},

{
    "location": "solver/misc.html#Base.show-Tuple{IO,ODLCommonTools.Boundary}",
    "page": "Assorted Function and Types",
    "title": "Base.show",
    "category": "method",
    "text": "Show method for Boundary objects\n\n\n\n"
},

{
    "location": "solver/misc.html#Base.show-Tuple{IO,ODLCommonTools.Interface}",
    "page": "Assorted Function and Types",
    "title": "Base.show",
    "category": "method",
    "text": "Show method for Interface objects\n\n\n\n"
},

{
    "location": "solver/misc.html#Boundaries-and-Interfaces-1",
    "page": "Assorted Function and Types",
    "title": "Boundaries and Interfaces",
    "category": "section",
    "text": "All physics modules need to apply boundary conditions and all DG schemes and some CG schemes need to do face integrals. The types and functions described here assist in identifying them. PDESolver does not track the faces of elements directly, instead it  tracks the element the face belongs to and the local face number, that is, the index of the face in the list of all faces that belong to the element. Using this representation, every interior face has two representations because it is part of two elements.Boundary\nInterface\ngetElementL\ngetFaceL\nshow(::IO, ::Boundary)\nshow(::IO, ::Interface)"
},

{
    "location": "solver/misc.html#ODLCommonTools.AbstractFunctional",
    "page": "Assorted Function and Types",
    "title": "ODLCommonTools.AbstractFunctional",
    "category": "type",
    "text": "ODLCommonTools.AbstractFunctional\n\nAbstract datatype for optimization related data structures. All data types corresponding to optimization problems should be a subtype of this.\n\nStatic Parameters\n\nTopt\n\n\n\n"
},

{
    "location": "solver/misc.html#ODLCommonTools.AbstractIntegralFunctional",
    "page": "Assorted Function and Types",
    "title": "ODLCommonTools.AbstractIntegralFunctional",
    "category": "type",
    "text": "Abstract type for integral objective functions.  Any integral objective   function should be a subtype of this, which is a subtype ofi   AbstractFunctional.\n\nStatic Parameters\n\nTopt\n\n\n\n"
},

{
    "location": "solver/misc.html#Functionals-1",
    "page": "Assorted Function and Types",
    "title": "Functionals",
    "category": "section",
    "text": "AbstractFunctional\nAbstractIntegralFunctional"
},

{
    "location": "solver/SolverCommon.html#",
    "page": "Solver Common",
    "title": "Solver Common",
    "category": "page",
    "text": ""
},

{
    "location": "solver/SolverCommon.html#SolverCommon.call_nlsolver",
    "page": "Solver Common",
    "title": "SolverCommon.call_nlsolver",
    "category": "function",
    "text": "This function takes in the 4 principle object, fully initialized, and calls a nonlinear solver on them, according to the options in the dictionary. The evalResidual function is passed to the nonlinear solver Inputs: mesh: a mesh object sbp: an SBP operator eqn: an equation object opts: options dictionary, used to determine which nonlinear solver to call pmesh: mesh used for calculating preconditioning jacobian in Newton\'s method, default to using mesh if not specified Outputs: none Aliasing restrictions: none (specificaly, mesh and pmesh can be the same object)\n\n\n\n"
},

{
    "location": "solver/SolverCommon.html#SolverCommon.createMeshAndOperator-Tuple{Any,Any}",
    "page": "Solver Common",
    "title": "SolverCommon.createMeshAndOperator",
    "category": "method",
    "text": "Create a SBP operator and a mesh.  This is used by all physics modules   to create the right type of operator and mesh based on the input options.   It is type unstable, but that is ok.\n\nIf the options dictionary specifies a second SBP operator type, a second   mesh and SBP operator will be created and stored in the mesh2 and sbp2\n\nInputs:     opts: options dictonary     dofpernode: number of degrees of freedom on each node\n\nOutputs     sbp : an AbstractSBP     mesh : an AbstractMesh     pmesh : an AbstractMesh, used for preconditioning, may be same object as             mesh     Tsol : DataType that should be used for eqn.q     Tres : DataType that should be used for eqn.res     Tmsh : DataType of mesh.dxidx and friends     mesh_time : time in seconds for creation of mesh (Float64)\n\n\n\n"
},

{
    "location": "solver/SolverCommon.html#SolverCommon.getDataTypes-Tuple{Dict}",
    "page": "Solver Common",
    "title": "SolverCommon.getDataTypes",
    "category": "method",
    "text": "This function determines the datatypes of the elements of the arrays of the   mesh quantities, sbp operator, solution variables and residual.\n\nIf the datatypes cannot be determined, an error is thrown.\n\nInputs:     opts: the options dictionary\n\nOutputs     Tmsh     Tsbp     Tsol     Tres\n\nOptions Keys\n\nrun_type\njac_method\nforce_solution_complex\nforce_mesh_complex\n\n\n\n"
},

{
    "location": "solver/SolverCommon.html#SolverCommon.loadRestartState",
    "page": "Solver Common",
    "title": "SolverCommon.loadRestartState",
    "category": "function",
    "text": "This function is used by all physics modules to load the most recently   saved state when restarting.\n\nInputs\n\nmesh: the mesh\nsbp: AbstractSBP\neqn: AbstractSolutionData, eqn.q_vec is overwritten with the saved state\nopts: options dictionary\n\nThe keys described in the Checkpointer    documentation are used to determine the most recent complete checkpoint.\n\nImplementation notes:      currently pmesh isn\'t used for anything because checkpointing does not      support mesh adaptation.  When this changes, this function will have to      be updated.\n\n\n\n"
},

{
    "location": "solver/SolverCommon.html#SolverCommon.createMesh",
    "page": "Solver Common",
    "title": "SolverCommon.createMesh",
    "category": "function",
    "text": "This function creates the mesh object and, optionally, a second mesh   used for preconditioning\n\nInputs:     opts: the options dictionary     sbp: an SBP operator     sbpface: an SBP face operator     topo: an ElementTopology describing the SBP reference element.  Only           needed for 3D DG, otherwise can be any value     Tmsh: the DataType of the elements of the mesh arrays (dxidx, jac, etc.)     dofpernode: number of degrees of freedom on every node     suffix: suffix added to options dictionary keys that describe the SBP             operator.  See createSBPOperator\n\nAll arguments except opts are typically provided by    createSBPOperator and getDataTypes\n\n\n\n"
},

{
    "location": "solver/SolverCommon.html#SolverCommon.createSBPOperator",
    "page": "Solver Common",
    "title": "SolverCommon.createSBPOperator",
    "category": "function",
    "text": "This function constructs the SBP operator and the associated SBP face   operator, as specified by the options dictionary.  It also determines   the shape_type that PumiInterface uses to describe the SBP operator to   Pumi.\n\nInputs\n\nopts: the options dictionary\nTsbp: the DataType specifying the Tsbp passed to the SBP operator      constructor\nsuffix: this suffix is added to all keys accessed in the options dictionary.        Usually the suffix is either the empty string or an integer.  This        provides a convenient way for the input file to specify several        different SBP operator and have this operator construct them.        Default value is the empty string.\n\nOutputs\n\nsbp: the SBP operator\nsbpface: the SBP face operator\nshape_type: an integer passed to the mesh constructor to describe the             operator\ntopo: in the 3D DG case, an ElementTopology describing the SBP reference       element, otherwise the integer 0.\n\nDG Operator Names\n\nSBPOmega: nodes on the interior of the element only\nSBPGamma: nodes on the faces of the element and the interior (similar           to Lagrange finite elements)\nSBPDiagonalE: operator with diagonal E matrix, with nodes on vertices,               faces, and interior (similar to Lagrange FE)\nSBPDiagonalE2: operator with diagonal E matrix, with nodes on faces                (but not vertices)\nSBPOmega2: attempt at optimized SBP Omega-type operators, probably not            working\nSBPOmega3: SBP Omega-type operator with degree 2p cubature rule for            all degree operators (p=1 and 2 are the same as SBPOmega),            unlike SBPOmega2, not optimized\n\nCG Operator Names\n\nSBPGamma: see above, this operator can be used to CG as well\n\n\n\n"
},

{
    "location": "solver/SolverCommon.html#Solver-Common-1",
    "page": "Solver Common",
    "title": "Solver Common",
    "category": "section",
    "text": "This page describes some functions that are used by all physics modules as part of initialization of a simulationModules = [SolverCommon]"
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
    "category": "type",
    "text": "Subtype of AbstractParamType.\n\nStatic Parameters:\n\nTsol\nTres\nTdim\n\nThis is a container passed to all low level function, useful for storing   miscellaneous parameters or constants\n\n\n\n"
},

{
    "location": "solver/advection/types.html#AdvectionEquationMod.ParamType2",
    "page": "Datatypes",
    "title": "AdvectionEquationMod.ParamType2",
    "category": "type",
    "text": "Convenient alias for all 2D ParamTypes\n\n\n\n"
},

{
    "location": "solver/advection/types.html#AdvectionEquationMod.ParamType3",
    "page": "Datatypes",
    "title": "AdvectionEquationMod.ParamType3",
    "category": "type",
    "text": "Convenient alias for all 3D ParamTypes\n\n\n\n"
},

{
    "location": "solver/advection/types.html#AdvectionEquationMod.AbstractAdvectionData",
    "page": "Datatypes",
    "title": "AdvectionEquationMod.AbstractAdvectionData",
    "category": "type",
    "text": "Direct subtype of AbstractSolutionData, inheriting Tsol and   Tres as static parameter\n\n\n\n"
},

{
    "location": "solver/advection/types.html#AdvectionEquationMod.AdvectionData",
    "page": "Datatypes",
    "title": "AdvectionEquationMod.AdvectionData",
    "category": "type",
    "text": "Subtype of AbstractAdvectionData, inheriting its static parameters   and adding Tdim.\n\n\n\n"
},

{
    "location": "solver/advection/types.html#AdvectionEquationMod.AdvectionData_",
    "page": "Datatypes",
    "title": "AdvectionEquationMod.AdvectionData_",
    "category": "type",
    "text": "AdvectionEquationMod.AdvectionData_\n\nThis type is an implementation of the abstract AdvectionData.  It is   parameterized by Tsol, the datatype of the solution variables and Tmsh,   the datatype of the mesh variables.   Tres is the \'maximum\' type of Tsol and Tmsh.   Tdim is the dimensionality of the equation being solve (2 or 3).\n\nThis type is (ultimately) a subtype of AbstractSolutionData    and contains all the required fields.\n\n\n\n"
},

{
    "location": "solver/advection/types.html#AdvectionEquationMod.PhysicsName",
    "page": "Datatypes",
    "title": "AdvectionEquationMod.PhysicsName",
    "category": "constant",
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
    "category": "function",
    "text": "AdvectionEquationMod.evalVolumeIntegrals\n\nEvaluates the volume integrals of the weak form.  eqn.res is updated   with the result.  Both the precompute_volume_flux and non    precompute_volume_flux versions are contained within this function\n\nInputs:\n\nmesh : mesh type    sbp  : Summation-by-parts operator    eqn  : Advection equation object    opts : options dictionary\n\nOutputs\n\nNone\n\n\n\n"
},

{
    "location": "solver/advection/volume.html#AdvectionEquationMod.calcAdvectionFlux",
    "page": "Volume Integrals",
    "title": "AdvectionEquationMod.calcAdvectionFlux",
    "category": "function",
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
    "category": "function",
    "text": "AdvectionEquationMod.evalFaceIntegrals\n\nThis function evaluates the interior face integrals for DG methods, using   the flux function from eqn.flux_func.  The solution variables are interpolated   to the faces, the flux computed, and then interpolated back to the   solution points.\n\nThis function also logs some quantities to disk (TODO: move this to   Utils/logging)\n\nInputs:     mesh:  an AbstractDGMesh     sbp     eqn     opts\n\n\n\n"
},

{
    "location": "solver/advection/flux.html#AdvectionEquationMod.calcFaceFlux",
    "page": "Face Integrals",
    "title": "AdvectionEquationMod.calcFaceFlux",
    "category": "function",
    "text": "This function calculates the DG flux between a specified set of faces,   using the solution data at the faces stored in eqn.q_face.   Note that the flux is negated because the face integrals have a    negative sign in the weak form.\n\nInputs:\n\nmesh\nsbp\neqn\nfunctor: the functor that calculates the flux at a node\ninterfaces: an array of type Interface that specifies which interfaces\n            to calculate the flux for\n\nInputs/Outputs:\n\nface_flux: array to store the flux in, numDofPerNode x nnodesPerFace\n           x length(interfaces)\n\nThe functor must have the signature:   func( uL, qR, alpha_x, alpha_y, dxidx, nrm, params)\n\nwhere uL and uR are the solution values for a node on the left and right   elements, alpha_x and alpha_y are the x and y advection velocities,   dxidx is the scaled mapping jacobian for elementL, and nrm is the face   normal in reference space.  params is eqn.params\n\n\n\n"
},

{
    "location": "solver/advection/flux.html#AdvectionEquationMod.calcFaceIntegrals_nopre",
    "page": "Face Integrals",
    "title": "AdvectionEquationMod.calcFaceIntegrals_nopre",
    "category": "function",
    "text": "Compute the face integrals without using eqn.q_face or eqn.flux_face.   The integral is computed directly and res is updated\n\nInputs:\n\nmesh\nsbp\neqn\nopts\nflux_func: the flux functor that computes the face flux at a node\n\n\n\n"
},

{
    "location": "solver/advection/flux.html#AdvectionEquationMod.calcSharedFaceIntegrals",
    "page": "Face Integrals",
    "title": "AdvectionEquationMod.calcSharedFaceIntegrals",
    "category": "function",
    "text": "Thin wrapper around calcSharedFaceIntegrals_inner.  This function is passed   to finishDataExchange, and internally calls calcSharedFaceIntegrals_inner.   See finishExchangeData for details on the interface and    calcSharedFaceIntegrals_inner for the integral that is computed.\n\n\n\n"
},

{
    "location": "solver/advection/flux.html#AdvectionEquationMod.calcSharedFaceIntegrals_inner",
    "page": "Face Integrals",
    "title": "AdvectionEquationMod.calcSharedFaceIntegrals_inner",
    "category": "function",
    "text": "AdvectionEquationMod.calcSharedFaceIntegrals\n\nThis function calculates the shared face integrals for the faces shared   with a single peer process.  This function is for   opts[\"parallel_type\"] == \"face\" and regular face integrals (ie. not the   entropy-stable face integrals) only.\n\nInputs:\n\nmesh\nsbp\neqn\nopts:\ndata: a SharedFaceData specifying which faces to compute\nfunctor: the FluxType to use for the face flux\n\n\n\n"
},

{
    "location": "solver/advection/flux.html#AdvectionEquationMod.calcSharedFaceIntegrals_inner_nopre",
    "page": "Face Integrals",
    "title": "AdvectionEquationMod.calcSharedFaceIntegrals_inner_nopre",
    "category": "function",
    "text": "Like calcSharedFaceIntegrals_inner_nopre, but it computes the integral one   face at a time rather than computing all the flux, storing it in   eqn.flux_sharedface and then doing the integral\n\nSee calcSharedFaceIntegrals_inner for a description of the arguments\n\n\n\n"
},

{
    "location": "solver/advection/flux.html#AdvectionEquationMod.calcSharedFaceIntegrals_element",
    "page": "Face Integrals",
    "title": "AdvectionEquationMod.calcSharedFaceIntegrals_element",
    "category": "function",
    "text": "Thin wrapper around calcSharedFaceIntegrals_inner.  This function is passed   to finishDataExchange, and internally calls calcSharedFaceIntegrals_inner.   See finishDataExchange for details on the interface and    calcSharedFaceIntegrals_inner for the integral that is computed.\n\n\n\n"
},

{
    "location": "solver/advection/flux.html#AdvectionEquationMod.calcSharedFaceIntegrals_element_inner",
    "page": "Face Integrals",
    "title": "AdvectionEquationMod.calcSharedFaceIntegrals_element_inner",
    "category": "function",
    "text": "Like calcSharedFaceIntegrals_inner, but for the case when   opts[\"parallel_data\"] == element.  This effectively means it has to   interpolate the solution from the elements to the faces and then do the   integral\n\n\n\n"
},

{
    "location": "solver/advection/flux.html#AdvectionEquationMod.calcSharedFaceIntegrals_element_inner_nopre",
    "page": "Face Integrals",
    "title": "AdvectionEquationMod.calcSharedFaceIntegrals_element_inner_nopre",
    "category": "function",
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
    "category": "function",
    "text": "Evaluate boundary integrals for advection equation, updating eqn.res with the result.\n\nInputs\n\nmesh : Abstract mesh type\nsbp  : Summation-by-parts operator\neqn  : Advection equation object\nopts : options dictionary\n\nOutputs\n\nNone\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.calcBoundaryFlux",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.calcBoundaryFlux",
    "category": "function",
    "text": "This function calculates the boundary flux for the portion of the boundary   with a particular boundary condition.  The eqn.q are converted to    conservative variables if needed.  For the DG version, eqn.q_bndry must   already be populated with the q variables interpolated to the boundary\n\nInputs:\n\nmesh : AbstractMesh   sbp : AbstractSBP   eqn : AdvectionEquation   functor : a callable object that calculates the boundary flux at a node   idx_range: the Range describing which Boundaries have the current BC   bndry_facenums:  An array with elements of type Boundary that tell which                    element faces have the boundary condition   Outputs:\n\nbndryflux : the array to store the boundary flux, corresponds to                bndry_facenums\n\nnote that bndry_facenums and bndryflux must be only the portion of the    their parent arrays that correspond to the Boundaries that have the    current boundary condition applied.\n\nThe functor must have the signature:   functor( q, aux_vars, x, dxidx, nrm, bndryflux_i, eqn.params)   where q are the conservative variables.   where all arguments (except params and nrm) are vectors of values at a node.\n\nparams is the ParamType associated with the the EulerEquation object   nrm = mesh.sbpface.normal[:, current_node]\n\nThis is a mid level function.\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.calcBoundaryFlux_nopre",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.calcBoundaryFlux_nopre",
    "category": "function",
    "text": "This function computes the boundary integrals (and should probably be renamed)   without using eqn.q_bndry of eqn.bndryflux.  eqn.res is updated with   the results.\n\nSee calcBoundaryFlux for the meaning of the arguments\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.BCDict",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.BCDict",
    "category": "constant",
    "text": "AdvectionEquationMod.BCDict\n\nIt stores all the possible boundary condition dictionary options. Whenever a  new boundary condition is created, it should get added to BCDict.\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.getBCFunctors",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.getBCFunctors",
    "category": "function",
    "text": "AdvectionEquationMod.getBCFunctors\n\nThis function uses the opts dictionary to populate mesh.bndry_funcs with the the functors\n\nThis is a high level function.\n\nInputs\n\nmesh : Abstract mesh type\nsbp  : Summation-by-parts operator\neqn  : Advection equation object\nopts : Input dictionary options\n\nOutputs\n\nNone\n\n\n\n"
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
    "category": "type",
    "text": "Default BC to calculate the boundary face integral (no numerical flux   functions)\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.exp2xplus2yBC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.exp2xplus2yBC",
    "category": "type",
    "text": "AdvectionEquationMod.exp2xplus2yBC\n\nUses the Roe solver to calculate the boundary flux using calc_exp2xplus2y to get the boundary state.\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.exp3xplusyBC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.exp3xplusyBC",
    "category": "type",
    "text": "AdvectionEquationMod.exp3xplusyBC\n\nUses the Roe solver to calculate the boundary flux using calc_exp3xplusy to get the boundary state.\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.exp5xplus4yplus2BC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.exp5xplus4yplus2BC",
    "category": "type",
    "text": "AdvectionEquationMod.exp5xplus4yplus2BC\n\nUses the Roe solver to calculate the boundary flux using calc_exp5xplus4yplus2  to get the boundary state.\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.exp5xplusyBC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.exp5xplusyBC",
    "category": "type",
    "text": "AdvectionEquationMod.exp5xplusyBC\n\nUses the Roe solver to calculate the boundary flux using calc_exp5xplusy to get the boundary state.\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.exp_xplusyBC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.exp_xplusyBC",
    "category": "type",
    "text": "AdvectionEquationMod.exp_xplusyBC\n\nCalculates q at the boundary which is equal to exp(x+y). It is a nodal  level function.\n\nInputs\n\nu : Advection variable (eqn.q)\nalpha_x & alpha_y : velocities in the X & Y directions\ncoords : Nodal coordinates\ndxidx  : Mapping Jacobian\nnrm    : SBP face-normal vectors\nbndryflux : Flux at the boundary\n\nOutputs\n\nNone\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.exp_xyBC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.exp_xyBC",
    "category": "type",
    "text": "AdvectionEquationMod.exp_xyBC\n\nUses the Roe solver to calculate the boundary flux using calc_exp_xy to get the boundary state\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.mms1BC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.mms1BC",
    "category": "type",
    "text": "AdvectionEquationMod.mms1BC\n\nUses the Roe solver to calculate the boundary flux using calc_mms1 to get   the boundary state.\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.p0BC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.p0BC",
    "category": "type",
    "text": "AdvectionEquationMod.p0BC\n\nUses the Roe solver to calculate the boundary flux using calc_p0 to   get the boundary state\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.p1BC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.p1BC",
    "category": "type",
    "text": "AdvectionEquationMod.p1BC\n\nUses the Roe solver to calculate the boundary flux using calc_p1 to   get the boundary state\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.p2BC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.p2BC",
    "category": "type",
    "text": "AdvectionEquationMod.p2BC\n\nUses the Roe solver to calculate the boundary flux using calc_p2 to   get the boundary state\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.p3BC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.p3BC",
    "category": "type",
    "text": "AdvectionEquationMod.p3BC\n\nUses the Roe solver to calculate the boundary flux using calc_p3 to   get the boundary state\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.p4BC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.p4BC",
    "category": "type",
    "text": "AdvectionEquationMod.p4BC\n\nUses the Roe solver to calculate the boundary flux using calc_p4 to   get the boundary state.\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.p5BC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.p5BC",
    "category": "type",
    "text": "AdvectionEquationMod.p5BC\n\nUses the Roe solver to calculate the boundary flux using calc_p5 to   get the boundary state.\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.sinwave_BC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.sinwave_BC",
    "category": "type",
    "text": "AdvectionEquationMod.sinwave_BC\n\nUses the Roe solver to calculate the boundary flux using calc_sinewave to   get the boundary state\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.sinwave_ampl_BC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.sinwave_ampl_BC",
    "category": "type",
    "text": "AdvectionEquationMod.sinwave_ampl_BC\n\nUses the Roe solver to calculate the boundary flux using calc_sinewave_ampl to   get the boundary state\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.sinwavey_BC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.sinwavey_BC",
    "category": "type",
    "text": "AdvectionEquationMod.sinwavey_BC\n\nUses the Roe solver to calculate the boundary flux using calc_sinewavey to   get the boundary state\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.sinwavey_pertBC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.sinwavey_pertBC",
    "category": "type",
    "text": "AdvectionEquationMod.sinwavey_pertBC\n\nUses the Roe solver to calculate the boundary flux using calc_sinewave_pert to   get the boundary state\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.unsteadymmsBC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.unsteadymmsBC",
    "category": "type",
    "text": "BC for unsteadymms\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.unsteadypolyBC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.unsteadypolyBC",
    "category": "type",
    "text": "BC for unsteadypoly\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.x4BC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.x4BC",
    "category": "type",
    "text": "AdvectionEquationMod.x4BC\n\nUses the Roe solver to calculate the boundary flux using calc_x4 to   get the boundary state.\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.x5plusy5BC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.x5plusy5BC",
    "category": "type",
    "text": "AdvectionEquationMod.x5plusy5BC\n\nCalculates q at the boundary which is equal to x^5 + y^5. It is a nodal  level function.\n\nInputs\n\nu : Advection variable (eqn.q)\nparams: the equation ParamType\ncoords : Nodal coordinates\nnrm_scaled    : scaled face normal vector in x-y space\nt:  current time value\n\nOutputs\n\nbndryflux : Flux at the boundary\nNone\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.xplusyBC",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.xplusyBC",
    "category": "type",
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
    "location": "solver/advection/bc.html#AdvectionEquationMod.RoeSolver-Union{Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol,Tsol,Any,Any}, Tuple{Tsol}} where Tsol",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.RoeSolver",
    "category": "method",
    "text": "AdvectionEquationMod.RoeSolver\n\nRoe solver for the advection equations. It determines the boundary flux on  each boundary. It is called at the nodal level\n\nInputs\n\nu    : Solution of advection equation at a particular node\nu_bc : Prescribed solution value at the boundary\nparams: the equation ParamType object\nnrm  : scaled face normal vector in x-y space\n\nOutputs\n\nbndryflux : Boundary flux at the particular node\n\n\n\n"
},

{
    "location": "solver/advection/bc.html#AdvectionEquationMod.flux1-Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol,Any,Any,Any,Any}",
    "page": "Boundary Integrals",
    "title": "AdvectionEquationMod.flux1",
    "category": "method",
    "text": "flux1\n\nCalculates the boundary flux for the advection equation. It works at the nodal level.\n\nInputs\n\nu_sbp_: The entry from u_sbp for this node\ndxidx : The jacobian for this node\nnrm   : nrm is the normal vector\nnet_flux:\nparams: the equation ParamType\n\nOutputs\n\nNone\n\n\n\n"
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
    "category": "constant",
    "text": "Dictionary that maps IC names to functions.  Every new IC should be added   to the list\n\n\n\n"
},

{
    "location": "solver/advection/ic.html#AdvectionEquationMod.ICFile-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP{Tsbp},AdvectionEquationMod.AdvectionData{Tsol,Tres,Tdim} where Tdim where Tres,Any,AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsbp}, Tuple{Tsol}} where Tsol where Tsbp where Tmsh",
    "page": "Initial Condition",
    "title": "AdvectionEquationMod.ICFile",
    "category": "method",
    "text": "AdvectionEquationMod.ICFile\n\nThis function reads a vector from a file on disk and set the solution to it. The vector must contain the same number of entries as there are degrees of  freedom in the mesh. \n\nThis function is useful for things like restarting from a checkpoint. In this case, the file should be the output of writedlm(eqn.q).  The degree  of freedom number must be the same for both simulation for this to work (the  file contains no degree of freedom number information).\n\nInputs\n\nmesh\nsbp\neqn\nopts\n\nInputs/Outputs\n\nu0: vector to populate with the solution\nAliasing restrictions: none.\n\n\n\n"
},

{
    "location": "solver/advection/ic.html#AdvectionEquationMod.ICexp_xplusy-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP{Tsbp},AdvectionEquationMod.AdvectionData{Tsol,Tres,2},Any,AbstractArray{Tsol,N} where N}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsbp}, Tuple{Tsol}} where Tres where Tsol where Tsbp where Tmsh",
    "page": "Initial Condition",
    "title": "AdvectionEquationMod.ICexp_xplusy",
    "category": "method",
    "text": "AdvectionEquationMod.ICexp_xplusy\n\nComputes the initial conditions for the state variable     u = exp(x + y + z), z is omitted in 2d\n\nInputs\n\nmesh : AbstractMesh type\nsbp  : Summation-by-parts operator\neqn  : Advection equation object\nopts : Options dictionary\nu0   : Array that stores inital state variables\n\nOutputs\n\nNone\n\n\n\n"
},

{
    "location": "solver/advection/ic.html#AdvectionEquationMod.ICx5plusy5-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP{Tsbp},AdvectionEquationMod.AdvectionData{Tsol,Tres,2},Any,AbstractArray{Tsol,N} where N}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsbp}, Tuple{Tsol}} where Tres where Tsol where Tsbp where Tmsh",
    "page": "Initial Condition",
    "title": "AdvectionEquationMod.ICx5plusy5",
    "category": "method",
    "text": "AdvectionEquationMod.ICx5plusy5\n\nComputes the initial conditions for the state variable     u = x^5 + y^5 + z^5; z is omitted it 2d\n\nInputs\n\nmesh : AbstractMesh type\nsbp  : Summation-by-parts operator\neqn  : Advection equation object\nopts : Options dictionary\nu0   : Array that stores inital state variables\n\nOutputs\n\nNone\n\n\n\n"
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
    "location": "solver/advection/source.html#AdvectionEquationMod.getSRCFunctors-Tuple{ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,AdvectionEquationMod.AdvectionData,Any}",
    "page": "Source Term",
    "title": "AdvectionEquationMod.getSRCFunctors",
    "category": "method",
    "text": "This function gets the functor specified by opts[\"SRCname\"] and stores   it to the equation object.  Currently one 1 source functor is allowed.\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCDict",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCDict",
    "category": "constant",
    "text": "AdvectionEquationMod.SRCDict\n\nIt stores all the possible source term dictionary options. Whenever a    new source is created, it should get added to SRCDict.\n\nAll functors must have the signature:\n\nsrc_func(params, coords::ParamType, t)\n\nwhere coords is the vector of length 2 containing the x and y coordinates   of the node, params.alpha_x and params.alpha_y are the advection velocities in the x an y   directions, and t is the current time\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRC0",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRC0",
    "category": "type",
    "text": "AdvectionEquationMod.SRC0\n\nThis is the zero source term.  This is the default of source term   is specified\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRC1",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRC1",
    "category": "type",
    "text": "AdvectionEquationMod.SRC1\n\nThis source term returns 1 everywhere.\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRC2",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRC2",
    "category": "type",
    "text": "AdvectionEquationMod.SRC2\n\nThis source term returns 1 everywhere.\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCexp2xplus2y",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCexp2xplus2y",
    "category": "type",
    "text": "AdvectionEquationMod.SRCexp2xplus2y\n\nCalculates the source term for q = exp(2x + 2y)\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCexp3xplusy",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCexp3xplusy",
    "category": "type",
    "text": "AdvectionEquationMod.SRCexp3xplusy\n\nCalculates the source term for q = exp(3*x + y)\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCexp5xplus4yplus2",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCexp5xplus4yplus2",
    "category": "type",
    "text": "AdvectionEquationMod.SRCexp5xplus4yplus2\n\nCalculates the source term for q = exp(5x + 4y +2)\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCexp5xplusy",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCexp5xplusy",
    "category": "type",
    "text": "AdvectionEquationMod.SRCexp5xplusy\n\nCalculates the source term for q = exp(5*x + y)\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCexp_xplusy",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCexp_xplusy",
    "category": "type",
    "text": "AdvectionEquationMod.SRCexp_xplusy\n\nThis is a source term that returns a source term for e^(x+y)\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCexp_xy",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCexp_xy",
    "category": "type",
    "text": "AdvectionEquationMod.SRCexp_xy\n\nCalculates the source term for q = exp(x*y)\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCmms1",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCmms1",
    "category": "type",
    "text": "AdvectionEquationMod.SRC1\n\nThis source term that returns: the derivative of mms1\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCp1",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCp1",
    "category": "type",
    "text": "AdvectionEquationMod.SRCp1\n\nThis source term that returns: the source term for a manufactured solution   using a 1st order polynomial\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCp2",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCp2",
    "category": "type",
    "text": "AdvectionEquationMod.SRCp3\n\nThis source term that returns: the source term for a manufactured solution   using a 2nd order polynomial\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCp3",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCp3",
    "category": "type",
    "text": "AdvectionEquationMod.SRCp3\n\nThis source term that returns: the source term for a manufactured solution   using a 3rd order polynomial\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCp4",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCp4",
    "category": "type",
    "text": "AdvectionEquationMod.SRCp4\n\nThis source term that returns: the source term for a manufactured solution   using a 4th order polynomial\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCp5",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCp5",
    "category": "type",
    "text": "AdvectionEquationMod.SRCp5\n\nThis source term that returns: the source term for a manufactured solution   using a 5th order polynomial\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCunsteadymms",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCunsteadymms",
    "category": "type",
    "text": "Source term for unsteady mms\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCunsteadypoly",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCunsteadypoly",
    "category": "type",
    "text": "Source term for unsteady poly\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCx",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCx",
    "category": "type",
    "text": "AdvectionEquationMod.SRCx\n\nThis source term that returns: f = x\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCx4",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCx4",
    "category": "type",
    "text": "AdvectionEquationMod.SRCx4\n\nThis source term that returns: the source term for a manufactured solution   using a 4th order polynomial\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCx5plusy5",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCx5plusy5",
    "category": "type",
    "text": "AdvectionEquationMod.SRCx5plusy5\n\nThis is a source term that returns a source term for e^(x+y)\n\n\n\n"
},

{
    "location": "solver/advection/source.html#AdvectionEquationMod.SRCxplusy",
    "page": "Source Term",
    "title": "AdvectionEquationMod.SRCxplusy",
    "category": "type",
    "text": "AdvectionEquationMod.SRCxplusy\n\ncalculates the source term for q = x + y\n\n\n\n"
},

{
    "location": "solver/advection/source.html#Advection-Source-Term-1",
    "page": "Source Term",
    "title": "Advection Source Term",
    "category": "section",
    "text": "This pages describes the functions that apply source terms  CurrentModule = AdvectionEquationMod  Modules = [AdvectionEquationMod]\n  Pages = [\"advection/source.jl\"]\n  Order = [:function, :constant, :type, :module, :macro]"
},

{
    "location": "solver/advection/common.html#",
    "page": "Common Functions",
    "title": "Common Functions",
    "category": "page",
    "text": ""
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_exp2xplus2y-Union{Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol,AbstractArray{Tmsh,N} where N,Any}, Tuple{Tmsh}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_exp2xplus2y",
    "category": "method",
    "text": "AdvectionEquationMod.calc_exp2xplus2y\n\nCalculates and return the expression u = exp(2x + 2y)\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_exp3xplusy-Union{Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol,AbstractArray{Tmsh,N} where N,Any}, Tuple{Tmsh}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_exp3xplusy",
    "category": "method",
    "text": "AdvectionEquationMod.calc_exp3xplusy\n\nCalculates and return the expression u = exp(3*x + y)\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_exp5xplus4yplus2-Union{Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol,AbstractArray{Tmsh,N} where N,Any}, Tuple{Tmsh}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_exp5xplus4yplus2",
    "category": "method",
    "text": "AdvectionEquationMod.calc_exp5xplus4yplus2\n\nCalculates and returns the expression u = exp(5x + 4y +2)\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_exp5xplusy-Union{Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol,AbstractArray{Tmsh,N} where N,Any}, Tuple{Tmsh}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_exp5xplusy",
    "category": "method",
    "text": "AdvectionEquationMod.calc_exp5xplusy\n\nCalculates and return the expression u = exp(5*x + y)\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_exp_xplusy-Union{Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol,AbstractArray{Tmsh,N} where N,Any}, Tuple{Tmsh}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_exp_xplusy",
    "category": "method",
    "text": "AdvectionEquationMod.calc_exp_xplusy\n\nCalculates and returns e^(x + y)\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_exp_xy-Union{Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol,AbstractArray{Tmsh,N} where N,Any}, Tuple{Tmsh}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_exp_xy",
    "category": "method",
    "text": "AdvectionEquationMod.calc_exp_xy\n\nCalculates and returns u = exp(x*y)\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_mms1-Union{Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol,AbstractArray{Tmsh,N} where N,Any}, Tuple{Tmsh}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_mms1",
    "category": "method",
    "text": "AdvectionEquationMod.calc_mms1\n\nCalculates and returns the value of the solution for doing Method of    Manufactured solutions.  This is for debugging only, and could change   at any time.\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_mms1dx-Union{Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol,AbstractArray{Tmsh,N} where N,Any}, Tuple{Tmsh}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_mms1dx",
    "category": "method",
    "text": "AdvectionEquationMod.calc_mms1dx\n\nCalculates and returns the x derivative of calc_mms1\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p0-Union{Tuple{Tmsh}, Tuple{Union{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol, AdvectionEquationMod.ParamType{Tsol,Tres,3} where Tres where Tsol},AbstractArray{Tmsh,N} where N,Any}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p0",
    "category": "method",
    "text": "calculates and returns a degree 0 polynomial (ie. a constant)\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p1-Union{Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol,AbstractArray{Tmsh,N} where N,Any}, Tuple{Tmsh}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p1",
    "category": "method",
    "text": "AdvectionEquationMod.calc_p1\n\nCalculates and returns a 1st order polynomial of x and y\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p1dx-Union{Tuple{Tmsh}, Tuple{Union{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol, AdvectionEquationMod.ParamType{Tsol,Tres,3} where Tres where Tsol},AbstractArray{Tmsh,N} where N,Any}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p1dx",
    "category": "method",
    "text": "AdvectionEquationMod.calc_p1dx\n\nCalculates and returns the x derivative of calc_p1\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p1dy-Union{Tuple{Tmsh}, Tuple{Union{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol, AdvectionEquationMod.ParamType{Tsol,Tres,3} where Tres where Tsol},AbstractArray{Tmsh,N} where N,Any}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p1dy",
    "category": "method",
    "text": "AdvectionEquationMod.calc_p1dy\n\nCalculates and returns the y derivative of calc_p1\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p1dz-Union{Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,3} where Tres where Tsol,AbstractArray{Tmsh,N} where N,Any}, Tuple{Tmsh}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p1dz",
    "category": "method",
    "text": "AdvectionEquationMod.calc_p1dz\n\nCalculates and returns the z derivative of calc_p1\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p2-Union{Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol,AbstractArray{Tmsh,N} where N,Any}, Tuple{Tmsh}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p2",
    "category": "method",
    "text": "AdvectionEquationMod.calc_p2\n\nCalculates and returns a 2nd order polynomial in x and y\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p2dx-Union{Tuple{Tmsh}, Tuple{Union{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol, AdvectionEquationMod.ParamType{Tsol,Tres,3} where Tres where Tsol},AbstractArray{Tmsh,N} where N,Any}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p2dx",
    "category": "method",
    "text": "AdvectionEquationMod.calc_p2dx\n\nCalculates and returns a the x derivative of calc_p2\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p2dy-Union{Tuple{Tmsh}, Tuple{Union{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol, AdvectionEquationMod.ParamType{Tsol,Tres,3} where Tres where Tsol},AbstractArray{Tmsh,N} where N,Any}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p2dy",
    "category": "method",
    "text": "AdvectionEquationMod.calc_p2dy\n\nCalculates and returns the y derivative of calc_p2\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p2dz-Union{Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,3} where Tres where Tsol,AbstractArray{Tmsh,N} where N,Any}, Tuple{Tmsh}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p2dz",
    "category": "method",
    "text": "AdvectionEquationMod.calc_p2dz\n\nCalculates and returns the z derivative of calc_p2\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p3-Union{Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol,AbstractArray{Tmsh,N} where N,Any}, Tuple{Tmsh}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p3",
    "category": "method",
    "text": "AdvectionEquationMod.calc_p3\n\nCalculates and returns a 3rd order polynomial in x and y (and z in 3d)\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p3dx-Union{Tuple{Tmsh}, Tuple{Union{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol, AdvectionEquationMod.ParamType{Tsol,Tres,3} where Tres where Tsol},AbstractArray{Tmsh,N} where N,Any}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p3dx",
    "category": "method",
    "text": "AdvectionEquationMod.calc_p3dx\n\nCalculates and returns the x derivataive of calc_p3\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p3dy-Union{Tuple{Tmsh}, Tuple{Union{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol, AdvectionEquationMod.ParamType{Tsol,Tres,3} where Tres where Tsol},AbstractArray{Tmsh,N} where N,Any}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p3dy",
    "category": "method",
    "text": "AdvectionEquationMod.calc_p3dy\n\nCalculates and returns the y derivative of calc_p3\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p3dz-Union{Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,3} where Tres where Tsol,AbstractArray{Tmsh,N} where N,Any}, Tuple{Tmsh}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p3dz",
    "category": "method",
    "text": "AdvectionEquationMod.calc_p3dz\n\nCalculates and returns the z derivative of calc_p3\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p4-Union{Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol,AbstractArray{Tmsh,N} where N,Any}, Tuple{Tmsh}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p4",
    "category": "method",
    "text": "AdvectionEquationMod.calc_p4\n\nCalculates and returns a 4th order polynomial in x and y\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p4dx-Union{Tuple{Tmsh}, Tuple{Union{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol, AdvectionEquationMod.ParamType{Tsol,Tres,3} where Tres where Tsol},AbstractArray{Tmsh,N} where N,Any}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p4dx",
    "category": "method",
    "text": "AdvectionEquationMod.calc_p4x\n\nCalculates and returns the x derivative of calc_p4\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p4dy-Union{Tuple{Tmsh}, Tuple{Union{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol, AdvectionEquationMod.ParamType{Tsol,Tres,3} where Tres where Tsol},AbstractArray{Tmsh,N} where N,Any}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p4dy",
    "category": "method",
    "text": "AdvectionEquationMod.calc_p4dy\n\nCalculates and returns the y derivative of calc_p4\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p4dz-Union{Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,3} where Tres where Tsol,AbstractArray{Tmsh,N} where N,Any}, Tuple{Tmsh}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p4dz",
    "category": "method",
    "text": "AdvectionEquationMod.calc_p4dz\n\nCalculates and returns the z derivative of calc_p4\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p5-Union{Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol,AbstractArray{Tmsh,N} where N,Any}, Tuple{Tmsh}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p5",
    "category": "method",
    "text": "AdvectionEquationMod.calc_p5\n\nCalculates and returns a 5th order polynomial in x and y (and z in 3d)\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p5dx-Union{Tuple{Tmsh}, Tuple{Union{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol, AdvectionEquationMod.ParamType{Tsol,Tres,3} where Tres where Tsol},AbstractArray{Tmsh,N} where N,Any}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p5dx",
    "category": "method",
    "text": "AdvectionEquationMod.calc_p5dx\n\nCalculates and returns the x derivative of calc_p5\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p5dy-Union{Tuple{Tmsh}, Tuple{Union{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol, AdvectionEquationMod.ParamType{Tsol,Tres,3} where Tres where Tsol},AbstractArray{Tmsh,N} where N,Any}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p5dy",
    "category": "method",
    "text": "AdvectionEquationMod.calc_p5y\n\nCalculates and returns the y derivative of calc_p5\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_p5dz-Union{Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,3} where Tres where Tsol,AbstractArray{Tmsh,N} where N,Any}, Tuple{Tmsh}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_p5dz",
    "category": "method",
    "text": "AdvectionEquationMod.calc_p5z\n\nCalculates and returns the z derivative of calc_p5\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_sinwave-Union{Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol,AbstractArray{Tmsh,N} where N,Any}, Tuple{Tmsh}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_sinwave",
    "category": "method",
    "text": "AdvectionEquationMod.calc_sinwave\n\nCalculates and returns sin(-x + t)\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_sinwave_ampl-Union{Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol,AbstractArray{Tmsh,N} where N,Any}, Tuple{Tmsh}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_sinwave_ampl",
    "category": "method",
    "text": "AdvectionEquationMod.calc_sinwave_ampl\n\nCalculates and returns Asin(-x + omegat)\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_sinwavey-Union{Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol,AbstractArray{Tmsh,N} where N,Any}, Tuple{Tmsh}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_sinwavey",
    "category": "method",
    "text": "AdvectionEquationMod.calc_sinwavey\n\nCalculates and returns sin(y)^2 + 5*sin(y) + 3/sin(y)\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_sinwavey_pert-Union{Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol,AbstractArray{Tmsh,N} where N,Any}, Tuple{Tmsh}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_sinwavey_pert",
    "category": "method",
    "text": "AdvectionEquationMod.calc_sinwavey_pert\n\nCalculates and returns 1000sin(x)calc_sinwavey\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_unsteadymms-Union{Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol,AbstractArray{Tmsh,N} where N,Any}, Tuple{Tmsh}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_unsteadymms",
    "category": "method",
    "text": "u = exp(x + y + z + t) in 3d (z = 0 in 2d)\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_x4-Union{Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol,AbstractArray{Tmsh,N} where N,Any}, Tuple{Tmsh}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_x4",
    "category": "method",
    "text": "AdvectionEquationMod.calc_x4\n\nCalculates and returns a 4th order polynomial in x\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_x4der-Union{Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol,AbstractArray{Tmsh,N} where N,Any}, Tuple{Tmsh}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_x4der",
    "category": "method",
    "text": "AdvectionEquationMod.calc_x5plusy5\n\nCalculates and returns the x derivative of calc_x4\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_x5plusy5-Union{Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol,AbstractArray{Tmsh,N} where N,Any}, Tuple{Tmsh}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_x5plusy5",
    "category": "method",
    "text": "AdvectionEquationMod.calc_x5plusy5\n\nCalculates and returns x^5 + y^5\n\n\n\n"
},

{
    "location": "solver/advection/common.html#AdvectionEquationMod.calc_xplusy-Union{Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol,AbstractArray{Tmsh,N} where N,Any}, Tuple{Tmsh}} where Tmsh",
    "page": "Common Functions",
    "title": "AdvectionEquationMod.calc_xplusy",
    "category": "method",
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
    "location": "solver/advection/adjoint.html#Advection-Equation-Steady-Adjoint-1",
    "page": "Adjoint",
    "title": "Advection Equation Steady Adjoint",
    "category": "section",
    "text": "PDESolver currently has the capability to compute the steady adjoint of a boundary functional. Recall the adjoint equation asfracpartial mathcalLpartial q = fracpartial mathcalJpartial q + psi^T fracpartial mathcalRpartial q = 0where, mathcalL is the Lagrangian for functional mathcalJ and q is the solution variable. The adjoint can be computed by calling the function calcAdjoint, which has been described below.  Modules = [AdvectionEquationMod]\n  Pages = [\"solver/advection/adjoint.jl\"]"
},

{
    "location": "solver/advection/boundary_functional.html#",
    "page": "Boundary Functional",
    "title": "Boundary Functional",
    "category": "page",
    "text": ""
},

{
    "location": "solver/advection/boundary_functional.html#AdvectionEquationMod.calcBndryFunctional-Union{Tuple{ODLCommonTools.AbstractDGMesh{Tmsh},SummationByParts.AbstractSBP,AdvectionEquationMod.AdvectionData{Tsol,Tres,Tdim} where Tdim where Tres,Any,ODLCommonTools.AbstractIntegralFunctional{Topt}}, Tuple{Tmsh}, Tuple{Topt}, Tuple{Tsol}} where Topt where Tsol where Tmsh",
    "page": "Boundary Functional",
    "title": "AdvectionEquationMod.calcBndryFunctional",
    "category": "method",
    "text": "AdvectionEquationMod.calcBndryfunctional\n\nThis function calculates the functional on a geometric boundary of a the computational space. This is a mid level function that should not be called from outside the module. Depending on the functional being computed, it may be necessary to define another method for this function based on a different boundary functional type.\n\nArguments\n\nmesh :  Abstract mesh object\nsbp  : Summation-By-Parts operator\neqn  : Advection equation object\nopts : Options dictionary\nfunctionalData : Object of the functional being computed\n\n\n\n"
},

{
    "location": "solver/advection/boundary_functional.html#AdvectionEquationMod.calcBoundaryFunctionalIntegrand-Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol,Any,Any,Any,AdvectionEquationMod.IntegralQData}",
    "page": "Boundary Functional",
    "title": "AdvectionEquationMod.calcBoundaryFunctionalIntegrand",
    "category": "method",
    "text": "Method for IntegralQData functional.\n\nNote that q should be scaled by the length of the normal vector so the   integration works correctly\n\n\n\n"
},

{
    "location": "solver/advection/boundary_functional.html#AdvectionEquationMod.calcBoundaryFunctionalIntegrand-Tuple{AdvectionEquationMod.ParamType{Tsol,Tres,2} where Tres where Tsol,Any,Any,Any,AdvectionEquationMod.QfluxData}",
    "page": "Boundary Functional",
    "title": "AdvectionEquationMod.calcBoundaryFunctionalIntegrand",
    "category": "method",
    "text": "AdvectionEquationMod.calcBoundaryFunctionalIntegrand\n\nComputes the integrand for boundary functional at a surface SBP node. Every functional needs to have its own method and the functional type determines which method is called.\n\nInputs\n\nparams : eqn.params object\nnx : X component of face normal vector\nny : Y component of face normal vector\nq  : Nodal solution variable\nfunctionalData : Object of the functional being computed\n\nOutputs\n\nfunctional_integrand : Computed integrand at the surface node\n\n\n\n"
},

{
    "location": "solver/advection/boundary_functional.html#PDESolver.evalFunctional-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP,AdvectionEquationMod.AdvectionData{Tsol,Tres,Tdim} where Tdim where Tres,Any,ODLCommonTools.AbstractFunctional}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Boundary Functional",
    "title": "PDESolver.evalFunctional",
    "category": "method",
    "text": "AdvectionEquationMod.evalFunctional\n\nHight level function that evaluates functionals specified in the options dictionary. The user must call this function for functional evaluation.This function is agnostic which type of a functional is being computed and calls a mid level type specific function for the actual functional evaluation.\n\nArguments\n\nmesh :  Abstract mesh object\nsbp  : Summation-By-Parts operator\neqn  : Advection equation object\nopts : Options dictionary\nfunctionalData : Object of the functional being computed.\n\n\n\n"
},

{
    "location": "solver/advection/boundary_functional.html#Advection-Boundary-Functional-1",
    "page": "Boundary Functional",
    "title": "Advection Boundary Functional",
    "category": "section",
    "text": "  CurrentModule = AdvectionEquationModThis page consists of all the functions necessary for computing a boundary functional along the geometric edges of a mesh for the advection equation.  Modules = [AdvectionEquationMod]\n  Pages = [\"solver/advection/boundary_functional.jl\"]"
},

{
    "location": "solver/advection/boundary_functional.html#AdvectionEquationMod.QfluxData",
    "page": "Boundary Functional",
    "title": "AdvectionEquationMod.QfluxData",
    "category": "type",
    "text": "###AdvectionEquationMod.QfluxData\n\nData type for storing relevant information pertaining to an a functional or an objective function.\n\nMembers\n\nbcnums : boundary condition groups on which the functional is to be                            computed.\nval : Computed value of the functional\ntarget_qFlux : Target value for the functional qFlux\n\n\n\n"
},

{
    "location": "solver/advection/boundary_functional.html#AdvectionEquationMod.IntegralQData",
    "page": "Boundary Functional",
    "title": "AdvectionEquationMod.IntegralQData",
    "category": "type",
    "text": "Functional that integrates the solution q over the specified boundary(/ies)\n\nFields\n\nbcnums: the boundary condition groups the functional is computed over\nval: the value of the functional, initially 0.0\n\n\n\n"
},

{
    "location": "solver/advection/boundary_functional.html#Functional-Types-1",
    "page": "Boundary Functional",
    "title": "Functional Types",
    "category": "section",
    "text": "QfluxData\nIntegralQData"
},

{
    "location": "solver/euler/euler.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "solver/euler/euler.html#Euler-Physics-1",
    "page": "Introduction",
    "title": "Euler Physics",
    "category": "section",
    "text": "Describe the equation being solved herefracpartial qpartial t = - nabla cdot F(q) + SWhere q are the conservative variables and F is the Euler flux.q = beginbmatrix rho rho u rho v rho w e endbmatrixwhere rho is density, u is the x velocity, v is the y velocity, w is the z velocity and e is the energy.The calloricaly perfect ideal gas law is used to close the system (TODO ref calcPressure)The x-y-z components of the Euler flux are:beginbmatrix rho u  rho v  rho w  rho u^2 + p  rho u v  rho u w  rho u v  rho v^2 + p  rho v w   rho u w  rho v w  rho w^2 + p   (e + p)u  (e + p)v  (e + p)w endbmatrixTODO: describe the physical quantities (include units)"
},

{
    "location": "solver/euler/euler.html#Discretizations-1",
    "page": "Introduction",
    "title": "Discretizations",
    "category": "section",
    "text": "The code currently implements three kinds of discretizations.The first is a standard discontinuous-Galerkin scheme using SBP operators and uses the numerical flux functions in flux functions section.The second scheme is an entropy stable discretization that uses a Hadamard product between the SBP operator matrices and matrices of flux values. The face integrals for the entropy stable scheme are described on the face element integrals page.The third scheme is a staggerd grid approach based on the entropy stable scheme.  It uses all the same mechanics as the entropy stable scheme, but requires interpolating data back and forth between the solution and flux grids. The functions for doing the interpolation are listed on the relevent pages for the volume and face integrals.   Pages = [ \"advection.md\"\n            \"types.md\"\n            \"volume.md\"\n            \"flux.md\"\n            \"faceElementIntegrals.md\"\n            \"bc.md\"\n            \"ic.md\"\n            \"source.md\"\n            \"common.md\"\n            \"conversion.md\"\n            \"flux_functions.md\"\n            \"stabilization.md\"\n            \"adjoint.md\"\n            \"boundary_functional.md\"\n            \"misc.md\"\n          ]\n  Depth = 1"
},

{
    "location": "solver/euler/types.html#",
    "page": "Datatypes",
    "title": "Datatypes",
    "category": "page",
    "text": ""
},

{
    "location": "solver/euler/types.html#EulerEquationMod.EulerData_",
    "page": "Datatypes",
    "title": "EulerEquationMod.EulerData_",
    "category": "type",
    "text": "EulerEquationMod.EulerData_\n\nThis type is an implementation of the abstract EulerData.  It is   parameterized by the residual datatype Tres and the mesh datatype Tmsh   because it stores some arrays of those types.  Tres is the \'maximum\' type of   Tsol and Tmsh, where Tsol is the type of the solution variables.   It is also paramterized by var_type, which should be a symbol describing   the set of variables stored in eqn.q.  Currently supported values are   :conservative and :entropy, which indicate the conservative variables and   the entropy variables described in:\n\n\'A New Finite Element Formulation for   Computational Fluid Dynamics: Part I\' by Hughes et al.`\n\nNote: this constructor does not fully populate all fields.  The           [init])@ref) function must be called to finish initialization.\n\nStatic Parameters:\n\nTsol : datatype of variables solution variables, ie. the       q vector and array\nTres : datatype of residual. ie. eltype(res_vec)\nTdim : dimensionality of equation, integer, (2 or 3, currently only 2 is       supported).\nTmsh : datatype of mesh related quantities\nvar_type : symbol describing variables used in weak form, (:conservative           or :entropy)\n\nFields\n\nThis type has many fields, not all of them are documented here.  A few   of the most important ones are:\n\ncomm: MPI communicator\ncommsize: size of MPI communicator\nmyrank: MPI rank of this process\n\nWhen computing the jacobian explicitly (options key calc_jac_explicit),   Tsol and Tres are typically Float64, however node-level operations    sometime use complex numbers or dual numbers.  Also, some operations on   entropy variables require doing parts of the computation with   conservative variables.  To support these use-cases, the fields\n\nparams: ParamType object with Tsol, Tres, Tmsh, and var_type matching the equation object\nparams_conservative: ParamType object with Tsol, Tres, and Tmsh matching the EulerData_ object, but var_type = :conservative\nparams_entropy: similar to param_conservative, but var_type = :entropy\nparams_complex: ParamType object with Tmsh and var_type matching the EulerData_ object, but Tsol = Tres = Complex128 \n\nexist.\n\n\n\n"
},

{
    "location": "solver/euler/types.html#EulerEquationMod.AbstractEulerData",
    "page": "Datatypes",
    "title": "EulerEquationMod.AbstractEulerData",
    "category": "type",
    "text": "EulerEquationMod.AbstractEulerData{Tsol, Tres}\n\nThis abstract type should be the supertype of all solution data objects   that are related to the Euler equations.\n\nIt should be used for specifying the type of a function argument only when   the function does no operations on the solution data object itself, it just   passes it onto other functions that do the work (thus AbstractEulerData   should be used for only the highest level functions).\n\nAnother way of saying say it is that this type should only be used when   the function only needs to ensure that it is solving the Euler equations,   but does not care even a little bit about how.\n\n\n\n"
},

{
    "location": "solver/euler/types.html#EulerEquationMod.EulerData",
    "page": "Datatypes",
    "title": "EulerEquationMod.EulerData",
    "category": "type",
    "text": "EulerEquationMod.EulerData\n\nThis type, although abstract, is the type functions should use for their   input arguments if they do any operations on the solution data object.   It stores all data used in evaluting the Euler Equations.\n\nIt is paramaterized on the types Tsol, the type of the   conservative variables q, and Tdim, the dimension of the equation\n\nIt should have the following fields:\n\n* res_type : datatype of residual (depreciated)\n* q  : 3D array holding conservative variables\n* q_vec  : vector to assemble q into\n* aux_vars : 3D array holding auxiliary variables\n* flux_parametric : 4D array [ndof per node, nnodes per element, nelements, Tdim]\n         holding the Euler flux in the xi and eta directions\n* res  : 3D array holding residual\n* res_vec   : vector form of res\n* edgestab_alpha : paramater used for edge stabilization, 4d array\n* bndryflux : 3D array holding boundary flux data\n* stabscale : 2D array holding edge stabilization scale factor\n* M : vector holding the mass matrix\n* Minv :  vector holding inverse mass matrix\n# Minv3D :  3D array holding inverse mass matrix for application to res (not res_vec)\n\n\n\n"
},

{
    "location": "solver/euler/types.html#EulerEquationMod.ParamType",
    "page": "Datatypes",
    "title": "EulerEquationMod.ParamType",
    "category": "type",
    "text": "EulerEquationMod.ParamType\n\nThis type holds the values of any constants or parameters needed during the   computation.  These parameters can be specified in the opts dictionary or   have default values set here.  If there is no reasonable default, values   are initialized to -1\n\nThere are also a bunch of arrays that are used as temporaries by low   level functions (to avoid having to allocate arrays themselves, which is   a performance trap).  In general, this Type is used as a container to pass   around values.\n\ngamma and R are the independent themodynamic variables\n\nWhether this type should be immutable or not is an open question\n\nThis type is paramaterized on the dimension of the equation for purposes   of multiple dispatch\n\nStatic Parameters:\n\nTdim : dimensionality of the equation, integer, (used for dispatch)\nvar_type : type of variables used used in the weak form, symbol, (used for         dispatch), currently supported values: :conservative, :entropy\nTsol : datatype of solution variables q\nTres : datatype of residual\nTmsh : datatype of mesh related quantities (mapping jacobian etc.)\n\nFields (with default values):\n\ncv  : specific heat constant\nR : specific gas constant (J/(Kg*K))\ngamma : ratio of specific heats\ngamma_1 : gamma - 1\n\nFields (without default values):\n\nMa  : free stream Mach number\nRe  : free stream Reynolds number\naoa : angle of attack (radians)\n\n\n\n"
},

{
    "location": "solver/euler/types.html#EulerEquationMod.ParamType2",
    "page": "Datatypes",
    "title": "EulerEquationMod.ParamType2",
    "category": "type",
    "text": "Useful alias for 2D ParamType\n\n\n\n"
},

{
    "location": "solver/euler/types.html#EulerEquationMod.ParamType3",
    "page": "Datatypes",
    "title": "EulerEquationMod.ParamType3",
    "category": "type",
    "text": "Useful alias for 3D ParamType\n\n\n\n"
},

{
    "location": "solver/euler/types.html#EulerEquationMod.calcElemSurfaceArea-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tdim where Tres where Tsol where Tmsh",
    "page": "Datatypes",
    "title": "EulerEquationMod.calcElemSurfaceArea",
    "category": "method",
    "text": "EulerEquationMod.calcElemFurfaceArea\n\nThis function calculates the wet area of each element. A weight of 2 is given to faces with Dirichlet boundary conditions. Arguments: mesh: AbstractMesh sbp: SBP operator eqn: an implementation of EulerData. Does not have to be fully initialized.\n\n\n\n"
},

{
    "location": "solver/euler/types.html#EulerEquationMod.calcTraceInverseInequalityConst-Union{Tuple{SummationByParts.AbstractSBP{Tsbp},SummationByParts.AbstractFace{Tsbp}}, Tuple{Tsbp}} where Tsbp",
    "page": "Datatypes",
    "title": "EulerEquationMod.calcTraceInverseInequalityConst",
    "category": "method",
    "text": "Compute the constant coefficent in inverse trace ineqality, i.e., the largest eigenvalue of  B^{1/2} R H^{-1} R^{T} B^{1/2}\n\nInput:   sbp Output:   cont_tii\n\n\n\n"
},

{
    "location": "solver/euler/types.html#EulerEquationMod.cleanup-Tuple{ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,EulerEquationMod.EulerData,Any}",
    "page": "Datatypes",
    "title": "EulerEquationMod.cleanup",
    "category": "method",
    "text": "This function performs all cleanup activities before the run_physics()   function returns.  The mesh, sbp, eqn, opts are returned by run_physics()   so there is not much cleanup that needs to be done, mostly closing files.\n\nInputs/Outputs:\n\nmesh: an AbstractMesh object\nsbp: an SBP operator\neqn: the EulerData object\nopts: the options dictionary\n\n\n\n"
},

{
    "location": "solver/euler/types.html#EulerEquationMod.getTypeParameters-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tres where Tsol where Tmsh",
    "page": "Datatypes",
    "title": "EulerEquationMod.getTypeParameters",
    "category": "method",
    "text": "EulerEquationMod.getTypeParameters\n\nGets the type parameters for mesh and equation objects.\n\nInput\n\nmesh : Object of abstract meshing type.\neqn  : Euler Equation object.\n\nOutput\n\nTmsh : Type parameter of the mesh.\nTsol : Type parameter of the solution array.\nTres : Type parameter of the residual array.\n\n\n\n"
},

{
    "location": "solver/euler/types.html#EulerEquationMod.openLoggingFiles-Tuple{Any,Any}",
    "page": "Datatypes",
    "title": "EulerEquationMod.openLoggingFiles",
    "category": "method",
    "text": "This function opens all used for logging data.  In particular, every data   file that has data appended to it in majorIterationCallback should be   opened here.  Most files are of type BufferedIO, so they must be flushed   periodically.\n\nThis function requires each output to have two keys: \"write_outname\"   and \"write_outname_fname\", where the first has a boolean value that   controls whether or not to write the output, and the second is the   file name (including extension) to write.\n\nThis function contains a list of all possible log files.  Every new   log file must be added to the list\n\nInputs:\n\nmesh: an AbstractMesh (needed for MPI Communicator)\nopts: options dictionary\n\nOutputs:\n\nfile_dict: dictionary mapping names of files to the file object             ie. opts[\"write_entropy_fname\"] => f\n\nExceptions: this function will throw an exception if any two file names               are the same\n\nImplementation notes:     When restarting, all files must be appended to.  Currently, files     are appended to in all cases.\n\n\n\n"
},

{
    "location": "solver/euler/types.html#ODLCommonTools.getAllTypeParams-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},EulerEquationMod.EulerData_{Tsol,Tres,Tdim,Tmsh,var_type},Any}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}, Tuple{var_type}} where var_type where Tdim where Tres where Tsol where Tmsh",
    "page": "Datatypes",
    "title": "ODLCommonTools.getAllTypeParams",
    "category": "method",
    "text": "EulerEquationMod.getAllTypeParameters\n\nGets the type parameters for mesh and equation objects.\n\nInput\n\nmesh : Object of abstract meshing type.\neqn  : Euler Equation object.\nopts : Options dictionary\n\nOutput\n\ntuple : Tuple of type parameters. Ordering is same as that of the concrete eqn object within this physics module.\n\n\n\n"
},

{
    "location": "solver/euler/types.html#EulerEquationMod.PhysicsName",
    "page": "Datatypes",
    "title": "EulerEquationMod.PhysicsName",
    "category": "constant",
    "text": "This physics is named Euler\n\n\n\n"
},

{
    "location": "solver/euler/types.html#EulerEquationMod.FaceElementIntegralType",
    "page": "Datatypes",
    "title": "EulerEquationMod.FaceElementIntegralType",
    "category": "type",
    "text": "Functor type for faceElementIntegrals.  These integrals operate on a face,   but require data from the entirety of the elements that make up the   face, rather than data interpolated to the face\n\n\n\n"
},

{
    "location": "solver/euler/types.html#Euler-Types-1",
    "page": "Datatypes",
    "title": "Euler Types",
    "category": "section",
    "text": "  CurrentModule = EulerEquationModThis page provides documentations for DataTypes and and simple functions that are defined in the Euler module  Modules = [EulerEquationMod]\n  Pages = [\"euler/types.jl\", \"euler/EulerEquationMod.jl\"]"
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
    "text": "  CurrentModule = EulerEquationModTODO: finish this page (after euler_funcs.jl has been broken up"
},

{
    "location": "solver/euler/volume.html#EulerEquationMod.evalVolumeIntegrals",
    "page": "Volume Integrals",
    "title": "EulerEquationMod.evalVolumeIntegrals",
    "category": "function",
    "text": "EulerEquationMod.evalVolumeIntegrals\n\nThis function evaluates the volume integrals of the Euler equations by   calling the appropriate SBP function on the Euler flux stored in   eqn.flux_parametric.   This function knows the dimension of the equation and does the right   thing for all dimensions.  eqn.res is updated with the result\n\nThis is a mid level function.\n\n\n\n"
},

{
    "location": "solver/euler/volume.html#Entry-Point-1",
    "page": "Volume Integrals",
    "title": "Entry Point",
    "category": "section",
    "text": "evalVolumeIntegrals"
},

{
    "location": "solver/euler/volume_diff.html#",
    "page": "Volume Integrals Jacobian",
    "title": "Volume Integrals Jacobian",
    "category": "page",
    "text": ""
},

{
    "location": "solver/euler/volume_diff.html#Volume-Terms-Differentiated-1",
    "page": "Volume Integrals Jacobian",
    "title": "Volume Terms Differentiated",
    "category": "section",
    "text": "  CurrentModule = EulerEquationModThis page documents the differentiated version of functions used to evaluate the volume integrals"
},

{
    "location": "solver/euler/volume_diff.html#EulerEquationMod.evalVolumeIntegrals_diff",
    "page": "Volume Integrals Jacobian",
    "title": "EulerEquationMod.evalVolumeIntegrals_diff",
    "category": "function",
    "text": "Differentiated version of evalVolumeIntegrals.  Throws an error   for unsupported schemes.\n\nInputs\n\nmesh\nsbp\neqn\nopts\nassembler\n\n\n\n"
},

{
    "location": "solver/euler/volume_diff.html#Entry-Point-1",
    "page": "Volume Integrals Jacobian",
    "title": "Entry Point",
    "category": "section",
    "text": "evalVolumeIntegrals_diff"
},

{
    "location": "solver/euler/volume_diff.html#EulerEquationMod.calcEulerFlux_diff-Union{Tuple{EulerEquationMod.ParamType{2,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1},AbstractArray{Tres,1},AbstractArray{Tmsh,N} where N,AbstractArray{Tsol,2}}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tres where Tsol where Tmsh",
    "page": "Volume Integrals Jacobian",
    "title": "EulerEquationMod.calcEulerFlux_diff",
    "category": "method",
    "text": "Computes the jacobian of calcEulerFlux with respect to q.   Methods are available for 2D and 3D\n\nInputs\n\nparams: ParamType, conservative variables only\nq: vector of conservative variables at node\naux_vars: auxiliary variables at the node\ndir: direction vector (possibly scaled) to compute the flux jacobian in\n\nInputs/Outputs\n\nFjac: flux jacobian, numDofPerNode x numDofPerNode\n\nAliasing restrictions: params.p_dot is overwritten\n\n\n\n"
},

{
    "location": "solver/euler/volume_diff.html#EulerEquationMod.calcPressure_diff-Union{Tuple{EulerEquationMod.ParamType{2,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1},AbstractArray{Tsol,1}}, Tuple{Tsol}} where Tsol",
    "page": "Volume Integrals Jacobian",
    "title": "EulerEquationMod.calcPressure_diff",
    "category": "method",
    "text": "Computes the gradient of pressure with respect to q at a node.   Methods are available in 2D and 3D\n\nInputs\n\nparams: ParamType, conservative variables only\nq: vector of conservative variables at the node\n\nInputs/Outputs\n\npdot: vector of length numDofPerNode, overwritten with derivative of p        wrt q\n\n\n\n"
},

{
    "location": "solver/euler/volume_diff.html#EulerEquationMod.calcVolumeIntegralsStrong_nopre_diff-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type,Any,PDESolver.AssembleElementData}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tdim where Tres where Tsol where Tmsh",
    "page": "Volume Integrals Jacobian",
    "title": "EulerEquationMod.calcVolumeIntegralsStrong_nopre_diff",
    "category": "method",
    "text": "Computes the derivative of the strong form volume terms with   respect to q, ie.\n\nd/dq (-Q * f)\n\nbut only the mesh.numDofPerNode x mesh.numDofPerNode diagonal block\n\nInputs\n\nmesh\nsbp\neqn\nopts\nassembler: used to assemble the contribution of each element into the            Jacobian\n\n\n\n"
},

{
    "location": "solver/euler/volume_diff.html#EulerEquationMod.calcVolumeIntegrals_nopre_diff-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type,Any,PDESolver.AssembleElementData}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tdim where Tres where Tsol where Tmsh",
    "page": "Volume Integrals Jacobian",
    "title": "EulerEquationMod.calcVolumeIntegrals_nopre_diff",
    "category": "method",
    "text": "Computes the derivataive of the volume term of the Roe scheme with respect   to q, ie\n\nd/dq (Q^T * f)\n\nInputs\n\nmesh\nsbp\neqn\nopts\nassembler: used to assemble the contribution of each element into the            Jacobian\n\n\n\n"
},

{
    "location": "solver/euler/volume_diff.html#Functions-1",
    "page": "Volume Integrals Jacobian",
    "title": "Functions",
    "category": "section",
    "text": "  Modules = [EulerEquationMod]\n  Order = [:function]\n  Pages = [\"euler/euler_funcs_diff.jl\"]"
},

{
    "location": "solver/euler/flux.html#",
    "page": "Face Integrals",
    "title": "Face Integrals",
    "category": "page",
    "text": ""
},

{
    "location": "solver/euler/flux.html#sec:euler_face_integrals-1",
    "page": "Face Integrals",
    "title": "Face Integrals",
    "category": "section",
    "text": "  CurrentModule = EulerEquationModThis page describes the functions that evaluate the face and shared face integrals."
},

{
    "location": "solver/euler/flux.html#EulerEquationMod.evalFaceIntegrals",
    "page": "Face Integrals",
    "title": "EulerEquationMod.evalFaceIntegrals",
    "category": "function",
    "text": "EulerEquationMod.evalFaceIntegrals\n\nThis function evaluates the face integrals in a DG formulation and   updates the residual.  The array eqn.flux_face must already be populated   with the face flux.\n\nInputs:\n\nmesh: an AbstractDGMesh\nsbp: an SBP operator\neqn: an EulerData object\nopts: the options dictonary\n\nOutputs:\n\nnone\n\n\n\n"
},

{
    "location": "solver/euler/flux.html#EulerEquationMod.evalSharedFaceIntegrals",
    "page": "Face Integrals",
    "title": "EulerEquationMod.evalSharedFaceIntegrals",
    "category": "function",
    "text": "EulerEquationMod.evalSharedFaceIntegrals\n\nThis function does the computation that needs the parallel   communication to have finished already, namely the face integrals   for the shared faces\n\nInputs:\n\nmesh\nsbp\neqn\nopts\n\n\n\n"
},

{
    "location": "solver/euler/flux.html#Entry-Point-1",
    "page": "Face Integrals",
    "title": "Entry Point",
    "category": "section",
    "text": "evalFaceIntegrals\nevalSharedFaceIntegrals"
},

{
    "location": "solver/euler/flux.html#EulerEquationMod.calcFaceFlux-Union{Tuple{ODLCommonTools.AbstractDGMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,:conservative},ODLCommonTools.FluxType,AbstractArray{ODLCommonTools.Interface,1},AbstractArray{Tres,3}}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tdim where Tres where Tsol where Tmsh",
    "page": "Face Integrals",
    "title": "EulerEquationMod.calcFaceFlux",
    "category": "method",
    "text": "EulerEquationMod.calcFaceFlux\n\nThis function calculates the DG flux between a specified set of faces,   using the solution data at the faces stored in eqn.q_face.   Note that the flux is negated because the face integrals have a   negative sign in the weak form.\n\nConservative variables only!\n\nInputs:\n\nmesh\nsbp\neqn\nfunctor: the functor that calculates the flux at a node\ninterfaces: an array of type Interface that specifies which interfaces            to calculate the flux for\n\nInputs/Outputs:\n\nface_flux: array to store the flux in, numDofPerNode x nnodesPerFace           x length(interfaces)\n\nThe functor must have the signature:\n\nfunc( uL, qR, aux_vars, dxidx, nrm, flux_j, eqn.params)\n\nwhere uL and uR are the solution values for a node on the left and right   elements, aux_vars are the auxiliary variables for the node,   dxidx is the scaled mapping jacobian for elementL, and nrm is the face   normal in reference space. flux_j is the array of length numDofPerNode to be   populated with the flux. params is eqn.params.\n\n\n\n"
},

{
    "location": "solver/euler/flux.html#EulerEquationMod.calcFaceIntegral_nopre-Union{Tuple{ODLCommonTools.AbstractDGMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,:conservative},Any,ODLCommonTools.FluxType,AbstractArray{ODLCommonTools.Interface,1}}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tdim where Tres where Tsol where Tmsh",
    "page": "Face Integrals",
    "title": "EulerEquationMod.calcFaceIntegral_nopre",
    "category": "method",
    "text": "Like calcFaceFlux, but computes the flux for a single element and   then integrates it immediately, updating eqn.res\n\nInputs:\n\nmesh\nsbp\neqn\nopts\nfunctor: a FluxType that evalutes the flux\ninterfaces: the vector of Interfaces to compute the integrals for\n\n\n\n"
},

{
    "location": "solver/euler/flux.html#EulerEquationMod.calcSharedFaceElementIntegralsStaggered_element_inner-Union{Tuple{ODLCommonTools.AbstractDGMesh{Tmsh},ODLCommonTools.AbstractDGMesh{Tmsh},SummationByParts.AbstractSBP,SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim,Any,Utils.SharedFaceData,EulerEquationMod.FaceElementIntegralType,ODLCommonTools.FluxType}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tres where Tsol where Tmsh",
    "page": "Face Integrals",
    "title": "EulerEquationMod.calcSharedFaceElementIntegralsStaggered_element_inner",
    "category": "method",
    "text": "Like calcSharedFaceElementIntegrals_element_inner, but for staggered grid.\n\ndata.q_recv is the solution grid data.  This function interpolates it to the   flux grid on the fly.\n\nInputs:\n\nmesh_s: solution grid mesh\nmesh_f: flux grid mesh\nsbp_s: SBP operator for solution grid\nsbp_f: SBP operator for the flux grid\nopts: options dictionary\ndata: SharedFaceData\nface_integral_functor: FaceElementIntegralType\nflux_functor: FluxType passed to face_integral functor\n\nInputs/Outputs:\n\n* eqn: equation object (lives on solution grid)\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/flux.html#EulerEquationMod.calcSharedFaceElementIntegrals_element-Union{Tuple{ODLCommonTools.AbstractDGMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim,Any,Utils.SharedFaceData}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tres where Tsol where Tmsh",
    "page": "Face Integrals",
    "title": "EulerEquationMod.calcSharedFaceElementIntegrals_element",
    "category": "method",
    "text": "This function is a thin wrapper around   calcSharedFaceElementIntegrals_element_inner,   presenting the interface needed by finishExchangeData.   See that function for the interface details.\n\n\n\n"
},

{
    "location": "solver/euler/flux.html#EulerEquationMod.calcSharedFaceElementIntegrals_element_inner-Union{Tuple{ODLCommonTools.AbstractDGMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim,Any,Utils.SharedFaceData,EulerEquationMod.FaceElementIntegralType,ODLCommonTools.FluxType}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tres where Tsol where Tmsh",
    "page": "Face Integrals",
    "title": "EulerEquationMod.calcSharedFaceElementIntegrals_element_inner",
    "category": "method",
    "text": "This function loops over given set of shared faces and computes a face   integral that   uses data from all volume nodes.  See FaceElementIntegralType   for details on the integral performed.\n\nInputs:\n\nmesh\nsbp\neqn\nopts\ndata: a SharedFaceData specifying which shared faces to compute\nface_integral_functor\nflux_functor\n\n\n\n"
},

{
    "location": "solver/euler/flux.html#EulerEquationMod.calcSharedFaceIntegrals-Union{Tuple{ODLCommonTools.AbstractDGMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim where Tres,Any,Utils.SharedFaceData}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Face Integrals",
    "title": "EulerEquationMod.calcSharedFaceIntegrals",
    "category": "method",
    "text": "This function is a thin wrapper around calcSharedFaceIntegrals_inner.   It present the interface needed by finishExchangeData.\n\n\n\n"
},

{
    "location": "solver/euler/flux.html#EulerEquationMod.calcSharedFaceIntegrals_element-Union{Tuple{ODLCommonTools.AbstractDGMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim where Tres,Any,Utils.SharedFaceData}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Face Integrals",
    "title": "EulerEquationMod.calcSharedFaceIntegrals_element",
    "category": "method",
    "text": "This function is a thin wrapper around   calcSharedFaceIntegrals_element_inner.   It presents the interface required by finishExchangeData\n\n\n\n"
},

{
    "location": "solver/euler/flux.html#EulerEquationMod.calcSharedFaceIntegrals_element_inner-Union{Tuple{ODLCommonTools.AbstractDGMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim where Tres,Any,Utils.SharedFaceData,ODLCommonTools.FluxType}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Face Integrals",
    "title": "EulerEquationMod.calcSharedFaceIntegrals_element_inner",
    "category": "method",
    "text": "This function calculates the shared face integrals for a given set of faces.   It uses the MPI send and receive buffers that contain the solution for the   elements on the boundary (rather than the data on the faces).  This   enables calculating a sparse jacobian with minimal parallel communication.\n\nInputs:\n\nmesh\nsbp\neqn\nopts\ndata: a SharedFaceData specifying the faces to calculate\nfunctor: the flux functor\n\n\n\n"
},

{
    "location": "solver/euler/flux.html#EulerEquationMod.calcSharedFaceIntegrals_inner-Union{Tuple{ODLCommonTools.AbstractDGMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim where Tres,Any,Utils.SharedFaceData,ODLCommonTools.FluxType}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Face Integrals",
    "title": "EulerEquationMod.calcSharedFaceIntegrals_inner",
    "category": "method",
    "text": "EulerEquationMod.calcSharedFaceIntegrals\n\nThis function calculates the shared face integrals over a given set of   faces.\n\nInputs:\n\nmesh\nsbp\neqn\nopts\ndata: the SharedFaceData specifying which faces to compute\nfunctor: the FluxType to use for the face flux\n\n\n\n"
},

{
    "location": "solver/euler/flux.html#EulerEquationMod.calcSharedFaceIntegrals_nopre_element_inner-Union{Tuple{ODLCommonTools.AbstractDGMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim,Any,Utils.SharedFaceData,ODLCommonTools.FluxType}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tres where Tsol where Tmsh",
    "page": "Face Integrals",
    "title": "EulerEquationMod.calcSharedFaceIntegrals_nopre_element_inner",
    "category": "method",
    "text": "Like calcSharedFaceIntegrals_element_inner, but performs the integration and   updates eqn.res rather than computing the flux only and storing it in   eqn.flux_sharedface\n\n\n\n"
},

{
    "location": "solver/euler/flux.html#EulerEquationMod.calcSharedFaceIntegrals_nopre_inner-Union{Tuple{ODLCommonTools.AbstractDGMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim,Any,Utils.SharedFaceData,ODLCommonTools.FluxType}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tres where Tsol where Tmsh",
    "page": "Face Integrals",
    "title": "EulerEquationMod.calcSharedFaceIntegrals_nopre_inner",
    "category": "method",
    "text": "Like calcSharedFaceIntegrals_inner, but performs the integration and   updates eqn.res rather than computing the flux and storing it in   eqn.flux_sharedface\n\n\n\n"
},

{
    "location": "solver/euler/flux.html#EulerEquationMod.getFaceElementIntegral-Union{Tuple{ODLCommonTools.AbstractDGMesh{Tmsh},ODLCommonTools.AbstractDGMesh{Tmsh},SummationByParts.AbstractSBP,SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type,EulerEquationMod.FaceElementIntegralType,ODLCommonTools.FluxType,SummationByParts.AbstractFace,AbstractArray{ODLCommonTools.Interface,1}}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tdim where Tres where Tsol where Tmsh",
    "page": "Face Integrals",
    "title": "EulerEquationMod.getFaceElementIntegral",
    "category": "method",
    "text": "This function loops over interfaces and computes a face integral that   uses data from all volume nodes. See FaceElementIntegralType for details on   the integral performed. This function works on the staggered grid.\n\nInputs:     mesh_s: mesh object on solution mesh     mesh_f: mesh object on flux mesh     sbp_s: SBP operator on solution mesh     sbp_f: SBP operator on flux mesh     face_integral_functor: a FaceElementIntegralType functor     flux_functor: a [FluxType] functor that is passed to face_integral_functor                   to use as the numerical flux function     sbpface: an AbstractFace for the flux grid     interfaces: vector of Interfaces that the integral will be                 computed for.\n\nInputs/Outputs:     eqn: equation object (implicitly lives on solution grid).  eqn.res is           updated with the results.\n\n\n\n"
},

{
    "location": "solver/euler/flux.html#EulerEquationMod.getFaceElementIntegral-Union{Tuple{ODLCommonTools.AbstractDGMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type,EulerEquationMod.FaceElementIntegralType,ODLCommonTools.FluxType,SummationByParts.AbstractFace,AbstractArray{ODLCommonTools.Interface,1}}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tdim where Tres where Tsol where Tmsh",
    "page": "Face Integrals",
    "title": "EulerEquationMod.getFaceElementIntegral",
    "category": "method",
    "text": "This function loops over interfaces and computes a face integral that   uses data from all volume nodes. See FaceElementIntegralType   for details on the integral performed.\n\n\n\n"
},

{
    "location": "solver/euler/flux.html#EulerEquationMod.getFluxFunctors-Tuple{ODLCommonTools.AbstractDGMesh,Any,Any,Any}",
    "page": "Face Integrals",
    "title": "EulerEquationMod.getFluxFunctors",
    "category": "method",
    "text": "EulerEquationMod.getFluxFunctors\n\nThis function retrieves the flux functors from the dictonary and   stores them to eqn.flux_func.\n\nInputs:     mesh: an AbstractDGMesh     sbp     eqn     opts\n\n\n\n"
},

{
    "location": "solver/euler/flux.html#EulerEquationMod.interpolateFace-Union{Tuple{ODLCommonTools.AbstractDGMesh,Any,Any,Any,AbstractArray{T,3} where T,AbstractArray{Tsol,4}}, Tuple{Tsol}} where Tsol",
    "page": "Face Integrals",
    "title": "EulerEquationMod.interpolateFace",
    "category": "method",
    "text": "EulerEquationMod.interpolateFace\n\nThis function interpolates the solution values from the internal nodes   to the face flux points of the elements\n\nInputs:\n\nmesh: an AbstractDGMesh\nsbp\neqn\nopts\nq: a 3D array of solution values at the nodes, numDofPerNode x   numNodesPerElement x numEl\n\nInputs/Outputs:\n\nq_face: a 4D array of solution values at each interface,        numDofPerNode x 2 x numfacenodes x numInterface        q_face[:, 1, j, i] stores the q values for elementL of interface        i node j and q_face[:, 2, j, i] stores the values for elementR\n\neqn.aux_vars_face is also populated\n\n\n\n"
},

{
    "location": "solver/euler/flux.html#EulerEquationMod.FluxDict",
    "page": "Face Integrals",
    "title": "EulerEquationMod.FluxDict",
    "category": "constant",
    "text": "EulerEquationMod.FluxDict\n\nThis dictonary maps the names of the fluxes (Strings) to the   functor object itself.  All flux functors should be added to the dictionary.\n\nAll fluxes have one method that calculates the flux in a particular direction   at a node.  Some fluxes have an additional method that computes the flux   in several directions at a node in a single function call, which can be   more efficient.  See calcEulerFlux_standard for an example.\n\nIn general, these functors call similarly-named function in bc_solvers.jl.   It is recommened to look at the documentation for those functions.\n\nTODO: document signature of the functors here\n\n\n\n"
},

{
    "location": "solver/euler/flux.html#Functions-1",
    "page": "Face Integrals",
    "title": "Functions",
    "category": "section",
    "text": "  Modules = [EulerEquationMod]\n  Order = [:function, :constant, :macro]\n  Pages = [\"euler/flux.jl\"]"
},

{
    "location": "solver/euler/flux.html#EulerEquationMod.DucrosFlux",
    "page": "Face Integrals",
    "title": "EulerEquationMod.DucrosFlux",
    "category": "type",
    "text": "Calls calcEulerFlux_Ducros\n\n\n\n"
},

{
    "location": "solver/euler/flux.html#EulerEquationMod.ErrorFlux",
    "page": "Face Integrals",
    "title": "EulerEquationMod.ErrorFlux",
    "category": "type",
    "text": "This flux function throws an error. Useful for defaults.\n\n\n\n"
},

{
    "location": "solver/euler/flux.html#EulerEquationMod.IRFlux",
    "page": "Face Integrals",
    "title": "EulerEquationMod.IRFlux",
    "category": "type",
    "text": "Calls calcEulerFlux_IR\n\n\n\n"
},

{
    "location": "solver/euler/flux.html#EulerEquationMod.IRSLFFlux",
    "page": "Face Integrals",
    "title": "EulerEquationMod.IRSLFFlux",
    "category": "type",
    "text": "Calls calcEulerFlux_IRSLF\n\n\n\n"
},

{
    "location": "solver/euler/flux.html#EulerEquationMod.IdentityFlux",
    "page": "Face Integrals",
    "title": "EulerEquationMod.IdentityFlux",
    "category": "type",
    "text": "This flux function sets F = q.  Useful for testing\n\n\n\n"
},

{
    "location": "solver/euler/flux.html#EulerEquationMod.LFFlux",
    "page": "Face Integrals",
    "title": "EulerEquationMod.LFFlux",
    "category": "type",
    "text": "Calls calcEulerFlux_standard\n\n\n\n"
},

{
    "location": "solver/euler/flux.html#EulerEquationMod.RoeFlux",
    "page": "Face Integrals",
    "title": "EulerEquationMod.RoeFlux",
    "category": "type",
    "text": "Calls the RoeSolver\n\n\n\n"
},

{
    "location": "solver/euler/flux.html#EulerEquationMod.SIPGViscousFlux",
    "page": "Face Integrals",
    "title": "EulerEquationMod.SIPGViscousFlux",
    "category": "type",
    "text": "Calls the SIPG (viscous) flux\n\n\n\n"
},

{
    "location": "solver/euler/flux.html#sec:euler_flux_functors-1",
    "page": "Face Integrals",
    "title": "Flux Functors",
    "category": "section",
    "text": "TODO: move this to another file (both code and docs)  Modules = [EulerEquationMod]\n  Order = [:type]\n  Pages = [\"euler/flux.jl\"]"
},

{
    "location": "solver/euler/flux_diff.html#",
    "page": "Face Integrals Jacobian",
    "title": "Face Integrals Jacobian",
    "category": "page",
    "text": ""
},

{
    "location": "solver/euler/flux_diff.html#sec:euler_face_integrals_diff-1",
    "page": "Face Integrals Jacobian",
    "title": "Face Integrals Differentiated",
    "category": "section",
    "text": "  CurrentModule = EulerEquationModThis page describes the functions that evaluate the face and shared face integrals."
},

{
    "location": "solver/euler/flux_diff.html#EulerEquationMod.evalFaceIntegrals_diff",
    "page": "Face Integrals Jacobian",
    "title": "EulerEquationMod.evalFaceIntegrals_diff",
    "category": "function",
    "text": "Differentiated version of evalFaceIntegrals_diff.  Throws an   error for unsupported schemes.\n\nInputs\n\nmesh\nsbp\neqn\nopts\nassembler\n\n\n\n"
},

{
    "location": "solver/euler/flux_diff.html#EulerEquationMod.evalSharedFaceIntegrals_diff",
    "page": "Face Integrals Jacobian",
    "title": "EulerEquationMod.evalSharedFaceIntegrals_diff",
    "category": "function",
    "text": "Differentiated version of evalSharedFaceIntegrals.   Currently only supports opts[parallel_data] == element because it is   used for Jacobian calculation.\n\nInputs\n\nmesh\nsbp\neqn\nopts\nassembler\n\n\n\n"
},

{
    "location": "solver/euler/flux_diff.html#Entry-Point-1",
    "page": "Face Integrals Jacobian",
    "title": "Entry Point",
    "category": "section",
    "text": "evalFaceIntegrals_diff\nevalSharedFaceIntegrals_diff"
},

{
    "location": "solver/euler/flux_diff.html#EulerEquationMod.calcSharedFaceIntegrals_element_diff-Union{Tuple{ODLCommonTools.AbstractDGMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim where Tres,Any,Utils.SharedFaceData}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Face Integrals Jacobian",
    "title": "EulerEquationMod.calcSharedFaceIntegrals_element_diff",
    "category": "method",
    "text": "Differentiated version of   calcSharedFaceIntegrals_element_inner.   It presents the interface required by finishExchangeData\n\n\n\n"
},

{
    "location": "solver/euler/flux_diff.html#EulerEquationMod.calcSharedFaceIntegrals_nopre_element_inner_diff-Union{Tuple{ODLCommonTools.AbstractDGMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim,Any,Utils.SharedFaceData,ODLCommonTools.FluxType_diff,PDESolver.AssembleElementData}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tres where Tsol where Tmsh",
    "page": "Face Integrals Jacobian",
    "title": "EulerEquationMod.calcSharedFaceIntegrals_nopre_element_inner_diff",
    "category": "method",
    "text": "Differentiated version of calcSharedFaceIntegrals_nopre_element_inner,\n\n\n\n"
},

{
    "location": "solver/euler/flux_diff.html#EulerEquationMod.getFluxFunctors_diff-Tuple{ODLCommonTools.AbstractDGMesh,Any,Any,Any}",
    "page": "Face Integrals Jacobian",
    "title": "EulerEquationMod.getFluxFunctors_diff",
    "category": "method",
    "text": "Gets the functor for the differentiated version of the flux function and   stores is to eqn.flux_func_diff\n\nInputs\n\nmesh\nsbp\neqn\nopts\n\nOptions Keys\n\nUses the Flux_name key to get the functor (same as for the non-differentiated   version). \n\nThis function only gets the functor if the jacobian will be computed explicitly   as specified by calc_jac_explicit, otherwise it uses a functor that   will throw an error if used.\n\n\n\n"
},

{
    "location": "solver/euler/flux_diff.html#EulerEquationMod.FluxDict_diff",
    "page": "Face Integrals Jacobian",
    "title": "EulerEquationMod.FluxDict_diff",
    "category": "constant",
    "text": "Container for all differentiated flux functors.  Maps name to object.   The names are exactly the same as the non-differentiated functor.\n\n\n\n"
},

{
    "location": "solver/euler/flux_diff.html#Functions-1",
    "page": "Face Integrals Jacobian",
    "title": "Functions",
    "category": "section",
    "text": "  Modules = [EulerEquationMod]\n  Order = [:function, :constant, :macro]\n  Pages = [\"euler/flux_diff.jl\"]"
},

{
    "location": "solver/euler/flux_diff.html#EulerEquationMod.RoeFlux_diff",
    "page": "Face Integrals Jacobian",
    "title": "EulerEquationMod.RoeFlux_diff",
    "category": "type",
    "text": "Calls the RoeSolver_diff\n\n\n\n"
},

{
    "location": "solver/euler/flux_diff.html#Flux-Functors-1",
    "page": "Face Integrals Jacobian",
    "title": "Flux Functors",
    "category": "section",
    "text": "TODO: move this to another file (both code and docs)  Modules = [EulerEquationMod]\n  Order = [:type]\n  Pages = [\"euler/flux_diff.jl\"]"
},

{
    "location": "solver/euler/faceElementIntegrals.html#",
    "page": "Face Element Integrals",
    "title": "Face Element Integrals",
    "category": "page",
    "text": ""
},

{
    "location": "solver/euler/faceElementIntegrals.html#sec:euler_face_element_integrals-1",
    "page": "Face Element Integrals",
    "title": "Face Element Integrals",
    "category": "section",
    "text": "  CurrentModule = EulerEquationModThis page describes the functions that evaluate the face element integrals for a single interface.  The functions that loop over all the interfaces are located on the face integrals page. These integrals require data from all the nodes of the elements rather than the face nodes as with regular face integrals.These integrals are used by the entropy stable scheme, and some of them internally use a numerical flux function.  This flux function must satisfy an entropy property for the resulting scheme to be entropy stable!  The IR flux function is typically used."
},

{
    "location": "solver/euler/faceElementIntegrals.html#EulerEquationMod.calcECFaceIntegral-Union{Tuple{ODLCommonTools.AbstractParamType{Tdim},SummationByParts.DenseFace,ODLCommonTools.Interface,AbstractArray{Tsol,2},AbstractArray{Tsol,2},AbstractArray{Tres,2},AbstractArray{Tmsh,2},ODLCommonTools.FluxType,AbstractArray{Tres,2},AbstractArray{Tres,2}}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tmsh where Tres where Tsol where Tdim",
    "page": "Face Element Integrals",
    "title": "EulerEquationMod.calcECFaceIntegral",
    "category": "method",
    "text": "Calculate the face integrals in an entropy conservative manner for a given   interface.  Unlike standard face integrals, this requires data from   the entirety of both elements, not just data interpolated to the face\n\nresL and resR are updated with the results of the computation for the    left and right elements, respectively.\n\nNote that nrm_xy must contains the normal vector in x-y space at the   face nodes.\n\nThe flux function must be symmetric!\n\nAliasing restrictions: none, although its unclear what the meaning of this                          function would be if resL and resR alias\n\nPerformance note: the version in the tests is the same speed as this one                     for p=1 Omega elements and about 10% faster for                      p=4 elements, but would not be able to take advantage of                      the sparsity of R for SBP Gamma elements\n\n\n\n"
},

{
    "location": "solver/euler/faceElementIntegrals.html#EulerEquationMod.calcECFaceIntegral-Union{Tuple{ODLCommonTools.AbstractParamType{Tdim},SummationByParts.SparseFace,ODLCommonTools.Interface,AbstractArray{Tsol,2},AbstractArray{Tsol,2},AbstractArray{Tres,2},AbstractArray{Tmsh,2},ODLCommonTools.FluxType,AbstractArray{Tres,2},AbstractArray{Tres,2}}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tmsh where Tres where Tsol where Tdim",
    "page": "Face Element Integrals",
    "title": "EulerEquationMod.calcECFaceIntegral",
    "category": "method",
    "text": "Method for sparse faces.  See other method for details\n\nAliasing restrictions: params.flux_vals1 must not be in use\n\n\n\n"
},

{
    "location": "solver/euler/faceElementIntegrals.html#EulerEquationMod.calcESLFFaceIntegral-Union{Tuple{ODLCommonTools.AbstractParamType{Tdim},SummationByParts.AbstractFace,ODLCommonTools.Interface,AbstractArray{Tsol,2},AbstractArray{Tsol,2},AbstractArray{Tres,2},AbstractArray{Tmsh,2},ODLCommonTools.FluxType,AbstractArray{Tres,2},AbstractArray{Tres,2}}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tmsh where Tres where Tsol where Tdim",
    "page": "Face Element Integrals",
    "title": "EulerEquationMod.calcESLFFaceIntegral",
    "category": "method",
    "text": "Calculate the face integral in an entropy stable manner using Lax-Friedrich   type dissipation.     This uses calcECFaceIntegral and calcLFEntropyPenaltyIntegral internally,    see those functions for details.\n\n\n\n"
},

{
    "location": "solver/euler/faceElementIntegrals.html#EulerEquationMod.calcESLW2FaceIntegral-Union{Tuple{ODLCommonTools.AbstractParamType{Tdim},SummationByParts.AbstractFace,ODLCommonTools.Interface,AbstractArray{Tsol,2},AbstractArray{Tsol,2},AbstractArray{Tres,2},AbstractArray{Tmsh,2},ODLCommonTools.FluxType,AbstractArray{Tres,2},AbstractArray{Tres,2}}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tmsh where Tres where Tsol where Tdim",
    "page": "Face Element Integrals",
    "title": "EulerEquationMod.calcESLW2FaceIntegral",
    "category": "method",
    "text": "Calculate the face integral in an entropy stable manner using   Lax-Wendroff type dissipation.     This uses calcECFaceIntegral and calcLW2EntropyPenaltyIntegral internally,    see those functions for details.\n\n\n\n"
},

{
    "location": "solver/euler/faceElementIntegrals.html#EulerEquationMod.calcESLWFaceIntegral-Union{Tuple{ODLCommonTools.AbstractParamType{Tdim},SummationByParts.AbstractFace,ODLCommonTools.Interface,AbstractArray{Tsol,2},AbstractArray{Tsol,2},AbstractArray{Tres,2},AbstractArray{Tmsh,2},ODLCommonTools.FluxType,AbstractArray{Tres,2},AbstractArray{Tres,2}}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tmsh where Tres where Tsol where Tdim",
    "page": "Face Element Integrals",
    "title": "EulerEquationMod.calcESLWFaceIntegral",
    "category": "method",
    "text": "Calculate the face integral in an entropy stable manner using approximate   Lax-Wendroff type dissipation.     This uses calcECFaceIntegral and calcLWEntropyPenaltyIntegral internally,    see those functions for details.\n\n\n\n"
},

{
    "location": "solver/euler/faceElementIntegrals.html#EulerEquationMod.calcEntropyFix-Tuple{EulerEquationMod.ParamType{2,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{T,1} where T}",
    "page": "Face Element Integrals",
    "title": "EulerEquationMod.calcEntropyFix",
    "category": "method",
    "text": "This function modifies the eigenvalues of the euler flux jacobian such   that if any value is zero, a little dissipation is still added.  The   absolute values of the eigenvalues modified eigenvalues are calculated.\n\nMethods are available for 2 and 3 dimensions\n\nThis function depends on the ordering of the eigenvalues produced by   calcEvals.\n\nInputs:     params: ParamType, used to dispatch to 2 or 3D method\n\nInputs/Outputs:     Lambda: vector of eigenvalues to be modified\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/faceElementIntegrals.html#EulerEquationMod.calcLFEntropyPenaltyIntegral-Union{Tuple{EulerEquationMod.ParamType{Tdim,:conservative,Tsol,Tres,Tmsh},SummationByParts.DenseFace,ODLCommonTools.Interface,AbstractArray{Tsol,2},AbstractArray{Tsol,2},AbstractArray{Tres,2},AbstractArray{Tmsh,2},AbstractArray{Tres,2},AbstractArray{Tres,2}}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tmsh where Tres where Tsol where Tdim",
    "page": "Face Element Integrals",
    "title": "EulerEquationMod.calcLFEntropyPenaltyIntegral",
    "category": "method",
    "text": "Calculate a term that provably dissipates (mathematical) entropy using a    Lax-Friedrich type of dissipation.     This   requires data from the left and right element volume nodes, rather than   face nodes for a regular face integral.\n\nNote that nrm_face must contain the scaled face normal vector in x-y space   at the face nodes, and qL, qR, resL, and resR are the arrays for the   entire element, not just the face.\n\nAliasing restrictions: params.nrm2, params.A0, w_vals_stencil, w_vals2_stencil\n\n\n\n"
},

{
    "location": "solver/euler/faceElementIntegrals.html#EulerEquationMod.calcLFEntropyPenaltyIntegral-Union{Tuple{EulerEquationMod.ParamType{Tdim,:conservative,Tsol,Tres,Tmsh},SummationByParts.SparseFace,ODLCommonTools.Interface,AbstractArray{Tsol,2},AbstractArray{Tsol,2},AbstractArray{Tres,2},AbstractArray{Tmsh,2},AbstractArray{Tres,2},AbstractArray{Tres,2}}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tmsh where Tres where Tsol where Tdim",
    "page": "Face Element Integrals",
    "title": "EulerEquationMod.calcLFEntropyPenaltyIntegral",
    "category": "method",
    "text": "Method for sparse faces.  See other method for details\n\nAliasing restrictions: params: v_vals, v_vals2, q_vals, A0, res_vals1, res_vals2\n\n\n\n"
},

{
    "location": "solver/euler/faceElementIntegrals.html#EulerEquationMod.calcLW2EntropyPenaltyIntegral-Union{Tuple{EulerEquationMod.ParamType{Tdim,:conservative,Tsol,Tres,Tmsh},SummationByParts.DenseFace,ODLCommonTools.Interface,AbstractArray{Tsol,2},AbstractArray{Tsol,2},AbstractArray{Tres,2},AbstractArray{Tmsh,2},AbstractArray{Tres,2},AbstractArray{Tres,2}}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tmsh where Tres where Tsol where Tdim",
    "page": "Face Element Integrals",
    "title": "EulerEquationMod.calcLW2EntropyPenaltyIntegral",
    "category": "method",
    "text": "Calculate a term that provably dissipates (mathematical) entropy using a    Lax-Wendroff type of dissipation.     This requires data from the left and right element volume nodes, rather than   face nodes for a regular face integral.\n\nNote nrm_face must contain the scaled normal vector in x-y space   at the face nodes, and qL, qR, resL, and resR are the arrays for the   entire element, not just the face.\n\nImplementation Detail:     Because the scaling does not exist in arbitrary directions for 3D,      the function projects q into n-t coordinates, computes the     eigendecomposition there, and then rotates back\n\nAliasing restrictions: from params the following fields are used:     Y, S2, Lambda, res_vals1, res_vals2,  w_vals_stencil,      w_vals2_stencil, v_vals, v_vals2, q_vals, q_vals2, nrm2, P\n\n\n\n"
},

{
    "location": "solver/euler/faceElementIntegrals.html#EulerEquationMod.calcLW2EntropyPenaltyIntegral-Union{Tuple{EulerEquationMod.ParamType{Tdim,:conservative,Tsol,Tres,Tmsh},SummationByParts.SparseFace,ODLCommonTools.Interface,AbstractArray{Tsol,2},AbstractArray{Tsol,2},AbstractArray{Tres,2},AbstractArray{Tmsh,2},AbstractArray{Tres,2},AbstractArray{Tres,2}}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tmsh where Tres where Tsol where Tdim",
    "page": "Face Element Integrals",
    "title": "EulerEquationMod.calcLW2EntropyPenaltyIntegral",
    "category": "method",
    "text": "Method for sparse faces.  See other method for details\n\nAliasing restrictions: params: v_vals, v_vals2, q_vals, q_vals2, A0, res_vals1, res_vals2                          A0, S2, Lambda, nrm2, P\n\n\n\n"
},

{
    "location": "solver/euler/faceElementIntegrals.html#EulerEquationMod.calcLWEntropyPenaltyIntegral-Union{Tuple{EulerEquationMod.ParamType{Tdim,:conservative,Tsol,Tres,Tmsh},SummationByParts.DenseFace,ODLCommonTools.Interface,AbstractArray{Tsol,2},AbstractArray{Tsol,2},AbstractArray{Tres,2},AbstractArray{Tmsh,2},AbstractArray{Tres,2},AbstractArray{Tres,2}}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tmsh where Tres where Tsol where Tdim",
    "page": "Face Element Integrals",
    "title": "EulerEquationMod.calcLWEntropyPenaltyIntegral",
    "category": "method",
    "text": "Calculate a term that provably dissipates (mathematical) entropy using a    an approximation to Lax-Wendroff type of dissipation.     This requires data from the left and right element volume nodes, rather than   face nodes for a regular face integral.\n\nNote that nrm_face must contain the scaled normal vector in x-y space   at the face nodes, and qL, qR, resL, and resR are the arrays for the   entire element, not just the face.\n\nThe approximation to Lax-Wendroff is the computation of\n\nfor i=1:Tdim     abs(niY_iS2_iLambda_iY_i.\')   end\n\nrather than computing the flux jacobian in the normal direction.\n\nAliasing restrictions: from params the following fields are used:     Y, S2, Lambda, res_vals1, res_vals2, res_vals3,  w_vals_stencil,      w_vals2_stencil, v_vals, v_vals2, q_vals, q_vals2\n\n\n\n"
},

{
    "location": "solver/euler/faceElementIntegrals.html#EulerEquationMod.getFaceElementFunctors-Tuple{Any,Any,EulerEquationMod.AbstractEulerData,Any}",
    "page": "Face Element Integrals",
    "title": "EulerEquationMod.getFaceElementFunctors",
    "category": "method",
    "text": "Populates the field(s) of the EulerData object with   FaceElementIntegralType functors as specified by the options   dictionary\n\nInputs\n\nmesh: an AbstractMesh\nsbp: an SBP operator\nopts: the options dictionary\n\nInputs/Outputs\n\neqn: the EulerData object\n\n\n\n"
},

{
    "location": "solver/euler/faceElementIntegrals.html#EulerEquationMod.getEntropyLFStab-Union{Tuple{EulerEquationMod.ParamType{Tdim,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1},AbstractArray{Tsol,1},AbstractArray{Tres,N} where N,AbstractArray{Tmsh,N} where N,AbstractArray{Tres,1}}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tdim where Tres where Tsol where Tmsh",
    "page": "Face Element Integrals",
    "title": "EulerEquationMod.getEntropyLFStab",
    "category": "method",
    "text": "This function computes the entropy dissipation term using Lax-Friedrich   type dissipation.  The term is evaluated using simple averaging of   qL and qR.  The term is subtracted off of F.\n\nThis function is dimension agnostic, but only works for conservative   variables.\n\nAliasing restrictions: params.q_vals3, see also getEntropyLFStab_inner\n\n\n\n"
},

{
    "location": "solver/euler/faceElementIntegrals.html#EulerEquationMod.getEntropyLFStab_inner-Union{Tuple{EulerEquationMod.ParamType{Tdim,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1},AbstractArray{Tsol,1},AbstractArray{Tsol,N} where N,AbstractArray{Tres,N} where N,AbstractArray{Tmsh,N} where N,AbstractArray{Tres,1}}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tdim where Tres where Tsol where Tmsh",
    "page": "Face Element Integrals",
    "title": "EulerEquationMod.getEntropyLFStab_inner",
    "category": "method",
    "text": "Updates the vector F with the stabilization term from Carpenter, Fisher,   Nielsen, Frankel, Entrpoy stable spectral collocation schemes for the    Navier-Stokes equatiosn: Discontinuous interfaces.  The term is subtracted   off from F.\n\nThe q_avg vector should some average of qL and qR, but the type of    averaging is left up to the user.\n\nThis function is agnostic to dimension, but only works for conservative   variables.\n\nAliasing: from params the following arrays are used: A0, v_vals               v_vals2.\n\n\n\n"
},

{
    "location": "solver/euler/faceElementIntegrals.html#EulerEquationMod.getIRA0-Union{Tuple{EulerEquationMod.ParamType{2,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tsol,1},AbstractArray{Tsol,2}}, Tuple{Tsol}} where Tsol",
    "page": "Face Element Integrals",
    "title": "EulerEquationMod.getIRA0",
    "category": "method",
    "text": "Computes dq/dv, where q are the conservative variables and v are the   IR entropy variables.  This is equiavlent to calcA0 scaled by gamma_1,   but computed from the conservative variables, which is much less expensive.\n\nMethods are available for 2 and 3 dimensions.   A0 is overwritten with the result\n\n\n\n"
},

{
    "location": "solver/euler/faceElementIntegrals.html#EulerEquationMod.getLambdaMaxSimple-Union{Tuple{EulerEquationMod.ParamType{Tdim,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tsol,1},AbstractArray{Tsol,1},AbstractArray{Tmsh,1}}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tsol}} where Tdim where Tmsh where Tsol",
    "page": "Face Element Integrals",
    "title": "EulerEquationMod.getLambdaMaxSimple",
    "category": "method",
    "text": "Calculates the maximum magnitude eigenvalue of the Euler flux    jacobian at the arithmatic average of two states.\n\nThis functions works in both 2D and 3D   Inputs:     params:  ParamType, conservative variable     qL: left state     qR: right state     dir: direction vector (does not have to be unit vector)\n\nOutputs:     lambda_max: eigenvalue of maximum magnitude\n\nAliasing restrictions: params.q_vals3 must be unused\n\n\n\n"
},

{
    "location": "solver/euler/faceElementIntegrals.html#Functions-1",
    "page": "Face Element Integrals",
    "title": "Functions",
    "category": "section",
    "text": "  Modules = [EulerEquationMod]\n  Order = [:function]\n  Pages = [\"euler/faceElementIntegrals.jl\", \"euler/IR_stab.jl\"]"
},

{
    "location": "solver/euler/faceElementIntegrals.html#EulerEquationMod.ECFaceIntegral",
    "page": "Face Element Integrals",
    "title": "EulerEquationMod.ECFaceIntegral",
    "category": "type",
    "text": "Entropy conservative term only\n\n\n\n"
},

{
    "location": "solver/euler/faceElementIntegrals.html#EulerEquationMod.ELFPenaltyFaceIntegral",
    "page": "Face Element Integrals",
    "title": "EulerEquationMod.ELFPenaltyFaceIntegral",
    "category": "type",
    "text": "Lax-Friedrich entropy penalty term only\n\n\n\n"
},

{
    "location": "solver/euler/faceElementIntegrals.html#EulerEquationMod.ELW2PenaltyFaceIntegral",
    "page": "Face Element Integrals",
    "title": "EulerEquationMod.ELW2PenaltyFaceIntegral",
    "category": "type",
    "text": "Lax-Wendroff entropy penalty term only\n\n\n\n"
},

{
    "location": "solver/euler/faceElementIntegrals.html#EulerEquationMod.ELWPenaltyFaceIntegral",
    "page": "Face Element Integrals",
    "title": "EulerEquationMod.ELWPenaltyFaceIntegral",
    "category": "type",
    "text": "Approximate Lax-Wendroff entropy penalty term only\n\n\n\n"
},

{
    "location": "solver/euler/faceElementIntegrals.html#EulerEquationMod.ESLFFaceIntegral",
    "page": "Face Element Integrals",
    "title": "EulerEquationMod.ESLFFaceIntegral",
    "category": "type",
    "text": "Entropy conservative integral + Lax-Friedrich penalty\n\n\n\n"
},

{
    "location": "solver/euler/faceElementIntegrals.html#EulerEquationMod.ESLW2FaceIntegral",
    "page": "Face Element Integrals",
    "title": "EulerEquationMod.ESLW2FaceIntegral",
    "category": "type",
    "text": "Entropy conservative integral + Lax-Wendroff penalty\n\n\n\n"
},

{
    "location": "solver/euler/faceElementIntegrals.html#EulerEquationMod.ESLWFaceIntegral",
    "page": "Face Element Integrals",
    "title": "EulerEquationMod.ESLWFaceIntegral",
    "category": "type",
    "text": "Entropy conservative integral + approximate Lax-Wendroff penalty\n\n\n\n"
},

{
    "location": "solver/euler/faceElementIntegrals.html#Flux-Functors-1",
    "page": "Face Element Integrals",
    "title": "Flux Functors",
    "category": "section",
    "text": "  Modules = [EulerEquationMod]\n  Order = [:constant, :type]\n  Pages = [\"euler/faceElementIntegrals.jl\"]"
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
    "text": "  CurrentModule = EulerEquationModThis page describes the functions that impose boundarya conditions. The boundary conditions are imposed weakly, using a penalty between the desired state (the boundary condition value) and the current state (the solution)."
},

{
    "location": "solver/euler/bc.html#EulerEquationMod.evalBoundaryIntegrals",
    "page": "Boundary Integrals",
    "title": "EulerEquationMod.evalBoundaryIntegrals",
    "category": "function",
    "text": "This function evaluates the boundary integrals in the Euler equations by   calling the appropriate SBP function on eqn.bndryflux, which must be populated   before calling this function.  eqn.res is updated with the result\n\nThis is a mid level function\n\n\n\n"
},

{
    "location": "solver/euler/bc.html#Entry-Points-1",
    "page": "Boundary Integrals",
    "title": "Entry Points",
    "category": "section",
    "text": "evalBoundaryIntegrals"
},

{
    "location": "solver/euler/bc.html#EulerEquationMod.calcBoundaryFlux-Union{Tuple{ODLCommonTools.AbstractCGMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim where Tres,ODLCommonTools.BCType,UnitRange,AbstractArray{ODLCommonTools.Boundary,1},AbstractArray{Tres,3}}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tres where Tsol where Tmsh",
    "page": "Boundary Integrals",
    "title": "EulerEquationMod.calcBoundaryFlux",
    "category": "method",
    "text": "EulerEquationMod.calcBoundaryFlux\n\nThis function calculates the boundary flux for the portion of the boundary   with a particular boundary condition.  The eqn.q are converted to   conservative variables if needed\n\nInputs:   mesh : AbstractMesh   sbp : AbstractSBP   eqn : EulerEquation   functor : a callable object that calculates the boundary flux at a node   idx_range: the Range describing which Boundaries have the current BC   bndry_facenums:  An array with elements of type Boundary that tell which                    element faces have the boundary condition   Outputs:   bndryflux : the array to store the boundary flux, corresponds to               bndry_facenums\n\nThe functor must have the signature   functor( q, aux_vars, x, nrm_xy, bndryflux_i, eqn.params)   where q are the conservative variables.   where all arguments (except params) are vectors of values at a node.\n\nparams is the ParamType associated with the the EulerEquation object   nrm = mesh.sbpface.normal[:, current_node]\n\nThis is a mid level function.\n\n\n\n"
},

{
    "location": "solver/euler/bc.html#EulerEquationMod.calcBoundaryFlux_nopre-Union{Tuple{ODLCommonTools.AbstractDGMesh{Tmsh},ODLCommonTools.AbstractDGMesh{Tmsh},SummationByParts.AbstractSBP,SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim,ODLCommonTools.BCType,UnitRange,AbstractArray{ODLCommonTools.Boundary,1}}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tres where Tsol where Tmsh",
    "page": "Boundary Integrals",
    "title": "EulerEquationMod.calcBoundaryFlux_nopre",
    "category": "method",
    "text": "Staggered grid version\n\n\n\n"
},

{
    "location": "solver/euler/bc.html#EulerEquationMod.calcBoundaryFlux_nopre-Union{Tuple{ODLCommonTools.AbstractDGMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim,ODLCommonTools.BCType,UnitRange,AbstractArray{ODLCommonTools.Boundary,1}}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tres where Tsol where Tmsh",
    "page": "Boundary Integrals",
    "title": "EulerEquationMod.calcBoundaryFlux_nopre",
    "category": "method",
    "text": "Like calcBoundaryFlux, but performs the integration and updates res rather   than storing the flux.\n\n\n\n"
},

{
    "location": "solver/euler/bc.html#EulerEquationMod.getBCFluxes-Tuple{ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,EulerEquationMod.EulerData,Any}",
    "page": "Boundary Integrals",
    "title": "EulerEquationMod.getBCFluxes",
    "category": "method",
    "text": "EulerEquationMod.getBCFluxes\n\nThis function calls other functions to calculate the boundary fluxes, passing   them pieces of the array needed.  This populates eqn.bndryflux.  It also   calls writeBoundary() to do any requested output.  If the options dictionary   specifies not to precompute the boundary flux, this function will do the   integration as well and update eqn.res.\n\nThis is a mid level function.\n\n\n\n"
},

{
    "location": "solver/euler/bc.html#EulerEquationMod.getBCFunctors-Tuple{ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,EulerEquationMod.EulerData,Any}",
    "page": "Boundary Integrals",
    "title": "EulerEquationMod.getBCFunctors",
    "category": "method",
    "text": "EulerEquationMod.getBCFunctors\n\nThis function uses the opts dictionary to populate mesh.bndry_funcs with   the the functors\n\nfunc(params::ParamType,\n     q::AbstractArray{Tsol,1},\n     aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},\n     nrm_xy::AbstractArray{Tmsh,1},\n     bndryflux::AbstractArray{Tres, 1},\n     bndry::BoundaryNode=NullBoundaryNode)\n\nThis is a high level function.\n\n\n\n"
},

{
    "location": "solver/euler/bc.html#EulerEquationMod.interpolateBoundary-Tuple{ODLCommonTools.AbstractDGMesh,Any,Any,Any,AbstractArray{T,3} where T,AbstractArray{T,3} where T}",
    "page": "Boundary Integrals",
    "title": "EulerEquationMod.interpolateBoundary",
    "category": "method",
    "text": "EulerEquationMod.interpolateBoundary\n\nInterpolates the solution variables to the exterior boundary of the mesh   and calculates any additional quantities at the boundary of the mesh.   DG only\n\nInputs:     mesh: an AbstractDGMesh     sbp     eqn     opts     q : the 3D array of solution variables for all elements, numdofpernode x         numNodesPerElement x numEl\n\nInputs/Outputs:     q_bndry: the array to be populated with the solution interpolated to              the boundary, numdofpernode x numNodesPerFace x num boundary faces\n\neqn.aux_vars_bndry is also populated\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/bc.html#EulerEquationMod.writeBoundary-NTuple{4,Any}",
    "page": "Boundary Integrals",
    "title": "EulerEquationMod.writeBoundary",
    "category": "method",
    "text": "EulerEquationMod.writeBoundary\n\nThis function writes information about the boundary faces and fluxes to files.   It is controlled by the input argument writeboundary, of type Bool.\n\nIt generates the files:     * boundaryfaces.dat : writes mesh.bndryfaces, an array with eltype Boundary                           to a file, one element per line     * boundaryflux.dat  : writes the element, local node number and boundary                           flux to a line in a human readable format     * boundaryflux2.dat : writes the real part ofmesh.bndryflux to space                           delimited file\n\nThis is a high level function.\n\n\n\n"
},

{
    "location": "solver/euler/bc.html#Functions-1",
    "page": "Boundary Integrals",
    "title": "Functions",
    "category": "section",
    "text": "  Modules = [EulerEquationMod]\n  Order = [:function]\n  Pages = [\"euler/bc.jl\"]"
},

{
    "location": "solver/euler/bc.html#EulerEquationMod.BCDict",
    "page": "Boundary Integrals",
    "title": "EulerEquationMod.BCDict",
    "category": "constant",
    "text": "Maps boundary conditions names to the functor objects.   Each functor should be callable with the signature\n\n\n\n"
},

{
    "location": "solver/euler/bc.html#EulerEquationMod.NullBoundaryNode",
    "page": "Boundary Integrals",
    "title": "EulerEquationMod.NullBoundaryNode",
    "category": "constant",
    "text": "Null boundary node (all fields zero).  This is useful as a default value   for a functiona argument (if the argument is unused).\n\n\n\n"
},

{
    "location": "solver/euler/bc.html#EulerEquationMod.BoundaryNode",
    "page": "Boundary Integrals",
    "title": "EulerEquationMod.BoundaryNode",
    "category": "type",
    "text": "Type that identified a particular node on a particular face of   a particular boundary.  It also contains the index of the face within   the array of Boundaries with the same boundary condition\n\nFields\n\nelement: the element the face is part of\nface: the local face number\nfaceidx: the index within the array of Boundaries with this BC\nnode: the node on the face\n\n\n\n"
},

{
    "location": "solver/euler/bc.html#EulerEquationMod.FreeStreamBC_dAlpha",
    "page": "Boundary Integrals",
    "title": "EulerEquationMod.FreeStreamBC_dAlpha",
    "category": "type",
    "text": "EulerEquationMod.FreeStreamBC_dAlpha <: BCTypes\n\nThis functor uses the Roe solver to calculate the flux for a boundary   state corresponding to the free stream velocity, using rho_free, Ma, aoa, and E_free\n\nThis is a low level functor\n\nArguments\n\nobj : Object of type BCType used for multiple dispatch. Every new boundary          condition needs to have its own type and entered in BCDict\nq   : Solution variable\naux_vars : Auxiliary variables\nx        : physical coordinates of the SBP node\nnrm_xy      : scaled normal vector in x-y space\nbndryflux : Computed flux value at the boundary\n\n\n\n"
},

{
    "location": "solver/euler/bc.html#EulerEquationMod.FreeStreamBC_revm",
    "page": "Boundary Integrals",
    "title": "EulerEquationMod.FreeStreamBC_revm",
    "category": "type",
    "text": "###EulerEquationMod.FreeStreamBC_revm\n\nReverse mode for FreeStreamBC.\n\nInputs\n\nobj : Type of the Boundary condition being evaluated. Its a subtype of         BCType_revm\nq   : Solution variable\naux_vars : Auxiliary variables\nx     : Node coordinates\nnrm_xy   : scaled normal vector in x-y space\nbndryflux_bar : Input flux value seed that is used to compute the reverse                   mode derivative.\nparams        : equation object parameters\n\nOutput\n\nnrm_bar : Derivative of bndryflux_bar w.r.t the mapping jacobian\n\n\n\n"
},

{
    "location": "solver/euler/bc.html#EulerEquationMod.allOnesBC",
    "page": "Boundary Integrals",
    "title": "EulerEquationMod.allOnesBC",
    "category": "type",
    "text": "EulerEquationMod.allOnesBC <: BCTypes\n\nThis functor uses the Roe solver to calculate the flux for a boundary   state where all the conservative variables have a value 1.0\n\nThis is a low level functor\n\n\n\n"
},

{
    "location": "solver/euler/bc.html#EulerEquationMod.allZerosBC",
    "page": "Boundary Integrals",
    "title": "EulerEquationMod.allZerosBC",
    "category": "type",
    "text": "EulerEquationMod.allZerosBC <: BCTypes\n\nThis functor uses the Roe solver to calculate the flux for a boundary   state where all the conservative variables have a value 0.0\n\nThis is a low level functor\n\n\n\n"
},

{
    "location": "solver/euler/bc.html#EulerEquationMod.isentropicVortexBC",
    "page": "Boundary Integrals",
    "title": "EulerEquationMod.isentropicVortexBC",
    "category": "type",
    "text": "EulerEquationMod.isentropicVortexBC <: BCTypes\n\nThis type and the associated call method define a functor to calculate   the flux using the Roe Solver using the exact InsentropicVortex solution   as boundary state.  See calcBoundaryFlux for the arguments all functors   must support.\n\nThis is a low level functor.\n\nArguments\n\nobj : Object of type BCType used for multiple dispatch. Every new boundary          condition needs to have its own type and entered in BCDict\nq   : Solution variable\naux_vars : Auxiliary variables\nx        : physical coordinates of the SBP node\nnrm_xy      : sclaed face normal in physical space\nbndryflux : Computed flux value at the boundary\n\n\n\n"
},

{
    "location": "solver/euler/bc.html#EulerEquationMod.isentropicVortexBC_physical",
    "page": "Boundary Integrals",
    "title": "EulerEquationMod.isentropicVortexBC_physical",
    "category": "type",
    "text": "EulerEquationMod.isentropicVortexBC_physical <: BCTypes\n\nThis type and the associated call method define a functor to calculate   the actual Euler flux  using the exact InsentropicVortex solution   as boundary state.  See calcBoundaryFlux for the arguments all functors   must support.\n\nThis is a low level functor.\n\n\n\n"
},

{
    "location": "solver/euler/bc.html#EulerEquationMod.isentropicVortexBC_revm",
    "page": "Boundary Integrals",
    "title": "EulerEquationMod.isentropicVortexBC_revm",
    "category": "type",
    "text": "###EulerEquationMod.isentropicVortexBC_revm\n\nReverse mode for isentropicVortexBC.\n\nInputs\n\nobj : Type of the Boundary condition being evaluated. Its a subtype of         BCType_revm\nq   : Solution variable\naux_vars : Auxiliary variables\nx     : Node coordinates\nnrm   : scaled normal vector in x-y space\nbndryflux_bar : Input flux value seed that is used to compute the reverse                   mode derivative.\nparams        : equation object parameters\n\nOutput\n\nnrm_bar : Derivative of flux w.r.t the nrm\n\n\n\n"
},

{
    "location": "solver/euler/bc.html#EulerEquationMod.noPenetrationBC_revm",
    "page": "Boundary Integrals",
    "title": "EulerEquationMod.noPenetrationBC_revm",
    "category": "type",
    "text": "###EulerEquationMod.noPenetrationBC_revm\n\nReverse mode for noPenetrationBC.\n\nInput\n\nobj : Type of the Boundary condition being evaluated. Its a subtype of         BCType_revm\nq   : Solution variable\naux_vars : Auxiliary variables\nx     : Node coordinates\ndxidx : Mapping jacobian matrix for the SBP node\nnrm   : sbpface normal vector\nbndryflux_bar : Input flux value seed that is used to compute the reverse                   mode derivative.\nparams        : equation object parameters\n\n\n\n"
},

{
    "location": "solver/euler/bc.html#EulerEquationMod.reanalysisBC",
    "page": "Boundary Integrals",
    "title": "EulerEquationMod.reanalysisBC",
    "category": "type",
    "text": "This is a special boundary condition used for imposing a numerical    solution as a boundary condition.\n\nFields\n\n* bc_vals: an array of size numDofPerNode x numNodesPerFace x numFaces with\n           this boundary condition on it.  This array can be accessed\n           using the `faceidx` field of [`BoundaryNode`](@ref).\n\nThis BC is special because it has to store the boundary values it is imposing.   Whenever this boundary condition is needed, the user should construct a   new object, with the data inside it, and then register it just before   constructing the EulerData object.\n\n\n\n"
},

{
    "location": "solver/euler/bc.html#EulerEquationMod.unsteadyVortexBC",
    "page": "Boundary Integrals",
    "title": "EulerEquationMod.unsteadyVortexBC",
    "category": "type",
    "text": "EulerEquationMod.unsteadyVortexBC <: BCTypes\n\nThis type and the associated call method define a functor to calculate   the flux using the Roe Solver using the exact InsentropicVortex solution   as boundary state.  See calcBoundaryFlux for the arguments all functors   must support.\n\nThis is a low level functor.\n\nArguments\n\nobj : Object of type BCType used for multiple dispatch. Every new boundary          condition needs to have its own type and entered in BCDict\nq   : Solution variable\naux_vars : Auxiliary variables\nx        : physical coordinates of the SBP node\ndxidx    : Mapping jacobian matrix for the SBP node\nnrm_xy      : SBP face normal\nbndryflux : Computed flux value at the boundary\n\n\n\n"
},

{
    "location": "solver/euler/bc.html#Boundary-Condition-Functors-1",
    "page": "Boundary Integrals",
    "title": "Boundary Condition Functors",
    "category": "section",
    "text": "Each functor defines a different type of boundary condition. Because the equation is solved in the weak form, the functors compute the flux associated with the penalty that imposes the boundary condition.  Modules = [EulerEquationMod]\n  Order = [:constant, :type]\n  Pages = [\"euler/bc.jl\"]"
},

{
    "location": "solver/euler/bc_diff.html#",
    "page": "Boundary Integrals Jacobian",
    "title": "Boundary Integrals Jacobian",
    "category": "page",
    "text": ""
},

{
    "location": "solver/euler/bc_diff.html#Boundary-Integrals-Differentiated-1",
    "page": "Boundary Integrals Jacobian",
    "title": "Boundary Integrals Differentiated",
    "category": "section",
    "text": "  CurrentModule = EulerEquationModThis page describes how the boundary integral contribution to the Jacobian is calculated."
},

{
    "location": "solver/euler/bc_diff.html#EulerEquationMod.evalBoundaryIntegrals_diff",
    "page": "Boundary Integrals Jacobian",
    "title": "EulerEquationMod.evalBoundaryIntegrals_diff",
    "category": "function",
    "text": "Differentiated version of evalBoundaryIntegrals.\n\nInputs\n\nmesh\nsbp\neqn\nopts\nassembler.\n\n\n\n"
},

{
    "location": "solver/euler/bc_diff.html#Entry-Point-1",
    "page": "Boundary Integrals Jacobian",
    "title": "Entry Point",
    "category": "section",
    "text": "evalBoundaryIntegrals_diff"
},

{
    "location": "solver/euler/bc_diff.html#EulerEquationMod.calcBoundaryFlux_nopre_diff-Union{Tuple{ODLCommonTools.AbstractDGMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol1,Tres1,Tdim,var_type} where var_type where Tdim,ODLCommonTools.BCType,UnitRange,AbstractArray{ODLCommonTools.Boundary,1},PDESolver.AssembleElementData}, Tuple{Tmsh}, Tuple{Tres1}, Tuple{Tsol1}} where Tmsh where Tres1 where Tsol1",
    "page": "Boundary Integrals Jacobian",
    "title": "EulerEquationMod.calcBoundaryFlux_nopre_diff",
    "category": "method",
    "text": "Computes the jacobian of the boundary condition terms for a given   BC and assembles it into the Jacobian.  See calcBoundaryFlux_nopre.   Note that this function does not depend on any quatiites being precomputed.\n\nInputs\n\nmesh\nsbp\neqn\nfunctor: this is the regular BCType functor, not BCType_diff\nidx_range: range of indices in mesh.bndryfaces that bndry_facenums            came from\nbndry_facenums: array of Boundary objects to compute the BC for\nassembler: AssembleElementData needed to assemble the element level            jacobian into the system jacobian\n\nImplementation Notes\n\nCurrently this computes the jacobian of the flux with respect to the   face via complex step (using eqn.params_complex).  This should be replaced   with dual numbers eventually\n\n\n\n"
},

{
    "location": "solver/euler/bc_diff.html#EulerEquationMod.getBCFluxes_diff-Tuple{ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,EulerEquationMod.EulerData,Any,PDESolver.AssembleElementData}",
    "page": "Boundary Integrals Jacobian",
    "title": "EulerEquationMod.getBCFluxes_diff",
    "category": "method",
    "text": "Computes the jacobian of the boundary condition terms computed by   getBCFluxes and assembles them into the Jacobian\n\nInputs\n\nmesh\nsbp\neqn\nopts\nassembler: AssembleElementData\n\n\n\n"
},

{
    "location": "solver/euler/bc_diff.html#Functions-1",
    "page": "Boundary Integrals Jacobian",
    "title": "Functions",
    "category": "section",
    "text": "  Modules = [EulerEquationMod]\n  Order = [:function]\n  Pages = [\"euler/bc_diff.jl\"]"
},

{
    "location": "solver/euler/bc_diff.html#Boundary-Condition-Functors-1",
    "page": "Boundary Integrals Jacobian",
    "title": "Boundary Condition Functors",
    "category": "section",
    "text": "The same functors used to compute the boundary integral are also used to compute the Jacobian contribution."
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
    "text": "  CurrentModule = EulerEquationModThis page describes the functions that apply initial conditions. The initial condition is loaded into the solution vector (not the 3D array). Users should call the appropriate function if they want to transfer the solution from the vector to the 3D array.It is important to note that the initial conditions functions do depend on the current time t.  In general, an initial condition function is really a function that computes an exact solution at a given time t. Applying an initial condition is a special case of this (because t = 0). Having the initial condition functions depend on time is used to calculate errors at the end of a simulation (for cases where the exact solution is known).Often, the initial condition function will call a function in common funcs to evalute the exact solution at each node on the mesh.Unlike boundary conditions, which are called many times during the execution of the code, initial conditions are called infrequently, so we do not create functors for them."
},

{
    "location": "solver/euler/ic.html#EulerEquationMod.ICDict",
    "page": "Initial Conditions",
    "title": "EulerEquationMod.ICDict",
    "category": "constant",
    "text": "Map IC names to functions.  Generally the name is the same as the function   name\n\n\n\n"
},

{
    "location": "solver/euler/ic.html#Accessing-Initial-Conditions-1",
    "page": "Initial Conditions",
    "title": "Accessing Initial Conditions",
    "category": "section",
    "text": "  Modules = [EulerEquationMod]\n  Order = [:constant]\n  Pages = [\"euler/ic.jl\"]"
},

{
    "location": "solver/euler/ic.html#EulerEquationMod.ICChannelMMS-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},Any,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim where Tres,Any,AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Initial Conditions",
    "title": "EulerEquationMod.ICChannelMMS",
    "category": "method",
    "text": "Initial condition of channel MMS\n\n\n\n"
},

{
    "location": "solver/euler/ic.html#EulerEquationMod.ICExp-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},Any,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim where Tres,Any,AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Initial Conditions",
    "title": "EulerEquationMod.ICExp",
    "category": "method",
    "text": "Assigns exp(kxy*z) as the initial condition, of each node, where k is    the index of the degree of freedom of the node\n\n\n\n"
},

{
    "location": "solver/euler/ic.html#EulerEquationMod.ICFile-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP{Tsbp},EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim where Tres,Any,AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsbp}, Tuple{Tsol}} where Tsol where Tsbp where Tmsh",
    "page": "Initial Conditions",
    "title": "EulerEquationMod.ICFile",
    "category": "method",
    "text": "EulerEquationMod.ICFile\n\nThis function reads a vector from a file on disk and set the solution to it.   The vector must contain the same number of entries as there are degrees of    freedom in the mesh. \n\nThis function is useful for things like restarting from a checkpoint.   In this case, the file should be the output of writedlm(eqn.q).  The degree    of freedom number must be the same for both simulation for this to work (the    file contains no degree of freedom number information).\n\nInputs:     mesh     sbp     eqn     opts\n\nInputs/Outputs:      u0: vector to populate with the solution\n\nAliasing restrictions: none.\n\n\n\n"
},

{
    "location": "solver/euler/ic.html#EulerEquationMod.ICFreeStream-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP{Tsbp},EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim where Tres,Any,AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsbp}, Tuple{Tsol}} where Tsol where Tsbp where Tmsh",
    "page": "Initial Conditions",
    "title": "EulerEquationMod.ICFreeStream",
    "category": "method",
    "text": "EulerEquationMod.ICFreeStream\n\nSets all components of the solution to the free stream condition according   to the angle of attack and and Mach number.\n\nInputs:     mesh     sbp     eqn     opts\n\nInputs/Outputs:      u0: vector to populate with the solution\n\nAliasing restrictions: none.\n\n\n\n"
},

{
    "location": "solver/euler/ic.html#EulerEquationMod.ICFreeStream0-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP{Tsbp},EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim where Tres,Any,AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsbp}, Tuple{Tsol}} where Tsol where Tsbp where Tmsh",
    "page": "Initial Conditions",
    "title": "EulerEquationMod.ICFreeStream0",
    "category": "method",
    "text": "Like calcFreeStream, but uses calcFreeStream0 instead of   calcFreeStream.\n\n\n\n"
},

{
    "location": "solver/euler/ic.html#EulerEquationMod.ICInvChannel-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP{Tsbp},EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim where Tres,Any,AbstractArray{Tsol,N} where N}, Tuple{Tmsh}, Tuple{Tsbp}, Tuple{Tsol}} where Tsol where Tsbp where Tmsh",
    "page": "Initial Conditions",
    "title": "EulerEquationMod.ICInvChannel",
    "category": "method",
    "text": "Initial condition for SU2 bump in inviscid channel case.  See also   The subsonic inflow and subsonic outflow boundary conditions.\n\n\n\n"
},

{
    "location": "solver/euler/ic.html#EulerEquationMod.ICIsentropicVortex-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP{Tsbp},EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim where Tres,Any,AbstractArray{Tsol,N} where N}, Tuple{Tmsh}, Tuple{Tsbp}, Tuple{Tsol}} where Tsol where Tsbp where Tmsh",
    "page": "Initial Conditions",
    "title": "EulerEquationMod.ICIsentropicVortex",
    "category": "method",
    "text": "EulerEquationMod.ICIsentropicVortex\n\nSets the solution to the steady isentropic vortex solution.\n\nInputs: mesh     sbp     eqn     opts\n\nInputs/Outputs:      u0: vector to populate with the solution\n\nAliasing restrictions: none.\n\n\n\n"
},

{
    "location": "solver/euler/ic.html#EulerEquationMod.ICIsentropicVortexWithNoise-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP{Tsbp},EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim where Tres,Any,AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsbp}, Tuple{Tsol}} where Tsol where Tsbp where Tmsh",
    "page": "Initial Conditions",
    "title": "EulerEquationMod.ICIsentropicVortexWithNoise",
    "category": "method",
    "text": "EulerEquationMod.ICIsentropicVortexWithNoise\n\nSets the solution to the steady isentropic vortex solution plus    a small random noise component.\n\nInputs:     mesh     sbp     eqn     opts\n\nInputs/Outputs:      u0: vector to populate with the solution\n\nAliasing restrictions: none.\n\n\n\n"
},

{
    "location": "solver/euler/ic.html#EulerEquationMod.ICOnes-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP{Tsbp},EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim where Tres,Any,AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsbp}, Tuple{Tsol}} where Tsol where Tsbp where Tmsh",
    "page": "Initial Conditions",
    "title": "EulerEquationMod.ICOnes",
    "category": "method",
    "text": "EulerEquationMod.ICOnes\n\nSets all components of the solution to 1.0\n\nInputs:     mesh     sbp     eqn     opts\n\nInputs/Outputs:      u0: vector to populate with the solution\n\nAliasing restrictions: none.\n\n\n\n"
},

{
    "location": "solver/euler/ic.html#EulerEquationMod.ICPeriodicMMS-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},Any,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim where Tres,Any,AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Initial Conditions",
    "title": "EulerEquationMod.ICPeriodicMMS",
    "category": "method",
    "text": "Writes calcPeriodicMMS to the initial condition vector u0\n\n\n\n"
},

{
    "location": "solver/euler/ic.html#EulerEquationMod.ICRho1E2-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP{Tsbp},EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim where Tres,Any,AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsbp}, Tuple{Tsol}} where Tsol where Tsbp where Tmsh",
    "page": "Initial Conditions",
    "title": "EulerEquationMod.ICRho1E2",
    "category": "method",
    "text": "EulerEquationMod.ICRho1E2\n\nSets all density values to 1.0 and energy values to 2.0\n\nInputs:     mesh     sbp     eqn     opts\n\nInputs/Outputs:      u0: vector to populate with the solution\n\nAliasing restrictions: none.\n\n\n\n"
},

{
    "location": "solver/euler/ic.html#EulerEquationMod.ICRho1E2U1VW0-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP{Tsbp},EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim where Tres,Any,AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsbp}, Tuple{Tsol}} where Tsol where Tsbp where Tmsh",
    "page": "Initial Conditions",
    "title": "EulerEquationMod.ICRho1E2U1VW0",
    "category": "method",
    "text": "EulerEquationMod.ICRho1E2U1VW0\n\nSets the density values 1.0, x momentum to 1.0,    v & w momenta to 0.0, and energy to 2.0 at a node.\n\nIt should work for 2D and 3D meshes.\n\nInputs:     mesh     sbp     eqn     opts\n\nInputs/Outputs:      u0: vector to populate with the solution\n\nAliasing restrictions: none.\n\n\n\n"
},

{
    "location": "solver/euler/ic.html#EulerEquationMod.ICRho1E2U3-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP{Tsbp},EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim where Tres,Any,AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsbp}, Tuple{Tsol}} where Tsol where Tsbp where Tmsh",
    "page": "Initial Conditions",
    "title": "EulerEquationMod.ICRho1E2U3",
    "category": "method",
    "text": "EulerEquationMod.ICRho1E2U3\n\nSets all components density values to 1.0, x and y momenta to 0.35355, and   energy to 2.0\n\nInputs:     mesh     sbp     eqn     opts\n\nInputs/Outputs:      u0: vector to populate with the solution\n\nAliasing restrictions: none.\n\n\n\n"
},

{
    "location": "solver/euler/ic.html#EulerEquationMod.ICSedovExplosion-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},Any,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim where Tres,Any,AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Initial Conditions",
    "title": "EulerEquationMod.ICSedovExplosion",
    "category": "method",
    "text": "Initial for square wave in 2D\n\n\n\n"
},

{
    "location": "solver/euler/ic.html#EulerEquationMod.ICSquare1D-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},Any,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim where Tres,Any,AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Initial Conditions",
    "title": "EulerEquationMod.ICSquare1D",
    "category": "method",
    "text": "Initial for square wave in 1D\n\n\n\n"
},

{
    "location": "solver/euler/ic.html#EulerEquationMod.ICSquare2D-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},Any,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim where Tres,Any,AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Initial Conditions",
    "title": "EulerEquationMod.ICSquare2D",
    "category": "method",
    "text": "Initial for square wave in 2D\n\n\n\n"
},

{
    "location": "solver/euler/ic.html#EulerEquationMod.ICTaylorGreen-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},Any,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim where Tres,Any,AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Initial Conditions",
    "title": "EulerEquationMod.ICTaylorGreen",
    "category": "method",
    "text": "This function applies the initial condition for the Taylor Green vortex,   using the constants in Gassner, Winters, and Kopriva\'s Split form Nodal   DG paper\n\n\n\n"
},

{
    "location": "solver/euler/ic.html#EulerEquationMod.ICUnsteadyVortex-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP{Tsbp},EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim where Tres,Any,AbstractArray{Tsol,N} where N}, Tuple{Tmsh}, Tuple{Tsbp}, Tuple{Tsol}} where Tsol where Tsbp where Tmsh",
    "page": "Initial Conditions",
    "title": "EulerEquationMod.ICUnsteadyVortex",
    "category": "method",
    "text": "EulerEquationMod.ICUnsteadyVortex\n\nSets the solution to the unsteady vortex problem.  eqn.params.t is used to   determine what time to use for the solution.\n\nInputs:     mesh     sbp     eqn     opts\n\nInputs/Outputs:      u0: vector to populate with the solution\n\nAliasing restrictions: none.\n\n\n\n"
},

{
    "location": "solver/euler/ic.html#EulerEquationMod.ICUnsteadyVortex2-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP{Tsbp},EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim where Tres,Any,AbstractArray{Tsol,N} where N}, Tuple{Tmsh}, Tuple{Tsbp}, Tuple{Tsol}} where Tsol where Tsbp where Tmsh",
    "page": "Initial Conditions",
    "title": "EulerEquationMod.ICUnsteadyVortex2",
    "category": "method",
    "text": "Vortex travelling at an angle\n\n\n\n"
},

{
    "location": "solver/euler/ic.html#EulerEquationMod.ICZero-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP{Tsbp},EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim where Tres,Any,AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsbp}, Tuple{Tsol}} where Tsol where Tsbp where Tmsh",
    "page": "Initial Conditions",
    "title": "EulerEquationMod.ICZero",
    "category": "method",
    "text": "EulerEquationMod.ICZero\n\nSets all components of the solution to zero\n\nInputs:     mesh     sbp     eqn     opts\n\nInputs/Outputs:      u0: vector to populate with the solution\n\nAliasing restrictions: none.\n\n\n\n"
},

{
    "location": "solver/euler/ic.html#EulerEquationMod.ICsmoothHeaviside-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP{Tsbp},EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim where Tres,Any,AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsbp}, Tuple{Tsol}} where Tsol where Tsbp where Tmsh",
    "page": "Initial Conditions",
    "title": "EulerEquationMod.ICsmoothHeaviside",
    "category": "method",
    "text": "EulerEquationMod.ICZero\n\nSets the density to the smooth Heaviside function, all other components to   zero.\n\nInputs:     mesh     sbp     eqn     opts\n\nInputs/Outputs:      u0: vector to populate with the solution\n\nAliasing restrictions: none.\n\n\n\n"
},

{
    "location": "solver/euler/ic.html#EulerEquationMod.ICsmoothHeavisideder-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP{Tsbp},EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim where Tres,Any,AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsbp}, Tuple{Tsol}} where Tsol where Tsbp where Tmsh",
    "page": "Initial Conditions",
    "title": "EulerEquationMod.ICsmoothHeavisideder",
    "category": "method",
    "text": "EulerEquationMod.ICsmoothHeavisideder\n\nSets the density to the derivative of the smooth Heaviside function, all    other components to zero.\n\nInputs:     mesh     sbp     eqn     opts\n\nInputs/Outputs:      u0: vector to populate with the solution\n\nAliasing restrictions: none.\n\n\n\n"
},

{
    "location": "solver/euler/ic.html#Initial-Condition-Functions-1",
    "page": "Initial Conditions",
    "title": "Initial Condition Functions",
    "category": "section",
    "text": "TODO: write  a macro for calling common funcs repeatedly  Modules = [EulerEquationMod]\n  Order = [:function]\n  Pages = [\"euler/ic.jl\"]"
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
    "text": "  CurrentModule = EulerEquationModThis page describes how to apply source terms."
},

{
    "location": "solver/euler/source.html#EulerEquationMod.evalSourceTerm",
    "page": "Source Term",
    "title": "EulerEquationMod.evalSourceTerm",
    "category": "function",
    "text": "EulerEquationMod.evalSourceTerm\n\nThis function performs all the actions necessary to update eqn.res   with the source term.  The source term is stored in eqn.src_func.  It is   an abstract field, so it cannot be accessed (performantly) directly, so   it is passed to an inner function.\n\nInputs:\n\nmesh : Abstract mesh type\nsbp  : Summation-by-parts operator\neqn  : Euler equation object\nopts : options dictonary\n\nOutputs: none\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/source.html#Entry-Point-1",
    "page": "Source Term",
    "title": "Entry Point",
    "category": "section",
    "text": "evalSourceTerm"
},

{
    "location": "solver/euler/source.html#EulerEquationMod.applySourceTerm-Tuple{Any,Any,Any,Any,ODLCommonTools.SRCType}",
    "page": "Source Term",
    "title": "EulerEquationMod.applySourceTerm",
    "category": "method",
    "text": "EulerEquationMod.applySourceTerm\n\nThis function updates eqn.res with the source term.  \n\nInputs:      mesh     sbp     eqn     opts     src_func:  the functor that returns the value of the source term at a node                This functor must have the signature:                src_func(q, coords, params, t)                where coords is a vector of length 2 or 3 containing the x and y                 coordinates of the node, params is the ParamType, t is the                 current time, and q is the vector to be populated with the                 source term values.\n\nOutputs: none\n\nAliasing restrictions: params.q_vals cannot be in use\n\n\n\n"
},

{
    "location": "solver/euler/source.html#EulerEquationMod.getSRCFunctors-Tuple{ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,EulerEquationMod.EulerData,Any}",
    "page": "Source Term",
    "title": "EulerEquationMod.getSRCFunctors",
    "category": "method",
    "text": "EulerEquationMod.getSRCFunctors\n\nThis function gets the functor specified by opts[\"SRCname\"] and stores   it to the equation object.  Currently one 1 source functor is allowed.\n\n\n\n"
},

{
    "location": "solver/euler/source.html#EulerEquationMod.SRCDict",
    "page": "Source Term",
    "title": "EulerEquationMod.SRCDict",
    "category": "constant",
    "text": "EulerEquationMod.SRCDict\n\nStores the functors that evaluate source terms at a node.  Every new    functor should be added to this dictonary\n\nAll functors must have the signature:\n\nsrc_func(q, coords, params::ParamType, t)\n\nwhere coords is the vector of length 2 containing the x and y coordinates   of the node, t is the current time, and q is the vector to be populated with   the source term values.\n\n\n\n"
},

{
    "location": "solver/euler/source.html#Functions-1",
    "page": "Source Term",
    "title": "Functions",
    "category": "section",
    "text": "  Modules = [EulerEquationMod]\n  Order = [:function, :constant]\n  Pages = [\"euler/source.jl\"]"
},

{
    "location": "solver/euler/source.html#EulerEquationMod.SRCExp",
    "page": "Source Term",
    "title": "EulerEquationMod.SRCExp",
    "category": "type",
    "text": "AdvectionEquationMod.SRC0\n\nThis is the zero source term.  This is the default of source term   is specified\n\n\n\n"
},

{
    "location": "solver/euler/source.html#EulerEquationMod.SRCPeriodicMMS",
    "page": "Source Term",
    "title": "EulerEquationMod.SRCPeriodicMMS",
    "category": "type",
    "text": "Functor for source term corresponding to ICPeriodicMMS\n\n\n\n"
},

{
    "location": "solver/euler/source.html#Source-Term-Functors-1",
    "page": "Source Term",
    "title": "Source Term Functors",
    "category": "section",
    "text": "  Modules = [EulerEquationMod]\n  Order = [:type]\n  Pages = [\"euler/source.jl\"]"
},

{
    "location": "solver/euler/common.html#",
    "page": "Common Functions",
    "title": "Common Functions",
    "category": "page",
    "text": ""
},

{
    "location": "solver/euler/common.html#EulerEquationMod.calcChannelMMS-Union{Tuple{EulerEquationMod.ParamType{2,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tmsh,1},AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Common Functions",
    "title": "EulerEquationMod.calcChannelMMS",
    "category": "method",
    "text": "Manufactured condition for a channel with y in [0, 1].   No source term.  2D only\n\n\n\n"
},

{
    "location": "solver/euler/common.html#EulerEquationMod.calcExp-Union{Tuple{EulerEquationMod.ParamType{2,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tmsh,1},AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Common Functions",
    "title": "EulerEquationMod.calcExp",
    "category": "method",
    "text": "Calculates a manufactured solution based on exponentials\n\n\n\n"
},

{
    "location": "solver/euler/common.html#EulerEquationMod.calcFreeStream-Union{Tuple{EulerEquationMod.ParamType{2,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tmsh,1},AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Common Functions",
    "title": "EulerEquationMod.calcFreeStream",
    "category": "method",
    "text": "EulerEquationMod.calcFreeStream\n\nThis function calculates the free stream solution for an airfoil problem   based on the angle of attack and Mach number in nondimensionalized variables.\n\nDensity and energy are set to params.rho_free (usually 1.0) and params.E_free,   (usually 1/(gammagamma_1) + 0.5Ma*Ma), and the x and y momenta as\n\nrhoMacos(angle of attack)  and rhoMasin(angle of attack).\n\nThe angle of attack must be in radians.\n\nThis function uses conservative variables regardless of the static parameter   of params.\n\nInputs:     coords: a vector of length 2 containing the x and y coordinates of the point     params: the params object.\n\nInputs/Outputs:     sol: vector of length 4 to be populated with the solution\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/common.html#EulerEquationMod.calcFreeStream0-Union{Tuple{EulerEquationMod.ParamType{2,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tmsh,1},AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Common Functions",
    "title": "EulerEquationMod.calcFreeStream0",
    "category": "method",
    "text": "Like calcFreeStream, but assumes zero angle of attack, regardless   of what the options dictionary says.\n\n\n\n"
},

{
    "location": "solver/euler/common.html#EulerEquationMod.calcFreeStream_dAlpha-Union{Tuple{EulerEquationMod.ParamType{2,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tmsh,1},AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Common Functions",
    "title": "EulerEquationMod.calcFreeStream_dAlpha",
    "category": "method",
    "text": "EulerEquationMod.calcFreeStream_daoa\n\nThis function calculates the free stream solution for an airfoil problem   based on the angle of attack and Mach number in nondimensionalized variables.\n\nDensity and energy are set to params.rho_free (usually 1.0) and params.E_free,   (usually 1/(gammagamma_1) + 0.5Ma*Ma), and the x and y momenta as\n\nrhoMacos(angle of attack)  and rhoMasin(angle of attack).\n\nThe angle of attack must be in radians.\n\nThis function uses conservative variables regardless of the static parameter   of params.\n\nInputs:     coords: a vector of length 2 containing the x and y coordinates of the point     params: the params object.\n\nInputs/Outputs:     sol: vector of length 4 to be populated with the solution\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/common.html#EulerEquationMod.calcInvChannelFreeStream-Union{Tuple{EulerEquationMod.ParamType{2,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tmsh,1},AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Common Functions",
    "title": "EulerEquationMod.calcInvChannelFreeStream",
    "category": "method",
    "text": "Free stream conditions for SU2 inviscid channel\n\n\n\n"
},

{
    "location": "solver/euler/common.html#EulerEquationMod.calcInvChannelIC-Union{Tuple{EulerEquationMod.ParamType{2,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tmsh,1},AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Common Functions",
    "title": "EulerEquationMod.calcInvChannelIC",
    "category": "method",
    "text": "Initial condition for SU2 inviscid channel\n\n\n\n"
},

{
    "location": "solver/euler/common.html#EulerEquationMod.calcIsentropicVortex-Union{Tuple{EulerEquationMod.ParamType{2,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tmsh,N} where N,AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Common Functions",
    "title": "EulerEquationMod.calcIsentropicVortex",
    "category": "method",
    "text": "EulerEquationMod.calcIsentropicVortex\n\nThis function calculates the isentropic vortex solution to the Euler   equations at a node.  This function uses an inner radius of 1, an inner   radius density of 2, an inner radius Mach number of 0.95, and an inner radius   pressure of 1/gamma.  The denstiy as a function of radius r can be found in,   for example,   \"Output Error Estimates for Summation-by-parts Finite-difference Schemes\",   JE Hicken.\n\nThis function uses conservative variables regardless of the static parameter   of params.\n\nInputs:     coords: a vector of length 2 containing the x and y coordinates of the point     params: the params object.\n\nInputs/Outputs:     sol: vector of length 4 to be populated with the solution\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/common.html#EulerEquationMod.calcOnes-Union{Tuple{EulerEquationMod.ParamType{2,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tmsh,1},AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Common Functions",
    "title": "EulerEquationMod.calcOnes",
    "category": "method",
    "text": "EulerEquationMod.calcOnes\n\nThis function sets all the solution variables at a node to 1.0\n\nThis function uses conservative variables regardless of the static parameter   of params.\n\nInputs:     coords: a vector of length 2 containing the x and y coordinates of the point     params: the params object.\n\nInputs/Outputs:     sol: vector of length 4 to be populated with the solution\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/common.html#EulerEquationMod.calcPeriodicMMS-Union{Tuple{EulerEquationMod.ParamType{2,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tmsh,1},AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Common Functions",
    "title": "EulerEquationMod.calcPeriodicMMS",
    "category": "method",
    "text": "Calculates a manufactured solution from Gassner, Winters, Kopriva: Split   Form Nodal DG Schemes with SBP Propertiy for the Compressible Euler Equations.   This is typically used with a mesh that spans [-1, 1] in all directions\n\n\n\n"
},

{
    "location": "solver/euler/common.html#EulerEquationMod.calcRho1Energy2-Union{Tuple{EulerEquationMod.ParamType{2,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tmsh,1},AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Common Functions",
    "title": "EulerEquationMod.calcRho1Energy2",
    "category": "method",
    "text": "EulerEquationMod.calcRho1Energy2\n\nThis function sets the density at a node to 1, energy to 2 and the momenta to   zero.\n\nThis function uses conservative variables regardless of the static parameter   of params.\n\nInputs:     coords: a vector of length 2 containing the x and y coordinates of the point     params: the params object.\n\nInputs/Outputs:     sol: vector of length 4 to be populated with the solution\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/common.html#EulerEquationMod.calcRho1Energy2U1VW0-Union{Tuple{EulerEquationMod.ParamType{2,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tmsh,N} where N,AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Common Functions",
    "title": "EulerEquationMod.calcRho1Energy2U1VW0",
    "category": "method",
    "text": "EulerEquationMod.calcRho1Energy2U1VW0\n\nSets the density values 1.0, x momentum to 1.0,    v & w momenta to 0.0, and energy to 2.0 at a node.\n\nIt should work for 2D and 3D meshes.\n\nThis function uses conservative variables regardless of the static parameter   of params.\n\nInputs:     coords: a vector of length 2 containing the x and y coordinates of the point     params: the params object.\n\nInputs/Outputs:     sol: vector of length 4 to be populated with the solution\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/common.html#EulerEquationMod.calcRho1Energy2U3-Union{Tuple{EulerEquationMod.ParamType{2,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tmsh,N} where N,AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Common Functions",
    "title": "EulerEquationMod.calcRho1Energy2U3",
    "category": "method",
    "text": "EulerEquationMod.calcRho1Energy2U3\n\nSets the density values 1.0, momenta to 0.35355, and   energy to 2.0 at a node.\n\nThis function uses conservative variables regardless of the static parameter   of params.\n\nInputs:     coords: a vector of length 2 containing the x and y coordinates of the point     params: the params object.\n\nInputs/Outputs:     sol: vector of length 4 to be populated with the solution\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/common.html#EulerEquationMod.calcSedovExplosion-Union{Tuple{EulerEquationMod.ParamType{2,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tmsh,1},AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Common Functions",
    "title": "EulerEquationMod.calcSedovExplosion",
    "category": "method",
    "text": "Sedov explosion, 2D only\n\nUniform fluid at rest with energetic region in the circle of radius    0.05 centered at the origin\n\n\n\n"
},

{
    "location": "solver/euler/common.html#EulerEquationMod.calcSquare1D-Union{Tuple{EulerEquationMod.ParamType{2,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tmsh,1},AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Common Functions",
    "title": "EulerEquationMod.calcSquare1D",
    "category": "method",
    "text": "One dimensional square wave, 2D only\n\nThe base flow is uniform, and if  -0.05 < x < 0.05, then it is perturbed   by a small amount.\n\n\n\n"
},

{
    "location": "solver/euler/common.html#EulerEquationMod.calcSquare2D-Union{Tuple{EulerEquationMod.ParamType{2,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tmsh,1},AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Common Functions",
    "title": "EulerEquationMod.calcSquare2D",
    "category": "method",
    "text": "Two dimensional square wave, 2D only\n\nThe base flow is uniform, and if  -0.05 < x, y < 0.05, then it is perturbed   by a small amount.\n\n\n\n"
},

{
    "location": "solver/euler/common.html#EulerEquationMod.calcUnsteadyVortex-Union{Tuple{EulerEquationMod.ParamType{Tdim,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tmsh,1},AbstractArray{Tsol,1}}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tsol}} where Tdim where Tsol where Tmsh",
    "page": "Common Functions",
    "title": "EulerEquationMod.calcUnsteadyVortex",
    "category": "method",
    "text": "EulerEquationMod.calcUnsteadyVortex\n\nThis function calculates the unsteady vortex solution to the Euler equations   at time params.t, where the vortex was centered at x = params.vortex_x0 at   time t=0.  The analytical solution can be found in, for example,   K. Mattsson et al./ Computers & Fluxs 36 (2007) 636-649\n\nThis function uses conservative variables regardless of the static parameter   of params.\n\nInputs:     coords: a vector of length 2 containing the x and y coordinates of the point     params: the params object.\n\nInputs/Outputs:     sol: vector of length 4 to be populated with the solution\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/common.html#EulerEquationMod.calcUnsteadyVortex2-Union{Tuple{EulerEquationMod.ParamType{Tdim,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tmsh,1},AbstractArray{Tsol,1}}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tsol}} where Tdim where Tsol where Tmsh",
    "page": "Common Functions",
    "title": "EulerEquationMod.calcUnsteadyVortex2",
    "category": "method",
    "text": "Unsteady vortex for Carpenters paper \"Entropy Stable Staggered Grid Spectral   Collocation for hte Burgers and Compressible Navier-Stokes Equations\n\nIt is similar to the other unsteady vortex, except the vortex travels   at an angle\n\nDoes not use params.Ma or params.vortex_strength, or vortex.x0\n\n\n\n"
},

{
    "location": "solver/euler/common.html#EulerEquationMod.calcVortex-Union{Tuple{EulerEquationMod.ParamType{2,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tmsh,1},AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Common Functions",
    "title": "EulerEquationMod.calcVortex",
    "category": "method",
    "text": "EulerEquationMod.calcVortex\n\nSets the density 1.0, energy to 2.0 at a node.  The momenta are calculated   according to solid body rotation with an angular velocity of 0.5 centered   at x = 0.\n\nThis function uses conservative variables regardless of the static parameter   of params.\n\nInputs:     coords: a vector of length 2 containing the x and y coordinates of the point     params: the params object.\n\nInputs/Outputs:     sol: vector of length 4 to be populated with the solution\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/common.html#EulerEquationMod.calcZeros-Union{Tuple{EulerEquationMod.ParamType{2,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tmsh,1},AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Common Functions",
    "title": "EulerEquationMod.calcZeros",
    "category": "method",
    "text": "EulerEquationMod.calcZeros\n\nThis function sets all the solution variables at a node to 0.0\n\nThis function uses conservative variables regardless of the static parameter   of params.\n\nInputs:     coords: a vector of length 2 containing the x and y coordinates of the point     params: the params object.\n\nInputs/Outputs:     sol: vector of length 4 to be populated with the solution\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/common.html#sec:euler_common_funcs-1",
    "page": "Common Functions",
    "title": "Common Functions",
    "category": "section",
    "text": "  CurrentModule = EulerEquationModThe functions here evaluate some function at a single point.  This is useful for defining manufactured solutions, initial conditions, and boundary conditions.  In fact, is is generally recommended that every initial condition and Dirichlet boundary condition create a function here to compute the state.  For boundary conditions, all the functor need to is call this function to get the state and call a numerical flux function.  Modules = [EulerEquationMod]\n  Pages = [\"euler/common_funcs.jl\"]"
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
    "text": "  CurrentModule = EulerEquationModIt is possible to write the Euler equations in terms of different sets of variables. The default is to use conservative variables, but there are others as well, such as entropy variables. This page describes how to convert between the different sets of variables. The functions perform conversions back and forth between conservative variables and the \"other\" variables. The code generally does not support conversion directly between two different sets of \"other variables. For example, if we call the conservative variables the q variables, the  first set of \"other\" variables the w1 variables and the second set of \"other\" variables the w2 variables, it is possible to convert between q and w1 as well as q and w2 but not between w1 and w2.The basic set of functions perform the conversion, and these function\'s names and in an underscore. Layered on top of these functions are a safer method of conversion, which uses the static parameters of the ParamType to determine what variables are stored in eqn.q and performs the appropriate conversion. This is important because if, for example, we are solving the equation using entropy variables, and the function convertToEntropy is called, it will (correctly) not do the conversion, because the variables are already the entropy variables. If convertToConservative is called, it will do the proper conversion from entropy to conservative variables.The function describe here provide one additional bit of functionality. Lets say a function needs to convert it current variables to conservative variables, do some operation, then convert the result back to entropy variables. Converting to conservative variables is easy, that\'s what convertToConservative does.  The problem is converting back. The code doesn\'t necessarily know what variables are stored in eqn.q, so it doesn\'t know what function to call.  Should it call convertToW1 or convertToW2? The solution is convertFromNaturalToWorkingVars, which converts from the conservative variables to whatever the variables stored in eqn.q are (as determined by the static parameters of the ParamType object)."
},

{
    "location": "solver/euler/conversion.html#EulerEquationMod.convertFromNaturalToWorkingVars-Union{Tuple{EulerEquationMod.ParamType{Tdim,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1},AbstractArray{Tsol,1}}, Tuple{Tdim}, Tuple{Tsol}} where Tsol where Tdim",
    "page": "Conversion",
    "title": "EulerEquationMod.convertFromNaturalToWorkingVars",
    "category": "method",
    "text": "EulerEquationMod.convertFromNaturlaToWorkingVars\n\nThis function converts a vector qc from the \"natural\" set of variables to    write the equation to the set of variables the equation is being solved in,    defined by the static parameter var_type of the params.  For the Euler    equations, the \"natural\" variables are the conservative variables.   Every new set of variables must extend this function with a new method.\n\nInputs:     params:  ParamType{Tdim, var_type}     qc: vector of \"natural\" variables\n\nInputs/Outputs:     qe: vector to be populated with the new variables\n\nLow level function.\n\nAliasing restrictions: none (this requires that every method of this function                          support in-place conversion).\n\n\n\n"
},

{
    "location": "solver/euler/conversion.html#EulerEquationMod.convertToConservative-Union{Tuple{EulerEquationMod.ParamType{Tdim,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1},AbstractArray{Tsol,1}}, Tuple{Tdim}, Tuple{Tsol}} where Tsol where Tdim",
    "page": "Conversion",
    "title": "EulerEquationMod.convertToConservative",
    "category": "method",
    "text": "low level function\n\nConverts conservative variables to conservative variables (ie. it   copies the input to the output).  This method exists to values can be    converted without knowing whether they are conservative or entropy.\n\nAliasing restrictions: none (qc and qe can be the same vector)\n\n\n\n"
},

{
    "location": "solver/euler/conversion.html#EulerEquationMod.convertToConservative-Union{Tuple{EulerEquationMod.ParamType{Tdim,:entropy,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1},AbstractArray{Tsol,1}}, Tuple{Tdim}, Tuple{Tsol}} where Tsol where Tdim",
    "page": "Conversion",
    "title": "EulerEquationMod.convertToConservative",
    "category": "method",
    "text": "EulerEquationMod.convertToConservative\n\nConverts  the entropy variables at a node to the conservative variables.\n\nInputs:   params  : ParamType{2, :entropy} used to dispatch to the proper method   qe  : vector (of length 4) of entropy variables\n\nInputs/outputs   qc : vector (of length 4) of conservative variables.  Contents of vector are        overwritten\n\nAliasing: none (qc and qe can be the same vector)\n\nThs is a low level function.\n\n\n\n"
},

{
    "location": "solver/euler/conversion.html#EulerEquationMod.convertToConservative-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,:entropy},Any,AbstractArray{Tsol,3}}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tres where Tdim where Tsol where Tmsh",
    "page": "Conversion",
    "title": "EulerEquationMod.convertToConservative",
    "category": "method",
    "text": "mid level function\n\nConverts the array (3D form) of entropy variables to conservative variables    in place.  If the array is already in conservative variables this is a no-op.   Other methods exist if q_arr is a vector.\n\nInputs:     mesh     sbp     eqn     opts\n\nInputs/Outputs:     q_arr: array of values (3D) to be converted\n\nAliasing: no restrictions\n\n\n\n"
},

{
    "location": "solver/euler/conversion.html#EulerEquationMod.convertToConservativeFromIR_-Union{Tuple{EulerEquationMod.ParamType{2,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tsol,1},AbstractArray{Tsol,1}}, Tuple{Tsol}} where Tsol",
    "page": "Conversion",
    "title": "EulerEquationMod.convertToConservativeFromIR_",
    "category": "method",
    "text": "EulerEquationMod.convertToConservativeFromIR_\n\nConverts the IR entropy variables at a node to the conservative variables   regardless of the static parameter var_type.\n\nInputs:   params  : ParamType{Tdim} used to dispatch to the proper method   qe  : vector (of length 4 or 5) of entropy variables\n\nInputs/outputs   qc : vector (of length 4 or 5) of conservative variables.  Contents of         vector are overwritten\n\nAliasing: none (qc and qe can be the same vector)\n\nThs is a low level function.\n\n\n\n"
},

{
    "location": "solver/euler/conversion.html#EulerEquationMod.convertToConservative_-Union{Tuple{EulerEquationMod.ParamType{2,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tsol,1},AbstractArray{Tsol,1}}, Tuple{Tsol}} where Tsol",
    "page": "Conversion",
    "title": "EulerEquationMod.convertToConservative_",
    "category": "method",
    "text": "EulerEquationMod.convertToConservative_\n\nConverts  the entropy variables at a node to the conservative variables   regardless of the static parameter var_type.\n\nInputs:   params  : ParamType{Tdim} used to dispatch to the proper method   qe  : vector (of length 4 or 5) of entropy variables\n\nInputs/outputs   qc : vector (of length 4 or 5) of conservative variables.  Contents of vector are        overwritten\n\nAliasing: none (qc and qe can be the same vector)\n\nThs is a low level function.\n\n\n\n"
},

{
    "location": "solver/euler/conversion.html#EulerEquationMod.convertToEntropy-Union{Tuple{EulerEquationMod.ParamType{Tdim,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1},AbstractArray{Tsol,1}}, Tuple{Tdim}, Tuple{Tsol}} where Tsol where Tdim",
    "page": "Conversion",
    "title": "EulerEquationMod.convertToEntropy",
    "category": "method",
    "text": "EulerEquationMod.convertToEntropy\n\nConverts the conservative variables at a node to entropy variables\n\nInput:   params : ParamType{s, :conservative}   qc  : vector (of length 4 or 5) of conservative variables\n\nOutputs:   qe : vector (of length 4 or 5) of conservative variables.  Contents of         vector areoverwritten\n\nlow level function\n\nAliasing restrictions: none (qc and qe can be the same vector)\n\n\n\n"
},

{
    "location": "solver/euler/conversion.html#EulerEquationMod.convertToEntropy-Union{Tuple{EulerEquationMod.ParamType{Tdim,:entropy,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1},AbstractArray{Tsol,1}}, Tuple{Tdim}, Tuple{Tsol}} where Tsol where Tdim",
    "page": "Conversion",
    "title": "EulerEquationMod.convertToEntropy",
    "category": "method",
    "text": "EulerEquationMod.convertToEntropy\n\nConverts the entropy variables to entropy variables (ie. it copies the    input to the output).  This method exists so variables can be converted    to entropy variables without knowing what type they are.\n\nlow level function\n\nAliasing restrictions: none (qc and qe can be the same vector)\n\n\n\n"
},

{
    "location": "solver/euler/conversion.html#EulerEquationMod.convertToEntropy-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,:conservative},Any,AbstractArray{Tsol,3}}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tres where Tdim where Tsol where Tmsh",
    "page": "Conversion",
    "title": "EulerEquationMod.convertToEntropy",
    "category": "method",
    "text": "mid level function\n\nConverts the array (3D form) of conservative variables to entropy variables    in place.  If the array is already in entropy variables this is a no-op\n\nMethods also exist for the 1D form.\n\nInputs:     mesh     sbp     eqn     opts\n\nInputs/Outputs:     q_arr: array of values (3D) to be converted\n\nAliasing: no restrictions\n\n\n\n"
},

{
    "location": "solver/euler/conversion.html#EulerEquationMod.convertToEntropy_-Union{Tuple{EulerEquationMod.ParamType{2,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tsol,1},AbstractArray{Tsol,1}}, Tuple{Tsol}} where Tsol",
    "page": "Conversion",
    "title": "EulerEquationMod.convertToEntropy_",
    "category": "method",
    "text": "EulerEquationMod.convertToEntropy_\n\nConverts the conservative variables at a node to the entropy variables   regardless of the static parameter var_type.\n\nInputs:   params  : ParamType{Tdim} used to dispatch to the proper method   qe  : vector (of length 4 or 5) of conservative variables\n\nInputs/outputs   qc : vector (of length 4 or 5) of entropy variables.  Contents of vector are        overwritten\n\nAliasing: none (qc and qe can be the same vector)\n\nThs is a low level function.\n\n\n\n"
},

{
    "location": "solver/euler/conversion.html#EulerEquationMod.convertToIR_-Union{Tuple{EulerEquationMod.ParamType{2,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tsol,1},AbstractArray{Tsol,1}}, Tuple{Tsol}} where Tsol",
    "page": "Conversion",
    "title": "EulerEquationMod.convertToIR_",
    "category": "method",
    "text": "EulerEquationMod.convertToIR_\n\nConverts the conservative variables at a node to the entropy variables   of Ismail and Roe regardless of the static parameter var_type.\n\nInputs:   params  : ParamType{Tdim} used to dispatch to the proper method   qe  : vector (of length 4 or 5) of conservative variables\n\nInputs/outputs   qc : vector (of length 4 or 5) of entropy variables.  Contents of vector are        overwritten\n\nAliasing: none (qc and qe can be the same vector)\n\nThs is a low level function.\n\n\n\n"
},

{
    "location": "solver/euler/conversion.html#Functions-1",
    "page": "Conversion",
    "title": "Functions",
    "category": "section",
    "text": "  Modules = [EulerEquationMod]\n  Pages = [\"euler/conversion.jl\"]"
},

{
    "location": "solver/euler/flux_functions.html#",
    "page": "Numerical Flux Functions",
    "title": "Numerical Flux Functions",
    "category": "page",
    "text": ""
},

{
    "location": "solver/euler/flux_functions.html#EulerEquationMod.RoeSolver-Union{Tuple{EulerEquationMod.ParamType{2,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tsol,1},AbstractArray{Tsol,1},AbstractArray{Tres,1},AbstractArray{Tmsh,1},AbstractArray{Tres,1}}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tres where Tsol where Tmsh",
    "page": "Numerical Flux Functions",
    "title": "EulerEquationMod.RoeSolver",
    "category": "method",
    "text": "EulerEquationMod.RoeSolver\n\nThis calculates the Roe flux for boundary conditions at a node. The inputs   must be in conservative variables.\n\nInputs\n\nparams : ParamType\nq  : conservative variables of the fluid\nqg : conservative variables of the boundary\naux_vars : vector of all auxiliary variables at this node\nnrm : scaled face normal vector (x-y space, outward normal of the face q    lives on)\n\nOutputs\n\nflux : vector to populate with solution\n\nAliasing restrictions:  none of the inputs can alias params.res_vals1,                           params.res_vals2, params.q_vals, params.flux_vals1, or                           params.sat\n\n\n\n"
},

{
    "location": "solver/euler/flux_functions.html#EulerEquationMod.RoeSolver-Union{Tuple{EulerEquationMod.ParamType{3,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tsol,1},AbstractArray{Tsol,1},AbstractArray{Tres,1},AbstractArray{Tmsh,1},AbstractArray{Tres,1}}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tres where Tsol where Tmsh",
    "page": "Numerical Flux Functions",
    "title": "EulerEquationMod.RoeSolver",
    "category": "method",
    "text": "The main Roe solver.  Populates flux with the computed flux.\n\n\n\n"
},

{
    "location": "solver/euler/flux_functions.html#EulerEquationMod.RoeSolver_revm-Union{Tuple{EulerEquationMod.ParamType{2,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tsol,1},AbstractArray{Tsol,1},AbstractArray{Tres,1},AbstractArray{Tmsh,1},AbstractArray{Tres,1},AbstractArray{Tmsh,1}}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tres where Tsol where Tmsh",
    "page": "Numerical Flux Functions",
    "title": "EulerEquationMod.RoeSolver_revm",
    "category": "method",
    "text": "###EulerEquationMod.RoeSolver_revm\n\nReverse mode of EulerEquationMod.RoeSolver. This function computes the reverse mode of the Roe flux w.r.t the mesh metrics\n\nInputs\n\nparams : Parameter object\nq  : Conservative variable of the fluid\nqg : Conservative variable of the boundary or the adjacent element\naux_vars : Auxiliary variables\nnrm : Element face normal vector in the physical space\nflux_bar : Flux value which needs to get differentiated\n\nOutput\n\nnrm_bar : derivaitve of the flux_bar w.r.t the mesh metrics\n\nAliasing Restrictions: Same as the forward function\n\n\n\n"
},

{
    "location": "solver/euler/flux_functions.html#EulerEquationMod.calcEulerFlux_Ducros-Union{Tuple{EulerEquationMod.ParamType{2,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1},AbstractArray{Tsol,1},AbstractArray{Tres,N} where N,AbstractArray{Tmsh,N} where N,AbstractArray{Tres,1}}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tres where Tsol where Tmsh",
    "page": "Numerical Flux Functions",
    "title": "EulerEquationMod.calcEulerFlux_Ducros",
    "category": "method",
    "text": "Calculates the numerical flux function associated with the Ducros flux   splitting.  Methods are available for 2D and 3D.\n\nInputs:\n\nparams:\nqL: the left state\nqR: the right state\naux_vars: the aux vars for the left state\ndir: the direction vector\n\nInputs/Outputs:\n\nF: vector to be populated with the flux\n\n\n\n"
},

{
    "location": "solver/euler/flux_functions.html#EulerEquationMod.calcEulerFlux_IR-Union{Tuple{EulerEquationMod.ParamType{2,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1},AbstractArray{Tsol,1},AbstractArray{Tres,N} where N,AbstractArray{Tmsh,1},AbstractArray{Tres,1}}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tres where Tsol where Tmsh",
    "page": "Numerical Flux Functions",
    "title": "EulerEquationMod.calcEulerFlux_IR",
    "category": "method",
    "text": "This function calculates the Ismail-Roe numerical flux at a node in a   specified direction\n\nInputs:\n\nparams: ParamType\nqL: left state vector\nqR: right state vector\naux_vars: auxiliary variable vector for qL\ndir: a direction vector of length Tdim\n\nInputs/Outputs:\n\nF: a numDofPerNode x Tdim matrix where each column will be populated with   the flux in the direction specified by the corresponding column of nrm\n\nAliasing Restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/flux_functions.html#EulerEquationMod.calcEulerFlux_IR-Union{Tuple{EulerEquationMod.ParamType{2,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1},AbstractArray{Tsol,1},AbstractArray{Tres,N} where N,AbstractArray{Tmsh,2},AbstractArray{Tres,2}}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tres where Tsol where Tmsh",
    "page": "Numerical Flux Functions",
    "title": "EulerEquationMod.calcEulerFlux_IR",
    "category": "method",
    "text": "This function computes the Ismail-Roe   flux in Tdim directions at once for a given state.  This is more efficient   Than calling the single direction method Tdim times.  Methods are available   for 2 and 3 dimensions.\n\nInputs:\n\nparams: ParamType\nqL: left state vector\nqR: right state vector\naux_vars: auxiliary variable vector for qL\ndir: a Tdim x Tdim matrix with each column containing a normal vector\n\nInputs/Outputs:\n\nF: a numDofPerNode x Tdim matrix where each column will be populated with   the flux in the direction specified by the corresponding column of nrm\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/flux_functions.html#EulerEquationMod.calcEulerFlux_IRSLF-Union{Tuple{EulerEquationMod.ParamType{Tdim,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1},AbstractArray{Tsol,1},AbstractArray{Tres,1},AbstractArray{Tmsh,1},AbstractArray{Tres,1}}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tdim where Tres where Tsol where Tmsh",
    "page": "Numerical Flux Functions",
    "title": "EulerEquationMod.calcEulerFlux_IRSLF",
    "category": "method",
    "text": "This is the second method that takes in a normal vector directly.   See the first method for a description of what this function does.\n\nInputs     qL, qR: vectors conservative variables at left and right states     aux_vars: aux_vars for qL     nrm: a normal vector in xy space\n\nInputs/Outputs     F: vector to be updated with the result\n\nAlising restrictions:     See getEntropyLFStab\n\n\n\n"
},

{
    "location": "solver/euler/flux_functions.html#EulerEquationMod.calcEulerFlux_IRSLW-Union{Tuple{EulerEquationMod.ParamType,AbstractArray{Tsol,1},AbstractArray{Tsol,1},AbstractArray{Tres,N} where N,AbstractArray{Tmsh,2},AbstractArray{Tmsh,N} where N,AbstractArray{Tres,1}}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tres where Tsol where Tmsh",
    "page": "Numerical Flux Functions",
    "title": "EulerEquationMod.calcEulerFlux_IRSLW",
    "category": "method",
    "text": "This function is similar to calcEulerFlux_IRSLF, but uses Lax-Wendroff   dissipation rather than Lax-Friedrich.\n\nAliasing restrictions: see getEntropyLWStab\n\n\n\n"
},

{
    "location": "solver/euler/flux_functions.html#EulerEquationMod.calcEulerFlux_standard-Union{Tuple{EulerEquationMod.ParamType{2,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1},AbstractArray{Tsol,1},AbstractArray{Tsol,1},AbstractArray{Tmsh,2},AbstractArray{Tres,2}}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tres where Tsol where Tmsh",
    "page": "Numerical Flux Functions",
    "title": "EulerEquationMod.calcEulerFlux_standard",
    "category": "method",
    "text": "This function computes the split form flux corresponding to the standard   flux in Tdim directions at once for a given state.  This is more efficient   Than calling the single direction method Tdim times.  Methods are available   for 2 and 3 dimensions.\n\nInputs:     params: ParamType     qL: left state vector     qR: right state vector     aux_vars: auxiliary variable vector for qL     dir: a Tdim x Tdim matrix with each column containing a normal vector\n\nInputs/Outputs:     F: a numDofPerNode x Tdim matrix where each column will be populated with        the flux in the direction specified by the corresponding column of nrm\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/flux_functions.html#EulerEquationMod.calcLFFlux-Union{Tuple{EulerEquationMod.ParamType{Tdim,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1},AbstractArray{Tsol,1},AbstractArray{Tsol,1},AbstractArray{Tmsh,1},AbstractArray{Tres,1}}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tdim where Tres where Tsol where Tmsh",
    "page": "Numerical Flux Functions",
    "title": "EulerEquationMod.calcLFFlux",
    "category": "method",
    "text": "Calculates the Lax-Friedrich flux function on the conservative variables\n\n\n\n"
},

{
    "location": "solver/euler/flux_functions.html#EulerEquationMod.calcSAT-Union{Tuple{EulerEquationMod.ParamType{2,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tsol,1},AbstractArray{Tsol,1},AbstractArray{Tmsh,1},AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Numerical Flux Functions",
    "title": "EulerEquationMod.calcSAT",
    "category": "method",
    "text": "###EulerEquationMod.calcSAT\n\nComputes the simultaneous approximation term for use in computing the numerical flux, namely:\n\n0.5*(|A_hat| - A)(qL - qR)\n\nArguments\n\nparams : Parameter object of type ParamType\n`roe_vars: [u, v, H] at roe average state\ndq  : difference between left and right states\nnrm : Normal to face in the physical space\nsat : Simultaneous approximation Term\n\n\n\n"
},

{
    "location": "solver/euler/flux_functions.html#EulerEquationMod.calcSAT_revm-Union{Tuple{EulerEquationMod.ParamType{2,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tmsh,1},AbstractArray{Tsol,1},AbstractArray{Tsol,1},Tsol,AbstractArray{Tsol,1},AbstractArray{Tmsh,1}}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Numerical Flux Functions",
    "title": "EulerEquationMod.calcSAT_revm",
    "category": "method",
    "text": "###EulerEquationMod.calcSAT_revm\n\nReverse mode of calcSAT\n\nInputs\n\nparams : Parameter object of type ParamType\nnrm : Normal to face in the physical space\ndq  : Boundary condition penalty variable\nvel : Velocities along X & Y directions in the physical space\nH   : Total enthalpy\nsat_bar : Inpute seed for sat flux whose derivative needs to be computed\n\nOutput\n\nnrm_bar : derivative of sat_bar w.r.t physical normal vector\n\n\n\n"
},

{
    "location": "solver/euler/flux_functions.html#sec:euler_flux_funcs-1",
    "page": "Numerical Flux Functions",
    "title": "Numerical Flux Functions",
    "category": "section",
    "text": "This page documents the numerical flux functions available in the codebc_solvers.jl should be renamed to this  CurrentModule = EulerEquationModModules = [EulerEquationMod]\nPages = [\"euler/bc_solvers.jl\"]"
},

{
    "location": "solver/euler/flux_functions_diff.html#",
    "page": "Numerical Flux Functions Jacobian",
    "title": "Numerical Flux Functions Jacobian",
    "category": "page",
    "text": ""
},

{
    "location": "solver/euler/flux_functions_diff.html#EulerEquationMod.RoeSolver_diff-Union{Tuple{EulerEquationMod.ParamType{2,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1},AbstractArray{Tsol,1},AbstractArray{Tres,1},AbstractArray{Tmsh,1},AbstractArray{Tres,2},AbstractArray{Tres,2}}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tres where Tsol where Tmsh",
    "page": "Numerical Flux Functions Jacobian",
    "title": "EulerEquationMod.RoeSolver_diff",
    "category": "method",
    "text": "Differentiated version of the RoeSolver.  Computes the jacobian   of the flux with respect to q and qg.  Methods are available for 2D   and 3D\n\nInputs\n\nparams: ParamType, conservative variables only\nq: vector of conservative variables for the left state\nqg: vector of conservative variables for the right state\naux_vars: auxiliary variables for the left state\nnrm: scaled normal vector in x-y space (outward wrt the element q lives on\n\nInputs/Outputs\n\nfluxL_dot: flux jacobian wrt q, numDofPerNode x numDofPerNode\nfluxR_dot: flux jacobian wrt qg, numDofPerNode x numDofPerNode\n\nAliasing restrictions:\n\nparams: sat, roe_vars, roe_vars_dot, euler_fluxjac, p_dot must not be used\n\n\n\n"
},

{
    "location": "solver/euler/flux_functions_diff.html#sec:euler_flux_funcs_diff-1",
    "page": "Numerical Flux Functions Jacobian",
    "title": "Numerical Flux Functions Differentiated",
    "category": "section",
    "text": "This page documents the differentiated numerical flux functions available in the codebc_solvers_diff.jl should be renamed to this  CurrentModule = EulerEquationModModules = [EulerEquationMod]\nPages = [\"euler/bc_solvers_diff.jl\"]"
},

{
    "location": "solver/euler/stabilization.html#",
    "page": "Stabilization",
    "title": "Stabilization",
    "category": "page",
    "text": ""
},

{
    "location": "solver/euler/stabilization.html#Stabilization-1",
    "page": "Stabilization",
    "title": "Stabilization",
    "category": "section",
    "text": "  CurrentModule = EulerEquationModThis page describes the stabilization methods used for continuous Galerkin formulations. Not all of these effectively stabilized the discretization, and since the move to discontinuous Galerkin, they may have bitrotted."
},

{
    "location": "solver/euler/stabilization.html#EulerEquationMod.edgestabilize!-Union{Tuple{SummationByParts.AbstractSBP{T},Array{ODLCommonTools.Interface,N} where N,ODLCommonTools.AbstractMesh,AbstractArray{T,2},AbstractArray{T,3},AbstractArray{T,4},AbstractArray{T,2},AbstractArray{T,4},AbstractArray{Tres,2},AbstractArray{T,2}}, Tuple{Tres}, Tuple{T}} where Tres where T",
    "page": "Stabilization",
    "title": "EulerEquationMod.edgestabilize!",
    "category": "method",
    "text": "EulerEquationMod.edgestabilize!\n\nAdds edge-stabilization (see Burman and Hansbo, doi:10.1016/j.cma.2003.12.032) to a given residual.  Different methods are available depending on the rank of u:\n\nFor scalar fields, it is assumed that u is a rank-2 array, with the first\n\ndimension for the local-node index, and the second dimension for the element index.\n\nFor vector fields, u is a rank-3 array, with the first dimension for the\n\nindex of the vector field, the second dimension for the local-node index, and the third dimension for the element index.\n\nNaturally, the number of entries in the dimension of u (and res) corresponding to the nodes must be equal to the number of nodes in the SBP operator sbp.\n\nInputs\n\nsbp: an SBP operator type\nifaces: list of element interfaces stored as an array of Interfaces\nu: the array of solution data, numDofPerNode x numNodesPerElement x numEl.      for scalar fields (ie. numDofPerNode), this can be a 2D array.\nx: Cartesian coordinates stored in (coord,node,element) format\ndxidx: scaled Jacobian of the mapping (as output from mappingjacobian!)\njac: determinant of the Jacobian array, numNodesPerElement x numEl\nalpha: array of transformation terms (see below)\nstabscale: numfaces x numNodes per face array of scaling parameter\n\nIn/Outs\n\nres: where the result of the integration is stored, same shape as u\n\nDetails\n\nThe array alpha is used to compute the directional derivative normal to the faces. For a 2-dimensional problem, it can be computed as follows:\n\n  for k = 1:mesh.numelem\n    for i = 1:sbp.numnodes\n      for di1 = 1:2\n        for di2 = 1:2\n          alpha[di1,di2,i,k] = (dxidx[di1,1,i,k].*dxidx[di2,1,i,k] + \n                            dxidx[di1,2,i,k].*dxidx[di2,2,i,k])*jac[i,k]\n        end\n      end\n    end\n  end\n\n\n\n"
},

{
    "location": "solver/euler/stabilization.html#EulerEquationMod.stabscale-Union{Tuple{AbstractArray{Tsol,1},AbstractArray{Tmsh,2},AbstractArray{Tmsh,1},EulerEquationMod.ParamType{2,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Stabilization",
    "title": "EulerEquationMod.stabscale",
    "category": "method",
    "text": "EulerEquationMod.stabscale\n\nThis function calculates the edge stabilization scalaing parameter at a node   and returns it.\n\nInputs:     * q : vector of conservative variables     * dxidx : jacobian of xi wrt x coordinates at the node     * nrm : normal vector in xi space     * params : ParamType{2}\n\nThis is a low level function\n\n\n\n"
},

{
    "location": "solver/euler/stabilization.html#EulerEquationMod.stabscale-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim where Tres}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Stabilization",
    "title": "EulerEquationMod.stabscale",
    "category": "method",
    "text": "EulerEquationMod.stabscale\n\nThis function calculate the stabilization scaling parameter across the   entire mesh by calling the low level method.  This populates eqn.stabscale.\n\nThis is a mid level function\n\n\n\n"
},

{
    "location": "solver/euler/stabilization.html#EulerEquationMod.GLS-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tdim where Tres where Tsol where Tmsh",
    "page": "Stabilization",
    "title": "EulerEquationMod.GLS",
    "category": "method",
    "text": "EulerEquationMod.GLS\n\nAdd Galerkin Least-Squares stabilization to the weak residual. This implementation is only for steady problems and conservative variables\n\nInputs\n\nmesh: AbstractMesh type\nsbp : Summation-by-parts operator\neqn : Equation object used elsewhere\n\nOutputs\n\nNone\n\n\n\n"
},

{
    "location": "solver/euler/stabilization.html#EulerEquationMod.SUPG-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tdim where Tres where Tsol where Tmsh",
    "page": "Stabilization",
    "title": "EulerEquationMod.SUPG",
    "category": "method",
    "text": "EulerEquationMod.SUPG\n\nAdd Streamwise Upwind Petrov-Galerkin stabilization to the weak residual. This implementation is only for steady problems\n\nInputs\n\nmesh: AbstractMesh type\nsbp : Summation-by-parts operator\neqn : Equation object used elsewhere\n\nOutputs\n\nNone\n\n\n\n"
},

{
    "location": "solver/euler/stabilization.html#EulerEquationMod.calcAxidxi-Union{Tuple{AbstractArray{Tsol,2},AbstractArray{Tsol,3},AbstractArray{Tsol,3},AbstractArray{Tsol,3},Any,Any}, Tuple{Tsol}} where Tsol",
    "page": "Stabilization",
    "title": "EulerEquationMod.calcAxidxi",
    "category": "method",
    "text": "EulerEquationMod.calcAxidxi\n\nCalculates AxiDxi + AetaDeta at the element level\n\nInputs\n\nAxidxi : Product of flux jacobian with shape function derivative             AxiDxi + AetaDeta\nshapefuncderiv: Shape function derivative (Dxi and Deta above)\nAxi : Flux jacobian in the xi direction\nAeta : Flux jacobian in the eta direction\nndof : Number of degrees of freedom per node\nnnpe : Number of nodes per element\n\nOutputs\n\nNone\n\n\n\n"
},

{
    "location": "solver/euler/stabilization.html#EulerEquationMod.calcEigenFactorization-Union{Tuple{EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type,AbstractArray{Tsol,1},AbstractArray{Tsol,1},AbstractArray{Tsol,2},AbstractArray{Tsol,2},AbstractArray{Tsol,2}}, Tuple{Tdim}, Tuple{Tres}, Tuple{Tsol}} where Tdim where Tres where Tsol",
    "page": "Stabilization",
    "title": "EulerEquationMod.calcEigenFactorization",
    "category": "method",
    "text": "EulerEquationMod.calcEigenFactorization\n\nComputes the eigen value facotrization of the flux jacobian in either xi or eta direction. This is a nodal level function. It is only for 2D conservative  variables.\n\nReference: Efficient Solution Methods for the NavierStokes Equations, T.H,            Pulliam, NASA Ames Research Center\n\nInputs\n\neqn  : Euler equation object\nq    : Conservative variables at a node\ndxidx: mapping jacobian vector. It can be either [dξ/dx dξ/dy] or            [dη/dx dη/dy] depending on what Axi or Aeta being factorized.\nT    : Matrix of eigen vectors.\nTinv : Inverse of T\nLambda : Diagonal matrix of eigen values\n\nOutputs\n\nNone\n\n\n\n"
},

{
    "location": "solver/euler/stabilization.html#EulerEquationMod.calcElementArea-Union{Tuple{AbstractArray{Tmsh,2}}, Tuple{Tmsh}} where Tmsh",
    "page": "Stabilization",
    "title": "EulerEquationMod.calcElementArea",
    "category": "method",
    "text": "EulerEquationMod.calcElementArea\n\nCalculates the element area of a 2D triangular element\n\nInputs\n\ncoords : Vertex coordinates of the element\n\nOutputs\n\narea : Area of the triangular element\n\n\n\n"
},

{
    "location": "solver/euler/stabilization.html#EulerEquationMod.calcFluxJacobian-Union{Tuple{EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type,AbstractArray{Tsol,1},AbstractArray{Tsol,2},AbstractArray{Tsol,2}}, Tuple{Tdim}, Tuple{Tres}, Tuple{Tsol}} where Tdim where Tres where Tsol",
    "page": "Stabilization",
    "title": "EulerEquationMod.calcFluxJacobian",
    "category": "method",
    "text": "EulerEquationMod.calcFluxJacobian\n\nIt computes the Euler flux Jacobian at the nodal level in te physical space  (X & Y axes). It uses conservative variables.\n\nInputs\n\neqn  : Equation object\nq    : Conservative variables at the node\nAx   : Flux Jacobian at a node in the X-direction\nAy   : Flux Jacobian at a node in the Y-direction\n\nOutputs\n\nNone \n\n\n\n"
},

{
    "location": "solver/euler/stabilization.html#EulerEquationMod.calcStabilizationTerm-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type,AbstractArray{Tsol,2}}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tdim where Tres where Tsol where Tmsh",
    "page": "Stabilization",
    "title": "EulerEquationMod.calcStabilizationTerm",
    "category": "method",
    "text": "EulerEquationMod.calcStabilizationTerm\n\nCalculates the stabilization term tau for all the nodes in the mesh\n\nInputs\n\nmesh : Abstract mesh type\nsbp  : Summation-by-parts operator\neqn  : Equation object within EulerEquationMod\ntau  : Stabilization term being calculated\n\nOutputs\n\nNone\n\n\n\n"
},

{
    "location": "solver/euler/stabilization.html#EulerEquationMod.calcTau-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type,AbstractArray{Tsol,4},Integer}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tdim where Tres where Tsol where Tmsh",
    "page": "Stabilization",
    "title": "EulerEquationMod.calcTau",
    "category": "method",
    "text": "EulerEquationMod.calcTau\n\nCalculates the stabilization term tau for GLS\n\nInputs\n\nmesh : Abstract mesh type\nsbp  : Summation-by-parts operator\neqn  : Euler equation object\ntau  : Stabilization term\n\nOutputs\n\nNone\n\n\n\n"
},

{
    "location": "solver/euler/stabilization.html#EulerEquationMod.circumcircleDiameter-Union{Tuple{AbstractArray{Tmsh,2}}, Tuple{Tmsh}} where Tmsh",
    "page": "Stabilization",
    "title": "EulerEquationMod.circumcircleDiameter",
    "category": "method",
    "text": "EulerEquationMod.circumcircleDiameter\n\nCalculates the circumcircle diameter of a triangular element\n\nInputs\n\ncoords : Vertex coordinates of the element\n\nOutputs\n\ndia : Diameter of the circumcircle\n\n\n\n"
},

{
    "location": "solver/euler/stabilization.html#EulerEquationMod.parametricFluxJacobian-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tdim where Tres where Tsol where Tmsh",
    "page": "Stabilization",
    "title": "EulerEquationMod.parametricFluxJacobian",
    "category": "method",
    "text": "EulerEquationMod.parametricFluxJacobian\n\nClaculates the flux jacobian in the parametric space for all the nodes in the mesh. It uses conservative variables\n\nInputs\n\nmesh : Abstract mesh object\nsbp  : Summation-by-parts operator\neqn  : Equation object\n\nOutputs\n\nNone\n\n\n\n"
},

{
    "location": "solver/euler/stabilization.html#EulerEquationMod.applyFilter-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP,ODLCommonTools.AbstractSolutionData{Tsol,Tres} where Tres,Any,AbstractArray{T,3} where T}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Stabilization",
    "title": "EulerEquationMod.applyFilter",
    "category": "method",
    "text": "EulerEquationMod.applyFilter\n\nThis function multiplies a filter matrix by the variables at every   node in the mesh.  The filter matrix is stored in eqn.params.filter_mat\n\nInputs:     mesh     sbp     eqn     opts\n\nKeyword Args:     trans=false  : transpose the filter matrix or not.\n\nInputs/Outputs:     arr: 3D array to multiply the filter matrix by\n\n\n\n"
},

{
    "location": "solver/euler/stabilization.html#EulerEquationMod.calcFilter-Tuple{SummationByParts.AbstractSBP,String,Any}",
    "page": "Stabilization",
    "title": "EulerEquationMod.calcFilter",
    "category": "method",
    "text": "EulerEquationMod.calcfilter\n\nThis function calculates the filter used by applyFilter().   The filter is calculated as Vfilter_kernelinv(V), where V is a matrix   that converts the interpolating SBP basis to a modal basis and filter_kernal   is a matrix (typically diagonal) that performs the filtering of the modal   basis coefficients.\n\nInputs:     sbp: SBP operator used to compute the solutions     filter_name: and String that matches the name of a filter function.                  This string is used to retrieve the function from a dictionary                  of all supported filter function.     opts:  options dictonary\n\nOutputs:     F_ret: a numNodesPerElement x numNodesPerElement filter operator\n\n\n\n"
},

{
    "location": "solver/euler/stabilization.html#EulerEquationMod.calcLowPassFilter-Tuple{SummationByParts.AbstractSBP,Any}",
    "page": "Stabilization",
    "title": "EulerEquationMod.calcLowPassFilter",
    "category": "method",
    "text": "EulerEquationMod.calcLowPassFilter\n\nThis function calculates a low pass filter diagonal filter matrix.\n\nInputs:     sbp: SBP operator     opts: options dictonary\n\nOutputs:     filt: numNodesPerElement x numNodesPerElement diagonal filter matrix\n\n\n\n"
},

{
    "location": "solver/euler/stabilization.html#EulerEquationMod.calcModalTransformationOp-Tuple{SummationByParts.AbstractSBP}",
    "page": "Stabilization",
    "title": "EulerEquationMod.calcModalTransformationOp",
    "category": "method",
    "text": "EulerEquationMod.calcModalTransformationOp\n\nThis function calculates a matrix operator V that transforms from a modal   basis to an interpolating basis for SBP operators.  Because the transformation   itself is rank deficient, V is augmented with the last n vector of its Q   matrix (from the full QR decomposition), where numNodesPerElement - rank(V)   = n\n\nInputs:     sbp:  SBP operator\n\nOutputs:     V_full: a numNodesPerElement x numNodesPerElement generalized             Vandermonde matrix\n\n\n\n"
},

{
    "location": "solver/euler/stabilization.html#EulerEquationMod.calcRaisedCosineFilter-Tuple{SummationByParts.AbstractSBP,Any}",
    "page": "Stabilization",
    "title": "EulerEquationMod.calcRaisedCosineFilter",
    "category": "method",
    "text": "EulerEquationMod.calcRaisedCosineFilter\n\nThis function constructs a diagonal matrix filt for a 2 dimensional basis   filter, using the raised cosine filter described in Spectral Methods for   the Euler Equations: Part I by Hussaini, Kproiva, Salas, Zang\n\nInputs     sbp: SBP operator     opts: options dictonary\n\nOutputs:     filt: an numNodesPerElement x numNodesPerElement diagonal filter matrix\n\n\n\n"
},

{
    "location": "solver/euler/stabilization.html#EulerEquationMod.getPascalLevel-Tuple{Integer}",
    "page": "Stabilization",
    "title": "EulerEquationMod.getPascalLevel",
    "category": "method",
    "text": "EulerEquationMod.getPascalLevel\n\nThis function returns the level of Pascals Triangle a particular node   lies in.\n\nInputs:     node:  node index\n\nOutputs:     level: integer describing level of Pascals triangle.\n\n\n\n"
},

{
    "location": "solver/euler/stabilization.html#EulerEquationMod.applyDissipation-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP,ODLCommonTools.AbstractSolutionData{Tsol,Tres} where Tres,Any,AbstractArray{T,3}}, Tuple{Tmsh}, Tuple{Tsol}, Tuple{T}} where T where Tsol where Tmsh",
    "page": "Stabilization",
    "title": "EulerEquationMod.applyDissipation",
    "category": "method",
    "text": "EulerEquationMod.applyDissipation\n\nThis function multiplies the dissipation matrix operator for each element   by the values in a 3D array (typically eqn.q) for each node of the element    and stores the results in eqn.res.\n\nThe dissipation matrix operators must be stored in   eqn.dissipation_mat, a numNodesPerElement x numNodesPerElement x numEl array.\n\nInputs:     mesh     sbp     eqn     opts     arr: a 3D array to apply the dissipation matrix to\n\nAliasing restrictions: no guarantees what happens if arr is eqn.res\n\n\n\n"
},

{
    "location": "solver/euler/stabilization.html#EulerEquationMod.calcDissipationOperator-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.AbstractEulerData{Tsol,Tres} where Tres,Any,String}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Stabilization",
    "title": "EulerEquationMod.calcDissipationOperator",
    "category": "method",
    "text": "EulerEquationMod.calcDissipationOperator\n\nThis function calculates the numNodesPerElement x numNodesPerElement    dissipation matrix operator for each element and returns them in an array   numNodesPerElement x numNodesPerElement x numEl array.\n\nThe dissipation matrix operator is calculated as epsilonfilt.\'h_jac*filt,   where filt is a (usually diagonal) filter matrix, h_jac is a scaling term    that incorporates information about the shape of the element.\n\nInputs:     mesh     sbp     eqn     opts     dissipation_name: an String of the function name used to retrieve                       the function that generates the matrix filt from a                       dictonary\n\nOutputs:     dissipation_mat: a numNodesPerElement x numNodesPerElement x numEl array\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/stabilization.html#EulerEquationMod.getDissipationFilterOperator-Union{Tuple{SummationByParts.TriSBP{T},Function}, Tuple{T}} where T",
    "page": "Stabilization",
    "title": "EulerEquationMod.getDissipationFilterOperator",
    "category": "method",
    "text": "EulerEquatoinMod.getDissipatoinFilterOperator\n\nThis function gets the dissipation filter operator used to construction   the artificial dissipation matrix.\n\nThe filter is calculated as Vfilter_kernelinv(V), where V is a matrix   that converts the interpolating SBP basis to a modal basis and filter_kernal   is a matrix (typically diagonal) that performs the filtering of the modal   basis coefficients.\n\nInputs:     sbp: SBP operator     filter: function to call to calculate the filter kernal\n\nOutputs:     F: a numNodesPerElement x numNodesPerElement filter matrix\n\n\n\n"
},

{
    "location": "solver/euler/stabilization.html#Functions-1",
    "page": "Stabilization",
    "title": "Functions",
    "category": "section",
    "text": "  Modules = [EulerEquationMod]\n  Pages = [\"euler/stabilization.jl\", \"euler/GLS.jl\", \"euler/filtering.jl\", \"euler/artificial_dissipation.jl\", ]"
},

{
    "location": "solver/euler/adjoint.html#",
    "page": "Adjoint",
    "title": "Adjoint",
    "category": "page",
    "text": ""
},

{
    "location": "solver/euler/adjoint.html#Euler-Equation-Steady-Adjoint-1",
    "page": "Adjoint",
    "title": "Euler Equation Steady Adjoint",
    "category": "section",
    "text": "PDESolver currently has the capability to compute the steady adjoint of a boundary functional. Recall the adjoint equation asfracpartial mathcalLpartial q = fracpartial mathcalJpartial q + psi^T fracpartial mathcalRpartial q = 0where, mathcalL is the Lagrangian for functional mathcalJ and q is the solution variable. The adjoint can be computed by calling the function calcAdjoint, which has been described below.  Modules = [EulerEquationMod]\n  Pages = [\"solver/euler/adjoint.jl\"]"
},

{
    "location": "solver/euler/boundary_functional.html#",
    "page": "Boundary Functional",
    "title": "Boundary Functional",
    "category": "page",
    "text": ""
},

{
    "location": "solver/euler/boundary_functional.html#EulerEquationMod.eval_dJdaoa-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim where Tres,Any,EulerEquationMod.BoundaryForceData,String,AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Boundary Functional",
    "title": "EulerEquationMod.eval_dJdaoa",
    "category": "method",
    "text": "EulerEquationMod.eval_dJdaoa\n\nCompute the complete derivative of a functional w.r.t angle of attack\n\nInputs\n\nmesh : Abstract mesh object\nsbp  : Summation-By-Parts operator\neqn  : Euler equation object\nopts : Options dictionary\nfunctionalData : Object of type AbstractFunctional. This is type is associated                    with the functional being computed and holds all the                    relevant data.\nfunctionalName : Name of the functional being evaluated\nadjoint_vec : Local portion of the adjoint vector owned by an MPI rank\n\nOutput\n\ndJdaoa : Complete derivative of the functional w.r.t angle of attack            This is a scalar value that is the same across all MPI ranks\n\n\n\n"
},

{
    "location": "solver/euler/boundary_functional.html#EulerEquationMod.calcBndryFunctional-Union{Tuple{ODLCommonTools.AbstractDGMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type,Any,EulerEquationMod.BoundaryForceData}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tdim where Tres where Tsol where Tmsh",
    "page": "Boundary Functional",
    "title": "EulerEquationMod.calcBndryFunctional",
    "category": "method",
    "text": "EulerEquationMod.calcBndryFunctional\n\nThis function calculates a functional on a geometric boundary of a the computational space. This is a mid level function that should not be called from outside the module. Depending on the functional being computed, it may be necessary to define another method for this function based on a different boundary functional type or parameters.\n\nInputs\n\nmesh :  Abstract mesh object\nsbp  : Summation-By-Parts operator\neqn  : Euler equation object\nopts : Options dictionary\nfunctionalData : Object which is a subtype of Abstract OptimizationData.                     This is type is associated with the functional being                     computed and holds all the relevant data.\n\n\n\n"
},

{
    "location": "solver/euler/boundary_functional.html#EulerEquationMod.calcBoundaryFunctionalIntegrand-Union{Tuple{EulerEquationMod.ParamType{2,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tsol,1},AbstractArray{Tres,1},AbstractArray{Tmsh,N} where N,AbstractArray{Int64,N} where N,EulerEquationMod.BoundaryForceData,AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tmsh where Tres where Tsol",
    "page": "Boundary Functional",
    "title": "EulerEquationMod.calcBoundaryFunctionalIntegrand",
    "category": "method",
    "text": "EulerEquationMod.calcBoundaryFunctionalIntegrand\n\nComputes the integrand for boundary functional at a surface SBP node. Every functional of a different type may need a corresponding method to compute the integrand. The type of the functional object, which is a subtype of AbstractFunctional.\n\nArguments\n\nparams : eqn.params object\nq : Nodal solution\naux_vars : Auxiliary variables\nnrm : Face normal vector in the physical space\nnode_info : Information about the SBP node\nobjective : Functional data type\nval : Function output value\n\n\n\n"
},

{
    "location": "solver/euler/boundary_functional.html#EulerEquationMod.calcBoundaryFunctionalIntegrand_revm-Union{Tuple{EulerEquationMod.ParamType{2,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tsol,1},AbstractArray{Tres,1},AbstractArray{Tmsh,N} where N,AbstractArray{Int64,N} where N,EulerEquationMod.BoundaryForceData,AbstractArray{Tmsh,1},AbstractArray{Tres,1}}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tmsh where Tres where Tsol",
    "page": "Boundary Functional",
    "title": "EulerEquationMod.calcBoundaryFunctionalIntegrand_revm",
    "category": "method",
    "text": "EulerEquationMod. calcBoundaryFunctionalIntegrand_revm\n\nReverse mode for boundary functional integrand w.r.t. nrm. Takes in input val_bar and return nrm_bar for further reverse propagation.\n\nArguments\n\nparams : eqn.params object\nq : Nodal solution\naux_vars : Auxiliary variables\nnrm : Face normal vector in the physical space\nnode_info : Information about the SBP node\nobjective : Functional data type\nnrm_bar : Resulting vector\nval_bar : Nodal portion of the seeding vector\n\n\n\n"
},

{
    "location": "solver/euler/boundary_functional.html#EulerEquationMod.evalFunctional_revm-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim where Tres,Any,ODLCommonTools.AbstractFunctional,String}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Boundary Functional",
    "title": "EulerEquationMod.evalFunctional_revm",
    "category": "method",
    "text": "EulerEquationMod.evalFunctional_revm\n\nReverse mode of EulerEquationMod.evalFunctional, It takes in functional value and return mesh.nrm_bndry_bar. Different functionals will need to be added to the if statement to further extend this function.\n\nArguments\n\nmesh :  Abstract mesh object\nsbp  : Summation-By-Parts operator\neqn  : Euler equation object\nopts : Options dictionary\nfunctionalData : Object of type AbstractFunctional. This is type is associated                     with the functional being computed and holds all the                     relevant data.\nfunctionalName : Name of the functional being evaluated.\n\n\n\n"
},

{
    "location": "solver/euler/boundary_functional.html#PDESolver.evalFunctional-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim where Tres,Any,ODLCommonTools.AbstractFunctional}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Boundary Functional",
    "title": "PDESolver.evalFunctional",
    "category": "method",
    "text": "EulerEquationMod.evalFunctional\n\nHigh level function that evaluates all the functionals specified over various edges. This function is agnostic to the type of the functional being computed and calls a mid level functional-type specific function for the actual evaluation.\n\nArguments\n\nmesh :  Abstract mesh object\nsbp  : Summation-By-Parts operator\neqn  : Euler equation object\nopts : Options dictionary\nfunctionalData : Object of type AbstractFunctional. This is type is associated                     with the functional being computed and holds all the                     relevant data.\n\nOutputs\n\nval: functional value\n\n\n\n"
},

{
    "location": "solver/euler/boundary_functional.html#Euler-Boundary-Functional-1",
    "page": "Boundary Functional",
    "title": "Euler Boundary Functional",
    "category": "section",
    "text": "This page consists of all the functions necessary for computing a boundary functional along the geometric edges of a mesh for Euler equations. A boundary functional should ALWAYS be evaluated by calling evalFunctional which is the highest level function.  Modules = [EulerEquationMod]\n  Pages = [\"solver/euler/boundary_functional.jl\"]"
},

{
    "location": "solver/euler/misc.html#",
    "page": "Misc",
    "title": "Misc",
    "category": "page",
    "text": ""
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.calcA0-Union{Tuple{EulerEquationMod.ParamType{2,:entropy,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1},AbstractArray{Tsol,2}}, Tuple{Tsol}} where Tsol",
    "page": "Misc",
    "title": "EulerEquationMod.calcA0",
    "category": "method",
    "text": "EulerEquationMod.calcA0\n\nThis function calculates the A0 (ie. dq/dv, where q are the conservative   and v are the entropy variables) for a node, and stores it A0\n\nThe formation of A0 is given in Hughes\n\nInputs:     params : ParamType{Tdim, :entropy}     q  : vector of entropy variables, length Tdim + 2\n\nInputs/Outputs:   A0 : (Tdim+2)x(Tdim+2) matrix populated with A0.  Overwritten\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.calcA0-Union{Tuple{EulerEquationMod.ParamType{Tdim,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1},AbstractArray{Tsol,2}}, Tuple{Tdim}, Tuple{Tsol}} where Tsol where Tdim",
    "page": "Misc",
    "title": "EulerEquationMod.calcA0",
    "category": "method",
    "text": "EulerEquationMod.calcA0\n\nThis function calculates the A0 (ie. dq/dq, where q are the conservative   variables at a node), and stores it A0.  This function is provided for   compatability purposes\n\nInputs:     params : ParamType{2, :entropy}     q  : vector of conservative variables, length 4\n\nInputs/Outputs:   A0 : 4x4 (or 5x5 in 3D)  matrix populated with A0.  Overwritten\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.calcA0Inv-Union{Tuple{EulerEquationMod.ParamType{2,:entropy,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1},AbstractArray{Tsol,2}}, Tuple{Tsol}} where Tsol",
    "page": "Misc",
    "title": "EulerEquationMod.calcA0Inv",
    "category": "method",
    "text": "EulerEquationMod.calcA0Inv\n\nCalculates inv(A0), where A0 = dq/dv, where q are the conservative variables   at a node and v are the entropy varaibles at a node, using the entropy   variables.\n\nInputs:     params : ParamType{Tdim, :entropy}     q  : vector (length 4 or 5) of entropy variables at a node\n\nInputs/Outputs:     A0inv : matrix to be populated with inv(A0).  Overwritten.\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.calcA0Inv-Union{Tuple{EulerEquationMod.ParamType{Tdim,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1},AbstractArray{Tsol,2}}, Tuple{Tdim}, Tuple{Tsol}} where Tsol where Tdim",
    "page": "Misc",
    "title": "EulerEquationMod.calcA0Inv",
    "category": "method",
    "text": "EulerEquationMod.calcA0Inv\n\nCalculates inv(A0), where A0 = dq/dq, where q are the conservative variables   at a node.  This function is provided for compatability purposes\n\nInputs:     params : ParamType{2, :entropy}     q  : vector (length Tdim + 2) of conservative variables at a node\n\nInputs/Outputs:     A0inv : matrix to be populated with inv(A0).  Overwritten.\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.calcA1-Union{Tuple{EulerEquationMod.ParamType{2,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1},AbstractArray{Tsol,2}}, Tuple{Tsol}} where Tsol",
    "page": "Misc",
    "title": "EulerEquationMod.calcA1",
    "category": "method",
    "text": "EulerEquationMod.calcA1\n\nThis function calculates the A1 (ie. dF1/dq, where F1 is the first column of   the Euler flux) for a node, aka the flux   Jacobian of the Euler flux in the x direction.  Methods are available for   both conservative and entropy variables.\n\nThe formation of A1 is given in Hughes paper\n\nInputs:     params : ParamType{2, :entropy}     q  : vector of entropy variables, length 4\n\nInputs/Outputs:   A1 : 4x4 matrix to be populated.  Overwritten\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.calcA2-Union{Tuple{EulerEquationMod.ParamType{2,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1},AbstractArray{Tsol,2}}, Tuple{Tsol}} where Tsol",
    "page": "Misc",
    "title": "EulerEquationMod.calcA2",
    "category": "method",
    "text": "EulerEquationMod.calcA2\n\nThis function calculates A2 (ie. dF2/dq, where F2 is the second column of the   Euler flux, aka the flux jacobian in the y direction.   Methods are available for both conservative and entropy variables.\n\nThe formation of A2 is given in Hughes\n\nInputs:     params : ParamType{2, :entropy},     q  : vector of entropy variables, length 4   Inputs/Outputs:   A2 : 4x4 matrix to be populated.  Overwritten\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.calcEntropy-Union{Tuple{EulerEquationMod.ParamType{2,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1}}, Tuple{Tsol}} where Tsol",
    "page": "Misc",
    "title": "EulerEquationMod.calcEntropy",
    "category": "method",
    "text": "EulerEquationMod.calcEntropy\n\nThis function calculates the entropy at a node and returns it.  Method are   available for conservative and entropy variables, 2D or 3D\n\nInputs:     params: ParamType{Tdim, var_type}, used to dispatch to the right method.     q: vector of solution variables at a node.\n\nReturns: entropy\n\nThis is a low level function\n\nAliasing Restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.calcEntropyIR-Union{Tuple{EulerEquationMod.ParamType{Tdim,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1}}, Tuple{Tdim}, Tuple{Tsol}} where Tsol where Tdim",
    "page": "Misc",
    "title": "EulerEquationMod.calcEntropyIR",
    "category": "method",
    "text": "This function calculates the entropy function U used to define the IR   variablesat a node and returns it.   It does not return the physical entropy   s.  This function is agnostic to the dimension of the equation.\n\nInputs:     params: a ParamType     q: a vector of conservative variables at a node\n\nOutputs:     U: the value of the entropy function\n\nAliasing restrictions: none.\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.calcEulerFlux-Union{Tuple{EulerEquationMod.ParamType{2,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1},AbstractArray{Tres,1},AbstractArray{Tmsh,N} where N,AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tres where Tsol where Tmsh",
    "page": "Misc",
    "title": "EulerEquationMod.calcEulerFlux",
    "category": "method",
    "text": "EulerEquationMod.calcEulerFlux\n\nThis function calculates the Euler flux from the conservative variables at    a single node in a particular direction.  2D only.\n\nInputs:    params  : ParamaterType{2, :conservative}    q  : vector of conservative variables    aux_vars : vector of auxiliary variables    dir :  vector in direction to calculate the flux\n\nInputs/Outputs:    F  : vector to populate with the flux\n\nThe Tdim paramater of params determine whether this method or the 3D    version is called.\n\nThis is a low level function\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.calcEulerFlux-Union{Tuple{EulerEquationMod.ParamType{2,:entropy,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1},AbstractArray{Tres,1},AbstractArray{Tmsh,N} where N,AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tres where Tsol where Tmsh",
    "page": "Misc",
    "title": "EulerEquationMod.calcEulerFlux",
    "category": "method",
    "text": "low level function\n\nCalculates the Euler flux from entropy variables\n\nInputs:\nparams : ParameterType{Tdim, :entropy}\nq : vector of entropy variables\naux_vars : vector of auxiliary variables\ndir : vector specifying the direction to caculate the flux\n\nInputs/Outputs:\nF  : vector to populate with the flux\n\nThis is a low level function.  The static parameters of\nthe ParameterType are used to dispatch to the right method for any\ncombination of variable type or equation dimension.\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.calcEulerFlux-Union{Tuple{EulerEquationMod.ParamType{3,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1},AbstractArray{Tres,1},AbstractArray{Tmsh,N} where N,AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tres where Tsol where Tmsh",
    "page": "Misc",
    "title": "EulerEquationMod.calcEulerFlux",
    "category": "method",
    "text": "EulerEquationMod.calcEulerFlux\n\nThis is the 3D method.  All arguments are same as the 2D version.\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.calcEulerFlux_revm-Union{Tuple{EulerEquationMod.ParamType{2,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1},Any,AbstractArray{Tmsh,1},AbstractArray{Tres,1},AbstractArray{Tmsh,1}}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tres where Tsol where Tmsh",
    "page": "Misc",
    "title": "EulerEquationMod.calcEulerFlux_revm",
    "category": "method",
    "text": "###EulerEquationMod.calcEulerFlux_revm\n\nCompute the derivative of the euler flux in reverse mode w.r.t to unit vector flux direction.\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.calcEulerFlux_revq-Union{Tuple{EulerEquationMod.ParamType{2,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1},Any,AbstractArray{Tmsh,1},AbstractArray{Tres,1},AbstractArray{Tsol,1}}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tres where Tsol where Tmsh",
    "page": "Misc",
    "title": "EulerEquationMod.calcEulerFlux_revq",
    "category": "method",
    "text": "###EulerEquationMod.calcEulerFlux_revq\n\nCompute the derivative of the Euler flux in reverse mode w.r.t q\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.calcMaxWaveSpeed-Union{Tuple{Any,Any,EulerEquationMod.EulerData{Tsol,Tres,Tdim,:conservative},Any}, Tuple{Tdim}, Tuple{Tres}, Tuple{Tsol}} where Tres where Tdim where Tsol",
    "page": "Misc",
    "title": "EulerEquationMod.calcMaxWaveSpeed",
    "category": "method",
    "text": "EulerEquationMod.calcMaxWaveSpeed\n\nThis function calculates the maximum wave speed (ie. acoustic wave speed)   present in the domain and returns it.  Methods are available for conservative   and entropy variables.  This function uses eqn.q_vec, not eqn.q to do the   calculation.\n\nThis is a mid level function\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.calcMomentContribution_rev!-Union{Tuple{SummationByParts.AbstractFace{Tsbp},AbstractArray{Tmsh,3},AbstractArray{Tmsh,3},AbstractArray{Tsol,3},AbstractArray{Tsol,3},AbstractArray{Tmsh,1},AbstractArray{Tres,1}}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsbp}, Tuple{Tsol}} where Tres where Tsol where Tmsh where Tsbp",
    "page": "Misc",
    "title": "EulerEquationMod.calcMomentContribution_rev!",
    "category": "method",
    "text": "calcMomentContribution_rev!\n\nThis is the reverse differentiated version of calcMomentContribution!.  See docs of calcMomentContribution! for further details of the primal method.  This function is differentiated with respect to the primal version\'s xsbp and dforce variables.\n\nInputs\n\nsbpface: an SBP face operator type\nxsbp: SBP-face nodes in physical space; [coord, sbp node, face]\ndforce: scaled force at the sbpface nodes; [coord, sbp node, face]\nxyz_about: point about which the moment is taken\nmoment_bar: left multiplies d(moment)/d(xsbp) and d(moment)/d(dforce)\n\nInOuts\n\nxsbp_bar: result of vector Jacobian product; [coord, sbp node, face]\ndforce_bar: result of vector Jacobian product; [coord, sbp node, face]\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.calcPressure-Union{Tuple{EulerEquationMod.ParamType{2,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1}}, Tuple{Tsol}} where Tsol",
    "page": "Misc",
    "title": "EulerEquationMod.calcPressure",
    "category": "method",
    "text": "EulerEquationMod.calcPressure\n\nThis function calculates the pressure from the conservative variables at a   node in 2D.  It returns a single value.\n\nInputs:     params : ParamType{Tdim, var_type }     q  : vector of conservative variables\n\nThe parameter of params determines whether the 2D or 3D, conservative   or entropy method is dispatched.\n\nThis is a low level function.\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.calcPressure-Union{Tuple{EulerEquationMod.ParamType{2,:entropy,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1}}, Tuple{Tsol}} where Tsol",
    "page": "Misc",
    "title": "EulerEquationMod.calcPressure",
    "category": "method",
    "text": "EulerEquationMod.calcPressure\n\nThis function calculates pressure using the entropy variables.\n\nInputs:     params : ParamType{2, :entropy}     q  : vector of entropy varaibles\n\nreturns pressure\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.calcPressure-Union{Tuple{EulerEquationMod.ParamType{3,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1}}, Tuple{Tsol}} where Tsol",
    "page": "Misc",
    "title": "EulerEquationMod.calcPressure",
    "category": "method",
    "text": "EulerEquationMod.calcPressure\n\n3D method.  See 2D method documentation\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.calcPressure_revq-Union{Tuple{EulerEquationMod.ParamType{2,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1},Any,Any}, Tuple{Tsol}} where Tsol",
    "page": "Misc",
    "title": "EulerEquationMod.calcPressure_revq",
    "category": "method",
    "text": "###EulerEquationMod.calcPressure_revq\n\nCompute the gradient of pressure w.r.t q in the reverse mode\n\nArguments\n\nparams : Parameter object\nq : Forward sweep solution variable\npress_bar : Reverse pressure gradient\nq_bar : Reverse mode solution gradient\n\n`\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.calcSpeedofSound-Union{Tuple{EulerEquationMod.ParamType{Tdim,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{Tsol,1}}, Tuple{Tdim}, Tuple{Tsol}} where Tsol where Tdim",
    "page": "Misc",
    "title": "EulerEquationMod.calcSpeedofSound",
    "category": "method",
    "text": "EulerEquationMod.calcSpeedofSound\n\nThis function calculates the speed of sound at a node and returns it.   Methods are available for both conservative and entropy variables.\n\nInputs:     params:  ParamType{Tdim, var_type}     q  vector of solution variables at a node\n\nReturns: speed of sound\n\nThis is a low level function\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.calcSteadyFluxJacobian-NTuple{4,Any}",
    "page": "Misc",
    "title": "EulerEquationMod.calcSteadyFluxJacobian",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.calcVolumeIntegralsSplitForm-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type,Any,ODLCommonTools.FluxType}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tdim where Tres where Tsol where Tmsh",
    "page": "Misc",
    "title": "EulerEquationMod.calcVolumeIntegralsSplitForm",
    "category": "method",
    "text": "Calculate (S .*F)1 where S is the skew symmetric part of sbp.Q and F   is a symmetric numerical flux function.  eqn.res is updated with the result.   Methods are available for curvilinear and non-curvilinear meshes\n\nInputs:     mesh     sbp     eqn     opts     functor: the numerical flux function F, of type FluxType\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.calcVolumeIntegralsSplitFormCurvilinear-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP,SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type,Any,ODLCommonTools.FluxType}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tdim where Tres where Tsol where Tmsh",
    "page": "Misc",
    "title": "EulerEquationMod.calcVolumeIntegralsSplitFormCurvilinear",
    "category": "method",
    "text": "This function calculates I_S2F.\'(S . F( uk(I_S2F wk), uk(I_S2F wk)))1.   This is similar to the other method but for the staggered grid   algorithm.\n\nInputs:     mesh_s: the mesh object for the solution grid     mesh_f: the mesh object for the flux grid     sbp_s: SBP operator for solution grid     sbp_f: SBP operator for flux grid     functor: the numerical flux functor (FluxType) used to compute F\n\nInputs/Outputs:     eqn: the equation object (implicitly on the solution grid).  eqn.res is          updated with the result\n\nAliasing Restrictions: none (the meshes and the sbp operators could alias,                          in which case this algorithm reduces to the                           non-staggered version\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.calcVolumeIntegralsSplitFormCurvilinear-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type,Any,ODLCommonTools.FluxType}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tdim where Tres where Tsol where Tmsh",
    "page": "Misc",
    "title": "EulerEquationMod.calcVolumeIntegralsSplitFormCurvilinear",
    "category": "method",
    "text": "Calculate (S .* F)1, where S is the skew-symmetric part of sbp.Q   and F is a symmetric numerical flux function.  eqn.res is updated   with the result.  This function is used for curvilinear meshes.\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.calcVolumeIntegrals_nopre-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type,Any}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tdim where Tres where Tsol where Tmsh",
    "page": "Misc",
    "title": "EulerEquationMod.calcVolumeIntegrals_nopre",
    "category": "method",
    "text": "Calculates the volume integrals for the weak form, computing the Euler flux   as needed, rather than using eqn.flux_parametric\n\nInputs:     mesh     sbp     eqn: eqn.res is updated with the result     opts\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.calcVorticity-Union{Tuple{EulerEquationMod.ParamType{3,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,Any,AbstractArray{Tsol,2},AbstractArray{Tmsh,3},AbstractArray{Tmsh,1},AbstractArray{Tres,2}}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tres where Tmsh where Tsol",
    "page": "Misc",
    "title": "EulerEquationMod.calcVorticity",
    "category": "method",
    "text": "This function calculates the vorticity at every node of an element   (because vorticity requires derivatives of the velocity, it would be   awkward to compute it at a node level).\n\n3D, conservative variables only\n\nInputs:     params: a ParamType     q: numDofPerNode x numNodesPerElement array of conservative variables at        each node     dxidx: the 3 x 3 x numNodesPerElement  scaled  mapping jacobian at the node     jac: numNodesPerElement vector of the mapping jacobian determinant at          each node\n\nInput/Outputs:     vorticity: a 3 x numNodesPerElement array  containing the vorticity in the                x, y, and z directions at each node, overwritten\n\nAliasing restrictions: from params: dxidx_element, velocities, velocity_deriv,                                       velocity_deriv_xy\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.getAuxVars-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tdim where Tres where Tsol where Tmsh",
    "page": "Misc",
    "title": "EulerEquationMod.getAuxVars",
    "category": "method",
    "text": "EulerEquationMod.getAuxVars\n\nThis function calculates any extra variables that are stored across the mesh   using the conservative variables eqn.q.  Currently only calculates pressure.\n\nThi is a mid level function\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.getEulerFlux-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type,Any}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tdim where Tres where Tsol where Tmsh",
    "page": "Misc",
    "title": "EulerEquationMod.getEulerFlux",
    "category": "method",
    "text": "EulerEquationMod.getEulerFlux\n\nThis function calculates the Euler flux across the entire mesh by passing   pieces of the eqn.q, eqn.aux_vars, eqn.f_xi and eqn.params to a low level   function.  The flux is calculated in the xi and eta directions,   scaled (mulitiplied) by the mapping jacobian (so that when performing the   integral we don\'t have to explictly divide by the jacobian, it just cancels   out with the jacobian factor introduced here.\n\nCalls writeFlux to do any requested output.\n\nThis is a mid level function\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.getEulerFlux2-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim where Tres,Any}, Tuple{Tmsh}, Tuple{Tsol}} where Tsol where Tmsh",
    "page": "Misc",
    "title": "EulerEquationMod.getEulerFlux2",
    "category": "method",
    "text": "EulerEquationMod.getEulerFlux2\n\nThis function calcules the euler flux over the entire mesh directly (ie.   does not call a low level function.  This function is deprecated, although   useful for benchmarking purposes.  2D only.\n\nThis is a mid level function\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.interpolateElementStaggered-Union{Tuple{EulerEquationMod.ParamType{Tdim,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,Any,AbstractArray{T,2} where T,AbstractArray{T,2} where T,AbstractArray{T,2} where T}, Tuple{Tdim}} where Tdim",
    "page": "Misc",
    "title": "EulerEquationMod.interpolateElementStaggered",
    "category": "method",
    "text": "This function does the interpolation from the solution grid to the flux   grid for a single element\n\nInputs:\n\nparams: a ParamType.  The qs_el1 and q_el1 fields are overwritten.\nqs_el: the conservative variables for the element on the solution grid\n\nInputs/Outputs:\n\naux_vars: array to be populated with the aux vars on the staggered grid\nqf_el: conservative variables on the flux grid\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.matVecA0-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,:entropy},Any,AbstractArray{Tsol,3}}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tres where Tdim where Tsol where Tmsh",
    "page": "Misc",
    "title": "EulerEquationMod.matVecA0",
    "category": "method",
    "text": "EulerEquationMod.matVecA0\n\nThis function is the same as matVecA0inv, except it multiplies by A0 not   A0inv.  See its documention.\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.matVecA0inv-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},SummationByParts.AbstractSBP,EulerEquationMod.EulerData{Tsol,Tres,Tdim,:entropy},Any,AbstractArray{Tsol,3}}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tres where Tdim where Tsol where Tmsh",
    "page": "Misc",
    "title": "EulerEquationMod.matVecA0inv",
    "category": "method",
    "text": "EulerEquationMod.matVecA0inv\n\nThis function takes a 3D array and multiplies it in place by the inv(A0)   matrix (calculated at each node), inplace, (where A0 =dq/dv, where q are the   conservative variables and v are some other variables), inplace.   Methods are available for conservative and entropy variables.\n\nFor conservative variables, A0 is the identity matrix and this is a no-op\n\nInputs:     mesh     sbp     eqn     opts\n\nInputs/Outputs:     res_arr: the array to multiply\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.writeFlux-NTuple{4,Any}",
    "page": "Misc",
    "title": "EulerEquationMod.writeFlux",
    "category": "method",
    "text": "EulerEquationMod.writeFlux\n\nThis function writes the real part of Euler flux to a file named Fxi.dat,   space delimited, controlled by the input options \'writeflux\', of type Bool.\n\nThis is a high level function.\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.calcEnstrophy-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},Any,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim,Any,AbstractArray{Tsol,3}}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tmsh where Tres where Tsol",
    "page": "Misc",
    "title": "EulerEquationMod.calcEnstrophy",
    "category": "method",
    "text": "Calculates ( 1/(2V) )*integral(rho * omega dot omega dV), where V is the   volume of the mesh omega is the vorticity vector, and rho is the density.   3D, conservative variables only.  This should work for CG and DG, but   has only been tested for the latter\n\nInputs:     mesh: an AbstractMesh     sbp: an SBP operator     eqn: an EulerData object     opts: options dictionary     q_arr: a 3D array of conservative variables \n\nOutputs:     val: the value of the integral (over the entire domain, not just the          part owned by this process)\n\nAliasing restrictions: see calcVorticity\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.calcKineticEnergy-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},Any,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type,Any,AbstractArray{Tsol,1}}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tmsh where Tdim where Tres where Tsol",
    "page": "Misc",
    "title": "EulerEquationMod.calcKineticEnergy",
    "category": "method",
    "text": "This function calculates ( 1/(2V) )integral(rho * v dot v dV), where   V is the volume of the mesh, v is the velocity vector, and rho is the density.   This is the total kinetic energy normalized by the volume of the domain   Conservative variables only.\n\nThis function contains an MPI blocking collective operation.  It must be   called by all processes at the same time.\n\nThis function relies on the sequential numbering of dofs on the same node\n\nInputs:     mesh: a mesh     sbp: an SBP Operator     eqn: an EulerData object     opts: options dictionary     q_vec: the vector of conservative variables for the entire mesh\n\nOutputs:     val: the value of the integral (over the entire domain, not just teh part          owned by this procss)\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.calcKineticEnergydt-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},Any,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type,Any,AbstractArray{Tsol,1},AbstractArray{Tres,1}}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tmsh where Tdim where Tres where Tsol",
    "page": "Misc",
    "title": "EulerEquationMod.calcKineticEnergydt",
    "category": "method",
    "text": "This function calclates the time derivative of calcKineticEnergy.\n\nThe idea is to expand the left hand side of d rho*u/dt = res using the    product rule and solve for du/dt, then use it to compute the integral.   Inputs:     mesh: a mesh     sbp: a SBP operator     eqn: an EulerData     opts: options dictionary     q_vec: vector of conservative variables for the entire mesh     res_vec: residual vector (dq/dt) of entire mesh\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#EulerEquationMod.contractResEntropyVars2-Union{Tuple{ODLCommonTools.AbstractDGMesh{Tmsh},Any,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type,Any,AbstractArray{T,1} where T,AbstractArray{T,1} where T}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tdim where Tmsh where Tres where Tsol",
    "page": "Misc",
    "title": "EulerEquationMod.contractResEntropyVars2",
    "category": "method",
    "text": "Like contractResEntropyVars, but uses a special summation technique\n\n\n\n"
},

{
    "location": "solver/euler/misc.html#Miscellaneous-Function-1",
    "page": "Misc",
    "title": "Miscellaneous Function",
    "category": "section",
    "text": "  CurrentModule = EulerEquationModThis page documents a bunch of functions that are useful but don\'t fit into another catagory.Actually, some of these functions do fit into another catagory and new files should be created to move them into.A bunch of the things in euler_funcs.jl  Modules = [EulerEquationMod]\n  Pages = [\"euler/euler_funcs.jl\", \"euler/entropy_flux.jl\"]"
},

{
    "location": "solver/euler/homotopy.html#",
    "page": "Homotopy",
    "title": "Homotopy",
    "category": "page",
    "text": ""
},

{
    "location": "solver/euler/homotopy.html#EulerEquationMod.calcHomotopyDiss-Union{Tuple{ODLCommonTools.AbstractDGMesh{Tmsh},Any,EulerEquationMod.EulerData{Tsol,Tres,Tdim,var_type} where var_type where Tdim,Any,AbstractArray{Tres,3}}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tmsh where Tres where Tsol",
    "page": "Homotopy",
    "title": "EulerEquationMod.calcHomotopyDiss",
    "category": "method",
    "text": "Calculate a first order accurate dissipation to use as a homotopy function\n\nInputs:     mesh: a DG mesh     sbp: an SBP operator     eqn: an EulerData object     opts: options dictionary\n\nInputs/Outputs:     res: 3D array to store the homotopy function in\n\nNote eqn.res is not modified by this function.\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/homotopy.html#EulerEquationMod.getLambdaMax-Union{Tuple{EulerEquationMod.ParamType{Tdim,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tsol,1},AbstractArray{Tmsh,1}}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tsol}} where Tdim where Tmsh where Tsol",
    "page": "Homotopy",
    "title": "EulerEquationMod.getLambdaMax",
    "category": "method",
    "text": "Calculate the maximum wave speed for a given state\n\nInputs:     params: a ParamType     qL: a vector of conservative variables at a node     dir: a direction vector, length 2 in 2D and 3 in 3D\n\nOutputs:     lambda_max: the maximum wave speed\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/homotopy.html#PDESolver.evalHomotopy",
    "page": "Homotopy",
    "title": "PDESolver.evalHomotopy",
    "category": "function",
    "text": "This function calls the appropriate homotopy function for the Euler module.\n\n\n\n"
},

{
    "location": "solver/euler/homotopy.html#Homotopy-1",
    "page": "Homotopy",
    "title": "Homotopy",
    "category": "section",
    "text": "This page documents the functions that compute the dissipation function needed by homotopy methods  Modules = [EulerEquationMod]\n  Pages = [\"euler/homotopy.jl\"]"
},

{
    "location": "solver/euler/homotopy_diff.html#",
    "page": "Homotopy Jacobian",
    "title": "Homotopy Jacobian",
    "category": "page",
    "text": ""
},

{
    "location": "solver/euler/homotopy_diff.html#EulerEquationMod.getLambdaMaxSimple_diff-Union{Tuple{EulerEquationMod.ParamType{Tdim,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tsol,1},AbstractArray{Tsol,1},AbstractArray{Tmsh,1},AbstractArray{Tres,1},AbstractArray{Tres,1}}, Tuple{Tdim}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tdim where Tmsh where Tres where Tsol",
    "page": "Homotopy Jacobian",
    "title": "EulerEquationMod.getLambdaMaxSimple_diff",
    "category": "method",
    "text": "Differentiated version of getLambdaMaxSimple\n\nInputs\n\nparams\nqL\nqR\ndir\n\nInputs/Outputs\n\nlambda_dotL: derivative of lambda wrt. qL\nlambda_dotR: derivative of lambda wrt. qR\n\n\n\n"
},

{
    "location": "solver/euler/homotopy_diff.html#EulerEquationMod.getLambdaMax_diff-Union{Tuple{EulerEquationMod.ParamType{2,var_type,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol where var_type,AbstractArray{Tsol,1},AbstractArray{Tmsh,1},AbstractArray{Tres,1}}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tmsh where Tres where Tsol",
    "page": "Homotopy Jacobian",
    "title": "EulerEquationMod.getLambdaMax_diff",
    "category": "method",
    "text": "Differentiated version of getLambdaMax\n\nInputs\n\nparams: ParamType\nqL: vector of conservative variables at a node\ndir: direction vector (can be scaled)\n\nInputs/Outputs\n\nqL_dot: derivative of lambda max wrt qL\n\nOutputs\n\nlambda_max: maximum eigenvalue\n\n\n\n"
},

{
    "location": "solver/euler/homotopy_diff.html#Homotopy-Jacobian-1",
    "page": "Homotopy Jacobian",
    "title": "Homotopy Jacobian",
    "category": "section",
    "text": "This page documents the functions that compute the Jacobian of the  dissipation function needed by homotopy methods  Modules = [EulerEquationMod]\n  Pages = [\"euler/homotopy_diff.jl\"]"
},

{
    "location": "solver/euler/eigensystem.html#",
    "page": "Eigensystem",
    "title": "Eigensystem",
    "category": "page",
    "text": ""
},

{
    "location": "solver/euler/eigensystem.html#Eigensystem-1",
    "page": "Eigensystem",
    "title": "Eigensystem",
    "category": "section",
    "text": "  CurrentModule = EulerEquationModThe eigenvalues and eigenvectors are important for constructing upwinding schemes.  The funtions describe on this page calculate them"
},

{
    "location": "solver/euler/eigensystem.html#EulerEquationMod.calcA-Tuple{EulerEquationMod.ParamType{2,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{T,1} where T,AbstractArray{T,2} where T}",
    "page": "Eigensystem",
    "title": "EulerEquationMod.calcA",
    "category": "method",
    "text": "This function calculate the Euler flux jacobian in the x direction from    the conservative variables.  Methods are available for 2 and 3 dimensions\n\nNote that this code calculates pressure internally (ie. it does not call   calcPressure).\n\nInputs:     params: a ParamType     q: a vector of conservative variables\n\nInputs/Outputs:     A: matrix to be populated with the flux jacobian (overwritten)\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/eigensystem.html#EulerEquationMod.calcB-Tuple{EulerEquationMod.ParamType{2,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{T,1} where T,AbstractArray{T,2} where T}",
    "page": "Eigensystem",
    "title": "EulerEquationMod.calcB",
    "category": "method",
    "text": "Like calcA, but in the y directions.  See that function for details\n\n\n\n"
},

{
    "location": "solver/euler/eigensystem.html#EulerEquationMod.calcC-Tuple{EulerEquationMod.ParamType{3,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{T,1} where T,AbstractArray{T,2} where T}",
    "page": "Eigensystem",
    "title": "EulerEquationMod.calcC",
    "category": "method",
    "text": "Like calcA, but in the z direction.  See that function for details\n\n\n\n"
},

{
    "location": "solver/euler/eigensystem.html#EulerEquationMod.calcEScalingx-Tuple{EulerEquationMod.ParamType{2,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{T,1} where T,AbstractArray{T,1} where T}",
    "page": "Eigensystem",
    "title": "EulerEquationMod.calcEScalingx",
    "category": "method",
    "text": "This function calculates the diagonal scaling matrix S^2 such that:     df/dq = (YS)Lambdainv(YS)  and     dq/dw = (YS)(Y*S).\'\n\nwhere df/dx is the Euler flux jacobian in the x direction, Y and Lambda are   its eigenvectors and eigenvalues, respective (which can be calculated via    calcEvecsx and calcEvalsx),   and dq/dw is the jacobian of the conservative variables wrt the IR entropy   variables.  Recall that eigenvectors are defined up to a multiplicitive    constant, so Y*S is also a matrix of eigenvectors.\n\nThis scaling was originally reported by Merriam: An Entropy-Based Approach   to Nonlinear Stability (NASA TR).\n\nMethods are available for 2 and 3 dimensions.\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/eigensystem.html#EulerEquationMod.calcEScalingy-Tuple{EulerEquationMod.ParamType{2,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{T,1} where T,AbstractArray{T,1} where T}",
    "page": "Eigensystem",
    "title": "EulerEquationMod.calcEScalingy",
    "category": "method",
    "text": "Like calcEScalingx, but in the y direction.  See that function for details.\n\n\n\n"
},

{
    "location": "solver/euler/eigensystem.html#EulerEquationMod.calcEScalingz-Tuple{EulerEquationMod.ParamType{3,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{T,1} where T,AbstractArray{T,1} where T}",
    "page": "Eigensystem",
    "title": "EulerEquationMod.calcEScalingz",
    "category": "method",
    "text": "Like calcEScalingx, but in the y direction.  See that function for details.\n\n\n\n"
},

{
    "location": "solver/euler/eigensystem.html#EulerEquationMod.calcEvalsx-Tuple{EulerEquationMod.ParamType{2,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{T,1} where T,AbstractArray{T,1} where T}",
    "page": "Eigensystem",
    "title": "EulerEquationMod.calcEvalsx",
    "category": "method",
    "text": "This function calculates the eigenvalues for the Euler flux jacobian in the   x direction corresponding to the eigenvectors from calcEvecsx.\n\nInputs:     params: a ParamType     q: a vector of conservative variables\n\nInputs/Outputs:     Lambda: vector to be populated with the eigenvalues (overwritten)\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/eigensystem.html#EulerEquationMod.calcEvalsy-Tuple{EulerEquationMod.ParamType{2,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{T,1} where T,AbstractArray{T,1} where T}",
    "page": "Eigensystem",
    "title": "EulerEquationMod.calcEvalsy",
    "category": "method",
    "text": "Like calcEvalsx, but in the y direction.  See that function for details.\n\n\n\n"
},

{
    "location": "solver/euler/eigensystem.html#EulerEquationMod.calcEvalsz-Tuple{EulerEquationMod.ParamType{3,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{T,1} where T,AbstractArray{T,1} where T}",
    "page": "Eigensystem",
    "title": "EulerEquationMod.calcEvalsz",
    "category": "method",
    "text": "Like calcEvalsx, but in the z direction.  See that function for details\n\n\n\n"
},

{
    "location": "solver/euler/eigensystem.html#EulerEquationMod.calcEvecsx-Tuple{EulerEquationMod.ParamType{2,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{T,1} where T,AbstractArray{T,2} where T}",
    "page": "Eigensystem",
    "title": "EulerEquationMod.calcEvecsx",
    "category": "method",
    "text": "This function calculates the (right) eigenvectors of the Euler flux jacobian   in the x direction at a given state q.  Methods are available in 2 and 3 dimensions\n\nInputs:     params: a ParamType     q: a vector of conservative variables\n\nInputs/Outputs:     R: Matrix whose columns will be populated with the eigenvectors (overwritten)\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "solver/euler/eigensystem.html#EulerEquationMod.calcEvecsy-Tuple{EulerEquationMod.ParamType{2,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{T,1} where T,AbstractArray{T,2} where T}",
    "page": "Eigensystem",
    "title": "EulerEquationMod.calcEvecsy",
    "category": "method",
    "text": "Like calcEvecsx, but for the y direction.  See that function for details\n\n\n\n"
},

{
    "location": "solver/euler/eigensystem.html#EulerEquationMod.calcEvecsz-Tuple{EulerEquationMod.ParamType{3,:conservative,Tsol,Tres,Tmsh} where Tmsh where Tres where Tsol,AbstractArray{T,1} where T,AbstractArray{T,2} where T}",
    "page": "Eigensystem",
    "title": "EulerEquationMod.calcEvecsz",
    "category": "method",
    "text": "Like calcEvecx, but in the z direction.  See that function for details.\n\n\n\n"
},

{
    "location": "solver/euler/eigensystem.html#Functions-1",
    "page": "Eigensystem",
    "title": "Functions",
    "category": "section",
    "text": "  Modules = [EulerEquationMod]\n  Pages = [\"euler/eigensystem.jl\"]"
},

{
    "location": "solver/euler/startup.html#",
    "page": "Startup",
    "title": "Startup",
    "category": "page",
    "text": ""
},

{
    "location": "solver/euler/startup.html#Startup-1",
    "page": "Startup",
    "title": "Startup",
    "category": "section",
    "text": "  CurrentModule = EulerEquationModThis page documents the functions that start the solver."
},

{
    "location": "solver/euler/startup.html#EulerEquationMod.createObjects-Tuple{Dict}",
    "page": "Startup",
    "title": "EulerEquationMod.createObjects",
    "category": "method",
    "text": "This function creates and initializes the mesh, sbp, eqn, and opts objects\n\nInputs:     opts: options dictionary\n\nOutputs:     mesh: an AbstractMesh.  The concrete type is determined by the options           dictionary     sbp: an AbstractSBP.  The concrete type is determined by the options          dictionary     eqn: an EulerData object     opts: the options dictionary     pmesh: mesh used for preconditioning, can be same object as mesh\n\n\n\n"
},

{
    "location": "solver/euler/startup.html#EulerEquationMod.createObjects-Tuple{ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,Dict}",
    "page": "Startup",
    "title": "EulerEquationMod.createObjects",
    "category": "method",
    "text": "Constructs a the EulerData object given a mesh, sbp, and options dictionary.\n\nUsed for submesh solves.\n\nInputs\n\nmesh\nsbp\nopts\n\nOutputs\n\nmesh\nsbp\neqn\nopts\npmesh: currently, always the same as mesh\n\n\n\n"
},

{
    "location": "solver/euler/startup.html#EulerEquationMod.postproc-NTuple{4,Any}",
    "page": "Startup",
    "title": "EulerEquationMod.postproc",
    "category": "method",
    "text": "This function does post processing, if requested by the input options.   Typical post processing includes calculation of errors, norms of important   quantities, writing of files. etc.\n\nInputs:     mesh     sbp     eqn     opts\n\n\n\n"
},

{
    "location": "solver/euler/startup.html#PDESolver.solvePDE",
    "page": "Startup",
    "title": "PDESolver.solvePDE",
    "category": "function",
    "text": "Given fully initialized mesh, sbp, eqn, opts, this function solves   the Euler equations.  The 4 object should be obtained from createObjects().\n\nSpecifically, it applies an initial condition and invokes a nonlinear   solver according to the options dictionary.\n\nInputs:     mesh: an AbstractMesh     sbp: an AbstractSBP     eqn: an AbstractEulerData     opts: the options dictionary.  This must be the options dictionary returned           by createObjects().  Changing values in the options dictionary after           calling createObjects() results in undefined behavior.     pmesh: mesh used for preconditioning, can be same object as mesh.            default value of mesh\n\n\n\n"
},

{
    "location": "solver/euler/startup.html#EulerEquationMod.checkOptions-Tuple{Any}",
    "page": "Startup",
    "title": "EulerEquationMod.checkOptions",
    "category": "method",
    "text": "This function checks the options in the options dictionary after default   values have been supplied and throws exceptions if unsupported/incompatable   options are specified\n\nInputs:     opts: options dictionary\n\nOutputs:     none\n\n\n\n"
},

{
    "location": "solver/euler/startup.html#Functions-1",
    "page": "Startup",
    "title": "Functions",
    "category": "section",
    "text": "  Modules = [EulerEquationMod]\n  Pages = [\"euler/startup_func.jl\", \"euler/check_options.jl\"]"
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
    "text": "Physics modules generally use the majorIterationCallback function to log important quantities to files.  Such logging should be controlled by two keys, \"write_outname\" where outname is the name of the quantity, which has a boolean value, and \"write_outname_fname\" that has a string value containing the name of the file to write (including extension).  Examples of things that can be logged are entropy and kinetic energy.  Both these keys should have default values, and users should generally not need to modify the second one."
},

{
    "location": "input/input.html#Input.make_input-Tuple{Dict,AbstractString}",
    "page": "Introduction",
    "title": "Input.make_input",
    "category": "method",
    "text": "Make an input file from an options dictionary.\n\nInputs:     dict: the options dictionary     fname: the file name (without extension)\n\nOutputs:     fname_ex: file name with extension\n\n\n\n"
},

{
    "location": "input/input.html#Input.read_input-Tuple{AbstractString}",
    "page": "Introduction",
    "title": "Input.read_input",
    "category": "method",
    "text": "PDESolver.read_input\n\nThis function read an input file, supplies default values whenever possible,   and does sanity checks on the options.   A dictionary (referred to as the options dictionary) is returned.   See read_input_file for the description of the input file format.\n\nAfter default values are supplied, the dictionary is printed to   arg_dict_output.jl (by MPI rank 0) in the format of an input file.   This is useful for rerunning a simulation.\n\nThis function prints warnings if keys in the dictionary do not correspond   to known options.  The file known_keys.jl lists all known keys.   See input_vals.txt for the list of keys, possible values, and their meanings.\n\nThis function is idempotent; this is essential for using arg_dict_output.jl   to rerun a simulation.\n\nInputs:     * fname : name of file to read, can be relative or absolute path.\n\nOutputs:     arg_dict: a Dict{Any, Any} containing the option keywords and values\n\n\n\n"
},

{
    "location": "input/input.html#Input.read_input-Tuple{Dict}",
    "page": "Introduction",
    "title": "Input.read_input",
    "category": "method",
    "text": "This method supplies default values to a given dictionary.  See the other   method for details.\n\nInputs\n\narg_dict: a dictionary\n\nOutputs\n\narg_dict: the same dictionary, with default values supplied\n\n\n\n"
},

{
    "location": "input/input.html#Input.read_input_file-Tuple{AbstractString}",
    "page": "Introduction",
    "title": "Input.read_input_file",
    "category": "method",
    "text": "This function reads an input file and turns it into a dictionary.   Unlike read_input, this function returns the contents of the   input file exactly as a dictionary, without supplying default values.\n\nInputs\n\nfname: a string containing the file name.  If the file name does not        start with a . or .. or /, then fname is treated as relative to        the pwd, otherwise fname is interpreted as an absolute path\n\nOutputs\n\narg_dict: the dictionary containing all the options.\n\nExceptions: if the input file does not contain a dictionary, an exception               is thrown.\n\nNotes: currently, the input file format is a literal dictionary that is          processed using eval().  This is likely to change to change in the          future to support static compilation.\n\n\n\n"
},

{
    "location": "input/input.html#Input.checkBCOptions-Tuple{Any}",
    "page": "Introduction",
    "title": "Input.checkBCOptions",
    "category": "method",
    "text": "Check that a model entity does not have more than one boundary condition   applied to it\n\n\n\n"
},

{
    "location": "input/input.html#Input.checkForIllegalOptions_post-Tuple{Any}",
    "page": "Introduction",
    "title": "Input.checkForIllegalOptions_post",
    "category": "method",
    "text": "Check the user supplied options for errors after supplying default options.\n\n\n\n"
},

{
    "location": "input/input.html#Input.checkForIllegalOptions_pre-Tuple{Any}",
    "page": "Introduction",
    "title": "Input.checkForIllegalOptions_pre",
    "category": "method",
    "text": "Check the user supplied options for errors before supplying default values\n\n\n\n"
},

{
    "location": "input/input.html#Input.checkKeys-Tuple{Any,Any}",
    "page": "Introduction",
    "title": "Input.checkKeys",
    "category": "method",
    "text": "PDESolver.checkKeys\n\nThis function verifies all the keys in the first argument are also keys   of the second argument and prints a warning to STDERR if they are not.\n\nInputs     arg_dict: first dictonary     known_keys: second dictonary\n\nOutputs:     cnt: number of unrecognized keys\n\n\n\n"
},

{
    "location": "input/input.html#Key-Validation-1",
    "page": "Introduction",
    "title": "Key Validation",
    "category": "section",
    "text": "After supplying default values, PDESolver checks that all keys in the dictonary are recognized keys.  It does this by comparing against the list of keys documented in the input_vals.txt file.  A warning of is printed to STDERR if an unrecognized key is found. There is a second file, input_vals_internal.txt, which lists keys that the user should never specify. These keys are completely dependent on other keys from input_vals.txt, and specifying them inconsistently can make the solver seriously misbehave.The mechanics of the key validation are as follows.  The extract_keys.jl script reads the input_vals.txt file (and the input_vals_internal.txt), looking for any string that is surrounded by double quotes and starts in the first character of a line.  All keys are written to a dictonary in the file known_keys.jl.  This file is included by PDESolver.The script extract_keys.jl is run as part of PDESolver installation.  If new keys are added it must be run again to silence warnings during key validation.  Modules = [Input]\n  Pages = [\"input/extract_keys.jl\",\n           \"input/Input.jl\",\n           \"input/make_input.jl\",\n           \"input/read_input.jl\"]"
},

{
    "location": "input/keys.html#",
    "page": "Important Keys",
    "title": "Important Keys",
    "category": "page",
    "text": ""
},

{
    "location": "input/keys.html#Important-Keys-1",
    "page": "Important Keys",
    "title": "Important Keys",
    "category": "section",
    "text": "  CurrentModule = PDESolverThe input dictionary contains many keys.  All possible user supplied keys are listed in input_vals.txt. The purpose of this page is to describe the most important keys and use them to elucidate some of the structure of the code.The most important key is physics which specifies which physics to solve. The name must be the name the physics module registers with the front end described in Registration Functions.The run_type key specifies what kind of solver to use.  There are a variety of options, including explicit and implicit time marching methods for unsteady problems and inexact-Newton methods for steady problems. Many of the time marching methods can also be used for pseudo-time stepping of steady problems. See call_nlsolver for the complete list of values.  Physics modules should never invoke a nonlinear solver directly, they should always use call_nlsolver().The jac_method key is needed if the nonlinear solver computes a Jacobian. PDESolver supports multiple methods of computing derivatives, including finite differences and the complex step method.  The functions in the Nonlinear solvers module are required to support all methods. This key specifies which method to use. Complex step is generally recommended, finite difference are used primarily for verification.The order key specifies the degree of the operator used to discretize the spatial terms. This is analagous to the degree of the polynomial basis functions used for finite elements.smb_name is the name of the mesh file.  This is passed directly to the AbstractMesh constructor.  The mesh should already be partitioned into the correct number of parts (ie. the number of MPI processes PDESolver was launched with).  For Pumi, the partitioned mesh is comprised of many numbered files.  The smb_name should not have the numbers appended to it, and should be the same for all processes.dmg_name is the name of the geometry file associated with the mesh.IC_name is the name of the initial condition.  See Registration Functions for how to add new initial conditions.  Note that initial conditions are physics module specific and that getting an error about an unrecognized initial condition is often the result of specifying the incorrect physics.operator_type: the name of the SBP operator to use.  See createSBPOperator."
},

{
    "location": "linearsolvers/linearsolvers.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "linearsolvers/linearsolvers.html#Linear-Solvers-1",
    "page": "Introduction",
    "title": "Linear Solvers",
    "category": "section",
    "text": "  CurrentModule = LinearSolversThe purpose of the LinearSolvers module is to provide a consistant API for preconditioning and solving linear systems, using both direct and iterative methods, serial and parallel.  The API is highly composable, allowing new preconditioners and linear operators to be built on top of existing ones.  One important feature is the ability to explictly control when the linear operator and preconditioner are recalculated, which is important for efficiently solving nonlinear problems.The interface for using a linear solver (including is component preconditioner and linear operator) is described on the LinearSolvers page.The interface for defining a new preconditioner is described on the Preconditioners page, and, similarly, the interface for new linear operators is described on the  Linear Operators](@ref sec:linearoperators) page.  These interfaces must be defined by every new preconditioner and linear operator, respectively, but they should not be used directly.  The LinearSolvers interface should be used instead.  Pages = [\"pc.md\"\n           \"lo.md\"\n           \"ls.md\"\n          ]\n  Depth = 1"
},

{
    "location": "linearsolvers/pc.html#",
    "page": "Preconditioners",
    "title": "Preconditioners",
    "category": "page",
    "text": ""
},

{
    "location": "linearsolvers/pc.html#sec:preconditioners-1",
    "page": "Preconditioners",
    "title": "Preconditioners",
    "category": "section",
    "text": "  CurrentModule = LinearSolversThis section describes the different catagories of preconditioners, including the API and the requirements for creating new preconditioners. Note that whenever a linear solver is used, the Linear Solver interface should be used rather than the Linear Operator interface described here."
},

{
    "location": "linearsolvers/pc.html#LinearSolvers.AbstractPC",
    "page": "Preconditioners",
    "title": "LinearSolvers.AbstractPC",
    "category": "type",
    "text": "Abstract supertype of all preconditioner types.  Preconditioners can be used   with iterative methods only.  When using a direct method, PCNone   should be used.\n\nThe purpose of this type is to provide a consistent interface for different   types of preconditioners.  In particular, it provides control over when the   preconditioner is recomputed.\n\nUsers should implement new preconditioners using composition with   one of the preconditioners defined here, namely PetscMatPC or   PetscMatFreePC, ie. they should define a new subtype of   AbstractPC that has either PetscMatPC or PetscMatFrePC as   a field.  This allows calling the existing functions for these types   (which compute a preconditioner for the Jacobian of the physics) and   modifying the preconditioner as needed.   User defined preconditioners should subtype either AbstractPetscMatPC   of AbstractPetscMatFreePC.\n\nFields\n\npc_inner: another AbstractPC\n\nNote that arbitrarily deep nesting of preconditioners is allowed.   The pc_inner field can be one of PCNone, PetscMatPC,   or PetscMatFreePC for a non-nested preconditioner, or some   other AbstractPC for a nested preconditioner.\n\n\n\n"
},

{
    "location": "linearsolvers/pc.html#LinearSolvers.AbstractPetscMatPC",
    "page": "Preconditioners",
    "title": "LinearSolvers.AbstractPetscMatPC",
    "category": "type",
    "text": "Abstract supertype of all Petsc matrix-explicit preconditioners\n\n\n\n"
},

{
    "location": "linearsolvers/pc.html#LinearSolvers.AbstractPetscMatFreePC",
    "page": "Preconditioners",
    "title": "LinearSolvers.AbstractPetscMatFreePC",
    "category": "type",
    "text": "Abstract supertype of all Petsc matrix-free preconditioners.\n\n\n\n"
},

{
    "location": "linearsolvers/pc.html#LinearSolvers.AbstractPetscPC",
    "page": "Preconditioners",
    "title": "LinearSolvers.AbstractPetscPC",
    "category": "constant",
    "text": "Alias for any kind of Petsc PC (matrix-explicit or matrix-free)\n\n\n\n"
},

{
    "location": "linearsolvers/pc.html#Type-hierarchy-1",
    "page": "Preconditioners",
    "title": "Type hierarchy",
    "category": "section",
    "text": "When creating a new preconditioner, it is very important to make it inherit from the proper supertype.AbstractPC\nAbstractPetscMatPC\nAbstractPetscMatFreePC\nAbstractPetscPC"
},

{
    "location": "linearsolvers/pc.html#LinearSolvers.calcPC-Tuple{LinearSolvers.AbstractPC,ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,ODLCommonTools.AbstractSolutionData,Dict,Any,Any}",
    "page": "Preconditioners",
    "title": "LinearSolvers.calcPC",
    "category": "method",
    "text": "This function calculates the preconditioner.  Every implementation of   AbstractPC should extend this function with a new method\n\nWhen creating new preconditioners, this function should generally be called   first on pc_inner, and modifications to the Jacobian should be made    subsequently.\n\nInputs\n\npc: the AbstractPC implementation\nmesh\nsbp\neqn\nopts\nctx_residual: the ctx required by physicsRhs\nt: current time\n\nImplementation Notes:     For matrix-explicit preconditioners, this might not actually calculate the     preconditioner.  Rather, it calculates the matrix the preconditioner is     based on, and the solve function calcultes the preconditioner from it.     Nevertheless, it supports the proper semantics for when the PC is updated     even when the PC matrix and LinearOperator matrix are the same (as long as     the solve function is called regularly).\n\n\n\n"
},

{
    "location": "linearsolvers/pc.html#LinearSolvers.applyPC-Tuple{LinearSolvers.AbstractPC,ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,ODLCommonTools.AbstractSolutionData,Dict,Any,AbstractArray{T,1} where T,AbstractArray{T,1} where T}",
    "page": "Preconditioners",
    "title": "LinearSolvers.applyPC",
    "category": "method",
    "text": "Applies the preconditioner, ie. x = inv(Ap)*b, where Ap is the approximation   to the matrix A.  Note that Ap itself may not be available for some   preconditioners, hence there is only an API for applying inv(Ap),   not Ap itself.\n\nMatrix-free preconditioners need to extend this function with a new method,\nmatrix-explicit preconditioners do not.\n\nInputs\n\npc: the AbstractPC implementation.\nmesh\nsbp\neqn: this argument should generally not be used because all the solution      related data should be stored in the pc object by calcPC\nopts\nt: current time\nb: a AbstractVector representing the local part of the solution (ie    eqn.q_vec)\n\nInputs/Outputs\n\nx: AbstractVector updated with results (same size as b) (do not overwrite)\n\n\n\n"
},

{
    "location": "linearsolvers/pc.html#LinearSolvers.applyPCTranspose-Tuple{LinearSolvers.AbstractPC,ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,ODLCommonTools.AbstractSolutionData,Dict,Any,AbstractArray{T,1} where T,AbstractArray{T,1} where T}",
    "page": "Preconditioners",
    "title": "LinearSolvers.applyPCTranspose",
    "category": "method",
    "text": "Applies the transpose of the preconditioner, ie. x = inv(Ap).\'*b.   Similar to applyPC, see that function for details.\n\nNote that not every preconditioning method supports this.\n\n\n\n"
},

{
    "location": "linearsolvers/pc.html#LinearSolvers.getBasePC",
    "page": "Preconditioners",
    "title": "LinearSolvers.getBasePC",
    "category": "function",
    "text": "This function returns the underlying preconditioner object, ie.   PCNone, PetscMatPC, or PetscMatFreePC.\n\nNote that arbitrarily deep nesting of preconditioners is allowed.   Users do not have to implement as long as the nested preconditioner is   stored in a field called pc_inner.\n\nFor matrix-explicit preconditioners, this function is useful for getting   the PetscMatPC object, which contains the preconditioning   Jacobian matrix.\n\nInputs\n\npc: the users AbstractPC\n\n\n\n"
},

{
    "location": "linearsolvers/pc.html#Utils.free-Tuple{LinearSolvers.AbstractPC}",
    "page": "Preconditioners",
    "title": "Utils.free",
    "category": "method",
    "text": "This function frees any memory belonging to external libraries.  Users must   call this function when they are finished with an AbstractPC   object.\n\nUsers do not have to define this function for their   AbstractPC types.\n\nInputs\n\npc: the AbstractPC object\n\n\n\n"
},

{
    "location": "linearsolvers/pc.html#API-1",
    "page": "Preconditioners",
    "title": "API",
    "category": "section",
    "text": "Every preconditioner supports the following functions.  When defining a new preconditioner, some of the functions must be defined for the new  AbstractPC type, while others are defined automatically based on the supertype of the preconditioner.calcPC(::AbstractPC, ::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict, ::Any, ::Any)\napplyPC(::AbstractPC, ::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict, ::Any, ::AbstractVector, ::AbstractVector)\napplyPCTranspose(::AbstractPC, ::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict, ::Any, ::AbstractVector, ::AbstractVector)\ngetBasePC\nfree(::AbstractPC)"
},

{
    "location": "linearsolvers/pc.html#LinearSolvers.PCNone",
    "page": "Preconditioners",
    "title": "LinearSolvers.PCNone",
    "category": "type",
    "text": "Preconditioner type for direct solve.   Do not use with PetscLinearOperator.\n\nPublic Fields\n\nnone\n\n\n\n"
},

{
    "location": "linearsolvers/pc.html#LinearSolvers.PCNone-Tuple{ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,ODLCommonTools.AbstractSolutionData,Dict}",
    "page": "Preconditioners",
    "title": "LinearSolvers.PCNone",
    "category": "method",
    "text": "Outer constructor\n\nInputs\n\nmesh\nsbp\neqn\nopts\n\n\n\n"
},

{
    "location": "linearsolvers/pc.html#LinearSolvers.PetscMatPC",
    "page": "Preconditioners",
    "title": "LinearSolvers.PetscMatPC",
    "category": "type",
    "text": "AbstractPC implementation for Petsc matrix-explicit preconditioners.\n\nPublic Fields\n\nA: a PetscMat object used to calculate the preconditioner\n\n\n\n"
},

{
    "location": "linearsolvers/pc.html#LinearSolvers.PetscMatPC-Tuple{ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,ODLCommonTools.AbstractSolutionData,Dict}",
    "page": "Preconditioners",
    "title": "LinearSolvers.PetscMatPC",
    "category": "method",
    "text": "Outer constructor\n\nInputs\n\nmesh\nsbp\neqn\nopts\n\n\n\n"
},

{
    "location": "linearsolvers/pc.html#LinearSolvers.PetscMatFreePC",
    "page": "Preconditioners",
    "title": "LinearSolvers.PetscMatFreePC",
    "category": "type",
    "text": "Type for managing a Petsc matrix-free preconditioner.  Actual preconditioners   can be built on top of this.  Because matrix-free preconditioners do not have   a particular form, this type requires the preconditioner to implement the   entire AbstractPC interface.  The benefit of using this type to   build a preconditioner is that is automatically exposes the preconditioner   to Petsc.\n\nThe calcPC function the user defines must call setPCCtx.\n\nPublic Fields\n\nnone\n\n\n\n"
},

{
    "location": "linearsolvers/pc.html#LinearSolvers.PetscMatFreePC-Tuple{ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,ODLCommonTools.AbstractSolutionData,Dict}",
    "page": "Preconditioners",
    "title": "LinearSolvers.PetscMatFreePC",
    "category": "method",
    "text": "Outer constructor\n\nInputs\n\nmesh\nsbp\neqn\nopts\n\n\n\n"
},

{
    "location": "linearsolvers/pc.html#LinearSolvers.setPCCtx",
    "page": "Preconditioners",
    "title": "LinearSolvers.setPCCtx",
    "category": "function",
    "text": "This function should be called by the users calcPC function.   The arguments passed to this function are passed to applyPC   during the next preconditioner application.   See that function for a description of the arguments\n\nInputs\n\npc:\nmesh\nsbp\neqn\nopts\nctx_residual\nt\n\n\n\n"
},

{
    "location": "linearsolvers/pc.html#Concrete-PC-Types-1",
    "page": "Preconditioners",
    "title": "Concrete PC Types",
    "category": "section",
    "text": "The following types are the \"Base\" PC types refered to by getBasePC. Every preconditioner must contain one of these (either directly, in the pc_inner field described by AbstractPC, of inside another preconditioner).PCNone\nPCNone(::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict)\nPetscMatPC\nPetscMatPC(::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict)\nPetscMatFreePC\nPetscMatFreePC(::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict)PetscMatFreePC has an additional function not needed by other preconditioner types:setPCCtx"
},

{
    "location": "linearsolvers/pc.html#Implementing-a-New-Preconditioner-1",
    "page": "Preconditioners",
    "title": "Implementing a New Preconditioner",
    "category": "section",
    "text": "Many of the API functions are automatically generated for matrix-explicit preconditioners.  Matrix-free preconditioners need to define most of the PC interface themselves.  A summary of the functions each preconditioner must implements is:Matrix-explicitcalcPCMatrix-freecalcPC\napplyPC\napplyPCTranspose"
},

{
    "location": "linearsolvers/pc.html#First-Level-PC-1",
    "page": "Preconditioners",
    "title": "First Level PC",
    "category": "section",
    "text": "This example shows how to implement a new preconditioner that directly   contains one of the Concrete PC Types described above.type ExamplePetscMatPC <: AbstractPetscMatPC\n  pc_inner::PetscMatPC  # the supertype of pc_inner and ExamplePetscMatPC\n                        # must match\nend\n\nfunction ExamplePetscMatPC(mesh::AbstractMesh, sbp::AbstractSBP,\n                    eqn::AbstractSolutionData, opts::Dict)\n\n  pc_inner = PetscMatPC(mesh, sbp, eqn, opts)\n\n  return ExamplePetscMatPC(pc_inner)\nend\n\n\nfunction calcPC(pc::ExamplePetscMatPC, mesh::AbstractMesh, sbp::AbstractSBP,\n                eqn::AbstractSolutionData, opts::Dict, ctx_residual, t)\n\n\n  calcPC(pc.pc_inner)\n\n  pc2 = getBasePC(pc)  # because pc_inner is the base pc, this returns\n                       # pc.pc_inner\n  # call set_values1!(pc2.Ap, ...) to put values into the preconditioning matrix\n\n  return nothing\nendBecause calcPC is defined and the base PC is one of the concrete PC types, namely PetscMatPC, the functions applyPC and applyPC are automatically defined for ExamplePetscMatPC. The new preconditioner is now fully usable."
},

{
    "location": "linearsolvers/pc.html#Multiple-Levels-of-Composition-1",
    "page": "Preconditioners",
    "title": "Multiple Levels of Composition",
    "category": "section",
    "text": "The same structure as shown in the previous structure can be used to build a new preconditioner out of any existing preconditioner, not only the Concrete PC Types.type OuterPetscMatPC <: AbstractPetscMatPC\n  pc_inner::ExamplePetscMatPC  # use the PC define in the previous example\nend\n\nfunction OuterPetscMatPC(mesh::AbstractMesh, sbp::AbstractSBP,\n                         eqn::AbstractSolutionData, opts::Dict)\n\n  pc_inner = ExamplePetscMatPC(mesh, sbp, eqn, opts)\n\n  return OuterPetscMatPC(pc_inner)\nend\n\n\nfunction calcPC(pc::OuterPetscMatPC, mesh::AbstractMesh, sbp::AbstractSBP,\n                eqn::AbstractSolutionData, opts::Dict, ctx_residual, t)\n\n\n  calcPC(pc.pc_inner)  # always call the inner PC calc function first\n\n  pc2 = getBasePC(pc)  # in this case, pc2 is the PetscMatPC inside \n                       # the ExamplePetscMatPC\n  # call set_values1!(pc2.Ap, ...) to put values into the preconditioning matrix  # because calcPC(pc.pc_inner) was called already, the values here should\n  # be added to those already in pc2.Ap\n\n  return nothing\nendAs before, this preconditioner is now fully usable."
},

{
    "location": "linearsolvers/pc.html#Matrix-Free-1",
    "page": "Preconditioners",
    "title": "Matrix Free",
    "category": "section",
    "text": "An outline of how a matrix free preconditioner should be implemented is:type ExamplePetscMatFreePC <: AbstractPetscMatFreePC\n  pc_inner::PetscMatFreePC\n  # other fields...\nend\n\nfunction ExamplePetscMatFreePC(mesh::AbstractMesh, sbp::AbstractSBP,\n                         eqn::AbstractSolutionData, opts::Dict)\n\n  pc_inner = PetscMatFreePC(mesh, sbp, eqn, opts)\n\n  return ExamplePetscMatFreePC(pc_inner)\nend\n\nfunction calcPC(pc::ExamplePetscMatFreePC, mesh::AbstractMesh, sbp::AbstractSBP,\n                eqn::AbstractSolutionData, opts::Dict, ctx_residual, t)\n\n  # do any setup operations, storing data into the \"other fields...\" of the pc\n\n  # these arguments will be passed into applyPC() the next time it is called\n  # by Petsc\n  # always do this last\n  setPCCtx(pc, mesh, sbp, eqn, opts, ctx_residual, t)\n\n  return nothing\nend\n\nfunction applyPC(pc::PetscMatFreePC, mesh::AbstractMesh, sbp::AbstractSBP,\n                 eqn::AbstractSolutionData, opts::Dict, t, b::AbstractVector, \n                 x::AbstractVector)\n\n  # use the data in pc to multiply the preconditioner by b, storing result in x\n\n  return nothing\nend\n\nfunction applyPCTranspose(pc::PetscMatFreePC, mesh::AbstractMesh,\n                 sbp::AbstractSBP,\n                 eqn::AbstractSolutionData, opts::Dict, t, b::AbstractVector, \n                 x::AbstractVector)\n\n  # use the data in pc to multiply the transpose of the preconditioner by b,\n  # storing result in x\n\n  return nothing\nendMatrix-free preconditioners support the same kind of composition as matrix explicit preconditioners, so it is recommended to add the result into x, rather than overwrite x when pc_inner is not one of the Concrete PC Types."
},

{
    "location": "linearsolvers/lo.html#",
    "page": "Linear Operators",
    "title": "Linear Operators",
    "category": "page",
    "text": ""
},

{
    "location": "linearsolvers/lo.html#sec:linearoperators-1",
    "page": "Linear Operators",
    "title": "Linear Operators",
    "category": "section",
    "text": "  CurrentModule = LinearSolversThis section describes the different catagories of linear operators, including the API and the requirements for defining new linear operators. Note that whenever a linear solver is used, the Linear Solver interface should be used rather than the Linear Operator interface described here."
},

{
    "location": "linearsolvers/lo.html#LinearSolvers.AbstractLO",
    "page": "Linear Operators",
    "title": "LinearSolvers.AbstractLO",
    "category": "type",
    "text": "Abstract supertype of all linear operators used for A when solving Ax = b.\n\nThe purpose of this type is to provide a consistent interface for different   types of linear operators.  This type really combines two notions: what the   type of the linear operator is, and how it should be solved.\n\nAny implementation of this type should subtype the appropriate catagory   of: AbstractDenseLO, AbstractSparseDirectLO,   AbstractPetscMatLO, AbstractPetscMatFreeLO\n\nNote that matrix-explicit implementations can often write a single function   for all these cases if using an matrix interface functions that are defined   for all the matrix types.  See MatExplicitLO\n\nRequired Fields\n\nlo_inner: another AbstractPC.  Can be one of DenseLO,           SparseDirectLO, PetscMatLO, or           PetscMatFreeLO, or any other user defined linear           operator.\n\n\n\n"
},

{
    "location": "linearsolvers/lo.html#LinearSolvers.AbstractDenseLO",
    "page": "Linear Operators",
    "title": "LinearSolvers.AbstractDenseLO",
    "category": "type",
    "text": "Linear operator type for Dense matrices.  This is generally used only for   debugging.\n\n\n\n"
},

{
    "location": "linearsolvers/lo.html#LinearSolvers.AbstractSparseDirectLO",
    "page": "Linear Operators",
    "title": "LinearSolvers.AbstractSparseDirectLO",
    "category": "type",
    "text": "Linear operator type for SparseMatrixCSC matrices, which use a direct   solver.\n\n\n\n"
},

{
    "location": "linearsolvers/lo.html#LinearSolvers.AbstractPetscMatLO",
    "page": "Linear Operators",
    "title": "LinearSolvers.AbstractPetscMatLO",
    "category": "type",
    "text": "Linear operator type for Petsc matrix-explicit.\n\n\n\n"
},

{
    "location": "linearsolvers/lo.html#LinearSolvers.AbstractPetscMatFreeLO",
    "page": "Linear Operators",
    "title": "LinearSolvers.AbstractPetscMatFreeLO",
    "category": "type",
    "text": "Linear operator type for Petsc matrix-free.\n\n\n\n"
},

{
    "location": "linearsolvers/lo.html#LinearSolvers.MatExplicitLO",
    "page": "Linear Operators",
    "title": "LinearSolvers.MatExplicitLO",
    "category": "constant",
    "text": "Useful union for all the matrix-explicit linear operator types.   Because matrices have a small set of common interface functions, it is   often possible to write a single function that works on all the different   types of matrices.\n\n\n\n"
},

{
    "location": "linearsolvers/lo.html#LinearSolvers.PetscLO",
    "page": "Linear Operators",
    "title": "LinearSolvers.PetscLO",
    "category": "constant",
    "text": "Union of Petsc linear operator types\n\n\n\n"
},

{
    "location": "linearsolvers/lo.html#LinearSolvers.DirectLO",
    "page": "Linear Operators",
    "title": "LinearSolvers.DirectLO",
    "category": "constant",
    "text": "Union of linear operators that do direct solves\n\n\n\n"
},

{
    "location": "linearsolvers/lo.html#Type-Hierarchy-1",
    "page": "Linear Operators",
    "title": "Type Hierarchy",
    "category": "section",
    "text": "When creating a new linear operator, it is very important to make it inherit from the proper abstract type.AbstractLO\nAbstractDenseLO\nAbstractSparseDirectLO\nAbstractPetscMatLO\nAbstractPetscMatFreeLOSeveral typealiases are also available:MatExplicitLO\nPetscLO\nDirectLO"
},

{
    "location": "linearsolvers/lo.html#LinearSolvers.calcLinearOperator",
    "page": "Linear Operators",
    "title": "LinearSolvers.calcLinearOperator",
    "category": "function",
    "text": "This function calculates the linear operator.  Every implementation of   AbstractLO should extend this function with a new   method.  For matrix-free operators, this function must exist but need   not perform any actions.\n\nFor matrix-explicit implementations, this function should be called on   lo_inner first and modifications to the Jacobian made subsequently.\n\nInputs\n\nlo: the AbstractLO implementation (fields may be updated)\nmesh\nsbp\neqn\nopts\nctx_residual: the ctx required by physicsRhs\nt: current time\n\nImplementation Notes:     For matrix-free operation, this function sets the Petsc ctx for the     PetscMat, which contains a reference to the mesh, sbp, eqn, opts arguments.     This could lead to unexpected behavior if those arguments are modified and     this function is not called again before the next solve.\n\n\n\nCalculates the linear operator.  Use this function only if you want to   calculate the linear operator and not the preconditioner.   Prefer calcPCandLO, which avoids calculating the matrix twice if   the preconditioner and linear operator share the same matrix\n\nInputs\n\nls: StandardLinearSolver\nmesh\nsbp\neqn\nopts\nctx_residual: the ctx required by physicsRhs like functions\nt: current time\n\nKeyword Arguments\n\nstart_comm: start parallel communication (if required by the lo), default              false.  This means the user is generally required to make sure              parallel communication is started before calling this              function.\n\n\n\n"
},

{
    "location": "linearsolvers/lo.html#LinearSolvers.applyLinearOperator",
    "page": "Linear Operators",
    "title": "LinearSolvers.applyLinearOperator",
    "category": "function",
    "text": "Applies the linear operator, ie. , Ax = b\n\nMatrix-explicit implementations AbstractLO do not have   to implement this function, though matrix-free implementations must extend   it with a new method.\n\nInputs\n\nlo: the AbstractLO implementation.\nmesh\nsbp\neqn\nopts\nctx_residual: the ctx for physicsRhs or the another right hand               side function built on top of it\nt: current time\nx: an AbstractVector (although never a PetscVec)\n\nInputs/Outputs\n\nb: vector updated with results (do not overwrite)\n\n\n\n"
},

{
    "location": "linearsolvers/lo.html#LinearSolvers.applyLinearOperatorTranspose",
    "page": "Linear Operators",
    "title": "LinearSolvers.applyLinearOperatorTranspose",
    "category": "function",
    "text": "Applies the transpose of the linear operator, ie. A.\'*x = b   Similar to applyLinearOperator, see that function for details.\n\nNote that not every method supports this.  In particular, Petsc   matrix-free LinearOperators don\'t currently expose this    (although they could with enough reverse-mode)\n\n\n\n"
},

{
    "location": "linearsolvers/lo.html#LinearSolvers.getBaseLO",
    "page": "Linear Operators",
    "title": "LinearSolvers.getBaseLO",
    "category": "function",
    "text": "Similar to getBasePC except it gets the underlying linear operator,   ie. one of DenseLO, SparseDirectLO, PetscMatLO   or PetscMatFreeLO.\n\nFor matrix-explicit methods, this is a good way of getting the underlying   linear operator object, which contains the matrix in the A field (for   all matrix-explicit linear operators).\n\nInputs\n\nlo: an AbstractLO\n\n\n\n"
},

{
    "location": "linearsolvers/lo.html#LinearSolvers.getBaseObject",
    "page": "Linear Operators",
    "title": "LinearSolvers.getBaseObject",
    "category": "function",
    "text": "This function calls either getBasePC or getBaseLO   depending on the type of its argument.\n\nInputs\n\nlo: either an AbstractLO or AbstractPC object\n\nOutputs\n\nthe base PC or LO object\n\n\n\n"
},

{
    "location": "linearsolvers/lo.html#Utils.free-Tuple{LinearSolvers.AbstractLO}",
    "page": "Linear Operators",
    "title": "Utils.free",
    "category": "method",
    "text": "This function frees any memory belonging to external libraries.  Users must   call this function when they are finished with an AbstractLO   object.\n\nUsers do not have to define this function for their   AbstractLO types.\n\nInputs\n\nlo: the AbstractLO object\n\n\n\n"
},

{
    "location": "linearsolvers/lo.html#API-1",
    "page": "Linear Operators",
    "title": "API",
    "category": "section",
    "text": "Every linear operator supports the following functions.  When defining a new linear operator, some of the functions must be extended with new methods, while others will be created automatically based on the supertype of the linear operator.calcLinearOperator\napplyLinearOperator\napplyLinearOperatorTranspose\ngetBaseLO\ngetBaseObject\nfree(::AbstractLO)"
},

{
    "location": "linearsolvers/lo.html#LinearSolvers.DenseLO",
    "page": "Linear Operators",
    "title": "LinearSolvers.DenseLO",
    "category": "type",
    "text": "Dense array linear operator.  Serial only.\n\nPublic Fields\n\nA: an Array{Float64, 1}\n\n\n\n"
},

{
    "location": "linearsolvers/lo.html#LinearSolvers.DenseLO-Tuple{LinearSolvers.PCNone,ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,ODLCommonTools.AbstractSolutionData,Dict}",
    "page": "Linear Operators",
    "title": "LinearSolvers.DenseLO",
    "category": "method",
    "text": "Outer constructor for DenseLO\n\nInputs\n\npc: a PCNone\nmesh\nsbp\neqn\nopts\n\n\n\n"
},

{
    "location": "linearsolvers/lo.html#LinearSolvers.SparseDirectLO",
    "page": "Linear Operators",
    "title": "LinearSolvers.SparseDirectLO",
    "category": "type",
    "text": "LinearOperator type to be used with sparse direct solve.  Note that the   sparsity pattern of the matrix A must be constant.\n\nFields\n\nA: the matrix\n\n\n\n"
},

{
    "location": "linearsolvers/lo.html#LinearSolvers.SparseDirectLO-Tuple{LinearSolvers.PCNone,ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,ODLCommonTools.AbstractSolutionData,Dict}",
    "page": "Linear Operators",
    "title": "LinearSolvers.SparseDirectLO",
    "category": "method",
    "text": "Outer constructor for SparseDirectLO\n\nInputs\n\npc: a PCNone\nmesh\nsbp\neqn\nopts\n\n\n\n"
},

{
    "location": "linearsolvers/lo.html#LinearSolvers.PetscMatLO",
    "page": "Linear Operators",
    "title": "LinearSolvers.PetscMatLO",
    "category": "type",
    "text": "Petsc matrix-explicit linear operator.\n\nNote that if an interface common to both Petsc and regular matrices is used   when accessing the Ap field, it is possible to write functions that   operate on both regular and Petsc linear operators.\n\nPublic Fields\n\nA: a PetscMat\n\n\n\n"
},

{
    "location": "linearsolvers/lo.html#LinearSolvers.PetscMatLO-Tuple{LinearSolvers.PetscMatPC,ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,ODLCommonTools.AbstractSolutionData,Dict}",
    "page": "Linear Operators",
    "title": "LinearSolvers.PetscMatLO",
    "category": "method",
    "text": "Outer constructor for PetscMatLO\n\nInputs\n\npc: a Petsc preconditioner of some kind\nmesh\nsbp\neqn\nopts\n\n\n\n"
},

{
    "location": "linearsolvers/lo.html#LinearSolvers.PetscMatFreeLO",
    "page": "Linear Operators",
    "title": "LinearSolvers.PetscMatFreeLO",
    "category": "type",
    "text": "Petsc matrix-free linear operator.\n\nPublic Fields\n\nnone\n\n\n\n"
},

{
    "location": "linearsolvers/lo.html#LinearSolvers.PetscMatFreeLO-Tuple{LinearSolvers.PetscMatPC,ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,ODLCommonTools.AbstractSolutionData,Dict}",
    "page": "Linear Operators",
    "title": "LinearSolvers.PetscMatFreeLO",
    "category": "method",
    "text": "Outer constructor for PetscMatFreeLO\n\nInputs\n\npc: a Petsc preconditioner of some kind\nmesh\nsbp\neqn\nopts\n\n\n\n"
},

{
    "location": "linearsolvers/lo.html#LinearSolvers.setLOCtx",
    "page": "Linear Operators",
    "title": "LinearSolvers.setLOCtx",
    "category": "function",
    "text": "This function sets the ctx for the underlying Petsc matrix object.  The   end result is that the mesh, sbp, eqn, and opts that are passed into this   function will be passed to the user-defined applyLinearOperator   during the next linear solve.  See that function for a description of   the arguments.\n\nInputs\n\nlo: the user defined linear operator\nmesh\nsbp\neqn\nopts\nctx_residual\nt\n\n\n\n"
},

{
    "location": "linearsolvers/lo.html#Concrete-LO-Types-1",
    "page": "Linear Operators",
    "title": "Concrete LO Types",
    "category": "section",
    "text": "DenseLO\nDenseLO(::PCNone, ::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict)\nSparseDirectLO\nSparseDirectLO(::PCNone, ::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict)\nPetscMatLO\nPetscMatLO(::PetscMatPC, ::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict)\nPetscMatFreeLO\nPetscMatFreeLO(::PetscMatPC, ::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict)PetscMatFreePC has an additional function not needed by other linear operator types:setLOCtx"
},

{
    "location": "linearsolvers/lo.html#Implementing-a-New-Linear-Operator-1",
    "page": "Linear Operators",
    "title": "Implementing a New Linear Operator",
    "category": "section",
    "text": "Many of the API functions are automatically generated for matrix-explicit linear operators.  Matrix-free linear operators need to define most of the LO interface themselves.  A summary of the functions each preconditioner must implement is:Matrix-explicitcalcLinearOperatorMatrix-freecalcLinearOperator\napplyLinearOperator\napplyLinearOperatorTranspose"
},

{
    "location": "linearsolvers/lo.html#Example-1",
    "page": "Linear Operators",
    "title": "Example",
    "category": "section",
    "text": "Linear operators use the same composition structure as preconditioners, see the Implementing a New Preconditioner section before proceeding.# Assume the ExamplePetscMatLO type has already been defined\n\ntype OuterPetscMatLO <: AbstractPetscMatPC\n  lo_inner::ExamplePetscMatLO\nend\n\n\nfunction OuterPetscMatPC(pc::OuterPetscMatLO, mesh::AbstractMesh,\n                         sbp::AbstractSBP, eqn::AbstractSolutionData, opts::Dict)\n\n  lo_inner = ExamplePetscMatLO(pc, mesh, sbp, eqn, opts)\n  return OuterPetscMatPC(lo_inner)\nend\n\nfunction calcLinearOperator(lo::OuterPetscMatLO, mesh::AbstractMesh,\n                            sbp::AbstractSBP, eqn::AbstractSolutionData,\n                            opts::Dict, ctx_residual, t)\n\n  # have inner LO do its computation\n  calcLinearOperator(lo.lo_inner, mesh, sbp, eqn, opts, ctx_residual, t)\n\n  lo2 = getBaseLO(lo)\n  # use set_values1!(lo2.A, ...) to modify the matrix A\n\n  return nothing\nendapplyLinearOperator and applyLinearOperatorTranspose are defined automatically for all AbstractPetscMatLO types (using getBaseLO)."
},

{
    "location": "linearsolvers/ls.html#",
    "page": "Linear Solvers",
    "title": "Linear Solvers",
    "category": "page",
    "text": ""
},

{
    "location": "linearsolvers/ls.html#sec:linearsolvers-1",
    "page": "Linear Solvers",
    "title": "Linear Solvers",
    "category": "section",
    "text": "  CurrentModule = LinearSolversThis section describes the interface for a linear solver.  This also provides a means to access some of the functions implemented by AbstractPC and AbstractLO. Whenever a linear solver is used, the linear solver API should be used rather than accessing the AbstractPC and AbstractLO"
},

{
    "location": "linearsolvers/ls.html#Type-Hierarchy-1",
    "page": "Linear Solvers",
    "title": "Type Hierarchy",
    "category": "section",
    "text": "LinearSolver\nStandardLinearSolver\nStandardLinearSolver(::Any, ::Any, ::MPI.Comm)"
},

{
    "location": "linearsolvers/ls.html#LinearSolvers.calcPC-Tuple{LinearSolvers.StandardLinearSolver,ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,ODLCommonTools.AbstractSolutionData,Dict,Any,Any}",
    "page": "Linear Solvers",
    "title": "LinearSolvers.calcPC",
    "category": "method",
    "text": "Calculates the preconditioner for the linear solver.  Thsi preconditioner   will be used for all linear solves until this function is called again.\n\nFor direct solvers, this function calculates the linear operator itself.   Prefer calcPCandLO whenever possible.\n\nInputs\n\nls: StandardLinearSolver\nmesh\nsbp\neqn\nopts\nctx_residual: the ctx required by physicsRhs like functions\nt: current time\n\nKeyword Arguments\n\nstart_comm: start parallel communication (if required by the PC), default              false.  This means the user is generally required to make sure              parallel communication is started before calling this              function.\n\n\n\n"
},

{
    "location": "linearsolvers/ls.html#LinearSolvers.calcLinearOperator-Tuple{LinearSolvers.StandardLinearSolver,ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,ODLCommonTools.AbstractSolutionData,Dict,Any,Any}",
    "page": "Linear Solvers",
    "title": "LinearSolvers.calcLinearOperator",
    "category": "method",
    "text": "Calculates the linear operator.  Use this function only if you want to   calculate the linear operator and not the preconditioner.   Prefer calcPCandLO, which avoids calculating the matrix twice if   the preconditioner and linear operator share the same matrix\n\nInputs\n\nls: StandardLinearSolver\nmesh\nsbp\neqn\nopts\nctx_residual: the ctx required by physicsRhs like functions\nt: current time\n\nKeyword Arguments\n\nstart_comm: start parallel communication (if required by the lo), default              false.  This means the user is generally required to make sure              parallel communication is started before calling this              function.\n\n\n\n"
},

{
    "location": "linearsolvers/ls.html#LinearSolvers.calcPCandLO-Tuple{LinearSolvers.StandardLinearSolver,ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,ODLCommonTools.AbstractSolutionData,Dict,Any,Any}",
    "page": "Linear Solvers",
    "title": "LinearSolvers.calcPCandLO",
    "category": "method",
    "text": "Calculates both the preconditioner and linear operator.  In the case where   they share the matrix, the calculation is only performed once.  This function   should be preferred to calling calcPC and calcLinearOperator one   after the other\n\nInputs\n\nls: StandardLinearSolver\nmesh\nsbp\neqn\nopts\nctx_residual: the ctx required by physicsRhs like functions\nt: current time\n\nKeyword Arguments\n\nstart_comm: start parallel communication (if required by the lo), default              false.  This means the user is generally required to make sure              parallel communication is started before calling this              function.  Note that it is not possible to interleave              communication and computation in this case, so performing              communication will be very expensive.\n\n\n\n"
},

{
    "location": "linearsolvers/ls.html#LinearSolvers.applyPC-Tuple{LinearSolvers.StandardLinearSolver,ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,ODLCommonTools.AbstractSolutionData,Dict,Any,AbstractArray{T,1} where T,AbstractArray{T,1} where T}",
    "page": "Linear Solvers",
    "title": "LinearSolvers.applyPC",
    "category": "method",
    "text": "Apply the preconditioner, ie. x = inv(Ap)*b, where Ap is an approximation   to the linear operator used for preconditioner.  This function also works   for matrix-free methods.\n\nNote that calcPC or calcPCandLO must be called before   this function can be used.\n\nFor direct methods, a linear solve is performed because there is no   preconditioner.\n\n\n\n"
},

{
    "location": "linearsolvers/ls.html#LinearSolvers.applyPCTranspose-Tuple{LinearSolvers.StandardLinearSolver,ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,ODLCommonTools.AbstractSolutionData,Dict,Any,AbstractArray{T,1} where T,AbstractArray{T,1} where T}",
    "page": "Linear Solvers",
    "title": "LinearSolvers.applyPCTranspose",
    "category": "method",
    "text": "Like applyPC, but applies the transpose (if possible, otherwise   throws an error).\n\n\n\n"
},

{
    "location": "linearsolvers/ls.html#LinearSolvers.linearSolve",
    "page": "Linear Solvers",
    "title": "LinearSolvers.linearSolve",
    "category": "function",
    "text": "Solves the linear system Ax = b for x. This function does not recompute the   precondtioner or linear operator,  The preconditioner and linear operator   used are the ones calculated by the most recent call to   calcPC, calcLinearOperator (or calcPCandLO.\n\nFor Petsc matrices, this function does the final matrix assembly.\n\nInputs\n\nls: the StandardLinearSolver object.\nb: the right hand side vector (local process portion only)\nverbose: verbosity level, default 5, < 4 indicates no output\n\nInputs/Outputs\n\nx: vector overwritten with result (local process portion only)\n\nImplementation Notes:\n\nThe preconditioner and matrix factorization (for direct solves) might   be computed during this function, but they will be computed at the state   correspondong to the last call to the functions that calculate them.\n\n\n\n"
},

{
    "location": "linearsolvers/ls.html#LinearSolvers.linearSolveTranspose",
    "page": "Linear Solvers",
    "title": "LinearSolvers.linearSolveTranspose",
    "category": "function",
    "text": "Similar to linearSolve, but solves A.\'x = b.  See that function   for details.\n\n\n\n"
},

{
    "location": "linearsolvers/ls.html#LinearSolvers.isLOMatFree",
    "page": "Linear Solvers",
    "title": "LinearSolvers.isLOMatFree",
    "category": "function",
    "text": "Returns true if the linear operator is matrix free, false otherwise\n\n\n\n"
},

{
    "location": "linearsolvers/ls.html#LinearSolvers.isPCMatFree",
    "page": "Linear Solvers",
    "title": "LinearSolvers.isPCMatFree",
    "category": "function",
    "text": "Returns true if the preconditioner is matrix free, false otherwise\n\n\n\n"
},

{
    "location": "linearsolvers/ls.html#LinearSolvers.setTolerances",
    "page": "Linear Solvers",
    "title": "LinearSolvers.setTolerances",
    "category": "function",
    "text": "Set the tolerances for iterative solves.  Has no effect for direct solves.   Supplying a negative value results in retaining the original value.\n\nInputs\n\nls: StandardLinearSolver\nreltol: relative residual tolerance\nabstol: absolute residual tolerance\ndtol: divergence tolerance\nitermax: maximum number of iterations\n\n\n\n"
},

{
    "location": "linearsolvers/ls.html#Utils.free-Tuple{LinearSolvers.StandardLinearSolver}",
    "page": "Linear Solvers",
    "title": "Utils.free",
    "category": "method",
    "text": "This function frees any memory owned by external libraries, both in the    StandardLinearSolver object itself and in the pc and lo objects.   Therefore, at the end of a run, all you need to do is free the   StandardLinearSolver and everything will be taken care of.\n\nIt is safe to call this function multiple times.\n\n\n\n"
},

{
    "location": "linearsolvers/ls.html#API-1",
    "page": "Linear Solvers",
    "title": "API",
    "category": "section",
    "text": "calcPC(::StandardLinearSolver, ::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict, ::Any, ::Any)\ncalcLinearOperator(::StandardLinearSolver, ::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict, ::Any, ::Any)\ncalcPCandLO(::StandardLinearSolver, ::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict, ::Any, ::Any)\napplyPC(::StandardLinearSolver, ::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict, ::Any, ::AbstractVector, ::AbstractVector)\napplyPCTranspose(::StandardLinearSolver, ::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict, ::Any, ::AbstractVector, ::AbstractVector)\nlinearSolve\nlinearSolveTranspose\nisLOMatFree\nisPCMatFree\nsetTolerances\nfree(::StandardLinearSolver)"
},

{
    "location": "NonlinearSolvers/nonlinearsolvers.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "NonlinearSolvers/nonlinearsolvers.html#sec:nonlinearsolvers-1",
    "page": "Introduction",
    "title": "NonlinearSolvers Introduction",
    "category": "section",
    "text": "This module contains both the time-stepping methods and the methods for solving steady problems. When solving a CFD problem, these are the functions that drive the code. The purpose of the physics modules is to evaluate the spatial residual.  The purpose of the NonlinearSolvers is to solve the problem. The steady methods and implicit time-marching methods generally use some kind of Newton iteration to solve a nonlinear problem.The first two pages describe the steady and unsteady methods, respectively, and the final pages describe the Newton implementation used by them.  Pages = [ \"steady.md\"\n            \"unsteady.md\"\n            \"newton.md\"\n            \"jacobian.md\"\n            \"jac_recalc.md\"\n            \"matrix.md\"\n            \"newton_inner.md\"\n          ]\n  Depth = 1"
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
    "text": "CN is selected with \"run_type\" => 20 in your problem\'s input options file. As this is an unsteady problem, you will also need to specify your    time step size (\"delta_t\" option) and your final solution time (\"t_max\" option)"
},

{
    "location": "NonlinearSolvers/unsteady/cn.html#Forward-in-time-1",
    "page": "Crank-Nicolson",
    "title": "Forward-in-time",
    "category": "section",
    "text": "\\begin{aligned} \\underbrace{\\left( \\mathbf{I} - \\frac{\\Delta t}{2} \\underbrace{\\frac{\\partial R\\left(u^{i+1}\\right)}{\\partial u^{i+1}}}_{\\text{jac from physicsJac}} \n\\right)}_{\\text{cnJac}} \\Delta u\n= \\\n\\underbrace{- \\left( \n  \\underbrace{u^{i+1}}_{\\text{eqn_nextstep.q_vec}} - \\      \\frac{\\Delta t}{2} \\underbrace{R\\left(u^{i+1}\\right)}_{\\text{eqn_nextstep.res_vec}} - \\    \\underbrace{u^i}_{\\text{eqn.q_vec}} - \\      \\frac{\\Delta t}{2} \\underbrace{R\\left(u^i\\right)}_{\\text{eqn.res_vec}} \n\\right)}_{\\text{cnRhs}} \\end{aligned}This equation is solved with PDESolver\'s Newton\'s method."
},

{
    "location": "NonlinearSolvers/unsteady/cn_uadj.html#",
    "page": "Crank-Nicolson: Unsteady Adjoint",
    "title": "Crank-Nicolson: Unsteady Adjoint",
    "category": "page",
    "text": ""
},

{
    "location": "NonlinearSolvers/unsteady/cn_uadj.html#Crank-Nicolson-Unsteady-Adjoint-–-EXPERIMENTAL,-INCORRECT-CODE-1",
    "page": "Crank-Nicolson: Unsteady Adjoint",
    "title": "Crank-Nicolson Unsteady Adjoint – EXPERIMENTAL, INCORRECT CODE",
    "category": "section",
    "text": ""
},

{
    "location": "NonlinearSolvers/unsteady/cn_uadj.html#Current-status-1",
    "page": "Crank-Nicolson: Unsteady Adjoint",
    "title": "Current status",
    "category": "section",
    "text": "As of mid-September, 2017, development of PDESolver\'s unsteady adjoint has been tabled for the time being. The code is preserved, and is accessible with the run_flag of 660. It runs with no errors, and produces a solution that qualititatively demonstrates properties of the correct unsteady adjoint. However, a test sensitivity check does not pass.Notation for this section:R is the residual\nA is the design variable: the amplitude of a sin function that is an exact solution to the advection equation:    \\begin{equation}   u = A \\sin(-x + \\omega t)   \\end{equation}\nJ is the objective function:   \\begin{equation}   J = \\int_{\\Gamma_1} u^2 d\\Gamma_1   \\end{equation}\n\\Gamma_1 is the right domain boundary in a square domainSome notes about current status that may be of assistance to further development or debugging:dRdA makes physical sense, and FD and CS methods match\ndJdu FD, CS, and analytical derivative match\nJ is verified by setting a polynomial solution for which SBP should be exact.\nloaded checkpoints match forward solve\nnorm of Jacobian calculated from loaded checkpoint matches forward solve\'s\ntime stepping bookkeeping in forward and reverse appears correct \nAdjoint initial condition equation appears to match the below derivation, as far as code-reading can show\nWhile the test sensitivity check is incorrect at all time steps,  the fact that it is incorrect at the adjoint initial condition indicates that the  bug manifests itself before the main reverse-sweep time-stepping loop.The test sensitivity check being performed is the comparison between these two derivatives:\\begin{equation} \\frac{d J}{d A} &= \\frac{\\partial J}{\\partial u} \\frac{\\partial u}{\\partial A} \\\n\\end{equation}\\begin{equation} \\frac{d J}{d A} &= \\psi^T \\left( - \\frac{\\partial R}{\\partial A}\\right) \\end{equation}For the above, note that: \\begin{equation} \\frac{\\partial J}{\\partial A} = 0 \\end{equation}"
},

{
    "location": "NonlinearSolvers/unsteady/cn_uadj.html#Unsteady-adjoint-derivation-1",
    "page": "Crank-Nicolson: Unsteady Adjoint",
    "title": "Unsteady adjoint derivation",
    "category": "section",
    "text": "The unsteady adjoint derivation starts with the generic Lagrangian equation:\\begin{equation} \\mathcal{L}(u, \\psi) = \\psi^T R(u) + J(u) \\end{equation}In the discrete context of CN, all of these variables are global-in-time. That is, the adjoint vector contains the adjoint at time step 1 concatenated with    the adjoint at time step 2, and so on, until time step n. Therefore, in this document we will rewrite the Lagrangian using bolded symbols to indicate    that a vector or matrix is global-in-time, as there will also be corresponding variables   specific to a particular time step:\\begin{equation} \\boldsymbol{\\mathcal{L}}(\\boldsymbol{u}, \\boldsymbol{\\psi}) = \\boldsymbol{\\psi}^T \\boldsymbol{R}(\\boldsymbol{u}) + \\boldsymbol{J}(\\boldsymbol{u}) \\end{equation}The global-in-time residual discretized according to the Crank-Nicolson method is:boldsymbolR(boldsymbolu) = beginbmatrix u_1 - u_0 - fracDelta t2 R(u_1) - fracDelta t2 R(u_0)  u_2 - u_1 - fracDelta t2 R(u_2) - fracDelta t2 R(u_1)  vdots  u_i - u_i-1 - fracDelta t2 R(u_i) - fracDelta t2 R(u_i-1)  u_i+1 - u_i - fracDelta t2 R(u_i+1) - fracDelta t2 R(u_i)  vdots  u_n - u_n-1 - fracDelta t2 R(u_n) - fracDelta t2 R(u_n-1) endbmatrixThe global-in-time adjoint vector is:boldsymbolpsi^T = psi_1^T psi_2^T dots psi_i^T psi_i+1^T dots psi_n^TNote that each time step\'s adjoint variable is a vector of length equal to the number of degrees of freedom in the mesh. And finally, the global-in-time objective function vector is:boldsymbolJ^T = J_1 J_2 dots J_i J_i+1 dots J_nTherefore, the full discrete Lagrangian is:boldsymbolmathcalL(boldsymbolu boldsymbolpsi) = boldsymbolpsi^T boldsymbolR(boldsymbolu) + boldsymbolJ(boldsymbolu) = beginbmatrix psi_1^T left( u_1 - u_0 - fracDelta t2 R(u_1) - fracDelta t2 R(u_0) right)  psi_2^T left( u_2 - u_1 - fracDelta t2 R(u_2) - fracDelta t2 R(u_1) right)  vdots  psi_i^T left( u_i - u_i-1 - fracDelta t2 R(u_i) - fracDelta t2 R(u_i-1) right)  psi_i+1^T left( u_i+1 - u_i - fracDelta t2 R(u_i+1) - fracDelta t2 R(u_i) right)  vdots  psi_n^T left( u_n - u_n-1 - fracDelta t2 R(u_n) - fracDelta t2 R(u_n-1) right) endbmatrix + beginbmatrix J(u_1)  J(u_2)  vdots  J(u_i)  J(u_i+1)  vdots  J(u_n) endbmatrixTaking the derivative of the Lagrangian with respect to the state at step i yields, for values of i not equal to 0 or n:fracpartial boldsymbolLpartial u_i = underbracepsi_i^T - psi_i^T fracDelta t2 fracpartial R(u_i)partial u_i_textcontribution from boldsymbolR(u_i) - underbracepsi_i+1^T - psi_i+1^T fracDelta t2 fracpartial R(u_i)partial u_textcontribution from boldsymbolR(u_i+1) + fracpartial J(u_i)partial u_i= 0^TOr, rearranging:fracpartial boldsymbolLpartial u_i = (psi_i - psi_i+1) - fracDelta t2 left( fracpartial R(u_i)partial u_i right)^T (psi_i + psi_i+1) + fracpartial J(u_i)partial u_i = 0"
},

{
    "location": "NonlinearSolvers/unsteady/cn_uadj.html#Initial-Condition-1",
    "page": "Crank-Nicolson: Unsteady Adjoint",
    "title": "Initial Condition",
    "category": "section",
    "text": "The derivative of the Lagrangian with respect to the state at the final step i = n is:fracpartial boldsymbolLpartial u_n = psi_n - fracDelta t2 left( fracpartial R(u_n)partial u_n right)^T psi_n + fracpartial J(u_n)partial u_n = 0Therefore, the value of the adjoint at time step n, which is the initial condition for the reverse sweep, is:psi_n = left( left(I - fracDelta t2 fracpartial R(u_n)partial u_n right)^T right)^-1 left( - fracpartial J(u_n)partial u_n right)^T"
},

{
    "location": "NonlinearSolvers/unsteady/cn_uadj.html#Direct-Solve-1",
    "page": "Crank-Nicolson: Unsteady Adjoint",
    "title": "Direct Solve",
    "category": "section",
    "text": "The method of performing a direct solve to advance the CN reverse sweep (as opposed to using Newton\'s method to converge each time step) starts with the restatement of the derivative of the Lagrangian at time step i:fracpartial boldsymbolLpartial u_i = underbracepsi_i^T - psi_i^T fracDelta t2 fracpartial R(u_i)partial u_i_textcontribution from boldsymbolR(u_i) - underbracepsi_i+1^T - psi_i+1^T fracDelta t2 fracpartial R(u_i)partial u_textcontribution from boldsymbolR(u_i+1) + fracpartial J(u_i)partial u_i= 0^TRearranging:left psi_i - fracDelta t2 left( fracpartial R(u_i)partial u_i right)^T psi_i right - left psi_i+1 + fracDelta t2 left( fracpartial R(u_i)partial u_i right)^T psi_i+1 right + fracpartial J(u_i)partial u_i = 0Grouping terms to isolate psi_i:left I - fracDelta t2 left( fracpartial R(u_i)partial u_i right)^T right psi_i = left psi_i+1 + fracDelta t2 left( fracpartial R(u_i)partial u_i right)^T psi_i+1 right - fracpartial J(u_i)partial u_iSolving for psi_i:psi_i = left I - fracDelta t2 left( fracpartial R(u_i)partial u_i right)^T right^-1 left( left psi_i+1 + fracDelta t2 left( fracpartial R(u_i)partial u_i right)^T psi_i+1 right - fracpartial J(u_i)partial u_i right)Therefore, psi_i is a function of 1) the Jacobian of the primal solution at step i, which is loaded from checkpointed data, 2) the derivative of the objective function with respect to the state, at step i, and 3) the adjoint solution at time step i+1.  The adjoint solution sweep is thus stepped backwards in time, starting at time step n."
},

{
    "location": "NonlinearSolvers/unsteady/cn_uadj.html#Checkpointing-1",
    "page": "Crank-Nicolson: Unsteady Adjoint",
    "title": "Checkpointing",
    "category": "section",
    "text": "Currently, all time steps are checkpointed.  Eventually, Revolve will be implemented, for which a separate Julia package has been developed.  See here for the publication discussing the Revolve algorithm."
},

{
    "location": "NonlinearSolvers/unsteady/cn_uadj.html#Global-in-time-Jacobian-1",
    "page": "Crank-Nicolson: Unsteady Adjoint",
    "title": "Global-in-time Jacobian",
    "category": "section",
    "text": "For reference, the structure of the global-in-time Jacobian is shown here. It should never be formed except in the course of debugging very simple use cases,    but it can be helpful for visualizing the matrix form of CN for all space and time.(Work in progress)"
},

{
    "location": "NonlinearSolvers/newton.html#",
    "page": "Newton\'s Method",
    "title": "Newton\'s Method",
    "category": "page",
    "text": ""
},

{
    "location": "NonlinearSolvers/newton.html#Newton\'s-method-1",
    "page": "Newton\'s Method",
    "title": "Newton\'s method",
    "category": "section",
    "text": "CurrentModule = NonlinearSolversNewton\'s method is intended to compute updates to some q by solving the following equation.\\begin{equation} \\frac{\\partial R(q)}{\\partial q} \\Delta q = -R(q) \\end{equation}In the most basic implementation of Newton\'s method in PDESolver, q corresponds to the solution,    and f(q) corresponds to the residual evaluation of the currently selected physics module.An example of a more sophisticated use of Newton\'s method is within the Crank-Nicolson timestepper,    which adds another layer on top of the physics residual:\\begin{equation} \\frac{\\partial g(R(q))}{\\partial q} \\Delta q = -g(R(q)) \\end{equation}"
},

{
    "location": "NonlinearSolvers/newton.html#NonlinearSolvers.newtonInner",
    "page": "Newton\'s Method",
    "title": "NonlinearSolvers.newtonInner",
    "category": "function",
    "text": "This function contains the algorithm for Newton\'s method.  It is intended   to be used by other functions rather than invoked as a solver directly.   For example newton uses newtonInner to solve steady problems while   crank_nicolson uses it to solve the nonlinear problem arising   from the time discretization.\n\nFor reference, Newton\'s method solves the problem R(q) = 0 using     dR/dq delta q = -R(q)   where R is the rsidual and q is the solution\n\nOn entry, eqn.q_vec must contain the initial guess for q.  On exit, eqn.q_vec   will contain the solution to f(q) = 0.  eqn.q will also be consistent with   eqn.q_vec, as will the send and receive buffers in eqn.shared_data\n\nOn exit, rhs_vec will have the residual corresponding to eqn.q_vec in it,   with the imaginary part set to zero.\n\nThe AbstractLO and AbstractPC supplied   inside the LinearSolver object must match the jac_type.\n\nWhen doing inexact Newton-Krylov, newonInner modifies the tolerances   of the linear solver.  Users calling newtonInner repeatedly with the   same linear solver object should reset the initial tolerances as needed.\n\nInputs:\n\nnewton_data: NewtonData object, typically obtained from setupNewton\nmesh: a mesh object\nsbp: an SBP operator\neqn: a solution data object\nopts: options dictionary\nls: a LinearSolver with the preconditioner and linear operator fully     initialized.  \nrhs_vec: vector to store R(q) in\nctx_residual: extra data required by rhs_func\n\nThe user must supply two functions, one to calculate the residual vector   (referred to as rhs_vec), and another to compute the Jacobian.\n\nrhs_func should compute (eqn.q_vec) -> (rhs_vec) and have the signature\n\nrhs_func(mesh, sbp, eqn, opts, rhs_vec, ctx_residual, t=0.0)\n\nNote that the ctx_residual passed into newtonInner is passed directly   to rhs_func.  The contents of ctx_residual can be whatever is needed by   rhs_func to perform its computation.  See physicsRhs for   the specific requirements on rhs_func and for an example implementation.\n\nThe same ctx_residual passed into newtonInner is passed directly to   calcLinearOperator..\n\nThis function supports jacobian/preconditioner freezing using the   prefix \"newton\".\n\nAliasing restrictions: None.  In particular, rhs_vec can alias eqn.res_vec,                          and this leads so some efficiency because it avoids                          needlessly copying data.\n\n\n\n"
},

{
    "location": "NonlinearSolvers/newton.html#NonlinearSolvers.setupNewton",
    "page": "Newton\'s Method",
    "title": "NonlinearSolvers.setupNewton",
    "category": "function",
    "text": "NonlinearSolvers.setupNewton\n\nPerforms setup work for newtonInner, including creating a    NewtonData object.\n\nThis function also resets the implicit Euler globalization.\n\nalloc_rhs: keyword arg to allocate a new object or not for rhs_vec                 true (default) allocates a new vector                 false will use eqn.res_vec\n\nrhs_func: only used for Petsc in matrix-free mode to do Jac-vec products             should be the rhs_func passed into newtonInner   ctx_residual: ctx_residual passed into newtonInner\n\nAllocates Jac & RHS\n\nSee cleanupNewton to the cleanup function\n\n\n\n"
},

{
    "location": "NonlinearSolvers/newton.html#Utils.free-Tuple{NonlinearSolvers.NewtonData}",
    "page": "Newton\'s Method",
    "title": "Utils.free",
    "category": "method",
    "text": "Cleans up after running Newton\'s method.\n\nInputs\n\nnewton_data: the NewtonData object\n\n\n\n"
},

{
    "location": "NonlinearSolvers/newton.html#NonlinearSolvers.getNewtonPCandLO",
    "page": "Newton\'s Method",
    "title": "NonlinearSolvers.getNewtonPCandLO",
    "category": "function",
    "text": "Returns the Newton precondtioner and linear operator specified by the options   dictionary\n\nInputs\n\nmesh\nsbp\neqn\nopts\nrhs_func: rhs_func required by newtonInner\n\n\n\n"
},

{
    "location": "NonlinearSolvers/newton.html#Features-1",
    "page": "Newton\'s Method",
    "title": "Features",
    "category": "section",
    "text": "PDESolver\'s Newton\'s method has a wide variety of features.  It contains the Jacobian calculation routines, which can be performed currently using:finite-differencing \ncomplex-stepThe Jacobian functions can act upon any arbitrary residual.Additionally, the following matrix forms are supported:Julia dense\nJulia sparse\nPETSc sparse\nMatrix-freeThe function that performs the Newton iteration is newtonInnerThe private data required by Newtons method is stored in the NewtonData object.setupNewton\nfree(::NewtonData)\ngetNewtonPCandLO"
},

{
    "location": "NonlinearSolvers/newton.html#Newton-internals-1",
    "page": "Newton\'s Method",
    "title": "Newton internals",
    "category": "section",
    "text": ""
},

{
    "location": "NonlinearSolvers/newton.html#NonlinearSolvers.NewtonData",
    "page": "Newton\'s Method",
    "title": "NonlinearSolvers.NewtonData",
    "category": "type",
    "text": "This type holds the data required by newtonInner as well as   configuration settings.\n\nPublic Fields\n\nmyrank: MPI rank\ncommsize: MPI communicator size\nitr: number of iterations\nres_norm_i: current iteration residual norm\nres_norm_i_1: previous iteration residual norm\nstep_norm_i: current iteration newton step norm\nstep_norm_i_1: previous iteration newton step norm\nres_norm_rel: norm of residual used as the reference point when computing               relative residuals.  If this is -1 on entry to newtonInner,               then the norm of the initial residual is used.\nstep_fac: factor used in step size limiter\nres_reltol: nonlinear relative residual tolerance\nres_abstol: nonlinear residual absolute tolerance\nstep_tol: step norm tolerance\nitermax: maximum number of newton iterations\nuse_inexact_nk: true if inexact-NK should be used, false otherwise\nkrylov_gamma: parameter used by inexact newton-krylov\nrecalc_policy: a RecalculationPolicy.\nls: a LinearSolver\nfconv: convergence.dat file handle (or DevNull if not used)\nverbose: how much logging/output to do\n\nOptions Keys\n\nIf res_reltol0 is negative, the residual of the initial condition will be   used for res_norm_rel\n\nnewton_verbosity is used to the verbose field\n\n\n\n"
},

{
    "location": "NonlinearSolvers/newton.html#NonlinearSolvers.reinitNewtonData",
    "page": "Newton\'s Method",
    "title": "NonlinearSolvers.reinitNewtonData",
    "category": "function",
    "text": "Reinitialized the NewtonData object for a new solve.\n\nNote that this does not reset the linear solver, which might be a problem   if inexact newton-krylov was used for the previous solve\n\n\n\n"
},

{
    "location": "NonlinearSolvers/newton.html#NonlinearSolvers.recordResNorm",
    "page": "Newton\'s Method",
    "title": "NonlinearSolvers.recordResNorm",
    "category": "function",
    "text": "Records the most recent nonlinear residual norm in the NewtonData object.   Also updates the implicit Euler globalization\n\nInputs\n\nnewton_data: the NewtonData\nres_norm: the residual norm\n\n\n\n"
},

{
    "location": "NonlinearSolvers/newton.html#NonlinearSolvers.recordStepNorm",
    "page": "Newton\'s Method",
    "title": "NonlinearSolvers.recordStepNorm",
    "category": "function",
    "text": "Records norm of the most recent newton step (ie. the norm of delta q)   in the NewtonData object\n\nInputs\n\nnewton_data: the NewtonData object\nstep_norm: the norm of the step\n\n\n\n"
},

{
    "location": "NonlinearSolvers/newton.html#NewtonData-API-1",
    "page": "Newton\'s Method",
    "title": "NewtonData API",
    "category": "section",
    "text": "NewtonData has a small API used by Newton\'s method.NewtonData\nreinitNewtonData\nrecordResNorm\nrecordStepNorm"
},

{
    "location": "NonlinearSolvers/newton.html#NonlinearSolvers.NewtonMatPC",
    "page": "Newton\'s Method",
    "title": "NonlinearSolvers.NewtonMatPC",
    "category": "type",
    "text": "Matrix-based Petsc preconditioner for Newton\'s method\n\n\n\n"
},

{
    "location": "NonlinearSolvers/newton.html#NonlinearSolvers.NewtonMatPC-Tuple{ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,ODLCommonTools.AbstractSolutionData,Dict}",
    "page": "Newton\'s Method",
    "title": "NonlinearSolvers.NewtonMatPC",
    "category": "method",
    "text": "Outer constructor for NewtonMatPC\n\n\n\n"
},

{
    "location": "NonlinearSolvers/newton.html#NonlinearSolvers.NewtonVolumePC",
    "page": "Newton\'s Method",
    "title": "NonlinearSolvers.NewtonVolumePC",
    "category": "type",
    "text": "This type holds all the data needed to calculate a preconditioner based only   on the volume integrals.  For DG methods, this matrix is block diagonal and   thus easily invertible.  This preconditioner is applied matrix-free.\n\njac_size is numDofPerNode * numNodesPerElememnt (the total number of unknowns   on an element).\n\nThis preconditioner is not suitable for use as an inner preconditioner, but   the functions calcVolumePC, factorVolumePC, and   applyVolumPC can be easily used to make a new PC.\n\nFields\n\nvolume_jac: jacobian of the volume integrals of each element,              jac_size x jac_size x numEl\nipiv: permutation information, jac_size x numEl\nis_factored: if true, volume_jac has been factored (LU with partial              pivoting), false otherwise\n\n\n\n"
},

{
    "location": "NonlinearSolvers/newton.html#NonlinearSolvers.NewtonVolumePreconditioner-Tuple{ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,ODLCommonTools.AbstractSolutionData,Dict}",
    "page": "Newton\'s Method",
    "title": "NonlinearSolvers.NewtonVolumePreconditioner",
    "category": "method",
    "text": "Regular constructor\n\nInputs\n\nmesh: the mesh\nsbp: the SBP operator\neqn: AbstractSolutionData\nopts: options dictionary\n\n\n\n"
},

{
    "location": "NonlinearSolvers/newton.html#NonlinearSolvers.NewtonDenseLO",
    "page": "Newton\'s Method",
    "title": "NonlinearSolvers.NewtonDenseLO",
    "category": "type",
    "text": "Dense linear operator for Newton\'s method.\n\nSubtype of AbstractDenseLO\n\n\n\n"
},

{
    "location": "NonlinearSolvers/newton.html#NonlinearSolvers.NewtonDenseLO-Tuple{LinearSolvers.PCNone,ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,ODLCommonTools.AbstractSolutionData,Dict}",
    "page": "Newton\'s Method",
    "title": "NonlinearSolvers.NewtonDenseLO",
    "category": "method",
    "text": "Outer constructor for NewtonDenseLO\n\n\n\n"
},

{
    "location": "NonlinearSolvers/newton.html#NonlinearSolvers.NewtonSparseDirectLO",
    "page": "Newton\'s Method",
    "title": "NonlinearSolvers.NewtonSparseDirectLO",
    "category": "type",
    "text": "Sparse direct linear operator for Newton\'s method.  Subtype of   AbstractSparseDirectLO\n\n\n\n"
},

{
    "location": "NonlinearSolvers/newton.html#NonlinearSolvers.NewtonSparseDirectLO-Tuple{LinearSolvers.PCNone,ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,ODLCommonTools.AbstractSolutionData,Dict}",
    "page": "Newton\'s Method",
    "title": "NonlinearSolvers.NewtonSparseDirectLO",
    "category": "method",
    "text": "Outer constructor for NewtonSparseDirectLO\n\n\n\n"
},

{
    "location": "NonlinearSolvers/newton.html#NonlinearSolvers.NewtonPetscMatLO",
    "page": "Newton\'s Method",
    "title": "NonlinearSolvers.NewtonPetscMatLO",
    "category": "type",
    "text": "Petsc matrix based linear operator for Newton\'s method.\n\nSubtype of AbstractPetscMatLO\n\n\n\n"
},

{
    "location": "NonlinearSolvers/newton.html#NonlinearSolvers.NewtonPetscMatLO-Tuple{Union{LinearSolvers.AbstractPetscMatFreePC, LinearSolvers.AbstractPetscMatPC},ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,ODLCommonTools.AbstractSolutionData,Dict}",
    "page": "Newton\'s Method",
    "title": "NonlinearSolvers.NewtonPetscMatLO",
    "category": "method",
    "text": "Outer constructor for NewtonPetscMatLO\n\n\n\n"
},

{
    "location": "NonlinearSolvers/newton.html#NonlinearSolvers.NewtonPetscMatFreeLO",
    "page": "Newton\'s Method",
    "title": "NonlinearSolvers.NewtonPetscMatFreeLO",
    "category": "type",
    "text": "Petsc matrix-free linear operator for Newton\'s method.\n\nSubtype of AbstractPetscMatFreeLO\n\n\n\n"
},

{
    "location": "NonlinearSolvers/newton.html#NonlinearSolvers.NewtonPetscMatFreeLO-Tuple{Union{LinearSolvers.AbstractPetscMatFreePC, LinearSolvers.AbstractPetscMatPC},ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,ODLCommonTools.AbstractSolutionData,Dict}",
    "page": "Newton\'s Method",
    "title": "NonlinearSolvers.NewtonPetscMatFreeLO",
    "category": "method",
    "text": "Newton mat-free linear operator constructor\n\nInputs\n\npc\nmesh\nsbp\neqn\nopts\nrhs_func: rhs_func from newtonInner\n\n\n\n"
},

{
    "location": "NonlinearSolvers/newton.html#Linear-Operators-and-Preconditioners-1",
    "page": "Newton\'s Method",
    "title": "Linear Operators and Preconditioners",
    "category": "section",
    "text": "A full set of linear operators and preconditioners are provided for solving a linear system where the linear operator is the jacobian of the physics..  These often are a good start for constructing linear operators for unsteady problems or more advanced methods for solving steady problems.NewtonMatPC\nNewtonMatPC(::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict)\nNewtonVolumePC\nNewtonVolumePreconditioner(::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict)\nNewtonDenseLO\nNewtonDenseLO(::PCNone, ::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict)\nNewtonSparseDirectLO\nNewtonSparseDirectLO(::PCNone, ::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict)\nNewtonPetscMatLO\nNewtonPetscMatLO(::AbstractPetscPC, ::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict)\nNewtonPetscMatFreeLO\nNewtonPetscMatFreeLO(::AbstractPetscPC, ::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict)"
},

{
    "location": "NonlinearSolvers/jacobian.html#",
    "page": "Jacobian Calculation",
    "title": "Jacobian Calculation",
    "category": "page",
    "text": ""
},

{
    "location": "NonlinearSolvers/jacobian.html#NonlinearSolvers.physicsJac",
    "page": "Jacobian Calculation",
    "title": "NonlinearSolvers.physicsJac",
    "category": "function",
    "text": "This function computes the Jacobian of an evalResidual-like function.   Specifically, it computes \\partial (eqn.q_vec)/ \\partial (eqn.res_vec).\n\nOther users of newtonInner (for example, implicit time marching   methods) will need to implemnt their own version of this function.  The   documentation for this function describes the requirements for the other   implementations.\n\n`newtonInner guarantees that eqn.q and eqn.q_vec will be consistent   when this function is called.  When using finite differences, eqn.res and   eqn.res_vec must contain the residual of the physics.\n\nInputs:\n\nmesh: an AbstractMesh\nsbp: an SBP operator\neqn: an AbstractSolutionData, eqn.res and eqn.res_vec may be overwritten\nopts: options dictonary\njac: the Jacobian, can be an Array, SparseMatrixCSC, or PetscMat\nctx_residual: a tuple of values.  ctx_residual[1] must be an                evalResidual-like function (eqn.q -> eqn.res) function               with signature func(mesh, sbp, eqn, opts, t).               See physicsRhs for a more thorough description                of ctx_residual.\nt: simulation time\n\nOptions Keys:\n\njac_type\njac_method\nepsilon\ncalc_jac_explicit\n\nImplementation Notes:\n\nThis function should not have to do any parallel communication. newtonInner   ensures that the rhs_func is called before jac_func, and rhs_func    handles   the parallel communication.\n\nImplementations of this function may perform either (eqn.q -> jacobian) or   (eqn.q_vec -> jacobian).  The first may be more computationally efficient,   but the second can be simpler for debugging.\n\nThis function supportes several types of jacobians (dense arrays,   SparseMatrixCSC, PetscMat), and several methods for calculating them   (finite difference and complex step).  Any function calling this function   should support them as well.\n\nIt is strongly recommneded to   use this function to compute the spatial jacobian and them modify the   resulting matrix (this function zeros the Jacobian matrix)\n\nWhen using Petsc matrices, the function may do intermediate assemblies   (PETSC_FLUSH_ASSEMBLY), but does not need to do the final assembly.\n\n\n\n"
},

{
    "location": "NonlinearSolvers/jacobian.html#Jacobian-calculation-1",
    "page": "Jacobian Calculation",
    "title": "Jacobian calculation",
    "category": "section",
    "text": "CurrentModule = NonlinearSolversSpecial mention of calcJacobianComplex.One of the most involved parts of Newton\'s method is forming the Jacobian. The NonlinearSolvers module contains the function physicsJac to compute the Jacobian of the physics fracpartial R(q)partial q, where R is evalResidual.  This function can be used by other methods as a starting point for computing the residual of g(R(q)) described in  Newton\'s method page.physicsJacphysicsJac calls one of several functions to compute the Jacobian.  Users should not call these functions directly, they should use physicsJac. There are two approaches to computing the Jacobian: using finite differences/complex step and calculating the entries of the matrix explicitly. The latter approach is much faster than the former, however the former does not require the physics module differentiate evalResidual."
},

{
    "location": "NonlinearSolvers/jacobian.html#NonlinearSolvers.calcJacFD",
    "page": "Jacobian Calculation",
    "title": "NonlinearSolvers.calcJacFD",
    "category": "function",
    "text": "NonlinearSolvers.calcJacFD\n\nThis function calculates the Jacobian using finite differences, perturbing   one degree of freedom at a time.  This is slow and not very accurate.     The Jacobian is calculated about the point in eqn.q_vec.\n\nInputs:     mesh: AbstractMesh     sbp:  SBP operator     eqn:  AbstractEquation object     opts: options dictionary     pert: perturbation to use the finite differences.  Must be of type Tsol.     func: residual evaluation function     res_0: vector containing residual at the point the Jacobian is calculated\n\nInputs/Outputs:     jac:: Jacobian matrix to be populated.  Must be a dense matrix\n\nAliasing restrictions: res_0 must not alias eqn.res_vec\n\nAt the start, calcJacFD assumes:     The Jacobian will be calculated at the state that is specified in eqn.q_vec .     res_0 should have the residual at that state in it\n\nAt exit, eqn.q_vec will have the same values as at the start.\n\neqn.q and eqn.res will be overwritten in the course of this function.\n\n\n\n"
},

{
    "location": "NonlinearSolvers/jacobian.html#NonlinearSolvers.calcJacobianComplex",
    "page": "Jacobian Calculation",
    "title": "NonlinearSolvers.calcJacobianComplex",
    "category": "function",
    "text": "NonlinearSolvers.calcJacComplex\n\nThis function calculates the Jacobian (dense) using the complex step method,    perturbing one degree of freedom at a time.  This is very slow.  The jacobian   is calculated about the point in eqn.q_vec.\n\nInputs:     mesh: AbstractMesh     sbp:  SBP operator     eqn:  AbstractEquation object     opts: options dictionary     pert: perturbation to use.  Must be of type Tsol.     func: residual evaluation function\n\nInputs/Outputs:     jac:: Jacobian matrix to be populated.  Must be a dense matrix\n\nAliasing restrictions: res_0 must not alias eqn.res_vec\n\n\n\n"
},

{
    "location": "NonlinearSolvers/jacobian.html#NonlinearSolvers.calcJacobianSparse",
    "page": "Jacobian Calculation",
    "title": "NonlinearSolvers.calcJacobianSparse",
    "category": "function",
    "text": "NonlinearSolvers.calcJacobianSparse\n\nThis function calculate the Jacobian sparsely (only the entries      within the sparsity bounds), using either finite differences or algorithmic      differentiation.  The jacobian is calculated about the point stored in      eqn.q (not eqn.q_vec).  A mesh coloring approach is used to compute the     jacobian.  Both eqn.q and the MPI send and receive buffers are perturbed     during this process.\n\nInputs:     mesh: AbstractMesh     sbp:  SBP operator     eqn:  AbstractEquation object     opts: options dictionary     pert: perturbation to use for the algorithmic differentiation.  Currently,           only complex numbers are supported.     func: residual evaluation function (eqn.q -> eqn.res)     res_0: element-based (3 dimensional) array containing the residual evaluated           at the point where the Jacobian is being calculated.            This is only used for finite differences (can be a 0 x 0 x 0 array            otherwise).\n\nInputs/Outputs:     jac:  Jacobian matrix.  Must be a sparse matrix type of some kind,            (including PetscMat).\n\nAliasing restrictions: res_0 must not alias eqn.res\n\n\n\n"
},

{
    "location": "NonlinearSolvers/jacobian.html#NonlinearSolvers.applyPerturbation",
    "page": "Jacobian Calculation",
    "title": "NonlinearSolvers.applyPerturbation",
    "category": "function",
    "text": "NonlinearSolvers.applyPerturbation\n\nThis function applies a perturbation to a the specified degree of freedom   on each element according to a mask.\n\nBecause this is element based perturbation, opts[\"parallel_data\"] must   be \"element\".\n\nInputs:     mesh: an AbstractMesh     color: the color to perturb     pert: perturbation to apply.  Can be any datatype     i: local degree of freedom number (in range 1:numDofPerNode) to perturb     j: local node number (in range 1:numNodesPerElement) to perturb\n\nInputs/Outputs:     arr: element based (3D) array of values to perturb     shared_data: array of SharedFaceData for ghost elements to be perturbed                  The receive buffers are perturbed according to the masks                  in mesh.shared_element_colormasks, the send buffers are                  perturbed consistently with arr.\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "NonlinearSolvers/jacobian.html#NonlinearSolvers.assembleElement",
    "page": "Jacobian Calculation",
    "title": "NonlinearSolvers.assembleElement",
    "category": "function",
    "text": "Assembles the volume integral contribution.\n\n\n\nThis function assembles the jacobian contribution for a single element   into the matrix, when the jacobian is computed using finite differences.   This is used by coloring-based methods for computing the jacobian.\n\nInputs:\n\nhelper: a AssembleData object\nmesh:  AbstractMesh object\neqn:  AbstractEquation\nres_arr: element-based (3D) array of perturbed residual values\nres_0:  element-based (3D) array of non-perturbed residual values\nel_res: element number of the element we are observing the change in\nel_pert: element number of the element that was perturbed\ndof_pert: the degree of freedom number of the perturbed dof\nepsilon: magnitude of perturbation\n\nInputs/Outputs:\n\njac: any kind of matrix (dense, SparseMatrixCSC, PetscMat)\n\nAliasing restrictions: res_arr and res_0 must not alias each other.\n\n\n\nSame as other method, but for complex numbers.  See that method for   details.  res_0 is not used in this case\n\n\n\nInputs\n\nhelper: an _AssembleElementData\nmesh: a mesh\nelnum: element number\njac: 4 dimensional array containing the jacobian of the element\n\njac contains the data for the jacobian of the volume terms for a given   element.  jaci j p q = partial Ri p elnum  partial eqnqj q elnum.   Its size is numDofPerNode x numDofPerNode x numNodesPerElement x numNodesPerElement.\n\n\n\n"
},

{
    "location": "NonlinearSolvers/jacobian.html#NonlinearSolvers.calcJacCol",
    "page": "Jacobian Calculation",
    "title": "NonlinearSolvers.calcJacCol",
    "category": "function",
    "text": "NonlinearSolvers.calcJacCol\n\nThis function extracts the entries for one column of the Jacobian from two residual evaluates that come from finite differences.\n\nInputs:     res_0: vector of unperturbed residual values     res: vector of perturbed residual values     epsilon: magnitude of perturbation\n\nInputs/Outputs:     jac_row = vector to be populated with the Jacobian entries\n\nAliasing restrictions: res_0 and res cannot alias (obviously).\n\n\n\nNonlinearSolvers.calcJacCol\n\nThis function extracts the entries for one column of the Jacobian from a    complex step residual evaluation\n\nInputs:     res: vector of perturbed residual values     epsilon: magnitude of perturbation\n\nInputs/Outputs:     jac_row = vector to be populated with the Jacobian entries\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "NonlinearSolvers/jacobian.html#Finite-Differencing/Complex-Step-1",
    "page": "Jacobian Calculation",
    "title": "Finite Differencing/Complex Step",
    "category": "section",
    "text": "The functions calcJacFD and calcJacobianComplex compute the Jacobian as a dense matrix, one column at a time.  This is very slow and only useful for debugging small cases. A more efficient method is implemented in calcJacobianSparse, which uses a graph coloring approach to perturb the solution at several nodes simultaneously. calcJacFD\ncalcJacobianComplex\ncalcJacobianSparse\napplyPerturbation\nassembleElement\ncalcJacCol"
},

{
    "location": "NonlinearSolvers/jacobian.html#NonlinearSolvers._AssembleElementData",
    "page": "Jacobian Calculation",
    "title": "NonlinearSolvers._AssembleElementData",
    "category": "type",
    "text": "Helper object for assembling element and interface jacobians into the   system matrix.\n\nFields\n\nA: the matrix (can be Array, SparseMatrixCSC, or PetscMat)\nidx: temporary array for row indices, length numDofPerNode\nidy: temporary array for column indices, length numDofPerNode\nvals: temporary array for matrix entries, size numDofPerNode square\nidx_i: temporary array for row indices when assembling interface         jacobians, length 2 x numDofPerNode\nidy_i: like idx_i, but for column indices\nvals_i: temporary array for storing matrix entries when assembling         interface jacobians, size 2 x numDofPerNode square\n\n\n\n"
},

{
    "location": "NonlinearSolvers/jacobian.html#NonlinearSolvers._AssembleElementData-Tuple{AbstractArray{T,2} where T,Any,Any,Any,Any}",
    "page": "Jacobian Calculation",
    "title": "NonlinearSolvers._AssembleElementData",
    "category": "method",
    "text": "Outer constructor for _AssembleElementData\n\n\n\n"
},

{
    "location": "NonlinearSolvers/jacobian.html#NonlinearSolvers.NullAssembleElementData",
    "page": "Jacobian Calculation",
    "title": "NonlinearSolvers.NullAssembleElementData",
    "category": "constant",
    "text": "An empty _AssembleElementData.  Useful for giving a default value   to fields.\n\n\n\n"
},

{
    "location": "NonlinearSolvers/jacobian.html#NonlinearSolvers.assembleElement-Tuple{NonlinearSolvers._AssembleElementData,ODLCommonTools.AbstractMesh,Integer,Array{Float64,4}}",
    "page": "Jacobian Calculation",
    "title": "NonlinearSolvers.assembleElement",
    "category": "method",
    "text": "Inputs\n\nhelper: an _AssembleElementData\nmesh: a mesh\nelnum: element number\njac: 4 dimensional array containing the jacobian of the element\n\njac contains the data for the jacobian of the volume terms for a given   element.  jaci j p q = partial Ri p elnum  partial eqnqj q elnum.   Its size is numDofPerNode x numDofPerNode x numNodesPerElement x numNodesPerElement.\n\n\n\n"
},

{
    "location": "NonlinearSolvers/jacobian.html#NonlinearSolvers.assembleInterface",
    "page": "Jacobian Calculation",
    "title": "NonlinearSolvers.assembleInterface",
    "category": "function",
    "text": "Assembles contributions into an element-block matrix.  jacLR and jacRL   are not used.\n\n\n\nAssembles the jacobian of an interface into the matrix.  Specialized   versions take advantage of the sparsity of the sbpface.\n\nInputs\n\nhelper: _AssembleElementData\nsbpface: an SBP face object\nmesh: a mesh\niface: an Interface object identify the interface to be assembled\njacLL: see below\njacLR:\njacRL\njacRR\n\njacAB where A = L or R and B = L or R, is the jacobian of the residual of   element A with respect to the solution of element B.\n\njacAB has the same size/layout as jac in assembleElement.\n\n\n\n"
},

{
    "location": "NonlinearSolvers/jacobian.html#NonlinearSolvers.assembleSharedFace",
    "page": "Jacobian Calculation",
    "title": "NonlinearSolvers.assembleSharedFace",
    "category": "function",
    "text": "Assembles contribution of the shared face terms into the element-block   matrix.  jacLR is not used.\n\n\n\nAssemble one half of an interface, used by shared face integrals.   See assembleInterface.\n\nInputs\n\nhelper: _AssembleElementData\nsbpface: an SBP face object\nmesh: a mesh\njacLL\njacLR\n\n\n\n"
},

{
    "location": "NonlinearSolvers/jacobian.html#NonlinearSolvers.assembleBoundary",
    "page": "Jacobian Calculation",
    "title": "NonlinearSolvers.assembleBoundary",
    "category": "function",
    "text": "Assembles the boundary terms into an element-block matrix.\n\n\n\nAssembles the jacobian of a boundary integral into the matrix.   Specialized versions take advantage of the sparsity of the sbpface.\n\nInputs\n\nhelper: _AssembleElementData\nsbpface: an SBP face object\nmesh: a mesh\njac: a jac, same layout as assembleElement\n\n\n\n"
},

{
    "location": "NonlinearSolvers/jacobian.html#Explicit-Calculation-1",
    "page": "Jacobian Calculation",
    "title": "Explicit Calculation",
    "category": "section",
    "text": "The most efficient method to compute the Jacobian is to explicitly compute all its entries and assemble them into a matrix. The evalJacobian function must be extended by each physics module for this to work. This function is passed all the same arguments as evalResidual plus an additional _AssembleElementData object which is used to assemble the contribution of each element or interface into the matrix._AssembleElementData\n_AssembleElementData(::AbstractMatrix, ::Any, ::Any, ::Any, ::Any)\nNullAssembleElementData\nassembleElement(::_AssembleElementData, ::AbstractMesh, ::Integer, ::Array{Float64, 4})\nassembleInterface\nassembleSharedFace\nassembleBoundary"
},

{
    "location": "NonlinearSolvers/jac_recalc.html#",
    "page": "Jacobian Freezing",
    "title": "Jacobian Freezing",
    "category": "page",
    "text": ""
},

{
    "location": "NonlinearSolvers/jac_recalc.html#Jacobian-and-Preconditioner-Freezing-1",
    "page": "Jacobian Freezing",
    "title": "Jacobian and Preconditioner Freezing",
    "category": "section",
    "text": "  CurrentModule = NonlinearSolversIt is often faster to reuse the Jacobian and/or preconditioner than calculate a new one each major iteration.  The functions described here provide a consistant API for doing so.Each method in NonlinearSolvers should specify if it supports Jacobian and PC freezing.  If it does, it should specify the prefix (see Construction).Note that some methods use newtonInner internally.  It is sometimes beneficial to have newtonInner recompute the Jacobian and PC.  In other cases, using the same Jacobian and PC for several calls to newtonInner is beneficial.  Both these use cases are supported. Both newtonInner and the outer method can have their own RecalculationPolicy objects. When the newtonInner RecalculationPolicy is RecalculateNever, it will never recalculate the Jacobian and PC and the outer method will be the soley responsible for updaing the Jacobian and PC.  In the reverse case, the outer method can use RecalculateNever and let newtonInner recalculate the Jacobian and PC."
},

{
    "location": "NonlinearSolvers/jac_recalc.html#NonlinearSolvers.doRecalculation",
    "page": "Jacobian Freezing",
    "title": "NonlinearSolvers.doRecalculation",
    "category": "function",
    "text": "This function uses decideRecalculation and the LinearSolver   interface to recalculate the PC and/or LO as specified by the policy.\n\nThis is the main interface the different methods should use for   recalculating the PC and/or LO.  The only case when   decideRecalculation is better is when the method is not using   the LinearSolvers module to define the PC and LO.\n\nInputs\n\npolicy: a RecalculationPolicy\nitr: current iteration number\n\nThe following arguments exactly match the signature of calcPCandLO\n\nls\nmesh\nsbp\neqn\nopts\nctx_residual\nt\n\nWhen this function exits, the PC and/or LO will have been updated, if   specified by the RecalculationPolicy.\n\n\n\n"
},

{
    "location": "NonlinearSolvers/jac_recalc.html#NonlinearSolvers.decideRecalculation",
    "page": "Jacobian Freezing",
    "title": "NonlinearSolvers.decideRecalculation",
    "category": "function",
    "text": "This function returns the enum specifying whether the PC and/or LO should   be recalculated.\n\nEvery implementation of RecalculationPolicy should extend this   function with a new method.\n\nInputs\n\npolicy: a RecalculationPolicy\nitr: the current iteration\n\nOutputs\n\nthe num: see RECALC_BOTH\n\n\n\n"
},

{
    "location": "NonlinearSolvers/jac_recalc.html#NonlinearSolvers.resetRecalculationPolicy",
    "page": "Jacobian Freezing",
    "title": "NonlinearSolvers.resetRecalculationPolicy",
    "category": "function",
    "text": "This function should reset the RecalculationPolicy to its initial state.   This function is called at the beginning of Newton\'s method (making it   safe to reuse the RecalculationPolicy for repeated calls to Newton)\n\nEvery RecalculationPolicy should extend function with a new   method.\n\nInputs\n\npolicy: a RecalculationPolicy\n\nOutputs\n\nnone\n\n\n\n"
},

{
    "location": "NonlinearSolvers/jac_recalc.html#NonlinearSolvers.RECALC_BOTH",
    "page": "Jacobian Freezing",
    "title": "NonlinearSolvers.RECALC_BOTH",
    "category": "constant",
    "text": "Enums returned by decideRecalculation telling the caller what to   recalculate.\n\nEnum Names\n\nRECALC_BOTH: recalculate PC and LO\nRECALC_PC: recalculate PC only\nRECALC_LO: recalculate LO only\nRECALC_NONE: don\'t recalculate anything\n\n\n\n"
},

{
    "location": "NonlinearSolvers/jac_recalc.html#API-1",
    "page": "Jacobian Freezing",
    "title": "API",
    "category": "section",
    "text": "doRecalculation\ndecideRecalculation\nresetRecalculationPolicy\nRECALC_BOTH"
},

{
    "location": "NonlinearSolvers/jac_recalc.html#NonlinearSolvers.getRecalculationPolicy",
    "page": "Jacobian Freezing",
    "title": "NonlinearSolvers.getRecalculationPolicy",
    "category": "function",
    "text": "This function constructs and returns the recalculation policy specified by   the options dictionary.\n\nMethods that want to use a RecalculationPolicy should this this function   to get it (do not call the constructors directly).\n\nInputs\n\nopts: options dictionary\nprefix: prefix of options to use.  This is usually the name of the method         the returned object will be used in, eg. \"newton\"\n\nOutputs\n\npolicy: a RecalculationPolicy of some kind\n\nOptions Keys\n\n\"prefix_recalculation_policy\": name of recalculation policy to get\n\nCurrent Policies\n\nString[\"RecalculateFixedIntervals\", \"RecalculateNever\"]\n\n\n\n"
},

{
    "location": "NonlinearSolvers/jac_recalc.html#NonlinearSolvers.RecalculationPolicyDict",
    "page": "Jacobian Freezing",
    "title": "NonlinearSolvers.RecalculationPolicyDict",
    "category": "constant",
    "text": "Maps names to RecalculationPolicy constructors.  All new   recalculation policies must be added to this list.\n\nThe constructors must have the signature:\n\nMyPolicyName(opts::Dict, prefix::String)\n\nwhere opts is the options dictionary and prefix is the prefix passed into   getRecalculationPolicy.\n\n\n\n"
},

{
    "location": "NonlinearSolvers/jac_recalc.html#Construction-1",
    "page": "Jacobian Freezing",
    "title": "Construction",
    "category": "section",
    "text": "getRecalculationPolicy\nRecalculationPolicyDict"
},

{
    "location": "NonlinearSolvers/jac_recalc.html#NonlinearSolvers.RecalculateFixedIntervals",
    "page": "Jacobian Freezing",
    "title": "NonlinearSolvers.RecalculateFixedIntervals",
    "category": "type",
    "text": "This RecalculationPolicy recalculates the LO and PC every x    number if iterations.\n\nOptions Keys\n\n\"prefix_prec_recalc_freq\": frequency of PC recalculation\n\"prefix_jac_recalc_freq\": frequency of LO recalculation\n\"prefix_recalc_first\": if true, recalculate the PC and LO on the first                          iteration, even if other criteria not met.\n\nwhere \"prefix\" is the prefix passed into getRecalculationPolicy\n\n\n\n"
},

{
    "location": "NonlinearSolvers/jac_recalc.html#NonlinearSolvers.RecalculateNever",
    "page": "Jacobian Freezing",
    "title": "NonlinearSolvers.RecalculateNever",
    "category": "type",
    "text": "This RecalculationPolicy never recalculates the PC or LO.  This is useful when Newton\'s method is used inside other methods to solve  easy problems.\n\n\n\n"
},

{
    "location": "NonlinearSolvers/jac_recalc.html#NonlinearSolvers.RecalculationPolicy",
    "page": "Jacobian Freezing",
    "title": "NonlinearSolvers.RecalculationPolicy",
    "category": "type",
    "text": "This is the supertype of all types that tell a method when to recalculate   the Jacobian.  Each subtype defines a different policy for when to   recalculate.\n\n\n\n"
},

{
    "location": "NonlinearSolvers/jac_recalc.html#Recalculation-Policies-1",
    "page": "Jacobian Freezing",
    "title": "Recalculation Policies",
    "category": "section",
    "text": "Modules = [NonlinearSolvers]\nPages = [\"jac_recalc.jl\"]\nOrder = [:type]"
},

{
    "location": "NonlinearSolvers/residual_evaluation.html#",
    "page": "Residual Evalution",
    "title": "Residual Evalution",
    "category": "page",
    "text": ""
},

{
    "location": "NonlinearSolvers/residual_evaluation.html#NonlinearSolvers.physicsRhs",
    "page": "Residual Evalution",
    "title": "NonlinearSolvers.physicsRhs",
    "category": "function",
    "text": "This function computes the vector form of of the residual from the vector   form of the solution, ie. q_vec -> rhs_vec, for a given physics.   This is one of the two functions required by newtonInner.\n\nUsers of newtonInner that are not physics modules (for example, implicit   time marching schemes) will need to implement their own version of this   function.  See the Implimentation Notes section below.\n\nInputs\n\nmesh: an AbstractMesh\nsbp: an AbstractSBP\neqn: an AbstractSolutionData (may be modified during this function)\nopts: the options dictionary\nctx_residual: a tuple of values.  ctx_residual[1] must be a function               that computes (q -> res).  Typically this is evalResidual.               The other entries of the tuple (if any) are not used.\n          The purpose of this argument is to make the signature of\n          the function generic enough so that implict time marching\n          methods can use it.  (If you are confused about this\n          programming pattern, google how callback are implemented\n          in C, this ctx is like a void* in C).\nt: the time at which to evalute the residual\n\nInputs/Outputs\n\nrhs_vec: vector to put the residual in\n\nOutput\n\nnorm of the residual vector\n\nImplementation Notes:\n\nThis function is really a wrapper around an evalResidual-like function.   It has to do 5 things:\n\n1. scatter eqn.q_vec -> eqn.q\n2. start parallel communication if needed\n3. call the evalResidual-like function to compute eqn.q -> eqn.res\n4. assemble the residual into the output vector, ie. eqn.res -> rhs_vec\n5. compute the norm of the rhs_vec\n\nAny implementation of this function (used with newtonInner) must have the   following properties:\n\n1. on exit, eqn.q and eqn.q_vec must be consistent\n2. this function start parallel communication if needed\n3. allow for the possibility that rhs_vec and eqn.res_vec alias\n\nIs is recommended to use calcNorm to compute the norm.\n\nOther implementations of this function are encouraged to use this function   to help construct their rhs_vec.\n\n\n\n"
},

{
    "location": "NonlinearSolvers/residual_evaluation.html#NonlinearSolvers.assembleResidual",
    "page": "Residual Evalution",
    "title": "NonlinearSolvers.assembleResidual",
    "category": "function",
    "text": "NonlinearSolvers.assembleResidual\n\nThis function takes the residual in eqn.res and performs an additive reduction   into res_vec.   This function wraps array3DTo1D.\n\nInputs:     mesh: AbstractMesh     sbp:  SBP operator     eqn:  AbstractEquation object     opts: options dictionary     res_vec: residual vector to put the residual into\n\nOutputs:     none\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "NonlinearSolvers/residual_evaluation.html#Utils.array1DTo3D-Tuple{Any,Any,Any,Any,AbstractArray{Float64,1}}",
    "page": "Residual Evalution",
    "title": "Utils.array1DTo3D",
    "category": "method",
    "text": "NonlinearSolvers.array1DTo3D\n\nThis function performs the scatter q_vec -> eqn.q\n\nInputs:     mesh: AbstractMesh     sbp:  SBP operator     eqn:  AbstractEquation object     opts: options dictionary     q_vec: vector containing solution variables \n\nOutputs:     none\n\nAliasing Restrictions: none\n\n\n\n"
},

{
    "location": "NonlinearSolvers/residual_evaluation.html#Residual-Evaluation-1",
    "page": "Residual Evalution",
    "title": "Residual Evaluation",
    "category": "section",
    "text": "CurrentModule = NonlinearSolversThis page describes some helper functions used by various methods in the NonlinearSolvers module.physicsRhs\nassembleResidual\narray1DTo3D(::Any, ::Any, ::Any, ::Any, ::AbstractArray{Float64, 1})"
},

{
    "location": "NonlinearSolvers/matrix.html#",
    "page": "Matrix Interface",
    "title": "Matrix Interface",
    "category": "page",
    "text": ""
},

{
    "location": "NonlinearSolvers/matrix.html#Matrix-Interface-1",
    "page": "Matrix Interface",
    "title": "Matrix Interface",
    "category": "section",
    "text": "One of the challenged of NonlinearSolvers is supporting several matrix type. The currently supported types areArray (a Julia dense matrix, column major)\nSparseMatrixCSC (a Julia implementation of compressed sparse column format)\nPetscMat (Julia wrappers of the PETSc library)Unfortunatly, these matrix types support somewhat different APIs that have little overlap. This is a problem, because it becomes difficult to write efficient,generic code, ie. code that works on all three  matrix types (without writing specialized code for each one) and performs nearly as well as an implementation specialized for a single implementation.To solve this problem, the Petsc wrappers implement an API that is much smaller than that of Array, but efficiently maps onto the capabilities of the different array types. This API should be used in NonlinearSolvers to ensure all functions support the different matrix types"
},

{
    "location": "NonlinearSolvers/matrix.html#API-1",
    "page": "Matrix Interface",
    "title": "API",
    "category": "section",
    "text": "TODO: when the Petsc docs are finished, document them here"
},

{
    "location": "NonlinearSolvers/newton_inner.html#",
    "page": "Newton Inner",
    "title": "Newton Inner",
    "category": "page",
    "text": ""
},

{
    "location": "NonlinearSolvers/newton_inner.html#Newton_inner-1",
    "page": "Newton Inner",
    "title": "Newton_inner",
    "category": "section",
    "text": "  CurrentModule = NonlinearSolvers"
},

{
    "location": "NonlinearSolvers/newton_inner.html#NewtonData-type-1",
    "page": "Newton Inner",
    "title": "NewtonData type",
    "category": "section",
    "text": ""
},

{
    "location": "NonlinearSolvers/newton_inner.html#setupNewton-1",
    "page": "Newton Inner",
    "title": "setupNewton",
    "category": "section",
    "text": ""
},

{
    "location": "NonlinearSolvers/newton_inner.html#newtonInner-1",
    "page": "Newton Inner",
    "title": "newtonInner",
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
    "text": "This module contains functions and types that are useful for the solver but independent of the equation being solved. Additional utility functions are located in the ODLCommontools. The functions defined in the Utils module are useful in the context of PDESolver and depend on the functions and datatypes defined in the other parts of the solver. The functions defined in ODLCommonTools are more general in nature and usable independent of PDESolver.  Pages = [\"io.md\"\n           \"logging.md\"\n           \"projections.md\"\n           \"parallel.md\"\n           \"misc.md\"\n          ]\n  Depth = 1"
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
    "category": "type",
    "text": "This type holds all the data necessary to perform MPI communication with   a given peer process that shared mesh edges (2D) or faces (3D) with the   current process.\n\nFields:\n\npeernum: the MPI rank of the peer process\npeeridx: the index of this peer in mesh.peer_parts\nmyrank: MPI rank of the current process\ncomm: MPI communicator used to define the above\n\nq_send: the send buffer, a 3D array of n x m x d.  While these dimensions\n        are arbitrary, there are two commonly used case.  If\n        opts[\"parallel_type\"] == face, then m is mesh.numNodesPerFace and\n        d is the number of faces shared with peernum.\n        If opts[\"parallel_type\"] == element, then \n        m = mesh.numNodesPerElement and d is the number of elements that\n        share faces with peernum.\nq_recv: the receive buffer.  Similar to q_send, except the size needs to\n        to be the number of entities on the *remote* process.\n\nsend_waited: has someone called MPI.Wait() on send_req yet?  Some MPI\n             implementations complain if Wait() is called on a Request\n             more than once, so use this field to avoid doing so.\nsend_req: the MPI.Request object for the Send/Isend/whatever other type of\n          Send\nsend_status: the MPI.Status object returned by calling Wait() on send_req\n\nrecv_waited: like send_waited, but for the receive\nrecv_req: like send_req, but for the receive\nrecv_status: like send_status, but for the receive\n\nbndries_local: Vector of Boundaries describing the faces from the local\n               side of the interface\nbndries_remote: Vector of Boundaries describing the facaes from the remote\n                side (see the documentation for PdePumiInterface before\n                using this field)\ninterfaces: Vector of Interfaces describing the faces from both sides (see\n            the documentation for PdePumiInterfaces, particularly the\n            mesh.shared_interfaces field, before using this field\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Utils.SharedFaceData-Union{Tuple{ODLCommonTools.AbstractMesh,Int64,Array{T,3},Array{T,3}}, Tuple{T}} where T",
    "page": "Parallel Constructs",
    "title": "Utils.SharedFaceData",
    "category": "method",
    "text": "Outer constructor for SharedFaceData.\n\nInputs:\n\nmesh: a mesh object\npeeridx: the index of a peer in mesh.peer_parts\nq_send: the send buffer\nq_recv: the receive buffer\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Utils.getSharedFaceData-Union{Tuple{Tsol}, Tuple{Type{Tsol},ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,Any}} where Tsol",
    "page": "Parallel Constructs",
    "title": "Utils.getSharedFaceData",
    "category": "method",
    "text": "This function returns a vector of SharedFaceData objects, one for each   peer processes the current process shares mesh edge (2D) or face (3D) with.   This function is intended to be used by the AbstractSolutionData constructors,   although it can be used to create additional vectors of SharedFaceData   objects.\n\nif opts[\"parallel_data\"] == \"face\", then the send and receive buffers   are numDofPerNode x numNodesPerFace x number of shared faces.\n\nif opts[\"parallel_data\"] == \"element\", the send and receive buffers are     numDofPerNode x numNodesPerElement x number of elements that share the     faces.  Note that the number of elements that share the faces can be     different for the send and receive buffers.\n\nInputs:\n\nTsol: element type of the arrays\nmesh: an AbstractMesh object\nsbp: an SBP operator\nopts: the options dictonary\n\nOutputs:\n\ndata_vec: Vector{SharedFaceData}.  data_vec[i] corresponds to \n          mesh.peer_parts[i]\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Base.copy!-Union{Tuple{T}, Tuple{Utils.SharedFaceData{T},Utils.SharedFaceData{T}}} where T",
    "page": "Parallel Constructs",
    "title": "Base.copy!",
    "category": "method",
    "text": "In-place copy for SharedFaceData.  This copies the buffers, but does not   retain the state of the Request and Status fields.  Instead they are   initialized the same as the constructor.\n\nThis function may only be called after receiving is complete,   otherwise an exception is thrown.\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Base.copy-Union{Tuple{T}, Tuple{Utils.SharedFaceData{T}}} where T",
    "page": "Parallel Constructs",
    "title": "Base.copy",
    "category": "method",
    "text": "Copy function for SharedFaceData.  Note that this does not retain the   send_req/send_status (and similarly for the recceive) state   of the original object.  Instead, they are initialized the same as the   constructor.\n\nThis function may only be called after receiving is complete,   otherwise an exception is thrown.\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Utils.assertReceivesConsistent-Union{Tuple{Array{Utils.SharedFaceData{T},1}}, Tuple{T}} where T",
    "page": "Parallel Constructs",
    "title": "Utils.assertReceivesConsistent",
    "category": "method",
    "text": "Like assertSendsConsistent, but for the receives\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Utils.assertReceivesWaited-Union{Tuple{Array{Utils.SharedFaceData{T},1}}, Tuple{T}} where T",
    "page": "Parallel Constructs",
    "title": "Utils.assertReceivesWaited",
    "category": "method",
    "text": "This function verifies all the receives have been waited on for the    supplied SharedFaceData objects\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Utils.assertSendsConsistent-Union{Tuple{Array{Utils.SharedFaceData{T},1}}, Tuple{T}} where T",
    "page": "Parallel Constructs",
    "title": "Utils.assertSendsConsistent",
    "category": "method",
    "text": "Verify either all or none of the sends have been waited on.  Throw an   exception otherwise.\n\nInputs:\n\nshared_data: Vector of SharedFaceData objects\n\nOutput:\n\nval: number of receives that have been waited on\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Utils.waitAllReceives-Union{Tuple{Array{Utils.SharedFaceData{T},1}}, Tuple{T}} where T",
    "page": "Parallel Constructs",
    "title": "Utils.waitAllReceives",
    "category": "method",
    "text": "This function is like MPI.Waitall, operating on the recvs of a vector of    SharedFaceData objects\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Utils.waitAllSends-Union{Tuple{Array{Utils.SharedFaceData{T},1}}, Tuple{T}} where T",
    "page": "Parallel Constructs",
    "title": "Utils.waitAllSends",
    "category": "method",
    "text": "This function is like MPI.Waitall, operating on the sends of a vector of    SharedFaceData objects\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Utils.waitAnyReceive-Union{Tuple{Array{Utils.SharedFaceData{T},1}}, Tuple{T}} where T",
    "page": "Parallel Constructs",
    "title": "Utils.waitAnyReceive",
    "category": "method",
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
    "category": "function",
    "text": "This function is a thin wrapper around exchangeData().  It is used for the   common case of sending and receiving the solution variables to other processes.   It uses eqn.shared_data to do the parallel communication.   eqn.shared_data must be passed into the corresponding finishDataExchange   call.\n\nInputs:     mesh: an AbstractMesh     sbp: an SBP operator     eqn: an AbstractSolutionData     opts: options dictionary\n\nKeyword arguments:     tag: MPI tag to use for communication, defaults to TAG_DEFAULT     wait: wait for sends and receives to finish before exiting\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Utils.exchangeData",
    "page": "Parallel Constructs",
    "title": "Utils.exchangeData",
    "category": "function",
    "text": "This function posts the MPI sends and receives for a vector of SharedFaceData.  It works for both opts[\"parallel_data\"] == \"face\" or \"element\".  The only   difference between these two cases is the populate_buffer() function.\n\nThe previous receives using these SharedFaceData objects should have   completed by the time this function is called.  An exception is throw   if this is not the case.\n\nThe previous sends are likely to have completed by the time this function   is called, but they are waited on if not.  This function might not perform   well if the previous sends have not completed.   #TODO: fix this using WaitAny\n\nInputs:     mesh: an AbstractMesh     sbp: an SBPOperator     eqn: an AbstractSolutionData     opts: the options dictionary     populate_buffer: function with the signature:                      populate_buffer(mesh, sbp, eqn, opts, data::SharedFaceData)                      that populates data.q_send   Inputs/Outputs:     shared_data: vector of SharedFaceData objects representing the parallel                  communication to be done\n\nKeyword Arguments:     tag: MPI tag to use for this communication, defaults to TAG_DEFAULT          This tag is typically used by the communication of the solution          variables to other processes.  Other users of this function should          provide their own tag\n\nwait: wait for the sends and receives to finish before returning.  This\n      is a debugging option only.  It will kill parallel performance.\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Utils.finishExchangeData",
    "page": "Parallel Constructs",
    "title": "Utils.finishExchangeData",
    "category": "function",
    "text": "This is the counterpart of exchangeData.  This function finishes the   receives started in exchangeData.\n\nThis function (efficiently) waits for a receive to finish and calls   a function to do calculations for on that data. If opts[\"parallel_data\"]   == \"face\", it also permutes the data in the receive buffers to agree   with the ordering of elementL.  For opts[\"parallel_data\"] == \"element\",   users should call SummationByParts.interiorFaceInterpolate to interpolate   the data to the face while ensuring proper permutation.\n\nInputs:     mesh: an AbstractMesh     sbp: an SBPOperator     eqn: an AbstractSolutionData     opts: the options dictonary     calc_func: function that does calculations for a set of shared faces                described by a single SharedFaceData.  It must have the signature                calc_func(mesh, sbp, eqn, opts, data::SharedFaceData)\n\nInputs/Outputs:     shared_data: vector of SharedFaceData, one for each peer process that                  needs to be communicated with.  By the time calc_func is                  called, the SharedFaceData passed to it has its q_recv field                  populated.  See note above about data permutation.\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Utils.TAG_DEFAULT",
    "page": "Parallel Constructs",
    "title": "Utils.TAG_DEFAULT",
    "category": "constant",
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
    "category": "function",
    "text": "Utils.verifyCommunication\n\nThis function checks the data provided by the Status object to verify a    communication completed successfully.  The sender\'s rank and the number of   elements is checked agains the expected sender and the buffer size\n\nInputs:     data: a SharedFaceData\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Utils.getSendDataFace",
    "page": "Parallel Constructs",
    "title": "Utils.getSendDataFace",
    "category": "function",
    "text": "This function populates the send buffer from eqn.q for    opts[\"parallle_data\"]  == \"face\"\n\nInputs:     mesh: a mesh     sbp: an SBP operator     eqn: an AbstractSolutionData     opts: options dictonary\n\nInputs/Outputs:     data: a SharedFaceData.  data.q_send will be overwritten\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Utils.getSendDataElement",
    "page": "Parallel Constructs",
    "title": "Utils.getSendDataElement",
    "category": "function",
    "text": "This function populates the send buffer from eqn.q for    opts[\"parallle_data\"]  == \"element\"\n\nInputs:\n\nmesh: a mesh\nsbp: an SBP operator\neqn: an AbstractSolutionData\nopts: options dictonary\n\nInputs/Outputs:\n\ndata: a SharedFaceData.  data.q_send will be overwritten\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Utils.@mpi_master",
    "page": "Parallel Constructs",
    "title": "Utils.@mpi_master",
    "category": "macro",
    "text": "Utils.mpi_master\n\nThis macro introduces an if statement that causes the expression to be    executed only if the variable myrank is equal to zero.  myrank must exist   in the scope of the caller\n\n\n\n"
},

{
    "location": "Utils/parallel.html#Utils.@time_all",
    "page": "Parallel Constructs",
    "title": "Utils.@time_all",
    "category": "macro",
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
    "location": "Utils/projections.html#Utils.calcLength-Tuple{ODLCommonTools.AbstractParamType{2},AbstractArray{T,1} where T}",
    "page": "Projections",
    "title": "Utils.calcLength",
    "category": "method",
    "text": "Calculates the length of a vector, using a manually unrolled loop.   Methods are available for 2D and 3D.\n\nInputs:\n\nparams: an AbstractParamType{Tdim}, used to dispatch to the right method\nnrm: the vector to calculate the length of.  Must have length 2 in 2D\n     and 3 in 3D\n\nOutputs:\n\nlength of nrm\n\n\n\n"
},

{
    "location": "Utils/projections.html#Utils.getProjectionMatrix-Tuple{ODLCommonTools.AbstractParamType{2},AbstractArray{T,1} where T,AbstractArray{T,2} where T}",
    "page": "Projections",
    "title": "Utils.getProjectionMatrix",
    "category": "method",
    "text": "This function populates the matrix P that project from x-y-z to n-t-b.   Methods are available for 2D and 3D, determined from the AbstractParamType   object.  This is a somewhat specialized routine in that it only works for   vectors like the state vector of the Euler or Navier-Stokes equations,   where the first and last equations are coordinate system invarient and   only the middle 2 (in 2D) or 3 (in 3D) equations need to be rotated.   Specialized multiplication routines are provied as projectToXY and   projectToNT.\n\nInputs:\n\nparams:  An AbstractParamType{Tdim}\nnrm: the normal direction in x-y coordinate system.  Must be a unit vector\n\nInputs/Outputs\n\nP:  matrix to be populated with projection matrix\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "Utils/projections.html#Utils.projectToNT-Tuple{ODLCommonTools.AbstractParamType{2},AbstractArray{T,2} where T,AbstractArray{T,1} where T,AbstractArray{T,1} where T}",
    "page": "Projections",
    "title": "Utils.projectToNT",
    "category": "method",
    "text": "This function projects a vector x from x-y coordinates to normal-tangential   (n-t) coordinates, where A is a projection matrix from x-y to n-t, obtained   from getProjectionMatrix.  This function is a specialized matrix-vector   multiplication routine and only works for matrices with the particular   sparsity pattern created by getProjectionMatrix.   Methods are available for 2D and 3D\n\nInputs:\n\nparams: an AbstractParamType{{Tdim} used to dispatch to the 2D or 3D\n        method\nP: the projection matrix\nx: vector to be projected\n\nInputs/Outputs:\n\nb: the result of the projection\n\nAliasing restrictions: x and b cannot alias\n\n\n\n"
},

{
    "location": "Utils/projections.html#Utils.projectToXY-Tuple{ODLCommonTools.AbstractParamType{2},AbstractArray{T,2} where T,AbstractArray{T,1} where T,AbstractArray{T,1} where T}",
    "page": "Projections",
    "title": "Utils.projectToXY",
    "category": "method",
    "text": "This function is similar to projectToNT, except it project from n-t   to x-y.  Note that P is still the projection matrix from getProjectionMatrix   that projects from x-y to n-t.\n\n\n\n"
},

{
    "location": "Utils/projections.html#Utils.getBinormalVector-NTuple{6,Number}",
    "page": "Projections",
    "title": "Utils.getBinormalVector",
    "category": "method",
    "text": "This function computes a vector normal to the 2 supplied vectors and   returns the components.\n\nInputs:\n\nn1, n2, n3: the components of the first vector\nt1, t2, t3: the components of the second vector\n\nOutputs:\n\nb1, b2, b3: the components of the binormal vector\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "Utils/projections.html#Utils.getOrthogonalVector-Tuple{ODLCommonTools.AbstractParamType{2},AbstractArray{T,1} where T}",
    "page": "Projections",
    "title": "Utils.getOrthogonalVector",
    "category": "method",
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
    "location": "Utils/logging.html#Utils.sharedFaceLogging-Union{Tuple{Any,Any,ODLCommonTools.AbstractSolutionData{Tsol,Tres} where Tres,Any,Utils.SharedFaceData,Any,Any}, Tuple{Tsol}} where Tsol",
    "page": "Logging",
    "title": "Utils.sharedFaceLogging",
    "category": "method",
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
    "category": "constant",
    "text": "Buffered version of STDERR.  This should always be used instead of STDERR\n\n\n\n"
},

{
    "location": "Utils/io.html#Utils.BSTDOUT",
    "page": "Input/Output",
    "title": "Utils.BSTDOUT",
    "category": "constant",
    "text": "Buffered version of STDOUT.  This should always be used instead of STDOUT\n\n\n\n"
},

{
    "location": "Utils/io.html#Utils.BufferedIO",
    "page": "Input/Output",
    "title": "Utils.BufferedIO",
    "category": "type",
    "text": "Utils.BufferedIO\n\nThis type provides a means to buffer IO in memory before writing it to a file.   Data written to the object is stored in the IOBuffer until flush() is called,    when the buffer is (efficiently) dumped into the file\n\n\n\n"
},

{
    "location": "Utils/io.html#Utils.BufferedIO",
    "page": "Input/Output",
    "title": "Utils.BufferedIO",
    "category": "type",
    "text": "Alternative constructor for BufferedIO, emulating the open() function.   This function creates the underlying file using open() and then creates   a BufferedIO around it.\n\nInputs:\n\nfname: AbstractString, name of file to open\nmode: file open mode, see documentation of open(), defaults to append\n\nOutputs:\n\na BufferedIO object\n\n\n\n"
},

{
    "location": "Utils/io.html#Utils.BufferedIO",
    "page": "Input/Output",
    "title": "Utils.BufferedIO",
    "category": "type",
    "text": "Utils.BufferedIO\n\nConstructor for BufferedIO type.  Takes an IOStream and creates an IOBuffer.   If no IOStream is given, a dummy stream is created.  If a dummy stream is   used, flush must never be called on the returned BufferedIO object\n\nInputs:\n\nf: an IOStream object, defaults to a dummy stream\n\nOutputs:\n\na BufferedIO object\n\n\n\n"
},

{
    "location": "Utils/io.html#Base.close-Tuple{Utils.BufferedIO}",
    "page": "Input/Output",
    "title": "Base.close",
    "category": "method",
    "text": "Base function close extended for BufferedIO\n\n\n\n"
},

{
    "location": "Utils/io.html#Base.flush-Tuple{Utils.BufferedIO}",
    "page": "Input/Output",
    "title": "Base.flush",
    "category": "method",
    "text": "Base function flush extended for BufferedIO\n\n\n\n"
},

{
    "location": "Utils/io.html#Detailed-Documentation-1",
    "page": "Input/Output",
    "title": "Detailed Documentation",
    "category": "section",
    "text": "  Modules = [Utils]\n  Pages = [\"Utils/io.jl\"]"
},

{
    "location": "Utils/checkpoint.html#",
    "page": "Checkpointing",
    "title": "Checkpointing",
    "category": "page",
    "text": ""
},

{
    "location": "Utils/checkpoint.html#Checkpointing-1",
    "page": "Checkpointing",
    "title": "Checkpointing",
    "category": "section",
    "text": "  CurrentModule = UtilsCheckpointing has two use cases: saving the state of the solver to be loaded later (for example, in unsteady adjoint calculations), and to restart the solver after a crash. This checkpointing functionality described here is useful for both of these functions. Its purpose is to provide an interface for writing the current state to a file that can be read back later. As long as at least 2 checkpoints are saved, the implementation guarantees that at least one checkpoint is loadable at any time, even if the code is terminated while writing a checkpoint.Log files require special handling when restarting. Log files should be opened in append mode when restarting (in fact, the code usually opens log files in append mode even when not restarting). One of the consequences of restarting is that there may be duplicate lines in log files.  For example, if a code is doing a run of 10,000 timesteps, checkpointing every 1,000 timesteps, and gets terminated at timestep 6,500, it can be restarted from the checkpoint at timestep 6,000 and run to the end.  This means any output written between timesteps 6,000 and 6,500 will  appear twice in the log files. To deal with this problem, the script log_cleaner.jl (located in src/scripts) can be used to post-process delimited data files (space, comma, tab, etc.) and remove duplicate entries in favor of the entry appearing last in the log file. It also checks removed entries to make sure they are the same (to floating point tolerance) as the entries that are preserved."
},

{
    "location": "Utils/checkpoint.html#API-1",
    "page": "Checkpointing",
    "title": "API",
    "category": "section",
    "text": ""
},

{
    "location": "Utils/checkpoint.html#Utils.Checkpointer",
    "page": "Checkpointing",
    "title": "Utils.Checkpointer",
    "category": "type",
    "text": "This type keeps track of which checkpoints are in use and has an API   for loading and saving checkpoints\n\nThe fields of the type should never be directly accessed.  Use the API   instead.\n\nEvery checkpoint must be in its own directory.\n\nCheckpointing has two uses: loading a previous state and restarting.   The first one is rather easy, all that needs to happen is the solution   variables get loaded.  Restarting is more involved because all the local   data of the NonlinearSolver that was running needs to be stored and then   loaded again in a new session.\n\nAn input file called input_vals_restart is created in the current directory   that can be used to restart the code.\n\nFields\n\nncheckpoints: total number of checkpoints\npaths: absolute paths to the checkpoints (array of strings)\nstatus: status of each checkpoint (used or free), (array of Ints)\nhistory: list of checkpoints from most recently used to least recently          used, unused entries set to -1\n\n\n\n"
},

{
    "location": "Utils/checkpoint.html#Utils.Checkpointer-Tuple{Integer,Integer,String}",
    "page": "Checkpointing",
    "title": "Utils.Checkpointer",
    "category": "method",
    "text": "Constructs a Checkpointer with a given number of checkpoints.\n\nThe checkpoints are put into directories called \"checkpoint\", where the   \"\" is replaced by the index of the checkpoint (from 1 to ncheckpoints).   This constructor does not error if the directories already exist, however   any data in the checkpoints may be overwritten.\n\nInputs\n\nmyrank: MPI rank of this process\nncheckpoints: number of checkpoints, defaults to 2.  Using fewer than               2 is not recommended (there will be a moment in time during               which the old checkpoint has been partially overwritten but               the new checkpoint is not complete yet)\nprefix: added to directory names (with an underscore added in-between).          Defaults to empty string. This is useful if there are multiple          sets of checkpoints (and therefore multiple Checkpointer objects).\n\nOutputs\n\na Checkpointer object, fully initialized\n\n\n\n"
},

{
    "location": "Utils/checkpoint.html#Utils.Checkpointer-Tuple{Dict,Integer}",
    "page": "Checkpointing",
    "title": "Utils.Checkpointer",
    "category": "method",
    "text": "This constructor loads a Checkpointer object from a file.  This should   be used for restarting.  Do not use the other constructor for restarting.   It will load the most recently written checkpoint that is complete (an   important distinction if the code was killed in the middle of writing   a checkpoint)\n\nInputs\n\nopts: the options dictionary.  \n\nOutputs\n\na Checkpointer object, fully initialized and aware of which checkpoints are in use and which are free\n\nThis function only loads the Checkpointer from the checkpoint.   See loadLastCheckpoint and readCheckpointData to   load the rest of checkpoint data.   The checkpoint that was loaded can be accessed via getLastCheckpoint\n\nImplementation notes:\n\nUses 4 keys in the options dictionary:\n\n\"writing checkpoint\": index of checkpoint that might be complete\n\"writing_checkpoint_path\": absolute path of checkpoint\n\"most_recent_checkpoint\": index of checkpoint that is definitely                            complete, but might be older than the above\n\"most_recent_checkpoint_path\": absolute path of checkpoint\n\nThis system ensure it is possible to restart even if the code is   killed in the middle of writing a checkpoint.\n\nThe Checkpointer object is the same on all processes.\n\n\n\n"
},

{
    "location": "Utils/checkpoint.html#Base.copy-Tuple{Utils.Checkpointer}",
    "page": "Checkpointing",
    "title": "Base.copy",
    "category": "method",
    "text": "Recursively copies all fields to make a new Checkpointer object.   Note that changes to one Checkpointer object will not affect the other   (ie. the record of which checkpoints are used and which are not).   This could easily lead to corrupting a checkpoint.   For this reason, this function should rarely be used.\n\n\n\n"
},

{
    "location": "Utils/checkpoint.html#Base.copy!-Tuple{Utils.Checkpointer,Utils.Checkpointer}",
    "page": "Checkpointing",
    "title": "Base.copy!",
    "category": "method",
    "text": "2 argument version of copy().  See that function for details.\n\n\n\n"
},

{
    "location": "Utils/checkpoint.html#Utils.AbstractCheckpointData",
    "page": "Checkpointing",
    "title": "Utils.AbstractCheckpointData",
    "category": "type",
    "text": "Every function that wants to checkpoint must implement an   AbstractCheckpointData containing all the data that needs to be loaded   when the code restarts.  The values in the AbstractSolutionData can be   different for each MPI process.\n\nThis type must contain only \"julia\" data, ie. objects that are managed   by Julia and have deterministic values.  No pointers, MPI communicators etc.\n\nFile IO is expensive, so include only as much data as needed in the   AbstractCheckpointData.  Anything that can be recalculated should be.\n\nGenerally, only the outermost function (ie. the time stepper or the   nonlinear solver for steady problems) will need to checkpoint.\n\nThere are no required fields for this type\n\n\n\n"
},

{
    "location": "Utils/checkpoint.html#Types-1",
    "page": "Checkpointing",
    "title": "Types",
    "category": "section",
    "text": "Checkpointer\nCheckpointer(::Integer, ::Integer, ::String)\nCheckpointer(::Dict, ::Integer)\ncopy(::Checkpointer)\ncopy!(::Checkpointer, ::Checkpointer)\nAbstractCheckpointData"
},

{
    "location": "Utils/checkpoint.html#Utils.saveNextFreeCheckpoint",
    "page": "Checkpointing",
    "title": "Utils.saveNextFreeCheckpoint",
    "category": "function",
    "text": "This function saves to the next free checkpoint\n\nInputs\n\ncheckpointer: the Checkpointer object\nmesh: an AbstractMesh\nsbp: an SBP operator\neqn: an AbstractSolutionData object\ncheckpoint_data: the users AbstractCheckpointData implementation\n\nInputs/Outputs\n\nopts: the options dictionary\n\nOutputs\n\ncheckpoint: the index of the checkpoint saved.\n\n\n\n"
},

{
    "location": "Utils/checkpoint.html#Utils.loadLastCheckpoint",
    "page": "Checkpointing",
    "title": "Utils.loadLastCheckpoint",
    "category": "function",
    "text": "This function loads the most recently saved checkpoint.\n\nInputs\n\ncheckpointer: the CheckPointer object\nmesh: an AbstractMesh\nsbp: an SBP operator\nopts: the options dictionary\n\nInputs/Outputs\n\neqn: an AbstractSolutionData object\n\nOutputs\n\nthe checkpoint index that was loaded\n\n\n\n"
},

{
    "location": "Utils/checkpoint.html#Utils.readLastCheckpointData",
    "page": "Checkpointing",
    "title": "Utils.readLastCheckpointData",
    "category": "function",
    "text": "Simple wrapper around readCheckpointData to loading the most   recently saved checkpoint data\n\nInputs\n\nchkpointer: a Checkpointer\ncomm_rank: MPI rank of this process\n\n\n\n"
},

{
    "location": "Utils/checkpoint.html#Utils.readCheckpointData",
    "page": "Checkpointing",
    "title": "Utils.readCheckpointData",
    "category": "function",
    "text": "This function reads an object of type AbstractSolutionData from a file   and reuturns the object.  This is type unstable, so every AbstractCheckpointData   implementation should create a constructor that calls this function and   uses a type assertion to specify that the object must be of their concrete   type.\n\nInputs\n\nchkpointer: a Checkpointer\nchkpoint: the index of hte checkpoint\ncomm_rank: the MPI rank of this process\n\nOutputs\n\nthe AbstractCheckpointData loaded from the checkpoint\n\n\n\n"
},

{
    "location": "Utils/checkpoint.html#Utils.countFreeCheckpoints",
    "page": "Checkpointing",
    "title": "Utils.countFreeCheckpoints",
    "category": "function",
    "text": "This function returns the number of free checkpoints\n\nInputs\n\ncheckpointer: the Checkpoint object\n\n\n\n"
},

{
    "location": "Utils/checkpoint.html#Utils.getNextFreeCheckpoint",
    "page": "Checkpointing",
    "title": "Utils.getNextFreeCheckpoint",
    "category": "function",
    "text": "This function returns the index of the next free checkpoint,  or 0 if   there are no free checkpoints\n\nInput\n\ncheckpointer: a Checkpointer\n\n\n\n"
},

{
    "location": "Utils/checkpoint.html#Utils.getLastCheckpoint",
    "page": "Checkpointing",
    "title": "Utils.getLastCheckpoint",
    "category": "function",
    "text": "This function returns the index of the most recently written checkpoint\n\nInputs\n\ncheckpointer: the Checkpointer object\n\nOutputs\n\nthe index of the most recently saved checkpoint\n\n\n\n"
},

{
    "location": "Utils/checkpoint.html#Utils.getOldestCheckpoint",
    "page": "Checkpointing",
    "title": "Utils.getOldestCheckpoint",
    "category": "function",
    "text": "Get the least recently written checkpoint\n\nInputs\n\ncheckpointer: the Checkpointer\n\nOutputs\n\nthe index of the least recently written checkpoint\n\n\n\n"
},

{
    "location": "Utils/checkpoint.html#Utils.freeOldestCheckpoint",
    "page": "Checkpointing",
    "title": "Utils.freeOldestCheckpoint",
    "category": "function",
    "text": "Frees the least recently written checkpoint. Calling this function when   all checkpoints are free is allowed.\n\nInputs\n\ncheckpointer: the Checkpointer\n\nOutputs\n\nreturns the checkpoint freed (0 all checkpoints were free on entry)\n\n\n\n"
},

{
    "location": "Utils/checkpoint.html#Utils.freeCheckpoint",
    "page": "Checkpointing",
    "title": "Utils.freeCheckpoint",
    "category": "function",
    "text": "This function frees a checkpoint (marks it as available to be overwritten)   Unlike certain algorithms that use free lists (cough, cough, malloc)   freeing an already free checkpoint is allowed.\n\nInputs\n\ncheckpointer: the Checkpointer object\ncheckpoint: the index of the checkpoint to free\n\nThe user must explictly free checkpoints (loading a checkpoint does not   free it).   Note that a free checkpoint can still be loaded for restarting if it has   not been saved to yet.\n\n\n\n"
},

{
    "location": "Utils/checkpoint.html#Functions-1",
    "page": "Checkpointing",
    "title": "Functions",
    "category": "section",
    "text": "These function provide the basic operations required for checkpointingsaveNextFreeCheckpoint\nloadLastCheckpoint\nreadLastCheckpointData\nreadCheckpointData\ncountFreeCheckpoints\ngetNextFreeCheckpoint\ngetLastCheckpoint\ngetOldestCheckpoint\nfreeOldestCheckpoint\nfreeCheckpoint"
},

{
    "location": "Utils/checkpoint.html#Example-1",
    "page": "Checkpointing",
    "title": "Example",
    "category": "section",
    "text": "This section provides a simple example of how to use checkpointing for an explicit time marching scheme. The original scheme (without checkpointing) looks like:function explicit_timemarching_example(mesh, sbp, eqn, opts, func::Function)\n\n  nsteps = round(Int, opts[\"tmax\"]/opts[\"delta_t\"])\n  t = 0.0\n\n  for step=1:nsteps\n\n    t = (step - 1)*delta_t\n    # q_vec -> q\n    array1DTo3D(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)\n    func(mesh, sbp, eqn, opts, t)\n\n    # res -> res_vec\n    array3DTo1D(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)\n\n    for i=1:mesh.numDof\n      eqn.q_vec[i] = ... # some function of eqn.res_vec[i]\n    end\n  end  # end loop over timesteps\n\n  return nothing\nendWith checkpointing, it becomes:\ntype ExampleCheckpointData <: AbstractCheckpointData\n  step::Int  # for this method, the only thing that needs to be stored is the\n             # current time step, everything else can be recalculated\nend\n\n# this constructor is used when restarting, see below\nfunction ExampleCheckpointData(chkpointer::Checkpointer, comm_rank)\n\n  chkpoint_data = readLastCheckpointData(chkpointer, comm_rank)\n\n  return chkpoint_data::ExampleCheckpointData\nend\n\nfunction explicit_timemarching_example(mesh, sbp, eqn, opts, func::Function)\n\n  nsteps = round(Int, opts[\"tmax\"]/opts[\"delta_t\"])\n  t = 0.0\n\n  if !opts[\"is_restart\"]\n    # regular starting of a run\n    stepstart = 1\n    chkpointdata = ExampleCheckpointData(stepstart)\n    # this is valid even if we don\'t intend to checkpoint (ncheckpoints = 0)\n    chkpointer = Checkpointer(mesh.myrank, opts[\"ncheckpoints\"]\n    skip_checkpoint = false\n  else\n    chkpointer = Checkpointer(opts, myrank)  # read Checkpointer from most recent\n                                             # checkpoint\n    # now get the checkpoint data telling which timestep we are on\n    # eqn.q_vec already holds the state saved in the checkpoint.  This is handled\n    # during the initialization of the equation object\n    chkpointdata = ExampleCheckpointData(chkpoinnter, mesh.myrank)\n    stepstart = chkpointdata.step\n    skip_checkpoint = true  # skip writing the first checkpoint.  Without this, a\n                            # checkpoint would be written immediately.\n                            # Writing the checkpoint is not harmful, but not useful\n                            # either.\n   end\n    \n  for step=stepstart:nsteps  # start the loop from stepstart\n\n    t = (step - 1)*delta_t   # t is recalculated from step, which was set using\n                             # stepstart, thus only stepstart needs to be saved\n\n    # save checkpoint, if needed\n    if opts[\"use_checkpointing\"] && step % opts[\"checkpoint_freq\"] == 0 && !skip_checkpoint\n      skip_checkpoint = false\n\n      # save all needed variables to the ExampleCheckpointData object\n      # For simple time marching schemes, step is the only needed data\n      chkpointdata.step = step\n\n      if countFreeCheckpoints(chkpointer) == 0\n        freeOldestCheckpoint(chkpointer)  # make room for a new checkpoint\n      end\n\n      saveNextFreeCheckpoint(chkpointer, mesh, sbp, eqn, opts, chkpointdata)\n    end\n\n    # q_vec -> q\n    array1DTo3D(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)\n    func(mesh, sbp, eqn, opts, t)\n\n    # res -> res_vec\n    array3DTo1D(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)\n\n    for i=1:mesh.numDof\n      eqn.q_vec[i] = ... # some function of eqn.res_vec[i]\n    end\n  end  # end loop over timesteps\n\n  return nothing\nend"
},

{
    "location": "Utils/checkpoint.html#Utils.writeCheckpointer",
    "page": "Checkpointing",
    "title": "Utils.writeCheckpointer",
    "category": "function",
    "text": "Call on master process only\n\nDoes not check if the checkpoint is already used (for reading the   state back after restart, this checkpoint should already be marked as   used)\n\n\n\n"
},

{
    "location": "Utils/checkpoint.html#Utils.saveCheckpoint",
    "page": "Checkpointing",
    "title": "Utils.saveCheckpoint",
    "category": "function",
    "text": "Save to a specified checkpoint.  Throw error if checkpoint is not free.   Users should not generally call this function directly.  Instead, they   should prefer saveNextFreeCheckpoint.\n\nThis function automatically saves eqn.q_vec and the checkpointer to   a file.  Any additional data should be in checkpoint_data.\n\nNote: mesh adaptation is not compatable with checkpointing   #TODO; add a field to the mesh to record that it has been modified\n\nInputs\n\ncheckpointer: the CheckPointer\ncheckpoint: the index of the checkpoint\nmesh: an AbstractMesh object\nsbp: SBP operator\neqn: an AbstractSolutionData\ncheckpoint_data: an AbstractCheckpointData.  This is the random bag                  of data the user needs saved.\n\nInputs/Outputs\n\n* opts: options dictionary\n\nImplementation Notes:\n\nUses options dictionary keys described by Checkpointer   Note that the checkpoint is eagerly marked as used, before finishing writing   the checkpoint.  Upon restart the code needs to check if this checkpoint   is really finished.\n\n\n\n"
},

{
    "location": "Utils/checkpoint.html#Utils.loadCheckpoint",
    "page": "Checkpointing",
    "title": "Utils.loadCheckpoint",
    "category": "function",
    "text": "Loads a specified checkpoint.  Only loads eqn.q_vec.  Does not load   the checkpointer (because this is loading a checkpoint and not a restart).   Also does not load the AbstractCheckpointData.   See readCheckpointData for loading of the AbstractCheckpointData.   Users should generally not call this function directly.  Users should   prefer loadLastCheckpoint whenever possible.\n\nNote that loading a checkpoint does not mark it as free.  Users must   explictly call freeCheckpoint.\n\nInputs\n\ncheckpointer: the Checkpointer object\ncheckpoint: the index of the checkpoint\nmesh: an AbstractMesh object\nsbp: the SBP operator\nopts: the options dictionary\n\nInputs/Output\n\neqn: an AbstractSolutionData, eqn.q_vec is overwritten\n\n\n\n"
},

{
    "location": "Utils/checkpoint.html#Utils.writeFlagFile",
    "page": "Checkpointing",
    "title": "Utils.writeFlagFile",
    "category": "function",
    "text": "Write the flag file that signifies a checkpoint is complete   Call on master process only!\n\nInputs\n\ncheckpointer: the Checkpointer object\ncheckpoint: the checkpoint index\n\n\n\n"
},

{
    "location": "Utils/checkpoint.html#Utils.checkFlagFile",
    "page": "Checkpointing",
    "title": "Utils.checkFlagFile",
    "category": "function",
    "text": "Returns true if the flag file exists and is consistent, returns false   otherwise\n\nInputs\n\ncheckpointer: the Checkpointer object\ncheckpoint: the checkpoint index\n\n\n\nSometimes need to check the flag file before the Checkpointer is available.   See the other method ofr details.\n\nInputs\n\nfpath: path to the checkpoint directory\n\n\n\n"
},

{
    "location": "Utils/checkpoint.html#Utils.deleteFlagFile",
    "page": "Checkpointing",
    "title": "Utils.deleteFlagFile",
    "category": "function",
    "text": "Deletes the flag file if it exists.  Does not error if flag file does not   exist.   Call on master process only\n\nInputs\n\ncheckpointer: the Checkpointer object\ncheckpoint: the checkpoint index\n\n\n\n"
},

{
    "location": "Utils/checkpoint.html#Utils.writeCheckpointData",
    "page": "Checkpointing",
    "title": "Utils.writeCheckpointData",
    "category": "function",
    "text": "Writes the AbstractCheckpointData to a file.  The file can be read by   readCheckpointData.\n\nInputs\n\ncheckpoint: Checkpointer object\ncheckpoint: the checkpoint index\nobj: the AbstractCheckpointData object\ncomm_rank: the MPI rank of this process\n\n\n\n"
},

{
    "location": "Utils/checkpoint.html#Utils.markCheckpointUsed",
    "page": "Checkpointing",
    "title": "Utils.markCheckpointUsed",
    "category": "function",
    "text": "Marks a checkpoint as unused.\n\nInputs\n\ncheckpointer: the Checkpointer object\ncheckpoint: the checkpoint index\n\n\n\n"
},

{
    "location": "Utils/checkpoint.html#Internal-Functions-1",
    "page": "Checkpointing",
    "title": "Internal Functions",
    "category": "section",
    "text": "The internal functions used for checkpointing are documented here. Users should not call these functions directly.  Improper use can cause checkpoint corruption.  writeCheckpointer\n  saveCheckpoint\n  loadCheckpoint\n  writeFlagFile\n  checkFlagFile\n  deleteFlagFile\n  writeCheckpointData\n  markCheckpointUsed"
},

{
    "location": "Utils/misc.html#",
    "page": "Misccellaneous",
    "title": "Misccellaneous",
    "category": "page",
    "text": ""
},

{
    "location": "Utils/misc.html#Utils",
    "page": "Misccellaneous",
    "title": "Utils",
    "category": "module",
    "text": "Module Utils:   This module holds miscellaneous functions used throughout the code\n\n\n\n"
},

{
    "location": "Utils/misc.html#Utils.Timings",
    "page": "Misccellaneous",
    "title": "Utils.Timings",
    "category": "type",
    "text": "Utils.Timings\n\nThis type accumulates the time spent in each part of the code.\n\n\n\n"
},

{
    "location": "Utils/misc.html#Utils.absvalue_deriv-Union{Tuple{Tval}, Tuple{Tval}} where Tval",
    "page": "Misccellaneous",
    "title": "Utils.absvalue_deriv",
    "category": "method",
    "text": "###Utils.absvalue_deriv\n\nComputes the derivative of the absolute value of a variable w.r.t itself\n\nInputs\n\nval : The variable whose derivative needs to be computed e.r.t itself\n\n\n\n"
},

{
    "location": "Utils/misc.html#Utils.applyPermColumn-Tuple{AbstractArray{T,1} where T,AbstractArray{T,2} where T,AbstractArray{T,2} where T}",
    "page": "Misccellaneous",
    "title": "Utils.applyPermColumn",
    "category": "method",
    "text": "Permute the columns of A according to permvec.  See applyPermRow for the   definition of a permvec.  Note that for column permutation the   interpretation is that destination column permvec[i] comes from   source column i.  This is consistent with the notion that post-multiplying   by a permutation matrix (obtained from permMatrix) is a column   permutation, i.e B = A*P\n\nAliasing: no aliasing allowed\n\n\n\n"
},

{
    "location": "Utils/misc.html#Utils.applyPermColumnInplace-Tuple{AbstractArray{T,1} where T,AbstractArray{T,2} where T,AbstractArray{T,2} where T}",
    "page": "Misccellaneous",
    "title": "Utils.applyPermColumnInplace",
    "category": "method",
    "text": "Like applyPermColumn, but the result is returned in A.  Both A and B   are overwritten\n\nAliasing: no aliasing allowed\n\n\n\n"
},

{
    "location": "Utils/misc.html#Utils.applyPermRow-Tuple{AbstractArray{T,1} where T,AbstractArray{T,2} where T,AbstractArray{T,2} where T}",
    "page": "Misccellaneous",
    "title": "Utils.applyPermRow",
    "category": "method",
    "text": "Permute the rows of A according to the permvec, storing the result in B   The permvec contains the source indices for each entry in B, ie.   B[i] comes from A[permvec[i]].  This is consistent with the mathematical   definition of a permutation that pre-multiplication by a permutation    matrix (obtained from permMatrix) is a row permutation, ie.   B = P*A\n\nAliasing: no aliasing allowed\n\n\n\n"
},

{
    "location": "Utils/misc.html#Utils.applyPermRowInplace-Tuple{AbstractArray{T,1} where T,AbstractArray{T,2} where T,AbstractArray{T,2} where T}",
    "page": "Misccellaneous",
    "title": "Utils.applyPermRowInplace",
    "category": "method",
    "text": "Like applyPermRow, but the result is returned in A.  Both A and B   are overwritten.\n\nAliasing: no aliasing allowed\n\n\n\n"
},

{
    "location": "Utils/misc.html#Utils.array1DTo3D-Union{Tuple{ODLCommonTools.AbstractCGMesh,Any,ODLCommonTools.AbstractSolutionData,Any,AbstractArray{T,1},AbstractArray{T,3}}, Tuple{T}} where T",
    "page": "Misccellaneous",
    "title": "Utils.array1DTo3D",
    "category": "method",
    "text": "Utils.array1DTo3D\n\nThis takes eqn.q_vec (the initial state), and disassembles it into eqn.q, the   3 dimensional array.  This function uses mesh.dofs   to speed the process.\n\nThis function also calls writeQ to do any requested output.\n\nInputs:     mesh     sbp     eqn     opts\n\nThis is a mid level function, and does the right thing regardless of equation   dimension.\n\nThe DG method for array1DTo3D assumes that q and q_vec refer to the     same memory address, and therefore does no explicit writing/copying.\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "Utils/misc.html#Utils.array3DTo1D-Union{Tuple{ODLCommonTools.AbstractCGMesh{Tmsh},Any,ODLCommonTools.AbstractSolutionData{Tsol,Tres} where Tres,Any,AbstractArray{T,3} where T,AbstractArray{Tres,1},Any}, Tuple{ODLCommonTools.AbstractCGMesh{Tmsh},Any,ODLCommonTools.AbstractSolutionData{Tsol,Tres} where Tres,Any,AbstractArray{T,3} where T,AbstractArray{Tres,1}}, Tuple{Tmsh}, Tuple{Tres}, Tuple{Tsol}} where Tres where Tsol where Tmsh",
    "page": "Misccellaneous",
    "title": "Utils.array3DTo1D",
    "category": "method",
    "text": "Utils.array3DTo1D\n\nThis function takes the 3D array of variables in arr and   reassembles it into the vector res_vec.  Note that   This is a reduction operation and zeros res_vec before performing the   operation, unless zero_resvec is set to false\n\nThis is a mid level function, and does the right thing regardless of   equation dimension\n\n\n\n"
},

{
    "location": "Utils/misc.html#Utils.calcBCNormal-Tuple{ODLCommonTools.AbstractParamType{2},AbstractArray{T,2} where T,AbstractArray{T,1} where T,AbstractArray{T,1} where T}",
    "page": "Misccellaneous",
    "title": "Utils.calcBCNormal",
    "category": "method",
    "text": "Calculate the scaled normal vector in parametric coordinates from the   face normal and scaled mapping jacobian.  nrm2 is overwritten with   the result.\n\n\n\n"
},

{
    "location": "Utils/misc.html#Utils.calcEuclidianNorm-Union{Tuple{MPI.Comm,AbstractArray{T,1}}, Tuple{T}} where T",
    "page": "Misccellaneous",
    "title": "Utils.calcEuclidianNorm",
    "category": "method",
    "text": "This function computes the Euclidian norm of a vector where each MPI   process owns part of the vector\n\nInputs:     comm: an MPI communicator     vec: the local part of the vector\n\nOutputs:     val: the Euclidian norm of the entire vector\n\nNote that unlike calcNorm, the time spent in the Allreduce is not logged   for this function.\n\n\n\n"
},

{
    "location": "Utils/misc.html#Utils.calcL2InnerProduct-Union{Tuple{ODLCommonTools.AbstractSolutionData,AbstractArray{T,N} where N,AbstractArray{T2,N} where N}, Tuple{T2}, Tuple{T}} where T2 where T",
    "page": "Misccellaneous",
    "title": "Utils.calcL2InnerProduct",
    "category": "method",
    "text": "This function calculate the SBP approximation to the L2 inner product of   two vectors of length mesh.numDof.\n\nL2_product = u^T H conj(v)\n\nNote that this function does not take a square root like calcNorm   does.\n\nInputs\n\neqn: AbstractSolutionData\nu: first vector, of length mesh.numDof\nv: second vector, of length mesh.numDof.  If complex, this vector gets    conjugated\n\nOutputs\n\nval: the value of the inner product\n\nKeyword Arguments\n\nglobalnrm: if true, computes the norm of the entire (parallel) vector,            if false, computes the norm of only the local part\n\n\n\n"
},

{
    "location": "Utils/misc.html#Utils.calcMeshH-Union{Tuple{ODLCommonTools.AbstractMesh{Tmsh},Any,Any,Any}, Tuple{Tmsh}} where Tmsh",
    "page": "Misccellaneous",
    "title": "Utils.calcMeshH",
    "category": "method",
    "text": "Utils.calcMeshH\n\nThis function calculates the average distance between nodes over the entire   mesh.  This function allocates a bunch of temporary memory, so don\'t call   it too often.  This is, strictly speaking, not quite accurate in parallel   because the divison by length happens before the allreduce.\n\nInputs:     mesh     eqn     opts\n\n\n\n"
},

{
    "location": "Utils/misc.html#Utils.calcNorm-Union{Tuple{ODLCommonTools.AbstractSolutionData,AbstractArray{T,N} where N}, Tuple{T}} where T",
    "page": "Misccellaneous",
    "title": "Utils.calcNorm",
    "category": "method",
    "text": "Utils.calcNorm\n\nThis function calculates the norm of a vector (of length numDof) using the     SBP norm.\n\nInputs:\n  eqn:  an AbstractSolutionData\n  res_vec:  vector to calculate the norm of\n\nKeyword arguments:\n  strongres: if res_vec is the residual of the weak form, then\n             strongres=true computes (efficiently) the norm of the strong\n             form residual.  Default false\n  globalnrm: compute the norm over all processes or not. Default true\n\nReturns:\n  val:  norm of solution using SBP norm (Float64)\n\nThere are no restrctions on the datatype of res_vec (ie. it can be complex)\n\nAliasing restrictions: none\n\n\n\n"
},

{
    "location": "Utils/misc.html#Utils.fastscale!-Tuple{AbstractArray,Number}",
    "page": "Misccellaneous",
    "title": "Utils.fastscale!",
    "category": "method",
    "text": "This function scales an array by a constant, and should be faster than scale!   because it is branch free\n\n\n\n"
},

{
    "location": "Utils/misc.html#Utils.fastzero!-Tuple{AbstractArray}",
    "page": "Misccellaneous",
    "title": "Utils.fastzero!",
    "category": "method",
    "text": "This function zeros out an array, and should be faster than fill! (branch   free)\n\nInputs/Outputs:     x: an array\n\n\n\n"
},

{
    "location": "Utils/misc.html#Utils.free-Tuple{Any}",
    "page": "Misccellaneous",
    "title": "Utils.free",
    "category": "method",
    "text": "Generic function to free any memory belonging to other libraries\n\n\n\n"
},

{
    "location": "Utils/misc.html#Utils.inversePerm-Tuple{AbstractArray{T,1} where T,AbstractArray{T,1} where T}",
    "page": "Misccellaneous",
    "title": "Utils.inversePerm",
    "category": "method",
    "text": "Compute the permutation vector that corresponds to the inverse permutation\n\nInputs:     permvec: the original permutation vector\n\nInputs/Outputs:     invperm: the inverse permutation vector\n\nAliasing: no aliasing allowed\n\n\n\n"
},

{
    "location": "Utils/misc.html#Utils.permMatrix!-Tuple{AbstractArray{T,1} where T,AbstractArray{T,2} where T}",
    "page": "Misccellaneous",
    "title": "Utils.permMatrix!",
    "category": "method",
    "text": "Create a permutation matrix from a permutation vector.  Only select   entries of A are overwritten.\n\n\n\n"
},

{
    "location": "Utils/misc.html#Utils.permMatrix-Tuple{AbstractArray{T,1} where T}",
    "page": "Misccellaneous",
    "title": "Utils.permMatrix",
    "category": "method",
    "text": "Create a permutation matrix from a permutation vector.  The element type   of the returned matrix is Int.\n\n\n\n"
},

{
    "location": "Utils/misc.html#Utils.removeComplex-Tuple{ODLCommonTools.AbstractMesh,SummationByParts.AbstractSBP,ODLCommonTools.AbstractSolutionData,Any}",
    "page": "Misccellaneous",
    "title": "Utils.removeComplex",
    "category": "method",
    "text": "Set the complex part of the solution to zero, including eqn.q, eqn.q_vec, and   the send and receive buffers in eqn.shared_data\n\nInputs\n\nmesh\nsbp\neqn\nopts\n\n\n\n"
},

{
    "location": "Utils/misc.html#Utils.writeQ-NTuple{4,Any}",
    "page": "Misccellaneous",
    "title": "Utils.writeQ",
    "category": "method",
    "text": "Utils.writeQ\n\nThis function writes the real part of the solution variables eqn.q to a space   delimited file called q.dat, controlled by the input options \'writeq\', of type bool\n\nThis is a high level function.\n\n\n\n"
},

{
    "location": "Utils/misc.html#Utils.write_timings-Tuple{Utils.Timings,AbstractString}",
    "page": "Misccellaneous",
    "title": "Utils.write_timings",
    "category": "method",
    "text": "Utils.write_timings\n\nWrite the values in a Timings object to a file.  Also writes the names of   fields to a separate file.\n\nInputs:     t: a  Timings object     fname: the file name, without extension\n\nThe values are written to the file fname.dat, and the names are written to   fname_names.dat\n\n\n\n"
},

{
    "location": "Utils/misc.html#Utils.@verbose1-Tuple{Any}",
    "page": "Misccellaneous",
    "title": "Utils.@verbose1",
    "category": "macro",
    "text": "Utils.verbose1\n\nThis macro introduces an if statement that causes the expression to be    executed only if the variable verbose is greater than or equal to 1.     verbose must exist in the scope of the caller\n\n\n\n"
},

{
    "location": "Utils/misc.html#Utils.@verbose2-Tuple{Any}",
    "page": "Misccellaneous",
    "title": "Utils.@verbose2",
    "category": "macro",
    "text": "Utils.verbose2\n\nThis macro introduces an if statement that causes the expression to be    executed only if the variable verbose is greater than or equal to 2.     verbose must exist in the scope of the caller\n\n\n\n"
},

{
    "location": "Utils/misc.html#Utils.@verbose3-Tuple{Any}",
    "page": "Misccellaneous",
    "title": "Utils.@verbose3",
    "category": "macro",
    "text": "Utils.verbose3\n\nThis macro introduces an if statement that causes the expression to be    executed only if the variable verbose is greater than or equal to 3.     verbose must exist in the scope of the caller\n\n\n\n"
},

{
    "location": "Utils/misc.html#Utils.@verbose4-Tuple{Any}",
    "page": "Misccellaneous",
    "title": "Utils.@verbose4",
    "category": "macro",
    "text": "Utils.verbose4\n\nThis macro introduces an if statement that causes the expression to be    executed only if the variable verbose is greater than or equal to 4.     verbose must exist in the scope of the caller\n\n\n\n"
},

{
    "location": "Utils/misc.html#Utils.@verbose5-Tuple{Any}",
    "page": "Misccellaneous",
    "title": "Utils.@verbose5",
    "category": "macro",
    "text": "Utils.verbose5\n\nThis macro introduces an if statement that causes the expression to be    executed only if the variable verbose is greater than or equal to 5.     verbose must exist in the scope of the caller\n\n\n\n"
},

{
    "location": "Utils/misc.html#Miscellaneous-1",
    "page": "Misccellaneous",
    "title": "Miscellaneous",
    "category": "section",
    "text": "  CurrentModule = Utils  Modules = [Utils]\n  Pages = [\"Utils/Utils.jl\"]"
},

{
    "location": "test/Testing.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "test/Testing.html#Testing-1",
    "page": "Introduction",
    "title": "Testing",
    "category": "section",
    "text": "PDESolver has a test system to verify the code is working properly. Every feature of the code should have a test that verifies the correctness of the feature.  No exceptions!Tests can take the form of unit tests, which verify the correctness of an individual function, or system integration tests, which verify several pieces of code (which were unit tested individually) are working together correctly.Some examples of unit tests are:testing that a flux function returns a known value (for cases when the numerical scheme is exact)\ntesting that a two-point numerica flux function returns the analytical flux when the left and right state are the same (the consistency property)Some examples of system integration tests are:Testing a uniform flow has zero residual\nTesting convergence rates"
},

{
    "location": "test/Testing.html#Contents-1",
    "page": "Introduction",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"Readme.md\", \"Travis.md\", \"TestSystem.md\"]\nDepth = 1"
},

{
    "location": "test/Readme.html#",
    "page": "Local Testing",
    "title": "Local Testing",
    "category": "page",
    "text": ""
},

{
    "location": "test/Readme.html#Test-System-1",
    "page": "Local Testing",
    "title": "Test System",
    "category": "section",
    "text": "The Ticon test system is composed of functions that test different aspects of the code.  The test system allows tags to be associated with each test function, and selective running of tests based on their tags. Currently, each physics module has its own set of tags and test list. The serial tests for a single physics module  are run as follows:  cd test/physics_module\n  julia ./runtests.jl  TAG1 TAG2 ...The parallel tests are run with:  cd test/physics_module\n  mpirun -np 2 julia ./runtests_parallel.jl  TAG1 TAG2 ...  cd test/physics_module\n  mpirun -np 4 julia ./runtests_parallel4.jl  TAG1 TAG2 ...for the 2 and 4 processors tests, respectivelyThe parallel tests for all physics modules are run with:  ./runtests_parallel.sh TAG1 TAG2...and all tests (serial + parallel) for all physics modules are run with  ./runtests.sh TAG1 TAG2Running the tests in a single julia session can be faster because the code does not have to be recompiled for each physics.  The serial tests can be run with:  julia ./runtests.jl TAG1 TAG2 ...the 2 processor tests with  mpirun -np 2 julia ./runtests_parallel2.jl TAG1 TAG2 ...and the 4 processor tests with  mpirun -np 4 julia ./runtests_parallel4.jl TAG1 TAG2 ...All the tests (serial + parallel) can be run with  ./runtests_fast.sh TAG1 TAG2 ...For all scripts except runtests_fast.sh, if no tags are specified, all tests are run. If any tags are specified, only functions matching those tags are run. For runtest_fast.sh only the tests with TAG_SHORTTEST are run by default. If any tags are specified, only tests matching the tags (and not TAG_SHORTTEST are run). The list of currently defined tags can be found in ./tags.jl Note that the serial tests must  be run before the parallel tests. All test must have one of LengthTags associated with them. Any test that takes less than 60 seconds to run should be specified as a TAG_SHORTTEST, otherwise TAG_LONGTEST."
},

{
    "location": "test/Readme.html#Adding-Tests-1",
    "page": "Local Testing",
    "title": "Adding Tests",
    "category": "section",
    "text": "Each physics module maintains a TestSet object which contains a list of tests and their associated tags.  See the previous paragraph for a note about the required tags.  TestSet.jl provides an API for adding tests to the list.All tests must be enclosed in functions.  These functions must fit into one of three catagories:zero argument functions\nfunctions that take the arguments (mesh, sbp, eqn, opts) as produced by an input file which already exists\nfunctions that take the arguments (mesh, sbp, eqn, opts) as produced by modifying an existing input fileThese catagories correspond to add_funcs1!, add_funcs2!, and add_funcs3! in TestSet.jl.The function run_testlist runs a test list.  See the documentation of these functions for details on how to use them.Note that it is acceptable for a test function to load input files internally. The benefit to type 2 functions over type 1 functions is in the case when several functions share an input file.  The test system will only load the file once if possible, rather than for every test function that uses is. Consequently, each test function must not depend on modifications made to the the (mesh, sbp, eqn,  opts) objects by other test functions, because the objects may or may not be re-initialized between calls to the test functions.  Also, because tags can be used to selectively run tests, it cannot be guaranteed that a particular test will run before another one."
},

{
    "location": "test/Travis.html#",
    "page": "CI Testing",
    "title": "CI Testing",
    "category": "page",
    "text": ""
},

{
    "location": "test/Travis.html#Travis-CI-Testing-1",
    "page": "CI Testing",
    "title": "Travis CI Testing",
    "category": "section",
    "text": "We currently use the Travis CI service to run automated tests on all branches.  This page describes the build system"
},

{
    "location": "test/Travis.html#Build-Matrix-1",
    "page": "CI Testing",
    "title": "Build Matrix",
    "category": "section",
    "text": "A build matrix is used to create one CI build per physics module. Each build sets an environment variable to indicate to the test/runtests_travis.sh script which physics module to run tests on. For example, the environment variablesTEST_ADVECTION=1\nTEST_EULER=1\nTEST_SIMPLEODE=1are used to decide whether or not to run tests for the Advection, Euler, and SimpleODE physics modules. Note that the environment variables set in  the env section of .travis.yml must exactly match the variables used in test/runtests_travis.sh."
},

{
    "location": "test/Travis.html#Caching-Dependencies-1",
    "page": "CI Testing",
    "title": "Caching Dependencies",
    "category": "section",
    "text": "It is possible to have Travis cache some of PDESolvers dependencies to save build time.  In particular, Petsc takes some time to compile and does not change very often, so it would be beneficial to cache it. One downside of caching is that a manual process is required to update the cache if we change version of Petsc.The configuration process is described here.  Most of the work is done in the deps/move_petsccache.sh script.Each branch has its own cache, and the cache needs to be created during the first build on a new branch.  Therefore, Travis must be configured to detect if the cache has been created and build the cache if needed. One additional wrinkle is that the directories created by the cache interfere with the PDESolver build system.  The .travis.yml configuration file executes the following procedure:If path/to/PETSc2/deps/petsc-3.7.6/arch-linux2-c-debug/lib/libpetsc.so exists, move the /lib directory to $HOME/lib, set PETSC_DIR andPETSC_ARCH, and delete /path/to/PETSc2.  If the cache did not exist,    the directory will still exist, but will be empty.  Delete /path/to/PETSc2 so the PDESolver build system will install the PETSc2 package.Install PDESolver and its dependencies\nIf PETSC_DIR and PETSC_ARCH are defined (done in step 1 if cache was retrieved), move $HOME/lib back to original location, creating any directories in the path as needed.If PETSC_DIR and PETSC_ARCH are not set in step 1, step 2 will build PETSc.  At the end of the Travis run, the /path/to/PETSc2/deps/petsc-3.7.6/arch-linux2-c-debug/lib directory will be uploaded to wherever Travis stores cache files.  It will then be available for the next Travis build.If the cache has already been built for a given branch and a different version of PETSc is needed, the cache must be deleted so step 2 can build the new version. The web interface to PETSc can be used to delete the cache."
},

{
    "location": "test/TestSystem.html#",
    "page": "Test API",
    "title": "Test API",
    "category": "page",
    "text": ""
},

{
    "location": "test/TestSystem.html#TestSystem.TestList",
    "page": "Test API",
    "title": "TestSystem.TestList",
    "category": "type",
    "text": "This type stores all the data that describes a set of tests and their   associated tags.\n\nFields:     funcs: list of functions that contain the tests.  The order in which            the tests run is the same as the order of insertion\n\nfunc_tags: Vector, same length as funcs, where each element is a vector\n           of tags associated with each func.\n\nfunc_type:  Vector, same length as funcs, describing how to run each\n            function. 1 = no arguments, 2 = args = (mesh, sbp, eqn, opts),\n            3 = (mesh, sbp, eqn, opts), using a modified version of an\n            existing input file\ninput_name: Vector, same length as funcs, containing the name of the input\n            file to use with each function.  For func_type == 1 this value\n            is unused, for func type == 2 it is used directly, for \n            func_type == 3 it is modified and written to a new name\n            according to mod_input\n\nmod_input: Vector, same length as funcs, where each element is a dictionary.\n           For func_type == 3 the input file specified by input_name is\n           modified with these keys.  There must also be a key called\n           \"new_fname\" that specifies name to write the modified file\n           to (not including file extension.  \n           This modified file will be loaded before the function \n           is called.\n\ntag_list: collection of all known tags\n\n\n\n"
},

{
    "location": "test/TestSystem.html#TestSystem.add_func1!",
    "page": "Test API",
    "title": "TestSystem.add_func1!",
    "category": "function",
    "text": "This function adds a new test function of func_type == 1 to the list   (a function that takes no arguments).\n\nInputs     test_list: the list of tests to append the function to     func: the function     tags: an array of tags associated with this function.  The user is free           to modify this array afterwards, but not to mutate the strings within           the array. (optional)\n\n\n\n"
},

{
    "location": "test/TestSystem.html#TestSystem.add_func2!",
    "page": "Test API",
    "title": "TestSystem.add_func2!",
    "category": "function",
    "text": "This function adds a new test function of func_type == 2 to the list   (a function that takes the arguments (mesh, sbp, eqn, opts), as   obtained from the input file specified).\n\nInputs     test_list: the list of tests to append the function to     func: the function     input_name: name of input file to be used with this function     tags: an array of tags associated with this function.  The user is free           to modify this array afterwards, but not to mutate the strings within           the array. (optional)\n\n\n\n"
},

{
    "location": "test/TestSystem.html#TestSystem.add_func3!",
    "page": "Test API",
    "title": "TestSystem.add_func3!",
    "category": "function",
    "text": "This function adds a new test function of func_type == 2 to the list   (a function that takes the arguments (mesh, sbp, eqn, opts), as   obtained from modifying the input file specified).\n\nInputs     test_list: the list of tests to append the function to     func: the function     input_name: name of input file to be modified for use with this function     mod_dict: dictionary of keys to be added (or replaced) in the input file     tags: an array of tags associated with this function.  The user is free           to modify this array afterwards, but not to mutate the strings within           the array. (optional)\n\n\n\n"
},

{
    "location": "test/TestSystem.html#TestSystem.not_isapprox-Tuple",
    "page": "Test API",
    "title": "TestSystem.not_isapprox",
    "category": "method",
    "text": "Helper function to compute !isapprox(args..., kwargs) in @test macros because\n\n@test !isapprox(args..., kwargs...)\n\ndoesn\'t work.\n\n\n\n"
},

{
    "location": "test/TestSystem.html#TestSystem.run_testlist",
    "page": "Test API",
    "title": "TestSystem.run_testlist",
    "category": "function",
    "text": "This function runs a test list.  Tests are run in the order they were   loaded into the TestList object.  This implementation handles the case of   several tests sets sharing an input file efficiencly, ie. if several test    functions use the same input file and are placed in the test list    consecutively, the input file will be loaded only once.\n\nA list of tags can be optionally supplied.  In this case, only the tests   that have the specified tags will be run.  If no list is supplied, all tags   are run.\n\nInputs:     testlist: a TestList loaded with functions     prep_func = function used with test functions of type 2 or 3.  It must                 have signature prep_func(fname::String). fname is the                 name of hte input file associated with the test function     tags: an array of tags (optional)\n\n\n\n"
},

{
    "location": "test/TestSystem.html#TestSystem.check_tags-Tuple{Function,Array{String,N} where N}",
    "page": "Test API",
    "title": "TestSystem.check_tags",
    "category": "method",
    "text": "This function checks that all necessary tags are present and throws   an exception is they are not.\n\nCurrently this only verifies that one of LengthTags is present\n\nInputs:     func: the function (only needed to produce error message)     tags: the tags for this function\n\nOutputs: none\n\n\n\n"
},

{
    "location": "test/TestSystem.html#TestSystem-API-1",
    "page": "Test API",
    "title": "TestSystem API",
    "category": "section",
    "text": "CurrentModule = TestSystemModules = [TestSystem]"
},

]}
