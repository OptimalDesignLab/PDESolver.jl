The Mesh:
  File PdePumiInterface.jl (part of PUMI.jl), defines a the mesh type PumiMesh2.  It has several useful fields, including the number of vertices, edges, elements, degrees of freedom, nodes, number of degrees of freedom per node and the number of edges on the boundary of the mesh.  The rest of this file defines some useful functions that take in the PumiMesh2 type.  Use these as much as possible for mesh operations.

  If you need to do something that these functions don't do, look in the PumiInterface.jl and PumiInterface2.jl files, which provide access to many of the functions in Pumi's APF interface.

The Equation
  The file Equation.jl (in PDESolver.jl/equation) defines an EulerEquation type, which constains some useful constants needed to evaluate F1 and F2.  It also contains the vector element stiffness matricies bigQT_xi and bigQT_eta ( by vector element, I mean an element with more than one degree of freedom per node).  These are the QT_zi double under bar and QT_eta double under bar in the weak form.  Because these matrices are stored here, you do not need to use Summation by Parts to evaluate the volume integrals (except for the source term, which we are omitting for now).

Summation by Parts
  You will need to use summation by parts to do the boundary integrals.  the function boundaryintegrate!() is what you want, and requires as arguments: an SBP operator (the variable named operator in the code), an array with elements of type Boundary (an type defined by SBP that contains the element number and the local face number of all faces (edges) on the boundary of the mesh), u, the array of data that is being integrated over the boundary, dxideta, the output from mappingjacobian!(), bndryflux, a function that evalautes the flux at the boundary (you will have to ask Prof. Hicken about this), and res, and array that the result of the boundary integral is stored in.  

  A few more details on the arguments might be helpful.  for bndryfaces, the local face (edge) number means the number of the face (edge) within the element (first, second, third...).  PdePumiInterface.jl provides a way to get this.  u, the array being integrated over the boundary is an array with dimension number of dofs per node by number of nodes per element by number of elements being integrated (its a 3d array).  I recommend doing one element at a time for simplicity.  See the weak form form exactly what the data inside u needs to be.  res needs to be the same shape as u.  After calling boundaryintegrate!, you will need to assemble res into the global solution vector.

Some Extra Functions
  There are three functions at the bottom of euler.jl (in PDESolver.jl/solver/euler), that evaluate terms you will need in the weak form and help assemble them into the global solution vector.  The first two are getF1() and getF2(), which evalaute F1 and F2 (see the weak form derivation).  You have to pass these function the arguments: a mesh object (like PumiMesh2), and SBP operator, and equation (of type EulerEquation), u0, the solution at the previous timestep, and element number, and a vector of length 12, and it will populate the vector with either F1 or F2, for getF1() and getF2() respectively.  F1 and F2 are used many times in the weak form, so having a function to evaluate them should make things easier.  The elements in the vector get overwritten, so its ok if there is still data in there from a previous iteration or something.

  The third function is assembleU(), and it has two methods (different sets are arguments).  The first one takes a vector of length 12, and element number, and the global solution vector and assemble the vector into it.  The elements in the vector must be ordered as follows: [the solution for node 1 dofs 1 through 4, the solution for node 2 dofs 1 through 4, the solution for node 3 dofs 1 through 4].  getF1 and getF2 provide F1 and F2 in this order, so as along as you keep things in the original order, these assembleU() will work properly.  The element argument (an integer) specifies which element the nodes belong to.

  The second method of assembleU() takes in a vector of length 3, an element number, a component number, and the global solution vector.  The method assembles the solution for a particular dof of all the nodes of a particular element into the global solution vector.  The element argument specifies which element, the component argument specifies which degree of freedom on the node (1 through 4).  I'm not sure if this function will be needed or not, but its here just in case.

Startup.jl
  This script runs the simulation.  It loads a mesh (creates the PumiMesh2 object), creates the SBP operator, creates the vectors to hold the solution at the current and previous timestep, and calls a function to apply an initial condition to the solution at the previous timestep, and calls the rk4 timestepping function.

The layout of the files is
PDESolver/Equation : Equation.jl
PDESolver/rk4  : rk4.jl
PDESolver/solver : Startup.jl
PDESolver/euler : euler.jl

