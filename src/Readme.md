# PDESolver
The PDESolver module is the main module of the code, and it provides a
common interface for solving different physics.  In so doing, it requires
the physics modules to follow a similar organizational structure.  The
physics modules have great flexibility as to their internal structure, but
are required to implement certain functions in order to be usable through
the common interface.

## Physics Module Interface
The required interface for each physics module is fairly modest. It must
  * Define an subtype of `AbstractSolutionData` (see `doc/interfaces.md`)
  * Define and export a function `run_physics(fname::AbstractString)`, where `
    physics` is replaced with the name of the physics for readability
  * Extend the function `evalResidual(mesh::AbstractMesh, sbp::AbstractSBP,
    eqn::AbstractSolutionData, opts::Dict)` with a method that is more specific in
    the type of `eqn` and equally specific in all other arguments.
  * Create (and *not* export) an `Associative{ASCIIString, BCType}` container
    named BCDict that maps boundary conditions names to functors.
  * Create (and *not* export) an `Associative{ASCIIString, Function}` container
    named ICDict that maps initial condition names to functions
  * Call the function `register_physics(physics_name::ASCIIString, mod::Module,
    startup_func::Function)` during module initialization

### `AbstractSolutionData` subtype
The purpose of defining a subtype of `AbstractSolutionData` is two-fold.  First
it allows functions like `evalResidual` to dispatch to the right physics module
based on the type of the `eqn` argument.  Second, it provides a place for the
physics module to store all its data.  The `AbstractSolutionData` type is
required to have certain fields (see `interfaces.md`), and is allowed to have
additional fields that are used internally by the physics module.  In some cases
there may be multiple equation objects with different states and residuals, so
it is imperative that all data that depends on the solution be stored in the
`AbstractSolutionData` object, and never stored in the physics module itself.

It may be beneficial for a physics module to define an abstract subtype of
`AbstractSolutionData` and then to define one or more concrete subtypes.
For example:

```julia
abstract PhysicsSolutionData{Tsol, Tres} <: AbstractSolutionData{Tsol, Tres}

type PhysicsSolutionData1{Tsol, Tres} <: PhysicsSolutionData{Tsol Tres}
# ... fields declared here
end

type PhysicsSolutionData2{Tsol, Tres} <: PhysicsSolutionData{Tsol Tres}
# ... fields declared here
end
```

Functions in the physics module can use `PhysicsSolutionData` for their argument
type, unless there are different versions of the functions for different
implementations of `PhysicsSolutionData` (ie. `PhysicsSolutionData1` or
`PhysicsSolutionData2`), in which case the different methods will be dispatched
to based on which `PhysicsSolutionData` is used.

### `run_physics`
This function takes in the file name or full file path to an input file, reads
it, initializes the 4 main objects, `mesh, sbp, eqn, opts`, and calls a solver
to solve the equation.  It may also do some post-processing (see the Optional
  Structure section).  This is the function that users call to run this physics
 module, and should be exported.

 The PDESolver module provides the function `createMeshAndOperator` to create
 the `mesh` and `sbp` objects.  The function `read_input` in the Input module
 creates the options dictionary.  The real work that `run_physics` has to do
 is create the `eqn` object and call a solver.

 This function returns the 4 main objects, so the user can write scripts that
 automate data collection.

### `evalResidual`
This function evalutes the residual of a steady equation (R(q) = 0), or the
spatial discretization of an unsteady equation (M dq/dt = R(q)), where M is
a mass matrix, q is the solution, and t is time.  A generic fallback is defined
in the PDESolver module.  This method throws an error.  Each physics module must
extend this function with a new method for the `AbstractSolutionData` subtype
defined in the physics module. The fact that `evalResidual` is a single function
(rather than a different function with a different name in each module) is
important for multi-physics simulations.  For example, solving the Advection
and Euler physics in a single run would require constructing an
`AbstractSolutionData` for both modules, and calling `evalResidual` twice:

```
evalResidual(mesh, sbp, eqn_advection, opts)
evalResidual(mesh, sbp, eqn_euler, opts)
```

The fact that evalResidual is a single function that dispatches to the right
physics depending on the type of the `eqn` argument makes writing multi-physics
solvers in a generic way much easier.

### ICDict and BCDict
The input file allows the user to specify the names of the initial and boundary
conditions to use.  Each physics module must maintain a mapping from these
names to the function or functor that evalutes them.  See the documentation
for the functions (defined in the PDESolver module) `register_ic` and
`register_bc` for the requirements on the IC function and BC functor.
The requirement that all physics modules implement these is necessary to
allow users to register new initial conditions and boundary conditions without
modifying the source code of the solver.

### `register_physics`
Every physics module must register itself with the PDESolver module. See the
documentation for the `register_physics` function in the PDESolver module.
This is necessary to implement a unified way of invoking the solver for any
physics.  Despite the names, these do not have to be dictionaries, but
dictionaries are good choices in general.  Whatever their concrete type, they
must be `Associative` containers.

## Optional Structure
It is recommended to have some additional structure. This is not required by
any interface, but it is beneficial because it encourages all physics modules
to have a similar layout and functionality.  This is good both for programmer
productivity and sanity, as well as code maintainability.

The recommended features are
  * An `init(mesh, sbp, eqn, opts)` function that finishes initializing the 4
    main objects.
  * A `checkOptions(opts::Dict)` function that checks for incompatible or
    unsupported options.
  * A `const PhysicsName` ASCIIString which is the name of the physics module
  * A `postproc(mesh, sbp, eqn, opts)` function that does post-processing

### `init` function
The purpose of the `run_physics` function is to create the 4 main objects and
call a solver.  Within `run_physics` it is helpful to have some additional
structure.  Each of the 4 objects can be created and its data fields allocated,
but the final step of initialization usually requires updating the fields of
several objects simultaneously to ensure consistency.  The `init` function
exists to do this final step.  This typically involves getting BC functors out of
BCDict, getting the flux functors (for DG) etc.

### `checkOptions` function
The `read_input` function creates the options dictionary and supplies default
values for all keys that can have reasonable defaults, and is shared by all
physics modules.  The purpose of `checkOptions` is to do any physics-specific
checking of options and throw errors if unsupported or incompatible options
are detected.  The general principle is that the code should avoid silently
doing something that is different than what the user specified, and that it is
better to throw and error and alert the user that something unusual has been
requested than to attempt to continue with the run and possibly do run under
the wrong conditions.  This function is usually called immediately after
`read_input`, so the user gets alerted right away if a problem is detected.

It is recommended that the author of a new physics modules read the list of
input values located in `src/input/input_vals.txt` and throw errors for values
that are not supported, keeping in mind that `read_input` supplies default values
for many of the keys.

### `PhysicsName` constant
The `register_physics` function requires the user to provide a name for the
physics.  This name is used as a value in the input file to specify which physics
to run.  Although not required, putting this value in a global constant might
be useful, because then it will be programmatically accessible in the future

### `postproc` function
It is often helpful for physics modules to do some amount of post-processing.
In particular, writing results to file, calculating errors, and
writing visualization files are commonly used.  Putting all this capability
in a function, and creating input options to specific which types of
post-processing to do is good for code organization.
This function should be called at the end of the `run_physics` function.
