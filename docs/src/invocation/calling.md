# Running the Code
Before running the code, you must `source` the shell script
`PumiInterface/src/use_julialib.sh`.  This enables Julia to find and use Pumi.

The code takes an input file that defines all options for the solver, including
the which physics to solve, which mesh to use, initial conditions,
boundary conditions, and discretization.  The file `src/input/input_vals.txt`
describes the input file format and all valid options.

## Simple Mode
Once an input file is prepared, the solver can be invoked with

```julia
mpirun -np x julia /path/to/PDESolver/src/startup.jl "input_file_name"
```

This will run the solver in parallel with `x` processes, solving the physics
specified in the input file.

## Advanced Mode
The code also provides a scripting interface.  This enables users to supply
initial conditions and boundary conditions, as well as do advanced
post-processing if desired.

The template script is:

```julia
using PDESolver  # load the interface to the code
using AdvectionEquationMod  # load the desired physics module

# register ICs and BCs here

input_fname = ARGS[1]  # use the first command line argument as the input file
                       # name

# solve the equation defined in the input file
mesh, sbp, eqn, opts = run_solver(input_fname)

# do post-processing using mesh, sbp, eqn, opts here

```

`mesh, sbp, eqn, opts` are the `AbstractMesh` object, the SummationByParts operator, the `AbstractSolutionData`, and the options dictionary, respectively.

The `AbstractMesh` object contains all the data about the mesh, including
coordinates and mapping jacobian information.  The `sbp` object is used for
numerically approximating derivatives. The `AbstractSolutionData` object contains
all the data about the solution and the residual of the equation.  The input
file gets parsed into the options dictionary, which is used to by the rest of
the code.

The PDESolver module exports several functions for use by the scripting
interface.  The first two `registerIC` and `registerBC` allow users to supply
initial conditions and boundary conditions.  See the documentation of the
functions for details.  The functions `printICNames` and `printBCNames` print
the names of all currently registers initial conditions and boundary conditions.

For example a user can, either in a Julia script or interactive (REPL) session:
```julia
using PDESolver
using AdvectionEquationMod

printBCNames(AdvectionEquationMod)
```

And this will print the names of all the boundary conditions currently known
to the Advection Equation module.


