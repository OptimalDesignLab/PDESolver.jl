# Starting a Simulation

```@meta
  CurrentModule = PDESolver
```

This page describes functions located in the PDESolver module that tie
together the physics modules and the Nonlinear solvers.  
```@autodocs
  Modules = [PDESolver]
  Pages = ["src/startup_func.jl", "src/initialization.jl"]
```


## Physics Module Startup

TODO: does this belong here, or in the physics module section?

Each physics module is required to do some of the setup work needed to
start a simulation.
The functions above facilitate doing so.
In particular, the physics module must

  * read the input dictionary
  * create an `AbstractMesh` and `AbstractSBP`
  * create an `AbstractSolutionData`
  * Load an initial condition
  * Calculate various quantities
  * Invoke a NonlinearSolver
  * Do postprocessing

### Input Dictionary

### Creating Mesh and Operator

### Create an Equation Object

### Load an initial condition

### Various calculations

### Invoke a NonlinearSolver

### Do Postprocessing

The functions here should be used by all physics modules to assist in
creating the `AbstractMesh` and `AbstractSBP` objects and calling a nonlinear
solver.

TODO: document what the physics module needs to implement for startup


