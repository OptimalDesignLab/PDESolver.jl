# PDESolver User Interface

```@meta
CurrentModule = PDESolver
```


## Invoking the Code
This page describes the API for running the code, as well as utility scripts
that provide simple operation of the code.

The entry point for running the code is the function [`run_solver`](@ref)

```@docs
run_solver
```

Each physics module is required to register itself when `using PhysicsModule`
is run.
When using this API to invoke the code, users should be sure to have loaded the
required physics module before calling `run_solver`.

In order to facilitate running the code, the script `src/startup.jl` loads
all the physics modules and passes the first command line argument as the
input file name

## Registering ICs and BCs

These functions allow users to register new initial conditions and
boundary conditions without modifying the source code

```@docs
registerIC
registerBC
```

## Interactive Functions
These functions are useful for inspecting the state of the front-end.

```@docs
printICNames
printBCNames
printPhysicsModules
```
