# PDESolver Structure

```@meta
  CurrentModule = PDESolver
```

This page describes functions located in the PDESolver module that tie
together the physics modules and the Nonlinear solvers.  
The functions here should be used by all physics modules to assist in
creating the `AbstractMesh` and `AbstractSBP` objects and calling a nonlinear
solver.

TODO: document what the physics module needs to implement for startup

```@autodocs
  Modules = [PDESolver]
  Pages = ["src/startup_func.jl", "src/initialization.jl"]
```
