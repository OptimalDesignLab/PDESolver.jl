# Documentation of the PDESolver Physics Module Interface

```@meta
CurrentModule = PDESolver
```

The `PDESolver` Module provides some structure for physics modules to plug into.
See [Interfaces in PDESolver](@ref) for details.
Each physics module should extend [`evalResidual`](@ref) with a new method to
which takes the `AbstractSolutionData` defined by the physics module as
an argument.
This structure allows `evalResidual` to be visible to the `NonlinearSolver`
module while being defined by the physics modules.

## Evaluating the Physics

```@docs
evalResidual
evalHomotopy
evalJacobian
evalHomotopyJacobian
```

## Registration Functions

These function provide the basic API for registering and retrieving
physics modules with the `PDESolver` module.

```@docs
register_physics
retrieve_physics
```

See also [`registerIC`](@ref) and [`registerBC`](@ref).
