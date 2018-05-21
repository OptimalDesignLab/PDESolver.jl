# Documentation of the PDESolver Physics Module Interface

```@meta
CurrentModule = PDESolver
```

The `PDESolver` Module provides some structure for physics modules to plug into.
See [Interfaces in PDESolver](@ref) for details.
Each physics module should create a new subtype of `AbstractSolutionData`.
This physics module should then extend the required methods, specializing
one of the arguments with the new (physics-specific) type.
This structure allows the API for the physics module to be visible to
other parts of the code (such as `NonlinearSolvers`), while being defined 
within the physics module.

One exception to this pattern of defining generic functions and then
specializing one argument is the initialization functions.
These functions run before the `AbstractSolutionData` object is created,
so they can not dispatch to different methods.
Instead, each physics module has to register itself and provide the
required initialization functions. See the [Registration Functions](@ref)
section.


## Functions to be Extended

All the functions in this section should be extended by each physics
module with a new method, specializing the type of the `eqn` argument.

```@autodocs
Modules = [PDESolver]
Pages = ["src/interface.jl"]
```

## Additional Functions

These functions are available for all physics module, but do not need to
be extended with a new method. They usually boil down to calling one of
the functions in the previous section.

```@autodocs
Modules = [PDESolver]
Pages = ["src/interface2.jl"]
```


## Registration Functions

These function provide the basic API for registering and retrieving
physics modules with the `PDESolver` module.

```@docs
register_physics
retrieve_physics
```

See also [`registerIC`](@ref) and [`registerBC`](@ref).
