# Initial Conditions

```@meta
  CurrentModule = EulerEquationMod
```

This page describes the functions that apply initial conditions.
The initial condition is loaded into the solution vector (not the 3D array).
Users should call the appropriate function if they want to transfer the
solution from the vector to the 3D array.

It is important to note that the initial conditions functions *do* depend
on the current time `t`.  In general, an initial condition function is
really a function that computes an exact solution at a given time `t`.
Applying an initial condition is a special case of this (because `t = 0`).
Having the initial condition functions depend on time is used to calculate
errors at the end of a simulation (for cases where the exact solution is
known).

Often, the initial condition function will call a function in
[common funcs](@ref sec:euler_common_funcs) to evalute the exact solution
at each node on the mesh.

Unlike boundary conditions, which are called many times during the execution
of the code, initial conditions are called infrequently, so we do not
create functors for them.

## Accessing Initial Conditions

```@autodocs
  Modules = [EulerEquationMod]
  Order = [:constant]
  Pages = ["euler/ic.jl"]
```

## Initial Condition Functions

TODO: write  a macro for calling common funcs repeatedly

```@autodocs
  Modules = [EulerEquationMod]
  Order = [:function]
  Pages = ["euler/ic.jl"]
```


