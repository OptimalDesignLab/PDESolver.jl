# [Common Functions](@id sec:euler_common_funcs)

```@meta
  CurrentModule = EulerEquationMod
```

The functions here evaluate some function at a single point.  This is
useful for defining manufactured solutions, initial conditions, and boundary
conditions.  In fact, is is generally recommended that every initial condition
and Dirichlet boundary condition create a function here to compute
the state.  For boundary conditions, all the functor need to is call this
function to get the state and call a numerical flux function.

```@autodocs
  Modules = [EulerEquationMod]
  Pages = ["euler/common_funcs.jl"]
```
