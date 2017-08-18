# Boundary Integrals

```@meta
  CurrentModule = EulerEquationMod
```

This page describes the functions that impose boundarya conditions.
The boundary conditions are imposed weakly, using a penalty between
the desired state (the boundary condition value) and the current state
(the solution).


## Entry Points


```@docs
evalBoundaryIntegrals
```

## Functions

```@autodocs
  Modules = [EulerEquationMod]
  Order = [:function]
  Pages = ["euler/bc.jl"]
```

## Boundary Condition Functors

Each functor defines a different type of boundary condition.
Because the equation is solved in the weak form, the functors compute the
flux associated with the penalty that imposes the boundary condition.

```@autodocs
  Modules = [EulerEquationMod]
  Order = [:constant, :type]
  Pages = ["euler/bc.jl"]
```


