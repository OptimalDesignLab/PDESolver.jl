# Boundary Integrals Differentiated

```@meta
  CurrentModule = EulerEquationMod
```

This page describes how the boundary integral contribution to the Jacobian
is calculated.


## Entry Point

```@docs
evalBoundaryIntegrals_diff
```

## Functions

```@autodocs
  Modules = [EulerEquationMod]
  Order = [:function]
  Pages = ["euler/bc_diff.jl"]
```


## Boundary Condition Functors

The same functors used to compute the boundary integral are also used
to compute the Jacobian contribution.



