# Stabilization

```@meta
  CurrentModule = EulerEquationMod
```

This page describes the stabilization methods used for continuous Galerkin
formulations.
Not all of these effectively stabilized the discretization, and since the
move to discontinuous Galerkin, they may have bitrotted.

## Functions

```@autodocs
  Modules = [EulerEquationMod]
  Pages = ["euler/stabilization.jl", "euler/GLS.jl", "euler/filtering.jl", "euler/artificial_dissipation.jl", ]
```


