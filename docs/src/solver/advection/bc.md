# Advection Boundary Integrals

```@meta
  CurrentModule = AdvectionEquationMod
```

This page describes the functions that compute the boundary integrals.

```@docs
  evalBoundaryIntegrals
  calcBoundaryFlux
  calcBoundaryFlux_nopre
  BCDict
  getBCFunctors
```

## Boundary Conditions
This section describes all boundary conditions currently available

```@autodocs
  Modules = [AdvectionEquationMod]
  Pages = ["advection/bc.jl"]
  Order = [:type]
```

## Numerical Flux Functions

This section lists the numerical flux functions used to impose the boundary
conditions weakly.

```@autodocs
  Modules = [AdvectionEquationMod]
  Pages = ["advection/bc_solvers.jl"]
```
