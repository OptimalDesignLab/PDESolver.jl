# Assorted Function and Types

```@meta
  CurrentModule = ODLCommonTools
```

This page contains several types and functions that are used throughout the
physics modules.  Many of these are defined in ODLCommonTools

## Abstract Functor Types
Abstract types are provided for commonly used [Functors](@ref):

```@docs
BCType
BCType_revm
SRCType
FluxType
FluxType_revm
```


## Boundaries and Interfaces

All physics modules need to apply boundary conditions and all DG schemes
and some CG schemes need to do face integrals.
The types and functions described here assist in identifying them:

```@docs
Boundary
Interface
getElementL
getFaceL
show(::IO, ::Boundary)
show(::IO, ::Interface)
```
