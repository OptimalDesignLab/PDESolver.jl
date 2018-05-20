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
The types and functions described here assist in identifying them.
PDESolver does not track the faces of elements directly, instead it 
tracks the element the face belongs to and the local face number, that is,
the index of the face in the list of all faces that belong to the element.
Using this representation, every interior face has two representations because
it is part of two elements.

```@docs
Boundary
Interface
getElementL
getFaceL
show(::IO, ::Boundary)
show(::IO, ::Interface)
```

## Functionals

```@docs
AbstractFunctional
AbstractIntegralFunctional
```
