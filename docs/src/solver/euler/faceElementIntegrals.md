# [Face Element Integrals](@id sec:euler_face_element_integrals)

```@meta
  CurrentModule = EulerEquationMod
```

This page describes the functions that evaluate the face element integrals
for a single interface.  The functions that loop over all the interfaces
are located on the [face integrals](@ref sec:euler_face_integrals) page.
These integrals require data from all the nodes of the elements rather than
the face nodes as with regular face integrals.

These integrals are used by the entropy stable scheme, and some
of them internally use
a [numerical flux function](@ref sec:euler_flux_functors).  This
flux function must satisfy an entropy property for the resulting scheme to
be entropy stable!  The IR flux function is typically used.


## Functions

```@autodocs
  Modules = [EulerEquationMod]
  Order = [:function]
  Pages = ["euler/faceElementIntegrals.jl", "euler/IR_stab.jl"]
```

## Flux Functors

```@autodocs
  Modules = [EulerEquationMod]
  Order = [:constant, :type]
  Pages = ["euler/faceElementIntegrals.jl"]
```
