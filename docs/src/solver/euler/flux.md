# [Face Integrals](@id sec:euler_face_integrals)

```@meta
  CurrentModule = EulerEquationMod
```

This page describes the functions that evaluate the face and shared face
integrals.


## Entry Point

```@docs
evalFaceIntegrals
evalSharedFaceIntegrals
```

## Functions

```@autodocs
  Modules = [EulerEquationMod]
  Order = [:function, :constant, :macro]
  Pages = ["euler/flux.jl"]
```

## [Flux Functors](@id sec:euler_flux_functors)
TODO: move this to another file (both code and docs)

```@autodocs
  Modules = [EulerEquationMod]
  Order = [:type]
  Pages = ["euler/flux.jl"]
```
