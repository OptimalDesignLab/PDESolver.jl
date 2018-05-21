# [Face Integrals Differentiated](@id sec:euler_face_integrals_diff)

```@meta
  CurrentModule = EulerEquationMod
```

This page describes the functions that evaluate the face and shared face
integrals.


## Entry Point

```@docs
evalFaceIntegrals_diff
evalSharedFaceIntegrals_diff
```

## Functions

```@autodocs
  Modules = [EulerEquationMod]
  Order = [:function, :constant, :macro]
  Pages = ["euler/flux_diff.jl"]
```

## Flux Functors
TODO: move this to another file (both code and docs)

```@autodocs
  Modules = [EulerEquationMod]
  Order = [:type]
  Pages = ["euler/flux_diff.jl"]
```
