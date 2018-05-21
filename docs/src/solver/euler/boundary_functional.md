# Euler Boundary Functional

This page consists of all the functions necessary for computing a boundary
functional along the geometric edges of a mesh for Euler equations. A boundary
functional should ALWAYS be evaluated by calling `evalFunctional` which is the
highest level function.

```@autodocs
  Modules = [EulerEquationMod]
  Pages = ["solver/euler/boundary_functional.jl"]
```
