# Advection Physics Documentation

This module evaluates the residual for the constant-coefficient advection
equation

$\frac{\partial q}{\partial t} = - \boldsymbol{a} \cdot \frac{\partial q}{\partial \boldsymbol{x}} + S$

where $\boldsymbol{a}$ is the vector of advection velocities in the x, y, and z
directions, $\boldsymbol{x}$ is the vector and x, y, and z directions, and $S$
is the source term.

This weak form of the equation is discretized as described in the [Introduction](@ref index_sbp).

```@contents
  Pages = [ "advection.md"
            "types.md"
            "volume.md"
            "flux.md"
            "bc.md"
            "ic.md"
            "source.md"
            "common.md"
            "adjoint.md"
            "boundary_functional.md"
          ]
  Depth = 1
```
