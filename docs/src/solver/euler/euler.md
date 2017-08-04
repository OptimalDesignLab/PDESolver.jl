# Euler Physics Documentation

Describe the equation being solved here

$\frac{\partial q}{\partial t} = - \nabla \cdot F(q) + S$

Where $q$ are the conservative variables and $F$ is the Euler flux.

$q = \begin{bmatrix} \rho, \rho u, \rho v, \rho w, e \end{bmatrix}$

where $\rho$ is density, $u$ is the x velocity, $v$ is the y velocity, $w$
is the z velocity and $e$ is the energy.

The calloricaly perfect ideal gas law is used to close the system (TODO ref calcPressure)

The x-y-z components of the Euler flux are:

$\begin{bmatrix} \rho u & \rho v & \rho w \\ \rho u^{2} + p & \rho u v & \rho u w \\ \rho u v & \rho v^2 + p & \rho v w  \\ \rho u w & \rho v w & \rho w^2 + p  \\ (e + p)u & (e + p)v & (e + p)w \end{bmatrix}$

TODO: describe the physical quantities (include units)

```@contents
  Pages = [ "advection.md"
            "types.md"
            "volume.md"
            "flux.md"
            "bc.md"
            "ic.md"
            "source.md"
            "common.md"
            "conversion.md"
            "flux_functions.md"
            "stabilization.md"
            "adjoint.md"
            "boundary_functional.md"
            "misc.md"
          ]
  Depth = 1
```
