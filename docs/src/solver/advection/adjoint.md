# Advection Equation Steady Adjoint

PDESolver currently has the capability to compute the steady adjoint of a
boundary functional. Recall the adjoint equation as

$\frac{\partial \mathcal{L}}{\partial q} = \frac{\partial \mathcal{J}}{\partial q} + \psi^T \frac{\partial \mathcal{R}}{\partial q} = 0$

where, $\mathcal{L}$ is the Lagrangian for functional $\mathcal{J}$ and $q$ is
the solution variable. The adjoint can be computed by calling the function
`calcAdjoint`, which has been described below.

```@autodocs
  Modules = [AdvectionEquationMod]
  Pages = ["solver/advection/adjoint.jl"]
```
