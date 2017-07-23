# PDESolver Documentation
Welcome to the PDESolver documentation.  These documents provide an overview
of PDESolver, a Julia based solver for partial differential equations.
This page will describe the form of the equations and the Summation-By-Parts
operators used to discretize them.  The [Table of Contents](@ref) links to
the components of PDESolver that calculate each term.

## Form of the Equation

PDESolver discretizes equation in the form:

$\frac{\partial q}{\partial t} = \mathcal{R}(q, t)$

where $q$ is the variable being solved for, $\mathcal{R}(u, t)$ is called the
residual and $t$ is time.
The residual contains all the spatial derivatives.  For example, the residual
for the 1D advection equation is $ a\frac{\partial q}{\partial x}$, where
$a$ is the advection velocity.
The code is capable of solving both unsteady problems, in the form shown above,
and steady problem of the form

$\mathcal{R}(q) = 0.$

The code is structured such that the physics modules are responsible for 
evaluating the residual and the Nonlinear Solvers are responsible for the
time term in the unsteady case, or solving the nonlinear rootfinding problem
for steady cases.

## [Summation-by-Parts operators](@id index_sbp)

Summation-by-Parts (SBP) operators are used to discretize the spatial derivatives in
the residual. In particular, we use the multi-dimensional Summation-by-Parts
operators first defined in 

```
  Hicken, J.E., Del Rey Fernandez, D.C, and Zingg, D.W., "Multi-dimensional 
  Summation-by-Parts Operators: General Theory and Application to 
  Simplex Elements", SIAM Journal on Scientific Computing, Vol, 38, No. 4, 2016,
   pp. A1935-A1958
```

See the introductor PDF (to be posted shortly) for an introduction to the
operators.

## Table of Contents
```@contents
Pages = ["invocation/calling.md",
         "solver/Readme.md",
         "solver/advection/advection.md",
         "solver/euler/euler.md",
         "solver/simpleODE/simpleODE.md",
         "NonlinearSolvers/nonlinearsolvers.md",
         "input/input.md",
         "Utils/Utils.md"
        ]
Depth=1
```
