# NonlinearSolvers Introduction

This module contains both the time-stepping methods and the methods for solving
steady problems.
When solving a CFD problem, these are the functions that drive the code.
The purpose of the [physics modules](@ref sec:physics_modules) is to
*evaluate* the spatial residual.  The purpose of the NonlinearSolvers is to
*solve* the problem.
The steady methods and implicit time-marching methods generally use some kind
of Newton iteration to solve a nonlinear problem.

The first two pages describe the steady and unsteady methods, respectively, and
the final pages describe the Newton implementation used by them.

```@contents
  Pages = [ "steady.md"
            "unsteady.md"
            "newton.md"
            "matrix.md"
            "newton_inner.md"
          ]
  Depth = 1
```
