# Linear Solvers

```@meta
  CurrentModule = LinearSolvers
```

The purpose of the LinearSolvers module is to provide a consistant API for
preconditioning and solving linear systems, using both direct and iterative
methods, serial and parallel.  The API is highly composable,
allowing new preconditioners and linear operators to be built on top of
existing ones.  One important feature is the ability to explictly control
when the linear operator and preconditioner are recalculated, which is
important for efficiently solving nonlinear problems.

The interface for using a linear solver (including is component preconditioner
and linear operator) is described on the [LinearSolvers](@ref sec:linearsolvers)
page.

The interface for defining a new preconditioner is described on the
[Preconditioners](@ref sec:preconditioners) page, and, similarly, the
interface for new linear operators is described on the 
Linear Operators](@ref sec:linearoperators) page.  These interfaces must
be defined by every new preconditioner and linear operator, respectively, but
they should not be used directly.  The
[LinearSolvers](@ref sec:linearsolvers) interface should be used instead.

```@contents
  Pages = ["pc.md"
           "lo.md"
           "ls.md"
          ]
  Depth = 1
```
