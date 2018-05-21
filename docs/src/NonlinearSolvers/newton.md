# Newton's method

```@meta
CurrentModule = NonlinearSolvers
```

Newton's method is intended to compute updates to some q by solving the following equation.

\begin{equation}
\frac{\partial R(q)}{\partial q} \Delta q = -R(q)
\end{equation}

In the most basic implementation of Newton's method in PDESolver, q corresponds to the solution, 
  and f(q) corresponds to the residual evaluation of the currently selected physics module.

An example of a more sophisticated use of Newton's method is within the Crank-Nicolson timestepper, 
  which adds another layer on top of the physics residual:

\begin{equation}
\frac{\partial g(R(q))}{\partial q} \Delta q = -g(R(q))
\end{equation}

## Features

PDESolver's Newton's method has a wide variety of features. 
It contains the Jacobian calculation routines, which can be performed currently using:

* finite-differencing 
* complex-step

The Jacobian functions can act upon any arbitrary residual.

Additionally, the following matrix forms are supported:

* Julia dense
* Julia sparse
* PETSc sparse
* Matrix-free

The function that performs the Newton iteration is 

```@docs
newtonInner
```

The private data required by Newtons method is stored in the `NewtonData`
object.

```@docs
setupNewton
free(::NewtonData)
getNewtonPCandLO
```

## Newton internals


# NewtonData API

[`NewtonData`](@ref) has a small API used by Newton's method.

```@docs
NewtonData
reinitNewtonData
recordResNorm
recordStepNorm
```

### Linear Operators and Preconditioners

A full set of linear operators and preconditioners are provided for
solving a linear system where the linear operator is the jacobian of the
physics..  These often are a good start for constructing
linear operators for unsteady problems or more advanced methods for solving
steady problems.

```@docs
NewtonMatPC
NewtonMatPC(::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict)
NewtonVolumePC
NewtonVolumePreconditioner(::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict)
NewtonDenseLO
NewtonDenseLO(::PCNone, ::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict)
NewtonSparseDirectLO
NewtonSparseDirectLO(::PCNone, ::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict)
NewtonPetscMatLO
NewtonPetscMatLO(::AbstractPetscPC, ::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict)
NewtonPetscMatFreeLO
NewtonPetscMatFreeLO(::AbstractPetscPC, ::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict)
```

