# Jacobian calculation

```@meta
CurrentModule = NonlinearSolvers
```


Special mention of calcJacobianComplex.

One of the most involved parts of Newton's method is forming the Jacobian.
The NonlinearSolvers module contains the function [`physicsJac`](@ref) to
compute the Jacobian of the physics $\frac{\partial R(q)}{\partial q}$, where
$R$ is [`evalResidual`](@ref).  This function can be used by other methods
as a starting point for computing the residual of $g(R(q))$ described in
[ Newton's method](@ref) page.


```@docs
physicsJac
```

`physicsJac` calls one of several functions to compute the Jacobian.  Users
should not call these functions directly, they should use `physicsJac`.
There are two approaches to computing the Jacobian: using finite differences/complex step and calculating the entries of the matrix explicitly.
The latter approach is much faster than the former, however the former does
not require the physics module differentiate [`evalResidual`](@ref).

## Finite Differencing/Complex Step

The functions [`calcJacFD`](@ref) and [`calcJacobianComplex`](@ref) compute
the Jacobian as a dense matrix, one column at a time.  This is very slow
and only useful for debugging small cases.
A more efficient method is implemented in [`calcJacobianSparse`](@ref),
which uses a graph coloring approach to perturb the solution at several
nodes simultaneously. 

```@docs
calcJacFD
calcJacobianComplex
calcJacobianSparse
applyPerturbation
assembleElement
calcJacCol
```

## Explicit Calculation

The most efficient method to compute the Jacobian is to explicitly compute
all its entries and assemble them into a matrix.
The [`evalJacobian`](@ref) function must be extended by each physics module
for this to work.
This function is passed all the same arguments as [`evalResidual`](@ref) plus an additional [`_AssembleElementData`](@ref) object which is used to assemble
the contribution of each element or interface into the matrix.

```@docs
_AssembleElementData
_AssembleElementData(::AbstractMatrix, ::Any, ::Any, ::Any, ::Any)
NullAssembleElementData
assembleElement(::_AssembleElementData, ::AbstractMesh, ::Integer, ::Array{Float64, 4})
assembleInterface
assembleSharedFace
assembleBoundary
```
