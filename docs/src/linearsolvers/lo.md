# [Linear Operators](@id sec:linearoperators)

```@meta
  CurrentModule = LinearSolvers
```

This section describes the different catagories of linear operators, including
the API and the requirements for defining new linear operators.
Note that whenever a linear solver is used, the
[Linear Solver](@ref sec:linearsolvers) interface should be used rather than
the Linear Operator interface described here.

## Type Hierarchy

When creating a new linear operator, it is very important to make it inherit
from the proper abstract type.

```@docs
AbstractLO
AbstractDenseLO
AbstractSparseDirectLO
AbstractPetscMatLO
AbstractPetscMatFreeLO
```

Several typealiases are also available:

```@docs
MatExplicitLO
PetscLO
DirectLO
```

## API

Every linear operator supports the following functions.  When defining a new
linear operator, some of the functions must be extended with new methods,
while others will be created automatically based on the supertype of the
linear operator.

```@docs
calcLinearOperator
applyLinearOperator
applyLinearOperatorTranspose
getBaseLO
getBaseObject
free(::AbstractLO)
```

## Concrete LO Types

```@docs
DenseLO
DenseLO(::PCNone, ::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict)
SparseDirectLO
SparseDirectLO(::PCNone, ::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict)
PetscMatLO
PetscMatLO(::PetscMatPC, ::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict)
PetscMatFreeLO
PetscMatFreeLO(::PetscMatPC, ::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict)
```

`PetscMatFreePC` has an additional function not needed by other linear operator
types:

```@docs
setLOCtx
```

## Implementing a New Linear Operator

Many of the API functions are automatically generated for matrix-explicit
linear operators.  Matrix-free linear operators need to define most of the
LO interface themselves.  A summary of the functions each preconditioner
must implement is:

**Matrix-explicit**

 * [`calcLinearOperator`](@ref)
 
**Matrix-free**

 * [`calcLinearOperator`](@ref)
 * [`applyLinearOperator`](@ref)
 * [`applyLinearOperatorTranspose`](@ref)


### Example

Linear operators use the same composition structure as preconditioners, see
the [Implementing a New Preconditioner](@ref) section before proceeding.


```
# Assume the ExamplePetscMatLO type has already been defined

type OuterPetscMatLO <: AbstractPetscMatPC
  lo_inner::ExamplePetscMatLO
end


function OuterPetscMatPC(pc::OuterPetscMatLO, mesh::AbstractMesh,
                         sbp::AbstractSBP, eqn::AbstractSolutionData, opts::Dict)

  lo_inner = ExamplePetscMatLO(pc, mesh, sbp, eqn, opts)
  return OuterPetscMatPC(lo_inner)
end

function calcLinearOperator(lo::OuterPetscMatLO, mesh::AbstractMesh,
                            sbp::AbstractSBP, eqn::AbstractSolutionData,
                            opts::Dict, ctx_residual, t)

  # have inner LO do its computation
  calcLinearOperator(lo.lo_inner, mesh, sbp, eqn, opts, ctx_residual, t)

  lo2 = getBaseLO(lo)
  # use set_values1!(lo2.A, ...) to modify the matrix A

  return nothing
end
```

[`applyLinearOperator`](@ref) and [`applyLinearOperatorTranspose`](@ref)
are defined automatically for all [`AbstractPetscMatLO`](@ref) types
(using [`getBaseLO`](@ref)).
  
