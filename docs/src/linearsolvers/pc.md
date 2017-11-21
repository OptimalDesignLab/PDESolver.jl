# [Preconditioners](@id sec:preconditioners)

```@meta
  CurrentModule = LinearSolvers
```

This section describes the different catagories of preconditioners.

## Type hierarchy

When creating a new preconditioner, it is very important to make it inherit
from the proper supertype.

```@docs
AbstractPC
AbstractPetscMatPC
AbstractPetscMatFreePC
AbstractPetscPC
```


## API

Every preconditioner supports the following functions.  When defining a new
preconditioner, some of the functions must be defined for the new 
[`AbstractPC`](@ref) type, while others are defined automatically based on the
supertype of the preconditioner.

```@docs
calcPC
applyPC
applyPCTranspose
getBasePC
free(::AbstractPC)
```

## Concrete PC Types

The following types are the "Base" PC types refered to by [`getBasePC`](@ref).
Every preconditioner must contain one of these (either directly, in the
`pc_inner` field described by [`AbstractPC`](@ref), of inside another
preconditioner).


```@docs
PCNone
PCNone(::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict)
PetscMatPC
PetscMatPC(::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict)
PetscMatFreePC
PetscMatFreePC(::AbstractMesh, ::AbstractSBP, ::AbstractSolutionData, ::Dict)
```

PetscMatFreePC has an additional function not needed by other preconditioner
types:

```@docs
setPCCtx
```

## Implementing a New Preconditioner

Many of the API functions are automatically generated for matrix-explicit
preconditioners.  Matrix-free preconditioners need to define most of the PC
interface themselves.  A summary of the functions each preconditioner must
implements is:

**Matrix-explicit**

 * [`calcPC`](@ref)

**Matrix-free**

 * [`calcPC`](@ref)
 * [`applyPC`](@ref)
 * [`applyPCTranspose`](@ref)

### First Level PC
  This example shows how to implement a new preconditioner that directly
  contains one of the [Concrete PC Types](@ref) described above.

```
type ExamplePetscMatPC <: AbstractPetscMatPC
  pc_inner::PetscMatPC  # the supertype of pc_inner and ExamplePetscMatPC
                        # must match
end

function ExamplePetscMatPC(mesh::AbstractMesh, sbp::AbstractSBP,
                    eqn::AbstractSolutionData, opts::Dict)

  pc_inner = PetscMatPC(mesh, sbp, eqn, opts)

  return ExamplePetscMatPC(pc_inner)
end


function calcPC(pc::ExamplePetscMatPC, mesh::AbstractMesh, sbp::AbstractSBP,
                eqn::AbstractSolutionData, opts::Dict, ctx_residual, t)


  calcPC(pc.pc_inner)

  pc2 = getBasePC(pc)  # because pc_inner is the base pc, this returns
                       # pc.pc_inner
  # call set_values1!(pc2.Ap, ...) to put values into the preconditioning matrix

  return nothing
end
```

Because `calcPC` is defined and the base PC is one of the concrete PC types,
namely [`PetscMatPC`](@ref), the functions [`applyPC`](@ref) and
[`applyPC`](@ref) are automatically defined for `ExamplePetscMatPC`.
The new preconditioner is now fully usable.

### Multiple Levels of Composition

The same structure as shown in the previous structure can be used to build
a new preconditioner out of any existing preconditioner, not only the [Concrete PC Types](@ref).

```
type OuterPetscMatPC <: AbstractPetscMatPC
  pc_inner::ExamplePetscMatPC  # use the PC define in the previous example
end

function OuterPetscMatPC(mesh::AbstractMesh, sbp::AbstractSBP,
                         eqn::AbstractSolutionData, opts::Dict)

  pc_inner = ExamplePetscMatPC(mesh, sbp, eqn, opts)

  return OuterPetscMatPC(pc_inner)
end


function calcPC(pc::OuterPetscMatPC, mesh::AbstractMesh, sbp::AbstractSBP,
                eqn::AbstractSolutionData, opts::Dict, ctx_residual, t)


  calcPC(pc.pc_inner)  # always call the inner PC calc function first

  pc2 = getBasePC(pc)  # in this case, pc2 is the PetscMatPC inside 
                       # the ExamplePetscMatPC
  # call set_values1!(pc2.Ap, ...) to put values into the preconditioning matrix  # because calcPC(pc.pc_inner) was called already, the values here should
  # be added to those already in pc2.Ap

  return nothing
end
```

As before, this preconditioner is now fully usable.
