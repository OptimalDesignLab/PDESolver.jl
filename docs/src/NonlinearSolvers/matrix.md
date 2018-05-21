# Matrix Interface

One of the challenged of NonlinearSolvers is supporting several matrix type.
The currently supported types are

 * Array (a Julia dense matrix, column major)
 * SparseMatrixCSC (a Julia implementation of compressed sparse column format)
 * PetscMat (Julia wrappers of the PETSc library)

Unfortunatly, these matrix types support somewhat different APIs that have
little overlap.
This is a problem, because it becomes difficult to write efficient,generic code, ie. code that works on all three 
matrix types (without writing specialized code for each one) and performs
nearly as well as an implementation specialized for a single implementation.

To solve this problem, the Petsc wrappers implement an API that is much
smaller than that of `Array`, but efficiently maps onto the capabilities of
the different array types.
This API should be used in NonlinearSolvers to ensure all functions support
the different matrix types

## API

TODO: when the Petsc docs are finished, document them here
