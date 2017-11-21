# Important Keys

```@meta
  CurrentModule = PDESolver
```

The input dictionary contains many keys.  All possible user supplied keys
are listed in [input_vals.txt](https://github.com/OptimalDesignLab/PDESolver.jl/blob/work/src/input/input_vals.txt).
The purpose of this page is to describe the most important keys and use
them to elucidate some of the structure of the code.

The most important key is `physics` which specifies which physics to solve.
The name must be the name the physics module registers with the front end
described in [Registration Functions](@ref).

The `run_type` key specifies what kind of solver to use.  There are a
variety of options, including explicit and implicit time marching methods
for unsteady problems and inexact-Newton methods for steady problems.
Many of the time marching methods can also be used for pseudo-time stepping
of steady problems. See [`call_nlsolver`](@ref) for the complete list of
values.  Physics modules should never invoke a nonlinear solver directly, they
should always use `call_nlsolver()`.

The `jac_method` key is needed if the nonlinear solver computes a Jacobian.
PDESolver supports multiple methods of computing derivatives, including
finite differences and the complex step method.  The functions in the
Nonlinear solvers module are required to support all methods.
This key specifies which method to use.
Complex step is generally recommended, finite difference are used primarily
for verification.

The `order` key specifies the degree of the operator used to discretize
the spatial terms. This is analagous to the degree of the polynomial basis
functions used for finite elements.

`smb_name` is the name of the mesh file.  This is passed directly to
the `AbstractMesh` constructor.  The mesh should already be partitioned into
the correct number of parts (ie. the number of MPI processes PDESolver was
launched with).  For Pumi, the partitioned mesh is comprised of many
numbered files.  The `smb_name` should *not* have the numbers appended to
it, and should be the same for all processes.

`dmg_name` is the name of the geometry file associated with the mesh.

`IC_name` is the name of the initial condition.  See [Registration Functions](@ref) for how to add new initial conditions.  Note that initial conditions are
physics module specific and that getting an error about an unrecognized
initial condition is often the result of specifying the incorrect physics.

`operator_type`: the name of the SBP operator to use.  See [`createSBPOperator`](@ref).



