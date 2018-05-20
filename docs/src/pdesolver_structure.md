# Starting a Simulation

```@meta
  CurrentModule = PDESolver
```

This page describes functions located in the PDESolver module that tie
together the physics modules and the Nonlinear solvers.  
```@autodocs
  Modules = [PDESolver]
  Pages = ["src/startup_func.jl", "src/initialization.jl"]
```

To run a simulation, the following steps must be done

  * read the input dictionary
  * create an `AbstractMesh` and `AbstractSBP`, and `AbstractSolutionData`
  * Load an initial condition
  * Calculate various quantities
  * Invoke a NonlinearSolver
  * Do postprocessing

Some of these steps are handled by the `PDESolver` module, and a few are
handled by the physics module.



### Input Dictionary

The first step is to read the input file.  
Reading input files is split into two parts.  The first part is done by the
[Input](@ref) module, which loads the file from disk and supplies default
values.
The second part is done by the physics-specific registered with [`register_physics`](@ref).
this  function that verifies the physics
module supports the given options (especially checking for combinations of
```@meta
  CurrentModule = EulerEquationMod
```
options that might not be supported).  See, for example, [`checkOptions`](@ref EulerEquationMod.checkOptions).
```@meta
  CurrentModule = PDESolver
```


### Creating Objects


The next thing the physics module needs to do is create `AbstractSBP` and [`AbstractMesh`](@ref), and [`AbstractSolutionData`](@ref)  objects.
The function [`createMeshAndOperator`](@ref) should be used by all physics modules to create the first two.

The `AbstractSolutionData` is created by the physics module itself.
The details of how to do this are left up to the physics module, but the
return values of [`createMeshAndOperator`](@ref) should be used for static
parameter values.

These creation of all three objects are performed by the `_createObjects` functions provided to
[`register_physics`](@ref).

### Load an initial condition

The function [`solvePDE`](@ref) is extended by each physics module to is do the remaining operations.


The details of how to load an initial condition are left up to the physics
module, but the end result must be the initial condition is present in
`eqn.q_vec`.

Physics modules generally use a Dictionary to map IC names (which is how 
ICs are referred to in the input file) to the function that applies the
IC.  See [`registerIC`](@ref) for further description.


### Various calculations

After loading the IC, the options dictionary may require the calculation of
a few quantities.  See [`solvePDE`](@ref) for the list of options keys
that must be supported.

### Invoke a NonlinearSolver

```@meta
  CurrentModule = PDESolver
```


The next step is calling a [Nonlinear Solver](@ref sec:nonlinearsolvers).
The function [`call_nlsolver`](@ref) takes the objects already constructed and
calls the appropriate nonlinear solver.
Currently, there are no strong guarantees about where the solution is stored
(`eqn.q` or `eqn.q_vec`) when this function returns (TODO: fix that).

### Do Postprocessing

```@meta
  CurrentModule = EulerEquationMod
```


The options dictionary may require post-processing be done, for example
calculating the solution error when the analytical solution is known.
Each physics module usually defines a function to do this.
See [`postproc`](@ref) for an example.

```@meta
  CurrentModule = PDESolver
```
