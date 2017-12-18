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


## Physics Module Startup

Each physics module is required to do some of the setup work needed to
start a simulation.
The functions above facilitate doing so.
In particular, the physics module must

  * read the input dictionary
  * create an `AbstractMesh` and `AbstractSBP`
  * create an `AbstractSolutionData`
  * Load an initial condition
  * Calculate various quantities
  * Invoke a NonlinearSolver
  * Do postprocessing

```@meta
  CurrentModule = EulerEquationMod
```


Physics modules should define a function called `run_physics` (ex. [`run_euler`](@ref)) that does all these operations (by calling other functions within the
physics module) and returns the mesh, sbp, eqn, and opts objects.

### Input Dictionary

The first thing a physics module must do is read the input file.  
Reading input files is split into two parts.  The first part is done by the
[Input](@ref) module, which loads the file from disk and supplies default
values.
The section part is done by the physics function that verifies the physics
module supports the given options (especially checking for combinations of
options that might not be supported).  See, for example, [`checkOptions`](@ref EulerEquationMod.checkOptions).

### Creating Mesh and Operator

```@meta
  CurrentModule = PDESolver
```


The next thing the physics module needs to do is create `AbstractSBP` and [`AbstractMesh`](@ref) objects.
The function [`createMeshAndOperator`](@ref) should be used by all physics modules to
do this.

### Create an Equation Object

Next, the physics module must create its [`AbstractSolutionData`](@ref) object.
The details of how to do this are left up to the physics module, but the
return values of [`createMeshAndOperator`](@ref) should be used for static
parameter values.

```@meta
  CurrentModule = EulerEquationMod
```


The creation of the mesh, sbp, and equation object are usually combined into
a single function called [`createObjects`](@ref).

### Load an initial condition

A function called `solve_physics` (ex. [`solve_euler`](@ref)) is created by
the physics module do the operations described in this section and the next two.

```@meta
  CurrentModule = PDESolver
```


The details of how to load an initial condition are left up to the physics
module, but the end result must be the initial condition is present in
`eqn.q_vec`.

Physics modules generally use a Dictionary to map IC names (which is how 
ICs are referred to in the input file) to the function that applies the
IC.  See [`registerIC`](@ref) for further description.


### Various calculations

```@meta
  CurrentModule = EulerEquationMod
```


After loading the IC, the options dictionary may require the calculation of
a few quantities.  See [`solve_euler`](@ref) for the list of options keys
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
  CurrentModule = EulerEquationMod
```
