# [Overview of Physics Modules](@id sec:physics_modules)

## `AbstractSolutionData` and Physics Module Implementation

This document describes some best practices for implementing a physics module.
These practices are not required, but have proven to be useful for producing
organized, readable, and reusable code.


## Levels of Functions

It is useful to divide functions into 3 catagories, high, mid, and low level
functions.  The purpose of high level functions is to decide which method of
performing an operation should be used and call other functions to do it.
For example, the Euler physics modules has `evalVolumeIntegrals` and
`evalBoundaryIntegrals` as high level functions.  There are several different
ways of calculating both the volume and boundary integrals.  The options
dictionary is used to decide what mid level function to call.  Each mid
level function implements a different way of doing the calculation.

The purpose of mid level functions is to loop over the mesh and call a
low level function for each node.  For example, the function `getEulerFlux`
loops over the nodes of the mesh and calls a function to calculate the Euler
flux at each node.  Mid level function names usually start with `get` to indicate
that their purpose is to calculate some quantity but don't do the calculation
themselves.

Low level functions calculate a quantity at a node.  For example, `calcEulerFlux`
calculates the Euler flux at a single node.  Low level function names usually
start with `calc` to indicate that they perform a specific calculation.
Often, different discretizations use the same structure of loops, but do a
slightly different calculation at each node.  Low level functions are called
inside the innermost loop of the code, so it would be too expensive to have
if statements to select which low level function to call, so various tricks
involving Julia's multiple dispatch system are used to get the compiler to
decide which low level function to call.  These will be described later in
this document.

It is often useful to dispatch to low level functions based on `Tdim` and
`var_type`.  For this reason the Euler equation implementation of `AbstractParamType`
is
```
type ParamType{Tdim, var_type, Tsol, Tres, Tmsh} <: AbstractParamType{Tdim}
```

The other static parameters are necessary because `ParamType` has fields of
those datatypes.

## `AbstractSolutionData` implementation

Each physics module should define and export a subtype of `AbstractSolutionData{Tsol, Tres}`.
The implementation of `AbstractSolutionData{Tsol, Tres}` must inherit the `Tsol`
and `Tres` static parameters, and may have additional static parameters as well.
It may also be helpful to define additional abstract types within the physics
module to provide different levels of abstractions.
For example, the Euler physics module defines:

```
abstract AbstractEulerData{Tsol, Tres} <: AbstractSolutionData{Tsol, Tres}
abstract EulerData {Tsol, Tdim, Tres, var_type} <: AbstractEulerData{Tsol, Tres}
type EulerData_{Tsol, Tres, Tdim, Tmsh, var_type} <: EulerData{Tsol, Tdim, Tres, var_type}
```

The first line is effectively just a name change and may not be necessary.
The second line adds the static parameters `Tdim`, and `var_type` while
inheriting the `Tsol` and `Tres` types from `AbstractEulerData`.
`Tdim` is the dimensionality of the equation, `Tres` is the datatype of the
residual variables, and `var_type` is a symbol indicating whether the equation
is being solved with conservative or entropy variables.
The third line defines a concrete type that implements all the features required
of an `AbstractSolutionData`, and adds a static parameter `Tmsh`, the datatype
of the mesh variables.  
The additional static parameter is necessary because one field of `EulerData_`
has type `Tmsh`.
Note that there could be multiple implementations of `AbstractSolutionData` for
the Euler equations, perhaps with different fields to store certain data or not.
All these implementations will need to have the static parameters `Tsol`,
`Tdim`, `Tres`, and `var_type`, so `EulerData` is defined as an abstract type,
 allowing all implementations to inherit from it.
All high level functions involved in evaluating the residual will take in an
argument of type `EulerData`.
Only when low level functions need to dispatch based on which implementation is
 used would it take in an `EulerData_` or another implementation.



### Variable Conversion

Some equations can be written in different variables, and need to convert
between them.  To do this, it is
`function convertFromNaturalToWorkingVars{Tsol}(params::ParamType{2, :var_type},
               qc::AbstractArray{Tsol,1}, qe::AbstractArray{Tsol,1})`

that converts from the "natural" variables in which to write an equation to
some other set of variables at a node.  For the Euler equations, the "natural"
variables would be the conservative variables, and one example of "other"
variables would be the entropy variables.

It is also sometimes useful to define the opposite conversion, ie. from
the working variables to the natural variables.



## Input Options

Many of the components of PDESolver have different options that control how they
work and what they do.
In order to  provide a unified method of specifying these options, an dictionary
 of type `Dict{String, Any}` is read in from a disk file.
This dictionary (called `opts` in function signatures), is passed to all high
and mid level functions so they can use values in the dictionary to determine their
 control flow.
Low level functions need to be extremely efficient, so they cannot have
conditional logic, therefore they are not passed the dictionary.
Note that retrieving values from a dictionary is very slow compared to accessing
the fields of a type, so all values that are accessed repeatedly should be stored
 as the field of a type.

## Functors

Functors are a trick used to get Julia's dispatch system to make decisions at
compile time rather than runtime.  This is particularly useful for boundary
conditions, where the list of mesh faces that have boundary conditions applied
is determined at runtime, but having conditional statements that execute for
every node on the mesh boundary would be slow.  Instead a construct is used
as follows:

```julia
mutable struct myBC <: BCType  # create a singleton type
end

function call(obj::myBC, q::AbstractVector, bndryflux::AbstractVector)
  # calculate boundary flux here
end
```

This defines a datatype and adds a method to the `call` function for that type.
The call function is what makes a datatype callable like a function.  This
method is called as follows:

```julia
functor = myBC()  # construct and object of type myBC
q = rand(4)
bndryflux = zeros(4)
functor(q, bndryflux)  # the Julia compiler turns this into call(functor, q, bndryflux)  
```

The way this is used for boundary conditions is through a two level construct
where an outer function passes a functor to an inner function.  Julia's JIT will
generate a method of the inner function that is specialized to the functor (this
is why it is important that the functor is a datatype).  For example:

```
function getBCFluxes(mesh, sbp, eqn, opts)

  for i=1:mesh.numBC  # loop over different boundary conditions
    functor_i = mesh.bndry_functor[i]  # get the functor for this boundary condition
    start_index = mesh.bndry_offsets[i]
    end_index = mesh.bndry_offsets[i+1] - 1
    # get data for boundary faces start_index:end_index

    calcBoundaryFlux(functor_i, data for boundary faces start_index:end_index)
  end
end  # end function

  function calcBoundaryFlux(functor_i::BCType, data for boundary faces start_index:end_index)
    for i=1:length(start_index:end_index)
      for j=1:num_nodes_on_face
        # get data for this boundary face node
        functor_i(data for this boundary face node)
      end
    end

  end  # end function
```

The benefit of this arrangement is that `mesh.numBC` different version of
calcBoundaryFlux get compiled, one for each functor, and each version knows
about the `call` method that was defined for the functor it is passed.  This two level
scheme allows the compiler to make all the decisions about what function to call
(ie. the `call` method of the functor), avoiding any conditional logic at runtime

This idea is also applicable to the flux functions used by DG methods.

## Creating a Manufactured Solution

Using the Method of Manufactured Solutions is an effective way to verify the
correctness of the code.  A guide to deriving the solution can be found [here](http://prod.sandia.gov/techlib/access-control.cgi/2000/001444.pdf).
### Solution and Derivative Functions

For simple equations such as linear advection, the general approach is to define
a function that evalutes the manufactured solution and its derivatives at a
given point, and from that create the required source term, initial condition,
and boundary condition.
Most physics modules have a file called `common_funcs.jl` where the solution
and derivative evaluation functions go.
Make sure the new function you create has the same signature as the other
functions in the file.
It is often useful to create different methods of the same function when
creating related manufactured solutions for 2D and 3D.
If the signature requires an [`AbstractParamType`](@ref), its static parameter
can be used to distinguish the methods.
Creating two methods of the same function, rather than two different functions,
will make it possible to define source terms, boundary conditions, and initial
conditions that work for both 2D and 3D.

For complicated equations such as Euler, it is tedious to construct the 
source term, boundary condition, and initial condition from the solution and
its derivatives.
In this case, it is better to use a symbolic math program to generate expressions for the source term, boundary condition, and initial condition directly.
Some symbolic math programs have the option to generate C or Fortran code, which
can be easily converted to Julia code.

### Source Term

To create the source term functor, locate the file where the source terms
are defined for the physics modules, usually called `source.jl` and create a
new functor and associated `call()` method (see description of functors above).
Make sure the functor object is a subtype of [`SRCType`](@ref) and the `call()`
method has the same signature (except for the first argument) as the other
call methods in the file.
The name of the functor should be `SRCfoo`, where `foo` is the name of the 
source term.

For simple equations such as linear advection, the body of the `call()`
function should construct the source term from the functions in `common_funcs.jl`.
For more complicated equations, the code that evalutes the source term at a
given point should be placed in the body of the `call()` function directly.

Note that the purpose of this function is to evalute the value of the source
term, not to do integration of any kind.

Once the functor is created, it should be added to the list of source terms
(usually a Dictionary located at the bottom of the file where the source terms
are defined).
Consult the physics module documentation for details.


### Boundary Condition

Construction of a boundary term is similar to construction of a source term.
Locate the file where the boundary conditions are defined, usually `bc.jl`, and
add a new functor.
Make sure the functor type is a subtype of [`BCType`](@ref) and the `call()`
method has the same signature (except for the first argument) as the other
`call()` methods in the file.
The naming convention for BC functors is `fooBC`, where `foo` is the name of the 
boundary condition.
The body of the `call()` method should evalute the flux caused by the imposition
of the boundary condition (because boundary conditions are imposed weakly).
This is typically accomplished by calculating the boundary condition state
and then calling a numerical flux function with both the current state and
the boundary state.

For simple equations, the boundary state should construct the boundary state
by calling the functions in `common_funcs.jl`.
For more complicated equations, the code to evalute the boundary state should
be contained in the `call()` method body.

Once the functor is created, it should be added to the list of boundary
conditions, usually a dictionary located at the bottom of the file where the
boundary conditions are defined.
Consults the physical module documentation for details.


### Initial condition

Initial conditions are a bit different than boundary conditions and source
terms because they do not use functors (functors are unnecessary because ICs
are evaluated infrequently).
Locate the file where initial conditions are defined, typically `ic.jl`, and
create a new function with the same signature as the the other functions in
the file.
This function should loop over all elements in the mesh, every node on the
element, and use `mesh.dofs` to assign the solution to proper indices of
the supplied vector.
The naming convention for IC functions is `ICfoo`, where `foo` is the name of
the initial condition.

For simple equations, the solution should be calculated using the functions
in `common_funcs.jl`, otherwise it should be calculated in the initial
condition function.

Initial condition functions are used to calculate errors during post-processing,
so it is important for the initial condition function to evaluate the solution
at the proper time for unsteady problems.

After the initial condition function is created, it should be added to the list
of initial conditions, usually a dictionary at the bottom of the file where
the initial conditions are defined.
See the physics module documentation for details.

# Initialization of a Simulation

This section lists an outline of how a simulation gets launched
After step 4, the procedure becomes a bit more complicated because there are optional steps.
Only the required steps are listed below.

1. The options dictionary is read in.  Default values are supplied for any key that is not specified, if a reasonable default value exists.
2. Second, the `sbp` operator is constructed.
3. The `mesh` object is constructed, using the options dictionary and the `sbp` operator.  Some of the options in the dictionary are used to determine how the mesh gets constructed.  For example, the options dictionary specifies what kind of mesh coloring to do.
4. The `eqn` object is constructed, using the `mesh`, `sbp`, and `opts` objects
5. The physics module `init` function is called, which initializes the physics module and finishes any initialization that `mesh` and `eqn` objects require.
6. The initial condition is applied to `eqn.q_vec`.
7. A nonlinear solver is called.  Which solver is called and what parameters it uses are determined by the options dictionary.
8. Post-processing is done, if required by the options dictionary.
