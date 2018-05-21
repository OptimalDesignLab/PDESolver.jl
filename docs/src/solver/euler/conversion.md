# Conversion Between Different Variables

```@meta
  CurrentModule = EulerEquationMod
```

It is possible to write the Euler equations in terms of different sets of
variables.
The default is to use conservative variables, but there are
others as well, such as entropy variables.
This page describes how to convert between the different sets of variables.
The functions perform conversions back and forth between conservative variables
and the "other" variables.
The code generally does not support conversion directly between two
different sets of "other variables.
For example, if we call the conservative variables the `q` variables, the 
first set of "other" variables the `w1` variables and the second set of
"other" variables the `w2` variables, it is possible to convert between
`q` and `w1` as well as `q` and `w2` but not between `w1` and `w2`.

The basic set of functions perform the conversion, and these function's names
and in an underscore.
Layered on top of these functions are a safer method of conversion, which
uses the static parameters of the [`ParamType`](@ref) to determine what
variables are stored in `eqn.q` and performs the appropriate conversion.
This is important because if, for example, we are solving the equation using
entropy variables, and the function [`convertToEntropy`](@ref) is called,
it will (correctly) not do the conversion, because the variables are already
the entropy variables.
If [`convertToConservative`](@ref) is called, it will do the proper
conversion from entropy to conservative variables.

The function describe here provide one additional bit of functionality.
Lets say a function needs to convert it current variables to conservative
variables, do some operation, then convert the result back to entropy variables.
Converting to conservative variables is easy, that's what [`convertToConservative`](@ref) does.  The problem is converting back.
The code doesn't necessarily know what variables are stored in `eqn.q`, so
it doesn't know what function to call.  Should it call `convertToW1` or
`convertToW2`?
The solution is [`convertFromNaturalToWorkingVars`](@ref), which converts from
the conservative variables to whatever the variables stored in `eqn.q` are
(as determined by the static parameters of the `ParamType` object).

## Functions

```@autodocs
  Modules = [EulerEquationMod]
  Pages = ["euler/conversion.jl"]
```


