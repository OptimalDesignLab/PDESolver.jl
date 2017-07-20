# Utilties

This module contains functions and types that are useful for the solver but
independent the equation being solved.
Additional utility functions are located in the
[ODLCommontools](https://github.com/OptimalDesignLab/ODLCommonTools.jl).
The functions defined in the `Utils` module are useful in the context
of `PDESolver` and depend on the functions and datatypes defined in the
other parts of the solver.
The functions defined in `ODLCommonTools` are more general in nature and
usable independent of `PDESolver`.


```@contents
  Pages = ["io.md"
           "logging.md"
           "projections.md"
           "parallel.md"
          ]
  Depth = 1
```
