# Testing

PDESolver has a test system to verify the code is working properly.
Every feature of the code should have a test that verifies the correctness
of the feature.  No exceptions!

Tests can take the form of unit tests, which verify the correctness of an
individual function, or system integration tests, which verify several
pieces of code (which were unit tested individually) are working together
correctly.

Some examples of unit tests are:

 * testing that a flux function returns a known value (for cases when the numerical scheme is exact)
 * testing that a two-point numerica flux function returns the analytical flux when the left and right state are the same (the consistency property)

Some examples of system integration tests are:

 * Testing a uniform flow has zero residual
 * Testing convergence rates


## Contents

```@contents
Pages = ["Readme.md", "Travis.md", "TestSystem.md"]
Depth = 1
```




