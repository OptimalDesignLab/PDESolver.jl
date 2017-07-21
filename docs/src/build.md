# Building PDESolver

PDESolver is not listed in Julia's `METADATA`, so you will have to clone the
repository and then build the code.

Before installing PDESolver, you should already have the following installed:
  * Julia v0.4
  * C/C++/Fortran compilers
  * An MPI implementation (including MPI compiler wrappers)
  * CMake v3.0.0 or later

If you are on a Debian based system, the following packages will install 
To do so, run the following Julia commands:

```julia
  Pkg.clone("https://github.com/OptimalDesignLab/PDESolver.jl.git")
  Pkg.resolve()  # install all packages in PDESolvers REQUIRE file
  Pkg.build("PDESolver")  # install all packages not in REQUIRE file
```


