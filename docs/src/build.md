# Building PDESolver

This page describes how to build the PDESolver package.
Julia has two kinds of dependencies, system packages and Julia packages.
System packages can either be installed using a Linux package manager
(such as `apt-get` on Debian-based system) or built from source.
Julia packages are install using the Julia package manager.

## System Dependencies

Before installing PDESolver, you should already have the following installed:

  * Julia v0.4
  * C/C++/Fortran compilers
  * An MPI implementation (including MPI compiler wrappers)
  * CMake v3.0.0 or later
  * BLAS and LAPACK

The build process for Julia itself can be found [here](https://github.com/JuliaLang/julia/tree/release-0.4).

If you are on a Debian based system, the following command will install
the remaining dependencies:

```
  sudo apt-get install build-essential gfortran cmake libblas-dev liblapack-dev mpich
```


## PDESolver Installation

PDESolver is not listed in Julia's `METADATA`, so you will have to clone the
repository and then build the code.

After installing the system [System Dependencies](@ref), run the following Julia commands to install the Julia dependencies and and PDESolver itself:


```julia
  Pkg.clone("https://github.com/OptimalDesignLab/PDESolver.jl.git")
  Pkg.resolve()  # install all packages in PDESolvers REQUIRE file
  Pkg.build("PDESolver")  # install all packages not in REQUIRE file and
                          # build PDESolver
```
This will install PDESolver and all Julia dependencies into the directory
specified by the `JULIA_PKGDIR` environmental variable.
If this variable is not specified, it will install to `~/.julia/v0.4/PDESolver`.

If there are errors building any Julia dependencies, see the
[Installation](@ref) page for methods of installing particular versions of the
dependencies.


After installation it is recommended to run the test suite.
To do so, run the following commands in the terminal:

```
  cd /path/to/pdesolver # (for example ~/.julia/v0.4/PDESolver)
  cd ./test
  ./runtests_fast.sh
```

If the tests complete without error, then the package is properly installed.

TODO: link to examples page
  

