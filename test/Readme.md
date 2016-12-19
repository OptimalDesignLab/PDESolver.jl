# Test System
The Ticon test system is composed of functions that test different aspects
of the code.  The test system allows tags to be associated with each test
function, and selective running of tests based on their tags.
Currently, each physics module has its own set of tags and test list.
The serial tests for a single physics module  are run as follows:

```julia
  julia ./runtest.jl  TAG1 TAG2 ...
```

The parallel tests for all physics modules are run with:

```julia
  ./runtests_parallel.sh TAG1 TAG2...
```

and all tests (serial + parallel) for all physics modules are run with

```julia
  ./runtests.sh TAG1 TAG 2
```

If no tags are specified, all tests are run.
The list of currently defined tags can be found at the top of runtests.jl.
Note that the serial tests must  be run before the parallel tests.


## Adding Tests
Each physics module maintains a `TestSet` object which contains a list of
tests and their associated tags.  TestSet.jl provides an API for adding
tests to the list.

All tests must be enclosed in functions.  These functions must fit into
one of three catagories:
  1. zero argument functions
  2. functions that take the arguments `(mesh, sbp, eqn, opts)` as produced by
     an input file which already exists
  3. functions that take the arguments `(mesh, sbp, eqn, opts)` as produced by
     modifying an existing input file

These catagories correspond to `add_funcs1!`, `add_funcs2!`, and `add_funcs3!`
in TestSet.jl.

The function `run_testlist` runs a test list.  See the documentation of
these functions for details on how to use them.

Note that it is acceptable for a test function to load input files internally.
The benefit to type 2 functions over type 1 functions is in the case when
several functions share an input file.  The test system will only load the
file once if possible, rather than for every test function that uses is.
Consequently, each test function must not depend on modifications made
to the the `(mesh, sbp, eqn,  opts)` objects by other test functions, because
the objects may or may not be re-initialized between calls to the test
functions.  Also, because tags can be used to selectively run tests, it
cannot be guaranteed that a particular test will run before another one.
