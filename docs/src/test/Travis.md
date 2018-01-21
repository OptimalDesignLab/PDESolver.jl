# Travis CI Testing

We currently use the Travis CI service to run automated tests on all
branches.  This page describes the build system

## Build Matrix

A build matrix is used to create one CI build per physics module.
Each build sets an environment variable to indicate to the `test/runtests_travis.sh`
script which physics module to run tests on.
For example, the environment variables

```
TEST_ADVECTION=1
TEST_EULER=1
TEST_SIMPLEODE=1
```

are used to decide whether or not to run tests for the Advection, Euler, and
SimpleODE physics modules. Note that the environment variables set in 
the `env` section of `.travis.yml` must *exactly* match the variables
used in `test/runtests_travis.sh`.

## Caching Dependencies

It is possible to have Travis cache some of PDESolvers dependencies to
save build time.  In particular, Petsc takes some time to compile and does
not change very often, so it would be beneficial to cache it.
One downside of caching is that a manual process is required to update the
cache if we change version of Petsc.


The configuration process is described here.  Most of the work is done in
the `deps/move_petsccache.sh` script.


Each branch has its own cache, and the cache needs to be created during
the first build on a new branch.  Therefore, Travis must be configured to
detect if the cache has been created and build the cache if needed.
One additional wrinkle is that the directories created by the cache
interfere with the PDESolver build system.  The `.travis.yml` configuration
file executes the following procedure:

  1. If `path/to/PETSc2/deps/petsc-3.7.6/arch-linux2-c-debug/lib/libpetsc.so` exists, move the `/lib` directory to `$HOME/lib`, set `PETSC_DIR` and
  `PETSC_ARCH`, and delete `/path/to/PETSc2`.  If the cache did not exist,
   the directory will still exist, but will be empty.  Delete `/path/to/PETSc2` so the PDESolver build system will install the `PETSc2` package.

  2. Install PDESolver and its dependencies

  3. If `PETSC_DIR` and `PETSC_ARCH` are defined (done in step 1 if cache was retrieved), move `$HOME/lib` back to original location, creating any directories in the path as needed.



If `PETSC_DIR` and `PETSC_ARCH` are not set in step 1, step 2 will build
PETSc.  At the end of the Travis run, the `/path/to/PETSc2/deps/petsc-3.7.6/arch-linux2-c-debug/lib` directory will be uploaded to wherever Travis stores cache files.  It will then be available for the next Travis build.


If the cache has already been built for a given branch and a different
version of PETSc is needed, the cache must be deleted so step 2 can build
the new version.
The web interface to PETSc can be used to delete the cache.
