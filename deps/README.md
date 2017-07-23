# Installation

PDESolver depends on several packages, some registered Julia packages, others 
not.

## Regular Packages

Like any Julia package, most of the dependencies that are registered packages
are listed in the `REQUIRE` file along with version numbers known to work.
The Julia package manager is used to resolve any dependency conflicts with 
other installed packages.

## Non-Standard Packages

The dependencies that are not registered packages are installed by manually 
cloning their repositories and building them.  Commit hashes are used to 
identify versions known to work.  The commits are checked out into a branch
called `pdesolver_version`.

If a package is already installed, by default it will left alone (ie. the 
specified commit will not be checked out).  This behavior can be overrided 
by declaring certain environmental variables in the shell before launching 
Julia (see below)


## Manual Installation Process

The packages listed in the `REQUIRE` file can be installed manually, by 
cloning the repository and checkout out a particular commit known to work.
This will force the checking out of he commit and re-building of the package 
even if it is already installed.

## Environmental Variables

`PDESOLVER_INSTALL_DEPS_MANUAL`: install the packages in the `REQUIRE` file
manually, forcefully.
`PDESOLVER_FORCE_DEP_INSTALL_ALL`: forces the checkout and re-installation 
 of the non-standard packages, even if already installed
`PDESOLVER_FORCE_DEP_INSTALL_pkg_name`: force the checkout and re-installation
of the package named `pkg_name`.

For all these environmental variables, the value is not important, only the 
existance of the variable.


## Offline Installation

If the internet is not accessible from the machine where the code is to be 
installed, it is possible to download the packages to a specified 
directory on another machine and then copy them to the machine where they
will be installed.  To do this, set the envinronmental variable 
`PDESOLVER_BUNDLE_DEPS` to tell the build script to perform bundling instead of 
installation, and set the path to the directory where the deps should be stored
in the variable `PDESOLVER_PKGDIR`.

After running the build script, copy the contents of the `PDESOLVER_PKGDIR` to
the package directory on the target machine, set `PDESOLVER_UNBUNDLE_DEPS`, and
run the `PDESolver/build.jl`.

## Logging

The file `deps/install.log` is created and written to with the progress of the
checkout process, including error information if a build fails.


