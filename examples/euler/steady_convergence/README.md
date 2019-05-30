# Stead convergence study example

This directory contains an sample convergence study.  The convergence study
uses the 2D isentropic vortex on meshes that are 5 x 5 elements to 20 x 20
elements, doubling the number of elements in each direction from one grid
to the next (3 meshes).

Although this problem is traditionally run on a quarter annulus domain, in
this case it is run on a square domain (x,y in the range [1, 3]) and
Dirichlet boundary conditions are imposed on all boundaries using a
characteristic approach.

This problem can also be run on the quarter annulus domain using curved
elements and the Euler wall boundary condition for the inner and outer
radii.  Dirichlet boundary conditions are applied to the inlet and outlet,
using a charactereistic based approach.  Because the mesh elements are
defined by quadratic curves, the convergence rate is sub-optimal for
degree 3 or 4 SBP operators.

The spatial discretization is an entropy stable formulation with diagonal-E
SBP operators.  The interface flux is the HLL Riemann solver, and the
Local-Projection Stabilization is used to stabilize the diagonal-E operators.
The nonlinear problem is solved using Newton's method.

# Organization:

  * `meshes`:
      contains the mesh files.  Their naming convention is
      `square_vortex_s_i_0.smb`, where `i` is the index of the mesh (1 through
      3) and the final `0` indicates the mesh file is used by the first MPI
      rank.  In this case the meshes are serial so there is only one
      MPI rank.  In general, the final number in the file name
      (before the extension) goes from `0` to `n-1` for a mesh that is
      partitioned for `n` MPI ranks.  The parallel meshes are named
      `square_vortex_p4_i_.smb`, 
      For the quarter annulus domain, the serial meshes are named
      `vortex_quadratic_s_i_0.smb` in serial and `vortex_quadratic_p4_i_.smb`
      in parallel.
      The parallel meshes are named `square_vortex_p4_i_*.smb`
  * `p1_serial`:
      convergence study files for p=1 SBP Diagonal-E operators on the square
      domain
  * `p1_serial_curve`:
     convergence study files for the quarter annulus domain in serial
  * `p1_parallel`:
      convergence study files for p=1 SBP Diagonal-E operators in
      parallel
  * `p1_parallel_curve`:
     convergence study files for the quarter annulus domain in parallel.
  * `calc_line.jl`:
      calculates the least squares slope of the error vs h on
      a log-log plot
  * `cleanup.sh`:
      delete all output files and directories, restoring directory
      to pristine state, should be run from the either the
      `p1_serial` or `p1_parallel` directory.

The `p1` directory contains 3 files. `input_vals1.jl` contains
a prototype input file.  This file specifies all the problem parameters that
the solver needs. `cptest.sh` creates the directories `m1` through `m3` and
copies the input file into each one.  It then modifies the mesh filename listed
in the input file so the solver will load the  `i`th mesh.
The file `runtests.jl` runs all 3 meshes, computes the error in each case,
and writes it to a file.
After all three runs have finished, running the `calc_line.jl` script located
in the top level directory will compute the convergence rate.

# Running the serial convergence study

To run the p=1 serial convergence study:

```
  cd p1_serial
  ./cptest
  julia ./runtest.jl
  julia ../calc_line.jl
```

These convergence tests take only a few seconds to run on my machine, however
compiling the code the first time takes several minutes

If you want to delete the results and run the study again, run the `cleanup.sh`
script from the `p1` directory.

To run the convergence study with p=2 operators, change the `order` option in
`p1_serial/input_vals1.jl` and run the above commands again.

The same steps can be followed for the `p1_serial_curve` directory`.

# Running the parallel convergence study

Similar to the serial case:

```
  cd p1_serial
  ./cptest
  mpirun -np 4 julia ./runtest.jl
  julia ../calc_line.jl
```

Because these meshes are so small, there is not much benefit to solving
in parallel, however this case does demonstrate the use of an iterative
linear solver in Newton's method.

Notice that `p1_parallel/input_vals1.jl` uses the parallel meshes
(`square_vortex_p4_i_.smb`) for the mesh file name.
Because the meshes are partition into 4 parts, the `mpirun` command above must
launch exactly 4 processes.


The same steps can be followed for the `p1_parallel_curve` directory`.
