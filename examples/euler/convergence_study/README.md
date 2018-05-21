# Convergence study example

This directory contains an sample convergence study.  The convergence study
uses the 2D isentropic vortex on meshes that are 20 x 20 elements to 60 x 60
elements, by 5s (69meshes).

# Organization:

  * `meshes`:
      contains the mesh files.  Their naming convention is
      `vortex_s_i_0.smb`, where `i` is the index of the mesh (1 through 9) and
      the final `0` indicates the mesh file is used by the first MPI
      rank.  In this case the meshes are serial so there is only one
      MPI rank.  In general, the final number in the file name
      (before the extension) goes from `0` to `n-1` for a mesh that is
      partitioned for `n` MPI ranks.
      The parallel meshes are named `vortex_p_i_*.smb`
      The directory also contains a `notes.txt` file describing the
      meshes.
  * `p1_serial`:
      convergence study files for p=1 SBP Omega operators
  * `p1_parallel`:
      convergence study files for p=1 SBP Omega operators in
      parallel
  * `calc_line.jl`:
      calculates the least squares slope of the error vs h on
      a log-log plot
  * `cleanup.sh`:
      delete all output files and directories, restoring directory
      to pristine state, should be run from the either the
      `p1_serial` or `p1_parallel` directory.

The `p1` directory contains 3 files. `input_vals1.jl` contains
a prototype input file.  This file specifies all the problem parameters that
the solver needs. `cptest.sh` creates the directories `m1` through `m6` and
copies the input file into each one.  It then modifies the mesh filename listed
in the input file so the solver will load the  `i`th mesh.
The file `runtests.jl` runs all six meshes, computes the error in each case,
and writes it to a file.
After all six runs have finished, running the `calc_line.jl` script located
in the top level directory will compute the convergence rate.

# Running the serial convergence study

To run the p=1 serial convergence study:

```
  cd p1_serial
  ./cptest
  julia ./runtest.jl
  julia ../calc_line.jl
```

These convergence tests take about 20 minutes on my machine.

If you want to delete the results and run the study again, run the `cleanup.sh`
script from the `p1` directory.

To run the convergence study with p=2 operators, change the `order` option in
`p1_serial/input_vals1.jl` and run the above commands again.

# Running the parallel convergence study

Similar to the serial case:

```
  cd p1_serial
  ./cptest
  mpirun -np 4 julia ./runtest.jl
  julia ../calc_line.jl
```

These convergence tests take about 6 minutes on my machine

Notice that `p1_parallel/input_vals1.jl` uses the parallel meshes
(`vortex_p_i_.smb`) for the mesh file name.
Because the meshes are partition into 4 parts, the `mpirun` command above must
launch exactly 4 processes.
