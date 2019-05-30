This examples solves a NACA 0012 airfoil at Mach number 0.8 and angle of attack
of 1.25 degrees.  The solution contains a strong shock wave on the upper surface
of the airfoil, as well as a weaker shock wave on the lower surface.
The free stream boundary is 40 chord length away from the airfoil.

The discretization is a p=1 SBP Diagonal-E operator with an entropy-stable
formulation and Local-Projection stabilization.  The shock sensor of
 Persons and Peraire is combined with a volume-based entropy-stable dissipation.
A coarse mesh with approximately 3,000 elements is used.  The elements
are curved to better conform to the shape of the boundary.

To run the case:

```
  julia ./run.jl input_vals1.jl
```

This will create a new directory called `work` and solve the case within it.
When you have finished inspecting the results, delete the `work` directory.

The case can also be run in parallel:

```
  mpirun -np 4 julia ./input_vals_parallel.jl
```
