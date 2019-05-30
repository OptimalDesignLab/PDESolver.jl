This examples solves a NACA 0012 airfoil at Mach number 0.5 and angle of attack
of 2 degrees.  It can be solved with either the Roe scheme or the entropy
stable scheme. To run the example with the Roe scheme:

```
  julia ./run.jl input_vals_roe.jl
```

This will create a new directory called `work` and solve the case within it.
When you have finished inspecting the results, delete the `work` directory.

To solve with the the entropy stable scheme

```
  julia ./run.jl input_vals_es.jl
```

