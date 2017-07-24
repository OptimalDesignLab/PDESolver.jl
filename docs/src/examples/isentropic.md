# Example 1: Steady Isentropic Vortex

Save the following as `input_vals.jl`.
```
arg_dict = Dict{Any, Any}(
# specify the physics and SBP operator
"physics" => "Euler",  # specify physics to run
"dimensions" => 2,  # this is a two dimensional mesh

# specify temporal and spatial discretization
"run_type" => 5,  # steady Newton (5 historically was complex-stepped Newton, as opposed to 4 being FD)
"jac_method" => 2,  # complex-step Newton Jacobian calculation
"jac_type" => 1,  # store the Jacobian as a Julia sparse matrix
"t_max" => 10.0,  # make time
"operator_type" => "SBPOmega",  # specify SBP operator
"order" => 1,  # p = 1 operator
"use_DG" => true,  # use discontinuous galerkin solver
"Flux_name" => "RoeFlux",  # numerical flux function used in face integrals
"CFL" => 0.10,  # CFL number
"itermax" => 20,

# specify the problem itself
"IC_name" => "ICIsentropicVortex",  # initial condtiion
"numBC" => 1,  # number of boundary conditions
"BC1" => [0, 1, 2, 3],  # geometric edges to apply the BC to
"BC1_name" => "isentropicVortexBC",   # name of boundary condition

# specify mesh
"smb_name" => "SRCMESHES/squarevortex_small.smb",

# misc options
"write_vis" => true,  # write paraview files
"do_postproc" => true,  # calculate error at end of run
"exact_soln_func" => "ICIsentropicVortex",  # use this function for the exact soltuion (to calculate error)
)
```

Run the case with `julia ~/.julia/v0.4/PDESolver/src/solver/euler/startup.jl input_vals.jl`.
