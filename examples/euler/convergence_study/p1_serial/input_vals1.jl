arg_dict = Dict{Any, Any}(
# specify the physics and SBP operator
"physics" => "Euler",  # specify physics to run
"dimensions" => 2,  # this is a two dimensional mesh

# specify temporal and spatial discretization
"run_type" => 30,  # use LSERK 54
"real_time" => true,  # this is an unsteady problem
"t_max" => 10.0,  # make time
"operator_type" => "SBPOmega",  # specify SBP operator
"order" => 1,  # p = 1 operator
"use_DG" => true,  # use discontinuous galerkin solver
"Flux_name" => "RoeFlux",  # numerical flux function used in face integrals
"CFL" => 0.10,  # CFL number

# specify the problem itself
"IC_name" => "ICUnsteadyVortex",  # initial condtiion
"vortex_x0" => 5.0,  # initial location of vortex
"numBC" => 1,  # number of boundary conditions
"BC1" => [0, 1, 2, 3],  # geometric edges to apply the BC to
"BC1_name" => "unsteadyVortexBC",   # name of boundary condition

# specify mesh
"smb_name" => "../../meshes/vortex_s_1_.smb",
"dmg_name" => ".null",

# misc options
"write_vis" => true,  # write paraview files
"output_freq" => 100,  # do output every this many time steps
"do_postproc" => true,  # calculate error at end of run
"exact_soln_func" => "ICUnsteadyVortex",  # use this function for the exact soltuion (to calculate error)
)
