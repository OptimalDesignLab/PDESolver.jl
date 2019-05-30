arg_dict = Dict{Any, Any}(
# specify the physics and SBP operator
"physics" => "Euler",  # specify physics to run
"operator_type" => "SBPDiagonalE",  # specify SBP operator
"dimensions" => 2,  # this is a two dimensional mesh

# specify temporal and spatial discretization
"run_type" => 5,  # use steady Newton's method
"jac_type" => 3,  # Petsc linear solver
"jac_method" => 2, # use complex step to compute Jacobian, or explicit
                   # Jacobian computation if possible
"real_time" => false,  # this is an steady problem
"order" => 1,  # p = 1 operator
"use_DG" => true,  # use discontinuous galerkin solver
"volume_integral_type" => 2,  # entropy-stable volume integral
"Volume_flux_name" => "IRFlux",  # flux function for volume integrals
"Flux_name" => "HLLFlux",  # numerical flux function for face integrals,
                             # HLL Riemann solver
# specify the problem itself
"IC_name" => "ICIsentropicVortex",  # initial condtiion
"numBC" => 2,  # number of boundary conditions
"BC1" => [0, 2],  # geometric edges to apply the BC to
"BC1_name" => "isentropicVortexBC",   # name of boundary condition
"BC2" => [1, 3],
"BC2_name" => "noPenetrationESBC",  # entropy stable version of Euler wall BC

# specify mesh
"smb_name" => "../../meshes/vortex_quadratic_p4_1_.smb",
"dmg_name" => ".null",

# misc options
"res_abstol" => 1e-10,
"res_reltol" => 1e-10,
"itermax" => 20,
"use_lps" => true,
"do_postproc" => true,  # calculate error at end of run
"exact_soln_func" => "ICIsentropicVortex",  # use this function for the exact soltuion (to calculate error)
)
