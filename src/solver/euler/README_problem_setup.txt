change:

== in startup.jl ==

  mesh file: smb_name
  initial conditions: needs to call a function in ic.jl 
    example: ICRho1E2(mesh, sbp, eqn, SL0)

== in euler.jl ==

  in evalBoundaryIntegrals
    BC function inside of boundaryintegrate! call, 6th argument
    ex: rho1Energy2BC, isentropicVortexBC
