# Adjoint Computation

function calcAdjoint{Tmsh, Tsol}(mesh::AbstractDGMesh{Tmsh},sbp::AbstractSBP,
                                 eqn::AdvectionData{Tsol}, opts, functional_number)

  # Get information corresponding to functional_number
  key = string("geom_edges_functional", functional_number)
  functional_edges = opts[key]

  # Calculate the Jacobian of the residual
  newton_data = NonlinearSolvers.NewtonData(mesh, sbp, eqn, opts)
  res_jac = zeros(Tres, mesh.numDof, mesh.numDof)
  pert = complex(0, opts["epsilon"])
  NonlinearSolvers.calcJacobianComplex(newton_data, mesh, sbp, eqn, opts, evalAdvection, pert, res_jac)


  # calculate the derivative of the function w.r.t q_vec
  func_deriv = zeros(Tsol, mesh.numDof)
  adjoint_vec = zeros(Tsol, mesh.numDof)
  # Loop over functioanl edges to calculate the derivative
  for i = 1:length(functional_edges)

  end  # end for i = 1:length(functional_edges)
  

  return nothing
end