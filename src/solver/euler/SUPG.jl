# SUPG implementation

function calcFluxJacobian{Tmsh, Tsol, Tdim}(mesh::AbstractMesh{Tmsh}, 
	                                        sbp::SBPOperator,  
	                                        eqn::EulerData{Tsol, Tdim}, 
	                                        q::AbstractArray{T,1},opts)
  # A1 & A2 represent the flux jacobian in this case. This code is only for 2D
  # This code works at the nodallevel
  # Source: http://www.theoretical-physics.net/dev/fluid-dynamics/euler.html
  
  # Create Flux Jacobian matrices at the nodal level
  A1 = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode)
  A2 = zeros(A1)
  gamma_1 = eqn.gamma_1
  gamma = eqn.gamma
  u = q[2]/q[1] # Get velocity in the x-direction 
  v = q[3]/q[1] # Get velocity in the x-direction 
  R = eqn.R     # Gas constant
  cv = eqn.cv   # Specific heat at constant volume
  intvar = (R/cv)*(q[4]/q[1] - 0.5*(u*u + v*v)) # intermediate variable

  # Populating A1
  A1[1,1] = 0
  A1[1,2] = 1
  A1[1,3] = 0
  A1[1,4] = 0
  A1[2,1] = -u*u + 0.5*R*(u*u + v*v)/cv 
  A1[2,2] = 2*u - R*u/cv
  A1[2,3] = -R*v/cv
  A1[2,4] = R/cv
  A1[3,1] = -uv
  A1[3,2] = v
  A1[3,3] = u
  A1[3,4] = 0
  A1[4,1] = -q[2]*q[4]/(q[1]*q[1]) - u*intvar + u*R*(u*u + v*v)/cv
  A1[4,2] = q[4]/q[1] + intvar - R*u*u/cv
  A1[4,3] = -R*u*v/cv
  A1[4,4] = u + R*u/cv

  # Populating A2
  A2[1,1] = 0
  A2[1,2] = 0
  A2[1,3] = 1
  A2[1,4] = 0
  A2[2,1] = -v*u 
  A2[2,2] = v
  A2[2,3] = u
  A2[2,4] = 0
  A2[3,1] = -v*v + 0.5*R*(u*u + v*v)/cv
  A2[3,2] = -R*u/cv
  A2[3,3] = 2*v - R*v/cv
  A2[3,4] = R/cv
  A2[4,1] = -q[3]*q[4]/(q[1]*q[1]) - v*intvar + v*(R/cv)*0.5*(u.u + v.v)
  A2[4,2] = -R*v*u/cv
  A2[4,3] = q[4]/q[1] + intvar - R*v*v/cv
  A2[4,4] = v + R*v/cv
  
  return nothing
end # end calcFluxJacobian


# Stabilization term
function calcStabilizationTerm(q::AbstractArray{T,1}, h,
                               eqn::EulerData{Tsol, Tdim})
  # Nodal level function. Calculates the stabilization term tau


end
