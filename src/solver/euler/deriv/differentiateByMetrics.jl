# This script holds differentiation by the mesh metrics. Write functions in here.
function getdFdm{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, 
                 sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim}, opts)

  dFluxdm = zeros(Tsol, mesh.numDofPerNode, Tdim, Tdim, sbp.numnodes, mesh.numEl)
  nrm = zeros(Tmsh, Tdim)
  for i = 1:mesh.numEl
  	for j = 1:mesh.numNodesPerElement
  	  q_vals = sview(eqn.q, :, j, i)
      aux_vars = sview(eqn.aux_vars, :, j, i)
      for k=1:Tdim  # loop over dimensions
        for p=1:Tdim
          nrm[p] = mesh.dxidx[k, p, j, i]
        end
        dfdm = sview(dFluxdm, :, :, k, j, i)
        calcdFluxdm(eqn.params, q_vals, aux_vars, nrm, dfdm)
      end
  	end
  end


  return dFluxdm
end

function calcdFluxdm{Tmsh, Tsol, Tres}(params::ParamType{2, :conservative}, 
                      q::AbstractArray{Tsol,1}, 
                      aux_vars::AbstractArray{Tres, 1}, 
                      dir::AbstractArray{Tmsh},  dF::AbstractArray{Tsol,2})


  # Function for computing the derivative of flux w.r.tthe mesh metrics

  press = calcPressure(params, q) # Calculate pressure

  # dF = zeros(Tsol, mesh.dim+2, mesh.dim)

  fac = 1/q[1]
  u = q[2]*fac
  v = q[3]*fac

  dF[1,1] = q[1]*u
  dF[2,1] = q[2]*u + press
  dF[3,1] = q[3]*u
  dF[4,1] = (q[4] + press)*u

  dF[1,2] = q[1]*v
  dF[2,2] = q[2]*v
  dF[3,2] = q[3]*v + press
  dF[4,2] = (q[4] + press)*v

  return nothing
end

function calcdFluxdm{Tmsh, Tsol, Tres}(params::ParamType{3, :conservative}, 
                      q::AbstractArray{Tsol,1}, 
                      aux_vars::AbstractArray{Tres, 1}, 
                      dir::AbstractArray{Tmsh},  dF::AbstractArray{Tsol,2})

  # Function for computing the derivative of flux w.r.tthe mesh metrics

  press = calcPressure(params, q) # Calculate pressure

# dF = zeros(Tsol, mesh.dim+2, mesh.dim)

  fac = 1/q[1]
  u = q[2]*fac
  v = q[3]*fac
  w = q[4]*fac

  dF[1,1] = q[1]*u
  dF[2,1] = q[2]*u + press
  dF[3,1] = q[3]*u
  dF[4,1] = q[4]*u
  dF[5,1] = (q[5] + press)*u

  dF[1,2] = q[1]*v
  dF[2,2] = q[2]*v
  dF[3,2] = q[3]*v + press
  dF[4,2] = q[4]*v
  dF[5,2] = (q[4] + press)*v

  dF[1,3] = q[1]*w
  dF[2,3] = q[2]*w
  dF[3,3] = q[3]*w
  dF[4,3] = q[4]*w + press
  dF[5,3] = (q[5] + press)*w

  return nothing
end
