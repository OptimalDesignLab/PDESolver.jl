# conversion.jl
# this file contains all functions related to converting from one set of
# variables to another

# organization:
#  The methods that blindly convert a vector of values from one set of 
#  variables to another have names ending with an underscore (convertToEntropy_)
#  The functions that look at the static parameters of the params object 
#  to determine what variables the input vector is in and then does the proper
#  conversion (or no-op, if that is the case) do not have an underscore (ie. 
#  (convertToConservative, convertToEntropy).  Note that it is safe to 
#  convert conservative variables to conservative variables, the proper
#  proper conversion is a no-op. Anyone writing physics code should call
#  the non-underscore methods, never the underscore methods.
#
#  There is also the concepts of the natural variables and the working variables
#  The natural variables are the ones the governing equation is usually 
#  written in (conservative for the Euler equations).  The working variables
#  are the ones the equation is being solved in, as determined by the 
#  static parameters of the params object.  The function 
#  convertFromNaturalToWorking vars does the transformation as indicated
#  It is ocasionally useful for writing generic code because it enables a 
#  function that doesn't know what type of variables it is using to convert 
#  to conservative and and then convert back (converting back is the problem
#  this function solves).

export convertFromNaturalToWorkingVars

#------------------------------------------------------------------------------
# convert variables without checking what type they currently are (unsafe)
#------------------------------------------------------------------------------
@doc """
### EulerEquationMod.convertToEntropy_

  Converts the conservative variables at a node to the entropy variables
  regardless of the static parameter var_type.

  Inputs:
  params  : ParamType{Tdim} used to dispatch to the proper method
  qe  : vector (of length 4 or 5) of conservative variables
  
  Inputs/outputs
  qc : vector (of length or 5) of entropy variables.  Contents of vector are
       overwritten

  Aliasing: none (qc and qe *can* be the same vector)

  Ths is a low level function.
"""->

function convertToEntropy_{Tsol}(params::ParamType{2}, 
                           qc::AbstractArray{Tsol,1}, 
                           qe::AbstractArray{Tsol, 1})

  gamma = params.gamma
  gamma_1 = params.gamma_1

  k1 = 0.5*(qc[2]^2 + qc[3]^2)/qc[1]
  rho_int = qc[4] -k1
  s = log(gamma_1*rho_int/(qc[1]^gamma))
  fac = 1.0/rho_int
  tmp1 = -qc[1]*fac
  qe[1] = (rho_int*(gamma + 1 -s) - qc[4])*fac
  qe[2] = qc[2]*fac
  qe[3] = qc[3]*fac
  qe[4] = tmp1

  return nothing
end

function convertToEntropy_{Tsol}(params::ParamType{3}, 
                           qc::AbstractArray{Tsol,1}, 
                           qe::AbstractArray{Tsol, 1})

  gamma = params.gamma
  gamma_1 = params.gamma_1

  k1 = 0.5*(qc[2]^2 + qc[3]^2 + qc[4]^2)/qc[1]
  rho_int = qc[5] - k1
  s = log(gamma_1*rho_int/(qc[1]^gamma))
  fac = 1.0/rho_int
  tmp1 = -qc[1]*fac
  qe[1] = (rho_int*(gamma + 1 -s) - qc[5])*fac
  qe[2] = qc[2]*fac
  qe[3] = qc[3]*fac
  qe[4] = qc[4]*fac
  qe[5] = tmp1

  return nothing
end


@doc """
### EulerEquationMod.convertToConservative_

  Converts  the entropy variables at a node to the conservative variables
  regardless of the static parameter var_type.

  Inputs:
  params  : ParamType{Tdim} used to dispatch to the proper method
  qe  : vector (of length 4 or 5) of entropy variables
  
  Inputs/outputs
  qc : vector (of length 4 or 5) of conservative variables.  Contents of vector are
       overwritten

  Aliasing: none (qc and qe *can* be the same vector)

  Ths is a low level function.
"""->
function convertToConservative_{Tsol}(params::ParamType{2}, 
                                qe::AbstractArray{Tsol,1}, 
                                qc::AbstractArray{Tsol, 1})

  gamma = params.gamma
  gamma_1 = params.gamma_1
  k1 = 0.5*(qe[2]^2 + qe[3]^2)/qe[4]
  s = gamma - qe[1] + k1
  rho_int = exp(-s/gamma_1)*(gamma_1/(-qe[4])^gamma)^(1/gamma_1)
  qc[1] = -qe[4]*rho_int
  qc[2] = qe[2]*rho_int
  qc[3] = qe[3]*rho_int
  qc[4] = (1.0 - k1)*rho_int
end

function convertToConservative_{Tsol}(params::ParamType{3}, 
                                qe::AbstractArray{Tsol,1}, 
                                qc::AbstractArray{Tsol, 1})

  gamma = params.gamma
  gamma_1 = params.gamma_1
  k1 = 0.5*(qe[2]^2 + qe[3]^2 + qe[4]^2)/qe[5]
  s = gamma - qe[1] + k1
  rho_int = exp(-s/gamma_1)*(gamma_1/(-qe[5])^gamma)^(1/gamma_1)
  qc[1] = -qe[5]*rho_int
  qc[2] = qe[2]*rho_int
  qc[3] = qe[3]*rho_int
  qc[4] = qe[4]*rho_int
  qc[5] = (1.0 - k1)*rho_int
end




#------------------------------------------------------------------------------
# converting to conservative variables safely
#------------------------------------------------------------------------------
@doc """
### EulerEquationMod.convertToConservative

  Converts  the entropy variables at a node to the conservative variables.

  Inputs:
  params  : ParamType{2, :entropy} used to dispatch to the proper method
  qe  : vector (of length 4) of entropy variables
  
  Inputs/outputs
  qc : vector (of length 4) of conservative variables.  Contents of vector are
       overwritten

  Aliasing: none (qc and qe *can* be the same vector)

  Ths is a low level function.
"""->
function convertToConservative{Tdim, Tsol}(params::ParamType{Tdim, :entropy}, 
                               qe::AbstractArray{Tsol,1}, 
                               qc::AbstractArray{Tsol, 1})

  convertToConservative_(params, qe, qc)

end

@doc """
# low level function

  Converts conservative variables to conservative variables (ie. it
  copies the input to the output).  This method exists to values can be 
  converted without knowing whether they are conservative or entropy.

  Aliasing restrictions: none (qc and qe *can* be the same vector)

"""->
function convertToConservative{Tdim, Tsol}(params::ParamType{Tdim, 
                              :conservative}, qe::AbstractArray{Tsol,1}, 
                              qc::AbstractArray{Tsol, 1})
# maybe the compiler will be smart enough to make this a no-op?
  for i=1:length(qe)
    qc[i] = qe[i]
  end

  return nothing
end


@doc """
### EulerEquationMod.convertFromNaturlaToWorkingVars

  This function converts a vector qc from the "natural" set of variables to 
  write the equation to the set of variables the equation is being solved in, 
  defined by the static parameter var_type of the params.  For the Euler 
  equations, the "natural" variables are the conservative variables.
  Every new set of variables must extend this function with a new method.

  Inputs:
    params:  ParamType{Tdim, var_type}
    qc: vector of "natural" variables

  Inputs/Outputs:
    qe: vector to be populated with the new variables

  Low level function.

  Aliasing restrictions: none (this requires that every method of this function
                         support in-place conversion).

"""->
function convertFromNaturalToWorkingVars{Tdim, Tsol}(
                                         params::ParamType{Tdim, :conservative},
                                         qc::AbstractArray{Tsol,1}, 
                                         qe::AbstractArray{Tsol,1})
  for i=1:length(qc)
    qe[i] = qc[i]
  end

end


@doc """
# mid level function

  Converts the array (3D form) of entropy variables to conservative variables 
  in place.  If the array is already in conservative variables this is a no-op.
  Other methods exist if q_arr is a vector.

  Inputs:
    mesh
    sbp
    eqn
    opts

  Inputs/Outputs:
    q_arr: array of values (3D) to be converted

  Aliasing: no restrictions
"""->
#3D array entropy -> conservative
function convertToConservative{Tmsh, Tsol, Tdim, Tres}(mesh::AbstractMesh{Tmsh},
                               sbp::AbstractSBP, 
                               eqn::EulerData{Tsol, Tres, Tdim, :entropy}, 
                               opts, q_arr::AbstractArray{Tsol, 3})

  for i=1:mesh.numEl  # loop over elements
    for j=1:mesh.numNodesPerElement
	q_view = sview(q_arr, :, j, i)  # reuse memory
	convertToConservative(eqn.params, q_view, q_view)
    end
  end

  return nothing
end

# 3D array conservative -> conservative
function convertToConservative{Tmsh, Tsol, Tdim, Tres}(mesh::AbstractMesh{Tmsh},
                               sbp::AbstractSBP, 
                               eqn::EulerData{Tsol, Tres, Tdim, :conservative},
                               opts, q_arr::AbstractArray{Tsol, 3})

  return nothing
end


# q_vec conversion entropy -> conservative
function convertToConservative{Tmsh, Tsol, Tdim, Tres}(mesh::AbstractMesh{Tmsh},
                               sbp::AbstractSBP, 
                               eqn::EulerData{Tsol, Tres, Tdim, :entropy}, 
                               opts, q_vec::AbstractArray{Tsol, 1})

  for i=1:mesh.numDofPerNode:mesh.numDof
    q_view = sview(q_vec, i:(i+mesh.numDofPerNode-1))
    convertToConservative(eqn.params, q_view, q_view)
  end

  return nothing
end

# q_vec conversion conservative -> conservative
function convertToConservative{Tmsh, Tsol, Tdim, Tres}(mesh::AbstractMesh{Tmsh},
                               sbp::AbstractSBP, 
                               eqn::EulerData{Tsol, Tres, Tdim, :conservative},
                               opts, q_arr::AbstractArray{Tsol, 1})

  return nothing
end



#------------------------------------------------------------------------------
# converting to entropy variables
#------------------------------------------------------------------------------

@doc """
### EulerEquationMod.convertToEntropy

  Converts the conservative variables at a node to entropy variables

  Input:
  params : ParamType{s, :conservative}
  qc  : vector (of length 4 or 5) of conservative variables

  Outputs:
  qe : vector (of length 4 or 5) of conservative variables.  Contents of 
       vector areoverwritten

  # low level function

  Aliasing restrictions: none (qc and qe *can* be the same vector)

"""->
function convertToEntropy{Tdim, Tsol}(params::ParamType{Tdim, :conservative}, 
               qc::AbstractArray{Tsol,1}, qe::AbstractArray{Tsol,1})

  convertToEntropy_(params, qc, qe)

end


@doc """
### EulerEquationMod.convertToEntropy

  Converts the entropy variables to entropy variables (ie. it copies the 
  input to the output).  This method exists so variables can be converted 
  to entropy variables without knowing what type they are.

  # low level function

  Aliasing restrictions: none (qc and qe *can* be the same vector)

"""->
function convertToEntropy{Tdim, Tsol}(params::ParamType{Tdim, :entropy}, 
               qc::AbstractArray{Tsol,1}, qe::AbstractArray{Tsol,1})

  for i=1:length(qc)
    qe[i] = qc[i]
  end

  return nothing
end

#doc """


function convertFromNaturalToWorkingVars{Tdim, Tsol}(
               params::ParamType{Tdim, :entropy}, 
               qc::AbstractArray{Tsol,1}, qe::AbstractArray{Tsol,1})

  convertToEntropy_(params, qc, qe)
end

@doc """
# mid level function

  Converts the array (3D form) of conservative variables to entropy variables 
  in place.  If the array is already in entropy variables this is a no-op

  Methods also exist for the 1D form.

  Inputs:
    mesh
    sbp
    eqn
    opts

  Inputs/Outputs:
    q_arr: array of values (3D) to be converted

  Aliasing: no restrictions
"""->
# 3D array conservative -> entropy
function convertToEntropy{Tmsh, Tsol, Tdim, Tres}(mesh::AbstractMesh{Tmsh}, 
                          sbp::AbstractSBP, 
                          eqn::EulerData{Tsol, Tres, Tdim, :conservative}, 
                          opts, q_arr::AbstractArray{Tsol, 3})



  for i=1:mesh.numEl  # loop over elements
    for j=1:mesh.numNodesPerElement
	q_view = sview(q_arr, :, j, i)  # reuse memory
	convertToEntropy(eqn.params, q_view, q_view)
    end
  end

  return nothing
end


# 3D array entropy -> entropy
function convertToEntropy{Tmsh, Tsol, Tdim, Tres}(mesh::AbstractMesh{Tmsh}, 
                          sbp::AbstractSBP, 
                          eqn::EulerData{Tsol, Tres, Tdim, :entropy}, opts, 
                          q_arr::AbstractArray{Tsol, 3})

  return nothing
end



# q_vec conversion conservative -> entropy
function convertToEntropy{Tmsh, Tsol, Tdim, Tres}(mesh::AbstractMesh{Tmsh}, 
                          sbp::AbstractSBP, 
                          eqn::EulerData{Tsol, Tres, Tdim, :conservative}, opts,
                          q_vec::AbstractArray{Tsol, 1})

  for i=1:mesh.numDofPerNode:mesh.numDof
    q_view = sview(q_vec, i:(i+mesh.numDofPerNode-1))
    convertToEntropy_(eqn.params, q_view, q_view)
  end

  return nothing
end

# q_vec conversion entropy -> entropy
function convertToEntropy{Tmsh, Tsol, Tdim, Tres}(mesh::AbstractMesh{Tmsh}, 
                          sbp::AbstractSBP, 
                          eqn::EulerData{Tsol, Tres, Tdim, :entropy}, opts, 
                          q_arr::AbstractArray{Tsol, 1})

  return nothing
end





