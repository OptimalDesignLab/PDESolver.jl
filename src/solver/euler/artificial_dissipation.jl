# artificial_dissipation.jl
# this file contains all the functions needede to construct and apply 
# artificial dissipation

@doc """
### EulerEquationMod.applyDissipation

  This function multiplies the dissipation matrix operator for each element
  by the values in a 3D array (typically eqn.q) for each node of the element 
  and stores the results in eqn.res.

  The dissipation matrix operators must be stored in
  eqn.dissipation_mat, a numNodesPerElement x numNodesPerElement x numEl array.

  Inputs:
    mesh
    sbp
    eqn
    opts
    arr: a 3D array to apply the dissipation matrix to

  Aliasing restrictions: no guarantees what happens if arr is eqn.res

"""->
function applyDissipation{Tmsh, Tsol, T}(mesh::AbstractMesh{Tmsh}, 
                          sbp::AbstractSBP, eqn::AbstractSolutionData{Tsol}, 
                          opts, arr::AbstractArray{T, 3})
# applies the artificial dissipation to the array arr
# arr must be mesh.numDofPerNode by sbp.numNodesPerElement by mesh.numEl

  # holds the filtered q variables
  q_filt = zeros(Tsol, mesh.numNodesPerElement, mesh.numDofPerNode)  
  len = mesh.numNodesPerElement*mesh.numDofPerNode
  for i=1:mesh.numEl
    filt_i = sview(eqn.dissipation_mat, :, :, i)
    q_vals = sview(eqn.q, :, :, i)  # mesh.numDofPerNode x sbp.numnodes
    # apply filter matrix to q_vals transposed, so it is applied to
    # all the rho components, then all the x momentum, etc.
    smallmatmatT!(filt_i, q_vals, q_filt)

    # update eqn.res, remembering to transpose q_filt
    for j=1:mesh.numNodesPerElement
      for k=1:mesh.numDofPerNode
	eqn.res[k, j, i] -= q_filt[j, k]
      end
    end

  end  # end loop over elements

  return nothing
end


@doc """
### EulerEquationMod.calcDissipationOperator

  This function calculates the numNodesPerElement x numNodesPerElement 
  dissipation matrix operator for each element and returns them in an array
  numNodesPerElement x numNodesPerElement x numEl array.

  The dissipation matrix operator is calculated as epsilon*filt.'*h_jac*filt,
  where filt is a (usually diagonal) filter matrix, h_jac is a scaling term 
  that incorporates information about the shape of the element.

  Inputs:
    mesh
    sbp
    eqn
    opts
    dissipation_name: an ASCIIString of the function name used to retrieve
                      the function that generates the matrix filt from a
                      dictonary

  Outputs:
    dissipation_mat: a numNodesPerElement x numNodesPerElement x numEl array

    Aliasing restrictions: none

"""->
function calcDissipationOperator{Tmsh, Tsol}(mesh::AbstractMesh{Tmsh}, 
                                 sbp::AbstractSBP, eqn::AbstractEulerData{Tsol},
                                 opts, dissipation_name::ASCIIString)
# calculates and returns the artificial dissipation operator array

  epsilon = eqn.params.dissipation_const  # dissipation constant

  # get the dissipation filter matrix for the reference element
  dissipation_func = dissipation_dict[dissipation_name]
  filt = getDissipationFilterOperator(sbp, dissipation_func)  

  # store a sbp.numnodes square matrix for each element
  dissipation_mat = zeros(Tmsh, sbp.numnodes, sbp.numnodes, mesh.numEl)

    
  # filter matrix for a non-reference element
  filt_i = Array(Tmsh, sbp.numnodes, sbp.numnodes)
  hinv = inv(diagm(sbp.w))
  h = diagm(sbp.w)
  for i=1:mesh.numEl
    # modify filter matrix to be in real (non reference) space
#=    
    for col = 1:sbp.numnodes
      for row = 1:sbp.numnodes
	filt_i[row, col] = filt[row, col]/mesh.jac[row, i]
      end
    end
=#
    # scale the mass (digonal) mass matrix by jacobian determinent
    # then multiply by h = (1/jac)^(1/p)
#    h_jac = sbp.w./(mesh.jac[:, i].^2.0)

    h_jac = sbp.w./(mesh.jac[:, i].^1.5)
    # JC modification 11/4
#    h_jac = sbp.w./(mesh.jac[:, i])

    h_jac_inv = 1./h_jac

    # this is the proper artificial dissipation
    dissipation_mat[:, :, i] = epsilon*filt.'*diagm(h_jac)*filt

    # this is the used for preconditioning the iterative solver
#    dissipation_mat[:, :, i] = epsilon*filt.'*diagm(sbp.w)*filt
  end  # end loop over elements


  return dissipation_mat

end  # end function



@doc """
### EulerEquatoinMod.getDissipatoinFilterOperator

  This function gets the dissipation filter operator used to construction
  the artificial dissipation matrix.

  The filter is calculated as V*filter_kernel*inv(V), where V is a matrix
  that converts the interpolating SBP basis to a modal basis and filter_kernal
  is a matrix (typically diagonal) that performs the filtering of the modal
  basis coefficients.


  Inputs:
    sbp: SBP operator
    filter: function to call to calculate the filter kernal

  Outputs:
    F: a numNodesPerElement x numNodesPerElement filter matrix
"""->
function getDissipationFilterOperator{T}(sbp::TriSBP{T}, filter::Function)
# calculate the filter operator (including the conversion to and from
# the modal basis) used for artificial dissipation
# the filter function defines the filter kernel

  vtx = [-1. -1.; 1. -1.; -1. 1.]
  x, y = SummationByParts.SymCubatures.calcnodes(sbp.cub, vtx)
  # loop over ortho polys up to degree d
  d = sbp.degree
  eta_c = (d-1)/(d+1)
  s = d
  alpha = 36.0

  V = zeros(T, (sbp.numnodes, div((d+1)*(d+2), 2)) )
#  V = zeros(T, (sbp.numnodes, convert(Int, (d+1)*(d+2)/2)) )
  lambda = zeros(T, (sbp.numnodes) )
  ptr = 0
  for r = 0:d
    for j = 0:r
      i = r-j
      V[:,ptr+1] = SummationByParts.OrthoPoly.proriolpoly(x, y, i, j)
      #TODO: generalize the filter interface
      lambda[ptr+1] = 1.0 - filter(r/(d+1), eta_c, alpha, s) 
      ptr += 1
    end
  end
  lambda[ptr+1:end] = 1.0
  #println("lambda = ",lambda)

  Z = nullspace(V')
  Vt = [V Z]
  F = Vt*diagm(lambda)*inv(Vt)
  return F
end


function damp1(eta, eta_c, alpha, s)
  if (eta <= eta_c)
    return 1.0
  else
    return exp(-alpha*((eta-eta_c)/(1-eta_c))^(2*s))
  end
end



global const dissipation_dict = Dict{ASCIIString, Function}(
"damp1" => damp1
)

