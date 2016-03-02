# filtering.jl
# this function contains all the functions to do filtering stabilization

@doc """
### EulerEquationMod.applyFilter

  This function multiplies a filter matrix by the variables at every
  node in the mesh.  The filter matrix is stored in eqn.params.filter_mat

  Inputs:
    mesh
    sbp
    eqn
    opts
   
  Keyword Args:
    trans=false  : transpose the filter matrix or not.

  Inputs/Outputs:
    arr: 3D array to multiply the filter matrix by


"""->
function applyFilter{Tmsh, Tsol}(mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP, 
                     eqn::AbstractSolutionData{Tsol}, opts, 
                     arr::Abstract3DArray; trans=false)
# applies filter to array arr
# trans determine whether or not to transpose the filter matrix

  if trans
    filter_mat = eqn.params.filter_mat.'
  else
    filter_mat = eqn.params.filter_mat
  end
  
  # holds the filtered results
  q_filt = zeros(Tsol, mesh.numNodesPerElement, mesh.numDofPerNode)  
  len = mesh.numNodesPerElement*mesh.numDofPerNode
  for i=1:mesh.numEl
    q_vals = view(arr, :, :, i)  # mesh.numDof x sbp.numnodes
    # apply filter matrix to q_vals transposed, so it is applied to
    # all the rho components, then all the x momentum, etc.
    smallmatmatT!(filter_mat, q_vals, q_filt)

    # copy values back into eqn.q, remembering to transpose q_filt
    for j=1:mesh.numNodesPerElement
      for k=1:mesh.numDofPerNode
	arr[k, j, i] = q_filt[j, k]
      end
    end

  end  # end loop over elements

  return nothing
end



@doc """
### EulerEquationMod.calcfilter

  This function calculates the filter used by applyFilter().
  The filter is calculated as V*filter_kernel*inv(V), where V is a matrix
  that converts the interpolating SBP basis to a modal basis and filter_kernal
  is a matrix (typically diagonal) that performs the filtering of the modal
  basis coefficients.

  Inputs:
    sbp: SBP operator used to compute the solutions
    filter_name: and ASCIIString that matches the name of a filter function.
                 This string is used to retrieve the function from a dictionary
                 of all supported filter function.
    opts:  options dictonary

  Outputs:
    F_ret: a numNodesPerElement x numNodesPerElement filter operator

"""->
function calcFilter(sbp::AbstractSBP, filter_name::ASCIIString, opts)
# calc the filter specified by filter_name


  filter_func = filter_dict[filter_name]
  filt_mat = filter_func(sbp, opts)

  # testing only
#  (m,n) = size(filt_mat)
#  filt_mat = eye(m)


  V = calcModalTransformationOp(sbp)
  println("filt_mat = \n", filt_mat)
  println("V = \n", V)
  println("cond(V) = ", cond(V))
  # calculate V*filter_mat*inv(V)
  # which converts from interpolating to modal, applies filter, then 
  # converts back to interpolating
  F_t = filt_mat.'
  V_t = V.'
  F_ret = (V_t\F_t).'
  F_ret = V*F_ret

  F_ret = V*filt_mat*inv(V)
  # for testing, return identity matrix
#  (m,n) = size(filt_mat)
#  F_ret = eye(m) 


  for i=1:length(F_ret)
    if abs(F_ret[i]) < 1e-15
      F_ret[i] = 0
    end
  end

  println("F_ret = \n", F_ret)
  return F_ret

end


@doc """
### EulerEquationMod.calcModalTransformationOp

  This function calculates a matrix operator V that transforms from a modal
  basis to an interpolating basis for SBP operators.  Because the transformation
  itself is rank deficient, V is augmented with the last n vector of its Q
  matrix (from the full QR decomposition), where numNodesPerElement - rank(V)
  = n

  Inputs:
    sbp:  SBP operator

  Outputs:
    V_full: a numNodesPerElement x numNodesPerElement generalized
            Vandermonde matrix

"""->
function calcModalTransformationOp(sbp::AbstractSBP)

  vtx = [-1. -1; 1 -1; -1 1]  # reference element
  x, y = SummationByParts.SymCubatures.calcnodes(sbp.cub, vtx)
  # loop over ortho polys up to degree d
  d = sbp.degree
  n = convert(Int, (d+1)*(d+2)/2)
  V = zeros(Float64, (sbp.numnodes, n))
  Vx = zeros(V)
  Vy = zeros(V)
  ptr = 0
  for r = 0:d
    for j = 0:r
      i = r-j
      V[:,ptr+1] = SummationByParts.OrthoPoly.proriolpoly(x, y, i, j)
      Vx[:,ptr+1], Vy[:,ptr+1] = SummationByParts.OrthoPoly.diffproriolpoly(x, y, i, j)
      ptr += 1
    end
  end


  Q, R = qr(V, thin=false)
  Qx, Rx = qr(Vx, thin=false)
  Qy, Ry = qr(Vy, thin=false)

  # make the V matrix square
  # for this to work, the QR decomposition *cannot* do column pivoting
  # if it does we will have to do a little more bookkepping
  V_full = zeros(sbp.numnodes, sbp.numnodes)
  V_full[:, 1:n] = V
  V_full[:, (n+1):end] = Q[:, (n+1):end]

  # make Vx, Vy square
  # also make them full rank by replacing linearly dependend columns with vectors from Q

  Vx_full = zeros(sbp.numnodes, sbp.numnodes)
  Vy_full = zeros(sbp.numnodes, sbp.numnodes)

  Vx_full[:, 1:n] = Vx
  Vy_full[:, 1:n] = Vy

  Vx_full[:, (n+1):end] = Qx[:, (n+1):end]
  Vy_full[:, (n+1):end] = Qy[:, (n+1):end]

  # also replace linearly dependent columns here?


  return V_full
end

@doc """
### EulerEquationMod.calcRaisedCosineFilter

  This function constructs a diagonal matrix filt for a 2 dimensional basis
  filter, using the raised cosine filter described in Spectral Methods for
  the Euler Equations: Part I by Hussaini, Kproiva, Salas, Zang

  Inputs
    sbp: SBP operator
    opts: options dictonary

  Outputs:
    filt: an numNodesPerElement x numNodesPerElement diagonal filter matrix

"""->
function calcRaisedCosineFilter(sbp::AbstractSBP, opts)
# calculates the 1D raised cosine filter
# from Spectral Methods for the Euler Equations: Part I - Fourier Methods and
# shock Capturing
# Hussaini, Koproiva, Salas, Zang

  filt = zeros(sbp.numnodes, sbp.numnodes)

  max_mode = getPascalLevel(sbp.numnodes)
  for i=1:sbp.numnodes
    mode_i = getPascalLevel(i)
    theta_i = (mode_i-1)/max_mode
    filt[i, i] = 0.5*(1 + cos(theta_i))
  end

  # get rid of nullspace component
  filt[sbp.numnodes, sbp.numnodes] = 0.0


  # increase the steepness of the filter
  for i=1:sbp.numnodes
    diff =  1 - filt[i, i]
    filt[i, i] -= 2*diff
  end

  for i=1:sbp.numnodes
    if filt[i,i] < 0
      filt[i,i] = 0
    end
  end

  
  return filt
end


@doc """
### EulerEquationMod.getPascalLevel

  This function returns the level of Pascals Triangle a particular node
  lies in.

  Inputs:
    node:  node index

  Outputs:
    level: integer describing level of Pascals triangle.
    
"""->
function getPascalLevel(node::Integer)
# get the current polynomial order of some entry node in 
# Pascals triangle
# this assumes the tree is being traversed in order, from 1 to n

  level = 1
  for i=1:(node+1) # loop over level of triangles
            # looping all the way to i is excessive
    # get the maximum node in current level of triangle
#    println("Level = ", level)
    max_n = div(level*(1 + level), 2)
#    println("max_n = ", max_n)
    # break when we are at the right level
    if  node <= max_n
#      println("breaking")
      break
    end

    level += 1  # increment level
  end


  return level
end


@doc """
### EulerEquationMod.calcLowPassFilter

  This function calculates a low pass filter diagonal filter matrix.

  Inputs:
    sbp: SBP operator
    opts: options dictonary

  Outputs:
    filt: numNodesPerElement x numNodesPerElement diagonal filter matrix

"""->
function calcLowPassFilter(sbp::AbstractSBP, opts)

  filt = zeros(sbp.numnodes, sbp.numnodes)
  for i=1:sbp.numnodes
    filt[i, i] = 1
  end

  filt[end, end] = 0.99

  
  return filt
end


global const filter_dict = Dict{ASCIIString, Function}(
"raisedCosineFilter" => calcRaisedCosineFilter,
"lowPassFilter" => calcLowPassFilter,
)



