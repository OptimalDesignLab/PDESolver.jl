# functions for projecting back and forth between x-y and n-t 
# (normal-tangential) coordinate systems
# in 3D, x-y-z to n-t-b (normal-tangential-binormal)
"""
  This function populates the matrix P that project from x-y-z to n-t-b.
  Methods are available for 2D and 3D, determined from the AbstractParamType
  object.  This is a somewhat specialized routine in that it only works for
  vectors like the state vector of the Euler or Navier-Stokes equations,
  where the first and last equations are coordinate system invarient and
  only the middle 2 (in 2D) or 3 (in 3D) equations need to be rotated.
  Specialized multiplication routines are provied as projectToXY and
  projectToNT.

  Inputs:

    params:  An AbstractParamType{Tdim}
    nrm: the normal direction in x-y coordinate system.  Must be a unit vector

  Inputs/Outputs

    P:  matrix to be populated with projection matrix

  Aliasing restrictions: none
"""
function getProjectionMatrix(params::AbstractParamType{2}, nrm::AbstractVector,
                             P::AbstractMatrix)
# test: P.' = inv(P)

  t1, t2 = getOrthogonalVector(params, nrm)

  # zero out P
  for i=1:length(P)
    P[i] = 0
  end
  P[1,1] = 1
  P[4,4] = 1
  P[2,2] = nrm[1]
  P[2,3] = nrm[2]
  P[3,2] = t1
  P[3,3] = t2

  return nothing
end

function getProjectionMatrix(params::AbstractParamType{3}, nrm::AbstractVector,
                             P::AbstractMatrix)
# test: P.' = inv(P)
  n1 = nrm[1]; n2 = nrm[2]; n3 = nrm[3]
  t1, t2, t3 = getOrthogonalVector(params, nrm)
  b1, b2, b3 = getBinormalVector(n1, n2, n3, t1, t2, t3)

  for i=1:length(P)
    P[i] = 0
  end

  P[1,1] = 1
  P[5,5] = 1
  P[2,2] = n1
  P[2,3] = n2
  P[2,4] = n3
  P[3,2] = t1
  P[3,3] = t2
  P[3,4] = t3
  P[4,2] = b1
  P[4,3] = b2
  P[4,4] = b3

  return nothing
end


"""
  This function generates a unit vector orthogonal to the input vector nrm and 
  returns its components.  Methods are availble for 2D and 3D

  Inputs:

    params: an AbstractParamType, used to dispatch to the 2D or 3D method
    nrm: the input vector

  Outputs:

    t1, t2, (and t3 in 3D): the components of the unit vector orthogonal to
                            nrm

  Aliasing restrictions: none
"""
function getOrthogonalVector(params::AbstractParamType{2}, nrm::AbstractVector)
# this code uses a clever trick to deal with the cases where one of the
# components of nrm is zero in a branch free manner

  # add this factor to the denominators to avoid NaNs
  # the round(nrm[2]) below ensure only the result with
  # denominator > 0.5 gets chosen, so adding 1e-50 doesn't
  # affect the results
  add_fac = 1e-50

  # use the fact that dot(nrm, tangent) = 0 to produce a single equation
  # in 2 unknowns
  # pick 2 solutions, one with the first component equal to 2, one with the
  # second component equal to 1.  Below the code figures out which one to use
  v1 = 1.0
  v2 = -nrm[1]/(nrm[2] + add_fac)
  w1 = -nrm[2]/(nrm[1] + add_fac)
  w2 = 1.0

  # if nrm[2] > 0.5, use v, otherwise use w
  fac = round(absvalue(nrm[2]))
  t1 = fac*v1 + (1-fac)*w1
  t2 = fac*v2 + (1-fac)*w2

  len = sqrt(t1*t1 + t2*t2)

  t1 = t1/len
  t2 = t2/len

  return t1, t2
end

function getOrthogonalVector(params::AbstractParamType{3}, nrm::AbstractVector)
# this code uses a clever trick to deal with the cases where one of the
# components of nrm is zero in a branch free manner

  # unpack normal vector
  n1 = nrm[1]
  n2 = nrm[2]
  n3 = nrm[3]


  # add this factor to the denominators to avoid NaNs
  # the round(nrm[3]) below ensure only the result with
  # denominator > 0.5 gets chosen, so adding 1e-50 doesn't
  # affect the results
  add_fac = 1e-50

  # compute the three possible vectors
  v1 = 1.0
  v2 = 1.0
  v3 = -(n1 + n2)/(n3 + add_fac)

  w1 = 1.0
  w2 = -(n1 + n3)/(n2 + add_fac)
  w3 = 1.0

  x1 = -(n2 + n3)/(n1 + add_fac)
  x2 = 1.0
  x3 = 1.0

  # pick either w or x
  fac = round(absvalue(nrm[3]))

  z1 = fac*x1 + (1-fac)*w1
  z2 = fac*x2 + (1-fac)*w2
  z3 = fac*x3 + (1-fac)*w3
  
  # pick either v or z
  fac = round(absvalue(nrm[1]))
  t1 = fac*v1 + (1-fac)*z1
  t2 = fac*v2 + (1-fac)*z2
  t3 = fac*v3 + (1-fac)*z3

  len = sqrt(t1*t1 + t2*t2 + t3*t3)

  t1 = t1/len
  t2 = t2/len
  t3 = t3/len

  return t1, t2, t3
end

"""
  This function computes a vector normal to the 2 supplied vectors and
  returns the components.

  Inputs:

    n1, n2, n3: the components of the first vector
    t1, t2, t3: the components of the second vector

  Outputs:

    b1, b2, b3: the components of the binormal vector

  Aliasing restrictions: none
"""
@inline function getBinormalVector(n1::Number, n2::Number, n3::Number,
                           t1::Number, t2::Number, t3::Number)

  b1 = n2*t3 - n3*t2
  b2 = -(n1*t3 - n3*t1)
  b3 = n1*t2 - n2*t1

  return b1, b2, b3
end

"""
  Calculates the length of a vector, using a manually unrolled loop.
  Methods are available for 2D and 3D.

  Inputs:

    params: an AbstractParamType{Tdim}, used to dispatch to the right method
    nrm: the vector to calculate the length of.  Must have length 2 in 2D
         and 3 in 3D

  Outputs:

    length of nrm
"""
@inline function calcLength(params::AbstractParamType{2}, nrm::AbstractVector)
  return sqrt(nrm[1]*nrm[1] + nrm[2]*nrm[2])
end

@inline function calcLength(params::AbstractParamType{3}, nrm::AbstractVector)
  return sqrt(nrm[1]*nrm[1] + nrm[2]*nrm[2] + nrm[3]*nrm[3]) 
end

function normalize_vec(params::AbstractParamType{T}, nrm::AbstractVector) where {T}
  len = calcLength(params, nrm)
  for d=1:T
    nrm[i] /= len
  end
end


"""
  This function projects a vector x from x-y coordinates to normal-tangential
  (n-t) coordinates, where A is a projection matrix from x-y to n-t, obtained
  from getProjectionMatrix.  This function is a specialized matrix-vector
  multiplication routine and only works for matrices with the particular
  sparsity pattern created by getProjectionMatrix.
  Methods are available for 2D and 3D

  Inputs:

    params: an AbstractParamType{{Tdim} used to dispatch to the 2D or 3D
            method
    P: the projection matrix
    x: vector to be projected

  Inputs/Outputs:

    b: the result of the projection

  Aliasing restrictions: x and b cannot alias
"""

function projectToNT(params::AbstractParamType{2}, P::AbstractMatrix, 
                     x::AbstractVector, b::AbstractVector)
# regular multiplication
  m, n = size(P)

  # first loop: overwrite b
  b[1] = x[1]*P[1, 1]
  b[2] = 0
  b[3] = 0
  b[4] = 0

  b[2] += P[2,2]*x[2]
  b[3] += P[3,2]*x[2]

  b[2] += P[2,3]*x[3]
  b[3] += P[3,3]*x[3]

  b[4] += P[4,4]*x[4]
  
  return nothing
end

function projectToNT(params::AbstractParamType{3}, P::AbstractMatrix, 
                     x::AbstractVector, b::AbstractVector)
# regular multiplication

  m, n = size(P)

  # first loop: overwrite b
  b[1] = x[1]*P[1, 1]
  b[2] = 0
  b[3] = 0
  b[4] = 0
  b[5] = 0

  b[2] += P[2,2]*x[2]
  b[3] += P[3,2]*x[2]
  b[4] += P[4,2]*x[2]

  b[2] += P[2,3]*x[3]
  b[3] += P[3,3]*x[3]
  b[4] += P[4,3]*x[3]

  b[2] += P[2,4]*x[4]
  b[3] += P[3,4]*x[4]
  b[4] += P[4,4]*x[4]

  b[5] += P[5,5]*x[5]
 
  return nothing
end

"""
  This function is similar to projectToNT, except it project from n-t
  to x-y.  Note that P is still the projection matrix from getProjectionMatrix
  that projects from x-y to n-t.
"""
function projectToXY(params::AbstractParamType{2}, P::AbstractMatrix, 
                     x::AbstractVector, b::AbstractVector)
# transpose multiplication

  m, n = size(P)


  b[1] = P[1,1]*x[1]

  b[2] = 0
  b[2] += P[2,2]*x[2]
  b[2] += P[3,2]*x[3]

  b[3] = 0
  b[3] += P[2,3]*x[2]
  b[3] += P[3,3]*x[3]

  b[4] = 0
  b[4] += P[4,4]*x[4]
  
  return nothing
end

function projectToXY(params::AbstractParamType{3}, P::AbstractMatrix, 
                     x::AbstractVector, b::AbstractVector)
# transposed multiplication

  m, n = size(P)


  b[1] = P[1,1]*x[1]

  b[2] = 0
  b[2] += P[2,2]*x[2]
  b[2] += P[3,2]*x[3]
  b[2] += P[4,2]*x[4]

  b[3] = 0
  b[3] += P[2,3]*x[2]
  b[3] += P[3,3]*x[3]
  b[3] += P[4,3]*x[4]

  b[4] = 0
  b[4] += P[2,4]*x[2]
  b[4] += P[3,4]*x[3]
  b[4] += P[4,4]*x[4]

  b[5] = 0
  b[5] += P[5,5]*x[5]
  
  return nothing
end
