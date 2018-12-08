# nan_tools.jl
# Helpful tools for debugging NaN's

#------------------------------------------------------------------------------
function hasnan(arr::AbstractArray{T, N}) where {T,N}

  len = length(find(isnan.(arr)))

  if len > 0
    return true
  else
    return false
  end
end
export hasnan

# Heres an example that was helpful in debugging a serial/parallel issue:
#   if hasnan(eqn.q_vec) println(BSTDOUT, " myrank: ", mesh.myrank, "  NaN in eqn.q_vec") end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Jared's nan tools
"""
  Returns true if an array contains a NaN, false otherwise
"""
function containsNaN(arr::AbstractArray)

  hasnan = false
  for i=1:length(arr)
    hasnan = hasnan || isnan(arr[i])
    if hasnan
      break
    end
  end

  return hasnan
end

"""
  Returns a tuple of indices of the first NaN of an array.  Returns
  a tuple of zeros if no NaNs are present
"""
function firstNaN(arr::AbstractArray{T, N}) where {T, N}

  hasnan = false
  idx = 0
  for i=1:length(arr)
    hasnan = hasnan || isnan(arr[i])
    if hasnan
      idx = i
      break
    end
  end

  if hasnan
    indices = ind2sub(size(arr), idx)
  else
    indices = ntuple(i -> 0, Val{N})
  end

  return indices
end

# Jared's examples
#=
A = rand(2,2)
A[2, 1] = NaN
println("containsNaN(A) = ", containsNaN(A))

B = rand(Complex128, 2,2)
#B[2, 1] = Complex128(NaN, NaN)
println("containsNaN(B) = ", containsNaN(B))

A[2, 1] = 0
println("firstNaN(A) = ", firstNaN(A))
=#
#------------------------------------------------------------------------------


