using ODLCommonTools
include("lapack.jl")

"""
  Holds temporary arrays for [`clipJac!`](@ref).  These arrays have element
  type that is real, not complex
"""
struct ClipJacData{T}
  jac_sym::Matrix{T}
  E_tmp::Matrix{T}
  Jac2::Matrix{T}
  lambda::Vector{T}
  work::Matrix{T}
  info::Vector{Int64}
end

"""
  Constructs a `ClipJacData{Float64}`.

  **Inputs**

   * m : the matrix passed to `clipJac!` is `m x m`
"""
function ClipJacData(m::Integer)

  jac_sym = zeros(Float64, m, m)
  E_tmp   = zeros(Float64, m, m)
  Jac2    = zeros(Float64, m, m)
  lambda  = zeros(Float64, m)
  work    = zeros(Float64, m, m)
  info    = zeros(Int, 1)

  return ClipJacData{Float64}(jac_sym, E_tmp, Jac2, lambda, work, info)
end


function clipJac2!(Jac::AbstractMatrix, data::ClipJacData)

  # compute the symmetric part of Jac
  m = size(Jac, 1)  # jac is square
  
  jac_sym = data.jac_sym
  @simd for i=1:m
    @simd for j=i:m  # lower triangle
      t1 = 0.5*(real(Jac[j, i]) + real(Jac[i, j]))
      jac_sym[j, i] = t1
      #jac_sym[i, j] = t1  # Lapack uses lower triangle
    end
  end

  # get eigenvalues
  lambda = data.lambda
  work = data.work
  info = data.info
  dsyev_fast('V', 'L', jac_sym, lambda, work, info)
  @assert info[1] == 0
  # now jac_sym contains the eigenvectors

  # perform clipping
  for i=1:m
    if lambda[i] > 0
      lambda[i] = 0
    end
  end

  # reconstitute Jac = E*diagm(lambda)*D.'
  #   Do tmp = E*diagm(lambda) first (column scaling)

  E_tmp = data.E_tmp
  @simd for i=1:m
    @simd for j=1:m
      E_tmp[j, i] = jac_sym[j, i]*lambda[i]
    end
  end

  Jac2 = data.Jac2
  LinAlg.BLAS.gemm!('N', 'T', 1.0, E_tmp, jac_sym, 0.0, Jac2)

  # subtract off of complex valued array
  @simd for i=1:length(Jac)
    Jac[i] -= Jac2[i]
  end

  return nothing
end

