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


function clipJacFast!(Jac::AbstractMatrix, data::ClipJacData, eigs_to_remove::String)
  
  numEigChgs = 0

  println(BSTDOUT, " in clipJacFast. eigs_to_remove: ", eigs_to_remove)

  if eigs_to_remove == "none"
    return numEigChgs
  end

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
  # dsyev_fast args:
  #   jobz, uplo, N, A, LDA, W, work, LWORK, info
  dsyev_fast('V', 'L', jac_sym, lambda, work, info)
  @assert info[1] == 0
  # now jac_sym contains the eigenvectors

  # perform clipping
    # if lambda[i] > 0      # default, origin was JEH's code - must be flipped for Ticon!
    # logic: We are subtracting this Jac from the original. We want to remove positive eigenvalues.
    # So, setting this matrix to only contain zero and positive eigenvalue entries, then subtracting it
    # from the original will leave a Jac with only negative and zero eigenvalue entries.
    # Ex:
    #   [4 -2 -1]   [4  0  0]   [0 -2 -1]
    #   [3 -1  1] - [3  0  1] = [0 -1  0]
    #   [0 -3  3]   [0  0  3]   [0 -3  0]
  if eigs_to_remove == "neg"              # putting eigs_to_remove outside of the for loop to maybe save comp time

    for i=1:m
      if lambda[i] < 0      
        lambda[i] = 0
        numEigChgs += 1
      end
    end

  elseif eigs_to_remove == "pos"

    for i=1:m
      if lambda[i] > 0      
        lambda[i] = 0
        numEigChgs += 1
      end
    end

  else
    error("eigs_to_remove specified improperly. Must be string(pos), string(neg), or string(none).")
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
    # Jac[i] -= Jac2[i]
    Jac[i] = -1.0*Jac2[i]
  end

  # return nothing
  return numEigChgs
end

function clipJac!(Jac::AbstractMatrix, eigs_to_remove::String)
                  # u::AbstractVector,
                  # A::AbstractVector{T}) where T

  numEigChgs = 0

  if eigs_to_remove == "none"
    return numEigChgs
  end

  # scale_u = 1e100   # NOTE: We do _not_ need to scale anything here.

  λ, E = eig(0.5*(Jac+Jac.'))
  n = length(λ)

  # println(BSTDOUT, "\nλ pre clip: ", λ)

  #------------------------------------------------------------------------------
  # Original clipping process
  if eigs_to_remove == "neg"
    for i = 1:n
      if λ[i] < 0.0
        λ[i] = 0.0
        numEigChgs += 1
      end
    end
  elseif eigs_to_remove == "pos"
    for i = 1:n
      if λ[i] > 0.0
        λ[i] = 0.0
        numEigChgs += 1
      end
    end
  else
    error("eigs_to_remove specified incorrectly.")
  end

  # println(BSTDOUT, "λ post clip: ", λ)

  D = diagm(λ)
  Jac[:,:] -= E*D*E.'

  return numEigChgs

end

