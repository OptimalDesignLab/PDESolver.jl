# wrappers for important Lapack functions

import Base.LinAlg.LAPACK.liblapack, Base.LinAlg.BlasInt, Base.LinAlg.BLAS.@blasfunc

@eval begin
function dsyev_fast(jobz::Char, uplo::Char, A::AbstractMatrix{Float64}, W::AbstractVector{Float64}, work::AbstractMatrix, info::Array{Int64})
  
  N = BlasInt(size(A, 1))
  LDA = BlasInt(size(A, 1))  # no strided matrices
  LWORK = BlasInt(length(work))

#  info = Ref{BlasInt}()  #TODO: this are still heap allocated? add an array argument instead

  @assert length(W) == N
  @assert length(work) > 3*N - 1
  #TODO: assert stride = 1

  ccall( (@blasfunc($(Symbol("dsyev_"))), liblapack), Void,
         (Ptr{UInt8}, Ptr{UInt8}, Ptr{BlasInt}, Ptr{Float64},
         Ptr{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt}),
         &jobz, &uplo, &N, A, &LDA, W, work, &LWORK, info)

  return info[1]
end


end  # end eval
#=
# test dsyev
n = 5
A = rand(n, n)
A = 0.5*(A + A.')
A_orig = copy(A)
W = zeros(n)
work = zeros(A)

info = dsyev_fast('V', 'L', A, W, work)
@assert info == 0

A2 = A*diagm(W)*A.'
println("A_orig = \n", A_orig)
println("A2 = \n", A2)
println("diff = ", A_orig - A2)
@assert maximum(abs.(A_orig - A2)) < 1e-13
=#


