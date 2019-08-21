# implementation of AbstractBinning for uniformly spaced bins

"""
  Implement uniformly spaced bins

  **Inputs**

  * vlow: the lowest value (the left endpoint of the first bin)
  * vhigh: the highest value (the right endpoint of the last bin)
  * N: the number of bins
  * comm: MPI communicator, defaults to MPI.COMM_WORLD
"""
mutable struct UniformBinning{K, V} <: AbstractBinning{K, V}
  vlow::V
  vhigh::V
  delta::V  # bin size
  N::Int    # number of bins
  bin_keys::Vector{Vector{K}}
  bin_vals::Vector{Vector{V}}

  # used for summation
  bin_keycount::Vector{Int}
  val_sum::Vector{V}

  comm::MPI.Comm

  function UniformBinning{K, V}(vlow, vhigh, N::Integer, comm::MPI.Comm=MPI.COMM_WORLD) where {K, V}

    @assert vhigh > vlow

    delta = (vhigh - vlow)/N
    bin_keys = Vector{Vector{K}}(N)
    bin_vals = Vector{Vector{V}}(N)
    for i=1:N
      bin_keys[i] = Vector{K}(0)
      bin_vals[i] = Vector{V}(0)
    end

    bin_keycount = zeros(Int, N)
    val_sum = zeros(Int, N)

    return new(vlow, vhigh, delta, N, bin_keys, bin_vals, bin_keycount, val_sum,
               comm)
  end
end


function calcBin(obj::UniformBinning, val)

  bin = Int(div(val - obj.vlow + obj.delta, obj.delta))

  # each bin is v_low <= v < v_high, so the value vhigh will be put in bin
  # N + 1.  Force it into bin N
  bin = min(bin, obj.N)

  return bin
end


function push!(obj::UniformBinning, key, val)

  bin = calcBin(obj, val)
  push!(obj.bin_keys[bin], key)
  push!(obj.bin_vals[bin], val)

  return nothing
end


function sumBins(obj::UniformBinning{K, V}) where {K, V}

  # the version of Allreduce with the same send and receive buffers is not
  # wrapped in MPI.jl, so use temporary arrays
  bin_keycount_local = zeros(Int, obj.N)
  val_sum_local = zeros(V, obj.N)

  for i=1:obj.N
    bin_keycount_local[i] = length(obj.bin_keys[i])
    for j=1:length(obj.bin_vals[i])
      val_sum_local[i] += obj.bin_vals[i][j]
    end
  end

  fill!(obj.bin_keycount, 0)
  fill!(obj.val_sum, zero(V))
  MPI.Allreduce!(bin_keycount_local, obj.bin_keycount, MPI.SUM, obj.comm)
  MPI.Allreduce!(val_sum_local, obj.val_sum, MPI.SUM, obj.comm)

  return ROView(obj.bin_keycount), ROView(obj.val_sum)
end
