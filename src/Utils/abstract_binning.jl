# abstract type for binning

import Base.push!

# K = key type, V = value type
"""
  Implements a binning algorithm.  The algorithm takes key-value pairs and
  sorts them into bins according to the values.

  **Static Parameters**

   * K: the type of the keys
   * V: the type of the values
"""
abstract type AbstractBinning{K, V} end


"""
  Compute the bin a given value falls in

  **Inputs**

   * obj: `AbstractBinning`
   * val: the value
"""
function calcBin(obj::AbstractBinning, val)

  error("abstract fallback for calcBin() reached")
end


"""
  Add a new key-value pair to the binning

  **Inputs**

   * obj: `AbstractBinning`
   * key: the key
   * val: the value
"""
function push!(obj::AbstractBinning, key, val)

  error("abstract fallback for calcBin() reached")
end


"""
  Computes the total number of keys and the sum of the values in each bin.
  This is a global collective operation.

  **Inputs**

   * obj: `AbstractBinning`

  **Outputs**

   * nkeys: vector of length `nBins`, containing the total number of keys in
            each bin
   * valsum: vector of length `nBins` the sum of the values in each bin
"""
function sumBins(obj::AbstractBinning)

  error("abstract fallback for sumBins() reached")
end


