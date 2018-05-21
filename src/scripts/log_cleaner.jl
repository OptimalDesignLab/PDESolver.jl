# this script removes duplicate sections from log files, in favor of the lower
#      one (in the file)
# ARGS[1] = file name
# ARGS[2] = column of file to inspect
# ARGS[3] = output file name (optional), if not specified, "_cleaned" is 
#           appended to the file name (before the extension)

using ODLCommonTools

function clean_log(fname::AbstractString, col::Int, output_fname::AbstractString)
 
  print("loading data...")
  data = readdlm(fname)
  println("done")

  @assert eltype(data) <: Number  # no arrays of strings allowed

  if size(data, 1) == 1  # handle special case
    writedlm(output_fname, data)
    return nothing
  end

  rnges = find_src_ranges(data, col)
  dest_rnges = find_dest_ranges(data, rnges, col)

  println("checking overlap")
  check_overlap(data, rnges, dest_rnges, col)

  print("copying data...")
  data2 = copy_blocks(data, rnges, dest_rnges, col)
  println("done")

  print("writing data...")
  writedlm(output_fname, data2)
  println("done")

  return nothing
end

"""
  Returns an array of Range objects that describe the monotonic blocks.
  Monotonicity is measured using the data in column col.

  **Inputs**

   * data: the array of data
   * col: the index of the column

  **Outputs**

   * rnges: array of Range objects containing the indices of each monotonic
            block
"""
function find_src_ranges(data::AbstractArray, col::Integer)

  rnges = Array(UnitRange{Int}, 0)

  start_idx = 1
  for i=2:size(data, 1)
    if data[i, col] <= data[i-1, col]
      rng = start_idx:(i-1)
      push!(rnges, rng)
      start_idx = i
    end
  end

  # add the last range
  rng = start_idx:size(data, 1)
  push!(rnges, rng)

  return rnges
end

"""
  This function figures out the ranges on the destination array corresponding
  to each block in the source array.  Non-used blocks have the range -1:-1.
  Not all unused blocks are guaranteed to have this sentinel value

  **Inputs**

   * data: the data array
   * rnges: source ranges, each must have montonically increasing indices in the
            data arrange
  * col: column of data array to look at

  **Outputs*

   * dest_rnges: ranges in the destination array
"""
function find_dest_ranges(data::AbstractArray, rnges::Array{UnitRange{Int}}, col::Integer)

  # figure out the destination ranges for each block
  # this assumes the blocks will be copies first to last, overwriting previous
  # data

  nblocks = length(rnges)
  dest_rnges = Array(UnitRange{Int}, nblocks)
  dest_rng_idx = 1

  # first range starts at 1
  dest_rng = 1:length(rnges[1])
  dest_rnges[dest_rng_idx] = dest_rng
  dest_rng_idx += 1

  # find where range i should start overwriting range i - 1
  block_offset = 0  # global offset where block - 1 starts
  block_offsets = zeros(Int, nblocks)  # store all of them, so we can
                                       # search older blocks
  blocks_used = Array(Bool, nblocks)
  fill!(blocks_used, true)
  for i=2:length(rnges)
    rng = rnges[i]
    block_offsets[i-1] = block_offset

    # search for first index in next block in the current block
    first_val = data[ first(rng), col]

    # search older blocks until match found
    found_block = 0
    found_idx = 0
    for block=(i-1):-1:1
      vals = sub(data, rnges[block], col)

      # TODO: this could be faster with a log(N) search
      first_idx = findfirst(vals, first_val)
      # if the index is not found, that means the input is not monotonic.
      # We would have to search the earlier blocks for the matching point
      if first_idx != 0
        found_block = block
        found_idx = first_idx
        break
      end
    end

    # construct the offset in the frame of reference of older block
    dest_rng = found_idx:(found_idx + length(rng) - 1)

    # shift it to the global frame of reference
    dest_rng += block_offsets[found_block]
    dest_rnges[dest_rng_idx] = dest_rng
    dest_rng_idx += 1

    # update the offset for the next iteration
    block_offset = block_offsets[found_block] + found_idx - 1
  end

  # TODO: check that the highest value is in the last block
  # I think that might break the algorithm

  vals = sub(data, :, col)
  max_val = maximum(vals)

  vals2 = sub(data, rnges[end], col)
  max_val2 = maximum(vals)

  if max_val != max_val2
    throw(ErrorException("Last interval does not contain the maximum index"))
  end

  # remove any blocks that exceed the length of the destination array
  nrows = last(dest_rnges[end])
  for i=1:nblocks
    rng = dest_rnges[i]
    if last(rng) > nrows
      dest_rnges[i] = -1:-1
    end
  end
  
  return dest_rnges
end

"""
  Copies the blocks from the source array into the destination array, which
  is returned

  **Inputs**

   * data: the data
   * src_rnges: monotonic ranges of data
   * dest_rnges: image of each src range in the destination array
   * col: column of data that is monotonic

  **Outputs**

   * data2: data extracted from data such that it is monotonic in the column
            col
"""
function copy_blocks(data::AbstractArray,
                     src_rnges::Array{UnitRange{Int}},
                     dest_rnges::Array{UnitRange{Int}}, col::Integer)

  # copy the blocks to the new array
  nrows = last(dest_rnges[end])
  data2 = zeros(eltype(data), nrows, size(data, 2))

  nblocks = length(src_rnges)

  for i=1:nblocks
    if dest_rnges[i] == -1:-1
      continue
    end
    src_block = sub(data, src_rnges[i], :)
    dest_block = sub(data2, dest_rnges[i], :)
    copy!(dest_block, src_block)
  end

  return data2
end

"""
  Checks that sections of blocks that are mapped to the same rows in the
  destination array are the same, throws an exception otherwise

  **Inputs**

   * data: the data
   * src_rnges: monotonic ranges of data
   * dest_rnges: image of each src range in the destination array
   * col: column of data that is monotonic
"""
function check_overlap(data::AbstractArray,
                       src_rnges::Array{UnitRange{Int}},
                       dest_rnges::Array{UnitRange{Int}}, col::Integer)

  nblocks = length(src_rnges)
  for i=1:nblocks
#    println("i = ", i)
    src_rng_i = src_rnges[i]
    dest_rng_i = dest_rnges[i]

    if dest_rng_i == -1:-1
      continue
    end
    for j=1:nblocks
      src_rng_j = src_rnges[j]
      dest_rng_j = dest_rnges[j]

      if dest_rng_j == -1:-1
        continue
      end

      overlap = intersect(dest_rng_i, dest_rng_j)
      overlap_start = first(overlap)
      overlap_last = last(overlap)

      if length(overlap) == 0
        continue
      end

      # figure out image of overlap in the two ranges
      begin_offset_i = overlap_start - first(dest_rng_i)
      end_offset_i = last(dest_rng_i) - overlap_last
      image_i = (first(src_rng_i) + begin_offset_i):(last(src_rng_i) - end_offset_i)
      @assert length(image_i) == length(overlap)

      begin_offset_j = overlap_start - first(dest_rng_j)
      end_offset_j = last(dest_rng_j) - overlap_last
      image_j = (first(src_rng_j) + begin_offset_j):(last(src_rng_j) - end_offset_j)

      @assert length(image_j) == length(overlap)
      compare_blocks(data, image_i, image_j, col, 1e-13)

    end
  end

  return nothing
end

"""
  Verifies that the specified rows in data are identical.

  **Inputs**

   * data: the data
   * image_i: first range of rows
   * image_j: second set of rows
   * col: column used to short
   * tol: tolerance for equality test

  **Outputs**

   * none

  This function throws an exception if the rows are not identical.
"""
function compare_blocks(data::AbstractArray, image_i::UnitRange{Int},
                        image_j::UnitRange{Int}, col::Integer,
                        tol::Number)

  @assert length(image_i) == length(image_j)
  # compare col, then the rest of the block
#=
  failflag1 = false
  for i=1:length(image_i)
    diff = abs(data[image_i[i], col] - data[image_j[i], col])
    failflag1 =  (diff > tol) || failflag1
#    @assert diff < tol
  end

  if failflag1
    throw(ErrorException("overlap does not refer to same range"))
  end
=#
  failflag = false
  for j=1:size(data, 2)
    for i=1:length(image_i)
      val = abs(data[ image_i[i], j] - data[ image_j[i], j] ) > tol
      failflag = val || failflag
    end
  end

  if failflag
    println(STDERR, "lines ", image_i, " and ", image_j, " should be identical but are not")
    throw(ErrorException(""))
  end
      
  return nothing
end

#------------------------------------------------------------------------------
# Run the code

if length(ARGS) < 2 || length(ARGS) > 3
  println(STDERR, "Usage: julia clean_log.jl fname column [output_fname]")
end

if length(ARGS) == 2
  fstem, fext = split_fname(ARGS[1])
  output_fname = string(fstem, "_cleaned", fext)
else
  output_fname = ARGS[3]
end

clean_log(ARGS[1], parse(Int, ARGS[2]), output_fname)
