# functions for deciding which elements to refine and how much

"""
  getTargetSizes for adapt_opt.strategy = 1.  Attempts to size elements
  such that the error target is satisfied in a single iteration.

  **Inputs**

   * adapt_opts: `AdaptOpts`
   * mesh
   * sbp
   * eqn
   * opts
   * el_error: the estimated error for each element
   * err_target: the error target

  **Inputs/Outputs**

   * el_sizes: on entry, contains the size of every element.  This arrays
               should be updated with the new element sizes
"""
function getTargetSizes_1(adapt_opts::AdaptOpts, mesh::AbstractMesh,
                         sbp::AbstractOperator,
                         eqn::AbstractSolutionData, opts,
                         el_error::AbstractVector,
                         err_target::Number,
                         el_sizes::AbstractVector)


  # currently we refine everything to equi-distribute the error.  This
  # likely leads to over-refinement

  err_target_el = err_target/mesh.numGlobalEl
  rate::Int = convert(Int, opts["order"] + 1)  # theoretical convergence rate

  for i=1:mesh.numEl
    size_i = el_sizes[i]
    err_i = el_error[i]

    el_sizes[i] = size_i*(err_i/err_target_el)^(-1/rate)
  end

end


"""
  getTargetSizes() for adapt_opts.strategy = 2.  Attempts to refine a fixed
  fraction of elements.

  Because refining exactly the x % of elements sorted by error is difficult
  to do in parallel, this function actually uses a binning algorithm to
  get approximately the top x % of elements.
"""
function getTargetSizes_2(adapt_opts::AdaptOpts, mesh::AbstractMesh,
                         sbp::AbstractOperator,
                         eqn::AbstractSolutionData, opts,
                         el_error::AbstractVector,
                         err_target::Number,
                         el_sizes::AbstractVector)

  bin_elcount, bin_errsum, bins = getBins(el_error, eqn.comm)
  Nbins = length(bin_elcount)
  el_sizes_orig = copy(el_sizes)

  # figure out how many bins
  el_target = Int(round(mesh.numGlobalEl*adapt_opts.fixed_fraction))
  nel = 0
  for i=Nbins:-1:1
    nel += bin_elcount[i]

    # mark these elements for refinement
    for j=1:length(bins.bin_keys[i])
      el_j = bins.bin_keys[i][j]
      el_sizes[el_j] /= 2
    end


    if nel > el_target
      println("intended to refine ", el_target, " elements, actually refined ", nel)
      for i=1:length(el_sizes)
        if el_sizes[i] < 0.99*el_sizes_orig[i]
          println("element ", i, " refined by factor of ", el_sizes_orig[i]/el_sizes[i])
        end
      end
        
      break
    end
  end  # end i

  return nothing
end



#------------------------------------------------------------------------------
# helper functions

"""
  Compute the number of elements in each bin and the total error in each bin.

  **Inputs**

   * el_error: the estimated error for each element
   * comm: MPI communicator

  **Outputs**

   * bin_elcount: number of elements in each bin
   * bin_errsum: total error in each bin
   * bins: a UniformBinning object
"""
function getBins(el_error::AbstractVector, comm::MPI.Comm)

  v_min = minimum(el_error); v_max = maximum(el_error)
  v_min = MPI.Allreduce(v_min, MPI.MIN, comm)
  v_max = MPI.Allreduce(v_max, MPI.MAX, comm)

  Nbins = 100  # arbitrary constant, hopefully this is enough bins to get
               # a good approximation 
  bins = UniformBinning{Int, Float64}(v_min, v_max, Nbins, comm)

  for i=1:length(el_error)
    push!(bins, i, el_error[i])
  end

  bin_elcount, bin_errsum = sumBins(bins)

  return bin_elcount, bin_errsum, bins
end


