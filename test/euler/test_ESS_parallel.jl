# test the entropy stable interface calculation in parallel
"""
  Test the entropy stable integrals in parallel
"""
function test_ESS_parallel()
  facts("----- testing ESS parallel -----") do

    if !MPI.Initialized()
      MPI.Init()
    end

    if MPI.Comm_rank(MPI.COMM_WORLD) == 1
      rmfile("./entropy.dat")
    end

    ARGS[1] = "input_vals_ESS_parallel.jl"
    include(STARTUP_PATH)

    data = readdlm("entropy.dat")

    npts, ncols = size(data)

    s_start = data[1, 3]
    s_end = data[end, 3]

    # check total change in entropy function
    @fact abs( (s_start - s_end)/s_start) --> less_than(1e-12)

    # check w^T * res at each timestep
    for i=1:npts
      @fact abs(data[i, 4]) --> less_than(1e-13)
    end

  end

  return nothing
end

#test_ESS_parallel()
add_func1!(EulerTests, test_ESS_parallel)
