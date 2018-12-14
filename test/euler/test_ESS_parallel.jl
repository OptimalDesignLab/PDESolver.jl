# test the entropy stable interface calculation in parallel
"""
  Test the entropy stable integrals in parallel
"""
function test_ESS_parallel()
  @testset "----- testing ESS parallel -----" begin

    if MPI.Comm_rank(MPI.COMM_WORLD) == 1
      rmfile("./entropy.dat")
    end

    fname = "input_vals_ESS_parallel.jl"
    mesh, sbp, eqn, opts = solvePDE(fname)
    data = readdlm("entropy.dat")

    npts, ncols = size(data)

    s_start = data[1, 3]
    s_end = data[end, 3]

    # check total change in entropy function
    @test  abs( (s_start - s_end)/s_start)  < 1e-12

    # check w^T * res at each timestep
    for i=1:npts
      @test  abs(data[i, 4])  < 1e-13
    end

  end

  return nothing
end

#test_ESS_parallel()
add_func1!(EulerTests, test_ESS_parallel, [TAG_SHORTTEST])
