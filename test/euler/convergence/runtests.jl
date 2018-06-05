# convergence/runtests.jl
# this file runs all the convergence tests

include("calc_line.jl")
include("calc_line2.jl")

include("./p1/conservative/runtests.jl")
include("./p1/conservative_3dg/runtests.jl")
include("./p1/conservative_dg/runtests.jl")
include("./p1/entropy/runtests.jl")
include("./p1/source_term/runtests.jl")

include("./p2/conservative/runtests.jl")
include("./p2/conservative_3dg/runtests.jl")
include("./p2/conservative_dg/runtests.jl")
include("./p2/entropy/runtests.jl")

include("./p3/conservative/runtests.jl")
include("./p3/conservative_dg/runtests.jl")
include("./p3/entropy/runtests.jl")

include("./p4/conservative/runtests.jl")
include("./p4/conservative_dg/runtests.jl")
include("./p4/entropy/runtests.jl")





function run_convergence_tests()
  @testset "----- testing convergence rates -----" begin

    original_dir = pwd()
    cd(dirname(@__FILE__))


    start_dir = pwd()
    println("start_dir = ", start_dir)
  #=
    # test p1 elements
    cd("./p1/conservative")
    run_convergence_p1_conservative()
  =#
    cd(start_dir)
    cd("./p1/conservative_dg")
    test_convergence_p1_dg()

    cd(start_dir)
    cd("./p1/source_term")
    test_convergence_p1_source()

    cd(start_dir)
    cd("./p1/conservative_3dg")
    test_convergence_p1_3dg()

  #=
    cd(start_dir)
    cd("./p1/entropy")
    test_convergence_p1_entropy()
  =#

  #----------------------------------------------------------------------------
  # p=2 elements
  #=
    cd(start_dir)
    cd("./p2/conservative")a
    test_convergence_p2_conservative()
  =#
    cd(start_dir)
    cd("./p2/conservative_dg")
    test_convergence_p2_dg()

    cd(start_dir)
    cd("./p2/conservative_3dg")
    test_convergence_p2_3dg()


  #=
    cd(start_dir)
    cd("./p2/entropy")
    test_convergence_p2_entropy()
  =#

  #----------------------------------------------------------------------------
  # p=3 elements
  #=
    cd(start_dir)
    cd("./p3/conservative")
    test_convergence_p3_conservative()
  =#
    cd(start_dir)
    cd("./p3/conservative_dg")
    test_convergence_p3_dg()

  #=
    cd(start_dir)
    cd("./p3/entropy")
    test_convergence_p3_entropy()
  =#

  #----------------------------------------------------------------------------
  # p=4 elements
  #=
    # test p4 elements
    cd(start_dir)
    cd("./p4/conservative")
    test_convergence_p4_conservative()
  =#
    cd(start_dir)
    cd("./p4/conservative_dg")
    test_convergence_p4_dg()

  #=
    cd(start_dir)
    cd("./p4/entropy")
    test_convergence_p4_entropy()
  =#
    # return to original directory
    cd(start_dir)

    cd(original_dir)

  end  # end facts block

  return nothing

end

#run_convergence_tests()
add_func1!(EulerTests, run_convergence_tests, [TAG_CONVERGENCE, TAG_NLSOLVERS, TAG_SHORTTEST])
