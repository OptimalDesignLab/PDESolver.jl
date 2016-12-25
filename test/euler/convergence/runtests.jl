# convergence/runtests.jl
# this file runs all the convergence tests

function run_convergence_tests()
  facts("----- testing convergence rates -----") do

    original_dir = pwd()
    cd(dirname(@__FILE__))

    start_dir = pwd()
  #=
    # test p1 elements
    cd("./p1/conservative")
    include(joinpath(pwd(), "runtests.jl"))
  =#
    cd(start_dir)
    cd("./p1/conservative_dg")
    include(joinpath(pwd(), "runtests.jl"))

    cd(start_dir)
    cd("./p1/source_term")
    include(joinpath(pwd(), "runtests.jl"))

    cd(start_dir)
    cd("./p1/conservative_3dg")
    include(joinpath(pwd(), "runtests.jl"))



  #=
    cd(start_dir)
    cd("./p1/entropy")
    include(joinpath(pwd(), "runtests.jl"))
  =#

  #=
    # test p2 elements
    cd(start_dir)
    cd("./p2/conservative")
    include(joinpath(pwd(), "runtests.jl"))
  =#
    cd(start_dir)
    cd("./p2/conservative_dg")
    include(joinpath(pwd(), "runtests.jl"))

    cd(start_dir)
    cd("./p2/conservative_3dg")
    include(joinpath(pwd(), "runtests.jl"))


  #=
    cd(start_dir)
    cd("./p2/entropy")
    include(joinpath(pwd(), "runtests.jl"))
  =#

  #=
    # test p3 elements
    cd(start_dir)
    cd("./p3/conservative")
    include(joinpath(pwd(), "runtests.jl"))
  =#
    cd(start_dir)
    cd("./p3/conservative_dg")
    include(joinpath(pwd(), "runtests.jl"))

  #=
    cd(start_dir)
    cd("./p3/entropy")
    include(joinpath(pwd(), "runtests.jl"))
  =#
  #=
    # test p4 elements
    cd(start_dir)
    cd("./p4/conservative")
    include(joinpath(pwd(), "runtests.jl"))
  =#
    cd(start_dir)
    cd("./p4/conservative_dg")
    include(joinpath(pwd(), "runtests.jl"))

  #=
    cd(start_dir)
    cd("./p4/entropy")
    include(joinpath(pwd(), "runtests.jl"))
  =#
    # return to original directory
    cd(start_dir)

    cd(original_dir)

  end  # end facts block

  return nothing

end

#run_convergence_tests()
add_func1!(EulerTests, run_convergence_tests, [TAG_CONVERGENCE, TAG_NLSOLVERS])
