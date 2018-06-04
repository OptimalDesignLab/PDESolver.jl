# this file lists all the input and output file names for the performance
# tests and then runs them

include("test_funcs.jl")

tests = [ "input_vals_2d_rk4.jl"    "perf_history_2d_rk4.txt";
          "input_vals_2d_rk4_noprecompute.jl"    "perf_history_2d_rk4_noprecompute.txt";
          "input_vals_3d_rk4.jl"    "perf_history_3d_rk4.txt";

          "input_vals_2d_newton.jl" "perf_history_2d_newton.txt";
          "input_vals_2d_newton_coloring.jl" "perf_history_2d_newton_coloring.txt"
          "input_vals_2d_newton_diagE.jl" "perf_history_2d_newton_diagE.txt"
        ]


checkInput(tests)  # make sure no duplicate file names

if length(ARGS) == 0  # run all tests
  for i=1:size(tests, 1)
    runtest(tests[i, 1], tests[i, 2])
  end
else  # run only the test specified (by name)

  for j=1:length(ARGS)

    # find the test
    idx = 0
    for i=1:size(tests, 1)
      if tests[i, 1] == ARGS[j]
        idx = i
        break
      end
    end

    # run the test or print error message
    if idx == 0
      println(STDERR, "Unrecognized test name ", ARGS[j])
      println(STDERR, "Known test names are:")
      for i=1:size(tests, 1)
        println(STDERR, "  ", tests[i, 1])
      end
      error("Unrecognized test")

    else  # found the test
      runtest(tests[idx, 1], tests[idx, 2])
    end

  end  # end loop j
end



