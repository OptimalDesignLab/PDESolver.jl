# functions for running performance tests

using ODLCommonTools
using PDESolver
using EulerEquationMod
using AdvectionEquationMod
using Utils

"""
  Function that runs a test case and updates the performance history file.
  The case gets run twice, once for warmup, once for the real timing.

  The tests actually get run in a newly created directory called `tmpdir`.
  This directory is deleted if it already exists.  This directory should
  *not* be added to the git repo!

  **Inputs**

   * fname: the input file name
   * outname: the name of the performance history file
"""
function runtest(fname, outname)

  if isdir("./tmpdir")
    rm("./tmpdir", recursive=true)
  end

  mkdir("./tmpdir")
  cd("./tmpdir")
  fname = joinpath("..", fname)
  outname = joinpath("..", outname)

  # run the test in a try block so that even if something goes wrong, we
  # still cd back to the original directory (so we can run tests in a loop)
  try
    println("----- warming up -----")
    val, telapsed, tgc, gc_bytes = @time_all run_solver(fname)

    println("----- Final run -----")
    val, telapsed, tgc, gc_bytes = @time_all run_solver(fname)

    f = open(outname, "a+")
    tstr = getTimeString()
    branchname = getBranchName()
    println(f, tstr, " ", branchname, " ", telapsed, " ", gc_bytes)
    close(f)
  finally
    cd("..")
  end

  return nothing
end

"""
  Function to verify that none of the tests have the same input or output file
  name

  **Inputs**

   * vals: a n x 2 array of Strings, where each line contains the input file
           name followed by the output file name
"""
function checkInput(vals)

  println("size(vals) = ", size(vals))
  for i=1:size(vals, 1)
    println("i = ", i)
    for j=(i+1):size(vals, 1)
      println("j = ", j)
      @assert vals[j, 1] != vals[i, 1]
      @assert vals[j, 2] != vals[i, 2]

    end

    # check that files exist
    for d=1:2
      if !isfile(vals[i, d])
        error( "file $(vals[i, d]) does not exist")
      end
    end
  end



  return nothing
end

#runtest()

