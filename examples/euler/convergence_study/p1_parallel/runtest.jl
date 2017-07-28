using PDESolver  # load main module
using EulerEquationMod  # load physics module 

fname = "input_vals1.jl"  # name of input file in each directory

for i=1:9  # loop over directories

  println("\nrunning convergence study in directory m$i") 

  cd("./m$i")  

  run_solver(fname)  # invoke the solver on the input file in the current
                     # directory

  cd("../")
end


