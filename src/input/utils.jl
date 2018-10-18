# misc utils for reading input files

"""
  Prints an options dictionary to the specified IO, sorting the keys first

  **Inputs**

   * f: IO to print to
   * opts: options dictionary
"""
function printOpts(f::IO, opts::Dict)

  kvals = sort(collect(keys(opts)))

  println(f, "Options dictionary has ", length(kvals), " keys:")
  for k in kvals
    println(f, "  ", k, " => ", opts[k])
  end

  return nothing
end

"""
  Method that prints to STDOUT

  **Inputs**

   * opts: options dictionary
"""
function printOpts(opts::Dict)

  printOpts(STDOUT, opts)
end


