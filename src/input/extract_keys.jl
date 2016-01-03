using ODLCommonTools

@doc """
### EulerEquationMod.extractKeys

  This function reads a text file that describes the keys used for the PDESolver
  options dictonary and generates a julia source file declaring a dictonary 
  named `known_keys`.  This dictonary can be used to check the options dictonary
  to see if any unrecognized keys are present.

  The text file must be formatted as follows:
    if a line contains a new key, it must be inclosed in double quotes, and 
    the first character of the line must be a double quote.

    Lines not meeting this specification are ignored.

    Inputs:
      fin: an ASCIIString containing the name (or path) of the file to be read
      fout: name of file to write the dictonary to
    
    Keyword Args:
      header: true if this is the first text file being processed
      footer: true if this is the last text file being processed

    The header and footer options enable reading several text files and writing
    them to a single output file

"""->
function extractKeys(fin::ASCIIString, fout::ASCIIString; header=true, footer=true)

  println("reading file ", fin)
  f = open(fin, "r")
  
  if header 
    println("removing file ", fout)
    rmfile(fout) 
  end
  fw= open(fout, "a")

  if header
    # write the dictonary declaration to the file
    write(fw, "# this file declares a dictonary of all known keys in the
               # options dictonary\n")
    write(fw, "known_keys = Dict{Any, Any}(\n")
  end

  cntr = 1
  for ln in eachline(f)
#    print(cntr, " ", ln)

     if ln[1] == Char(34)  # a double quote
       # find the matching end quote
       end_idx = 0
       for i=2:length(ln)
         if ln[i] == Char(34)
           end_idx = i
           break
         end
       end  # end search for matching quote

       if end_idx == 0
         println(STDERR, "Warning: no matching close quote in line 
                          $cntr of file ", fin)
       end

       # assuming we actually found a matching quote
       key = ln[1:end_idx]
       write(fw, key)
       write(fw, " => true,")
       write(fw, "\n")

     end  # end if first character is a quote

    cntr += 1
  end  # end loop over lines

  if footer
    # close the dictonary declaration
    write(fw, ")")
   end

   close(f)
   close(fw)

end

extractKeys("input_vals.txt", "known_keys.jl", header=true, footer=false)
extractKeys("input_vals_internal.txt", "known_keys.jl", header=false, footer=true)
