# Input
A PDESolver execution is controlled by an options dictonary.  The user must
supply this dictonary in the form of a Julia source file that declares a 
`Dict{Any, Any}` called `arg_dict`. The file path (relative to the users `pwd` 
must be passed to the solver as the the first optional argument.  For example,
to launch a PDESolver run, execute

```
julia /path/to/startup.jl "path/to/dictonary"
```

Descriptions of the valid keys and values can be found in the file 
input_vals.txt.
PDESolver supplies default values for all keys for which there is a sane
default.  Values that might be unused are initialized to -1 or "none".

## Conventions
Physics modules generally use the `majorIterationCallback` function to
log important quantities to files.  Such logging should be controlled by
two keys, `"write_outname" where `outname` is the name of the quantity, which
has a boolean value, and `"write_outname_fname"` that has a string value
containing the name of the file to write (including extension).  Examples
out things that can be logged are entropy and kinetic energy.  Both these keys
should have default values, and users should generally not need to modify the
second one.

##Key Validation
After supplying default values, PDESolver checks that all keys in the dictonary
are recognized keys.  It does this by comparing against the list of keys
documented in the input_vals.txt file.  A warning of is printed to STDERR if
an unrecognized key is found.

The mechanics of the key validation are as follows.  The `extract_keys.jl`
script reads the input_vals.txt file (and the input_vals_internal.txt), looking for any string that is surrounded
by double quotes and starts in the first character of a line.  All keys are 
written to a dictonary in the file `known_keys.jl`.  This file is included by
PDESolver.

The script `extract_keys.jl` is run as part of PDESolver installation.  If 
new keys are added it must be run again to silence warnings during key 
validation.
