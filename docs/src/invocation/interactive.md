# Running PDESolver in an interactive session using Julia's REPL

Julia's REPL is a powerful tool for debugging. 
Being able to run commands interactively can reveal important information 
 that would otherwise be available only through cumbersome print statements and
 error logs.

This page describes the first steps involved for running PDESolver interactively.
Calling PDESolver from Julia's REPL in the fashion described here is an experimental 
  method, and is not tested consistently.
Additionally, not all parts of PDESolver are usable after the steps below.
Instead, this document is intended to provide a springing-off point for the developer 
  who wishes to run PDESolver interactively, and can adapt the commands below to their needs.

# Script

```
using PDESolver
using AdvectionEquationMod
using ArrayViews
using ODLCommonTools
using SummationByParts
using PdePumiInterface
using NonlinearSolvers
using ForwardDiff
using Utils
import ODLCommonTools.sview
using MPI
using Input

# after here, hopefully this is the absolute minimum commands needed
MPI.Init()

# example for input file; substitute your own
opts = read_input("sine_adv_input_CN_adjoint.jl")

Tdim = opts["dimensions"]
# note: funny character in opts, last entry when loaded in REPL
opts["Tdim"] = Tdim
dofpernode = 1

sbp, mesh, pmesh, Tsol, Tres, Tmsh, mesh_time = createMeshAndOperator(opts, dofpernode)

var_type = opts["variable_type"]

# initialize the eqn object with the desired physics
# eqn = EulerData_{Tsol, Tres, Tdim, Tmsh, var_type}(mesh, sbp, opts)
eqn = AdvectionData_{Tsol, Tres, Tdim, Tmsh}(mesh, sbp, opts)

init(mesh, sbp, eqn, opts)

# sample commands for populating some fields of eqn
rand!(eqn.q)
rand!(eqn.res)
# q_vec and res_vec are populated after this too
```

# REPL debugging commands

Here are some helpful Julia commands for debugging/investigation within the REPL at this stage.
Refer to the Julia documentation for details, or press '?' within the REPL 
  and type the command you wish to see help for.

```
size()
length()
typeof()
fieldnames()
isimmutable()
whos()
```

