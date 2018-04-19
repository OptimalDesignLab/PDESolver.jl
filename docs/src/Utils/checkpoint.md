# Checkpointing

```@meta
  CurrentModule = Utils
```

Checkpointing has two use cases: saving the state of the solver to be loaded
later (for example, in unsteady adjoint calculations), and to restart the
solver after a crash.
This checkpointing functionality described here is useful for both of these
functions.
Its purpose is to provide an interface for writing the current state to a file
that can be read back later.
As long as at least 2 checkpoints are saved, the implementation guarantees
that at least one checkpoint is loadable at any time, even if the code is
terminated while writing a checkpoint.

Log files require special handling when restarting.
Log files should be opened in append mode when restarting (in fact, the code
usually opens log files in append mode even when not restarting).
One of the consequences of restarting is that there may be
duplicate lines in log files.  For example, if a code is doing a run of 10,000
timesteps, checkpointing every 1,000 timesteps, and gets terminated at timestep
6,500, it can be restarted from the checkpoint at timestep 6,000 and run to the
end.  This means any output written between timesteps 6,000 and 6,500 will 
appear twice in the log files.
To deal with this problem, the script `log_cleaner.jl` (located in `src/scripts`)
can be used to post-process delimited data files (space, comma, tab, etc.) and remove
duplicate entries in favor of the entry appearing last in the log file.
It also checks removed entries to make sure they are the same (to floating point
tolerance) as the entries that are preserved.


To facilitate restarting, the checkpointing code creates a new input file called `input_vals_restart` in
the current working directory.
This file can be used to continue running the code from the most recent
checkpoint.  For example

```
  julia /path/to/startup.jl "input_vals_restart"
```

## API

### Types

```@docs
Checkpointer
Checkpointer(::Integer, ::Integer, ::ASCIIString)
Checkpointer(::Dict, ::Integer)
copy(::Checkpointer)
copy!(::Checkpointer, ::Checkpointer)
AbstractCheckpointData
```

### Functions
These function provide the basic operations required for checkpointing
```@docs
saveNextFreeCheckpoint
loadLastCheckpoint
readLastCheckpointData
readCheckpointData
countFreeCheckpoints
getNextFreeCheckpoint
getLastCheckpoint
getOldestCheckpoint
freeOldestCheckpoint
freeCheckpoint
```
## Example

This section provides a simple example of how to use checkpointing for an explicit
time marching scheme.
The original scheme (without checkpointing) looks like:

```julia
function explicit_timemarching_example(mesh, sbp, eqn, opts, func::Function)

  nsteps = round(Int, opts["tmax"]/opts["delta_t"])
  t = 0.0

  for step=1:nsteps

    t = (step - 1)*delta_t
    # q_vec -> q
    disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    func(mesh, sbp, eqn, opts, t)

    # res -> res_vec
    assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)

    for i=1:mesh.numDof
      eqn.q_vec[i] = ... # some function of eqn.res_vec[i]
    end
  end  # end loop over timesteps

  return nothing
end
```

With checkpointing, it becomes:

```julia

type ExampleCheckpointData <: AbstractCheckpointData
  step::Int  # for this method, the only thing that needs to be stored is the
             # current time step, everything else can be recalculated
end

# this constructor is used when restarting, see below
function ExampleCheckpointData(chkpointer::Checkpointer, comm_rank)

  chkpoint_data = readLastCheckpointData(chkpointer, comm_rank)

  return chkpoint_data::ExampleCheckpointData
end

function explicit_timemarching_example(mesh, sbp, eqn, opts, func::Function)

  nsteps = round(Int, opts["tmax"]/opts["delta_t"])
  t = 0.0

  if !opts["is_restart"]
    # regular starting of a run
    stepstart = 1
    chkpointdata = ExampleCheckpointData(stepstart)
    # this is valid even if we don't intend to checkpoint (ncheckpoints = 0)
    chkpointer = Checkpointer(mesh.myrank, opts["ncheckpoints"]
    skip_checkpoint = false
  else
    chkpointer = Checkpointer(opts, myrank)  # read Checkpointer from most recent
                                             # checkpoint
    # now get the checkpoint data telling which timestep we are on
    # eqn.q_vec already holds the state saved in the checkpoint.  This is handled
    # during the initialization of the equation object
    chkpointdata = ExampleCheckpointData(chkpoinnter, mesh.myrank)
    stepstart = chkpointdata.step
    skip_checkpoint = true  # skip writing the first checkpoint.  Without this, a
                            # checkpoint would be written immediately.
                            # Writing the checkpoint is not harmful, but not useful
                            # either.
   end
    
  for step=stepstart:nsteps  # start the loop from stepstart

    t = (step - 1)*delta_t   # t is recalculated from step, which was set using
                             # stepstart, thus only stepstart needs to be saved

    # save checkpoint, if needed
    if opts["use_checkpointing"] && step % opts["checkpoint_freq"] == 0 && !skip_checkpoint
      skip_checkpoint = false

      # save all needed variables to the ExampleCheckpointData object
      # For simple time marching schemes, step is the only needed data
      chkpointdata.step = step

      if countFreeCheckpoints(chkpointer) == 0
        freeOldestCheckpoint(chkpointer)  # make room for a new checkpoint
      end

      saveNextFreeCheckpoint(chkpointer, mesh, sbp, eqn, opts, chkpointdata)
    end

    # q_vec -> q
    disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    func(mesh, sbp, eqn, opts, t)

    # res -> res_vec
    assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)

    for i=1:mesh.numDof
      eqn.q_vec[i] = ... # some function of eqn.res_vec[i]
    end
  end  # end loop over timesteps

  return nothing
end
```
  


## Internal Functions
The internal functions used for checkpointing are documented here.
Users should not call these functions directly.  Improper use can cause
checkpoint corruption.

```@docs
  writeCheckpointer
  saveCheckpoint
  loadCheckpoint
  writeFlagFile
  checkFlagFile
  deleteFlagFile
  writeCheckpointData
  markCheckpointUsed
```
