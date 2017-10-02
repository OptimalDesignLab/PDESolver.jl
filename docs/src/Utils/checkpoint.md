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

TODO: write a little how-to


## API

### Types

```@docs
Checkpointer
Checkpointer(::Integer, ::ASCIIString)
Checkpointer(::Dict)
copy(::Checkpointer)
copy!(::Checkpointer, ::Checkpointer)
AbstractCheckpointData
```

### Functions
These function provide the basic operations required for checkpointing
```@docs
saveNextFreeCheckpoint
loadLastCheckpoint
readCheckpointData
countFreeCheckpoints
getNextFreeCheckpoint
getLastCheckpoint
getOldestCheckpoint
freeOldestCheckpoint
freeCheckpoint
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
