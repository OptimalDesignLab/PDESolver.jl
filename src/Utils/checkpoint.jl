# datatype and API for doing checkpointing

# Enum for checkpoint states
global const CheckpointFree = 1
global const CheckpointUsed = 2

# file name to write the checkpointer to
global const Checkpointer_fname = "checkpointer.dat"

# file name stem for solution files
global const Solution_fname = "qvec"

# file name for input_file
global const Input_fname = "input_vals_restart"

# file name step for user checkpoint data
global const CheckpointData_fname = "checkpointdata"

"""
  Every function that wants to checkpoint must implement an
  AbstractCheckpointData containing all the data that needs to be loaded
  when the code restarts.  The values in the AbstractSolutionData can be
  different for each MPI process.
  
  This type must contain only "julia" data, ie. objects that are managed
  by Julia and have deterministic values.  No pointers, MPI communicators etc.

  File IO is expensive, so include only as much data as needed in the
  AbstractCheckpointData.  Anything that can be recalculated should be.

  Generally, only the outermost function (ie. the time stepper or the
  nonlinear solver for steady problems) will need to checkpoint.

  There are no required fields for this type
"""
abstract AbstractCheckpointData


"""
  This type keeps track of which checkpoints are in use and has an API
  for loading and saving checkpoints

  The fields of the type should *never* be directly accessed.  Use the API
  instead.

  Every checkpoint must be in its own directory.

  Checkpointing has two uses: loading a previous state and restarting.
  The first one is rather easy, all that needs to happen is the solution
  variables get loaded.  Restarting is more involved because all the local
  data of the NonlinearSolver that was running needs to be stored and then
  loaded again in a new session.

  An input file called $(Input_fname) is created in the current directory
  that can be used to restart the code.

  **Fields**

   * ncheckpoints: total number of checkpoints
   * paths: absolute paths to the checkpoints (array of strings)
   * status: status of each checkpoint (used or free), (array of Ints)
   * history: list of checkpoints from most recently used to least recently
              used, unused entries set to -1

"""
type Checkpointer
  ncheckpoints::Int  # number of checkpoints
  paths::Array{String, 1}  # paths to the directories
  status::Array{Int, 1}  # is checkpoint free
  history::Array{Int, 1}  # list of checkpoints, from most recently used
                          # (first) to least recently used
                          # unused entries are set to -1
end

#------------------------------------------------------------------------------
# outer constructors

"""
  Constructs a Checkpointer with a given number of checkpoints.

  The checkpoints are put into directories called "checkpoint*", where the
  "*" is replaced by the index of the checkpoint (from 1 to ncheckpoints).
  This constructor does not error if the directories already exist, however
  any data in the checkpoints may be overwritten.

  **Inputs**

   * myrank: MPI rank of this process

   * ncheckpoints: number of checkpoints, defaults to 2.  Using fewer than
                   2 is not recommended (there will be a moment in time during
                   which the old checkpoint has been partially overwritten but
                   the new checkpoint is not complete yet)

   * prefix: added to directory names (with an underscore added in-between).
              Defaults to empty string. This is useful if there are multiple
              sets of checkpoints (and therefore multiple Checkpointer objects).

  **Outputs**

   * a Checkpointer object, fully initialized
"""
function Checkpointer(myrank::Integer, ncheckpoints::Integer=2,
                      prefix::String="")

  paths = Array(String, ncheckpoints)
  status = Array(Int, ncheckpoints)
  history = Array(Int, ncheckpoints)
  fill!(history, -1)

  if prefix != ""
    prefix = string(prefix, "_")
  end

  # give initial values: all checkpoints free
  for i=1:ncheckpoints

    # create an absolute path in case the user calls cd() somewhere else in
    # the code
    paths[i] = joinpath_ascii(pwd(), "$(prefix)checkpoint$i")
    status[i] = CheckpointFree

    if !isdir(paths[i]) && myrank == 0
      mkdir(paths[i])
    end

  end

  checkpointer = Checkpointer(ncheckpoints, paths, status, history)

  # make sure the flag files don't exist (perhaps left over from a previous
  # run)


  if myrank == 0
    for i=1:ncheckpoints
      deleteFlagFile(checkpointer, i)
    end
  end


  return checkpointer
end

"""

  This constructor loads a Checkpointer object from a file.  This should
  be used for restarting.  Do *not* use the other constructor for restarting.
  It will load the most recently written checkpoint that is complete (an
  important distinction if the code was killed in the middle of writing
  a checkpoint)

  **Inputs**

   * opts: the options dictionary.  

  **Outputs**

   * a Checkpointer object, fully initialized and aware of which checkpoints
     are in use and which are free

  This function only loads the Checkpointer from the checkpoint.
  See [`loadLastCheckpoint`](@ref) and [`readCheckpointData`](@ref) to
  load the rest of checkpoint data.
  The checkpoint that was loaded can be accessed via [`getLastCheckpoint`](@ref)


  **Implementation notes:**

  Uses 4 keys in the options dictionary:

   * "writing checkpoint": index of checkpoint that might be complete
   * "writing_checkpoint_path": absolute path of checkpoint
   * "most_recent_checkpoint": index of checkpoint that is definitely
                                complete, but might be older than the above
   * "most_recent_checkpoint_path": absolute path of checkpoint

  This system ensure it is possible to restart even if the code is
  killed in the middle of writing a checkpoint.

  The Checkpointer object is the same on all processes.
"""
function Checkpointer(opts::Dict, myrank::Integer)

  # figure out whether checkpoint that was in the process of being written
  # finished.
  # if not, use the checkpoint before that.
  checkpoint1_path = opts["writing_checkpoint_path"]
  checkpoint2_path = opts["most_recent_checkpoint_path"]

  if checkFlagFile(checkpoint1_path)
    if myrank == 0
      println(BSTDOUT, "Loading checkpoint most recent checkpoint")
    end
    checkpoint = opts["writing_checkpoint"]
    checkpoint_path = checkpoint1_path
  else
    @assert checkFlagFile(checkpoint2_path)
    if myrank == 0
      println(BSTDOUT, "Most recent checkpoint incomplete, loading older checkpoint")
    end
    checkpoint = opts["most_recent_checkpoint"]
    checkpoint_path = checkpoint2_path
  end

  if myrank == 0
    println(BSTDOUT, "checkpoint path = ", checkpoint_path)
  end
  fname = joinpath_ascii(checkpoint_path, Checkpointer_fname)
  f = open(fname, "r")
  checkpointer = deserialize(f)
  close(f)

  return checkpointer::Checkpointer
end

import Base.copy, Base.copy!

"""
  Recursively copies all fields to make a new Checkpointer object.
  Note that changes to one Checkpointer object will not affect the other
  (ie. the record of which checkpoints are used and which are not).
  This could easily lead to corrupting a checkpoint.
  For this reason, this function should rarely be used.
"""
function copy(chkpointer::Checkpointer)

  ncheckpoints = chkpointer.ncheckpoints
  paths = copy(chkpointer.paths)
  for i=1:ncheckpoints
    paths[i] = copy(chkpointer.paths[i])
  end
  status = copy(checkpointer.status)
  history = copy(checkpointer.history)

  return Checkpointer(ncheckpoints, paths, status, history)
end

"""
  2 argument version of copy().  See that function for details.
"""
function copy!(dest::Checkpointer, src::Checkpointer)

  dest.ncheckpoints = src.ncheckpoints
  dest.paths = copy(src.paths)
  for i=1:dest.ncheckpoints
    dest.paths[i] = copy(src.paths[i])
  end
  dest.status = copy(src.status)
  dest.history = copy(src.history)

  return nothing
end
#------------------------------------------------------------------------------
# Internal functions

"""
  Call on master process only

  Does *not* check if the checkpoint is already used (for reading the
  state back after restart, this checkpoint should already be marked as
  used)

"""
function writeCheckpointer(chkpointer::Checkpointer, checkpoint::Int)

  fname = joinpath_ascii(chkpointer.paths[checkpoint], Checkpointer_fname)
  f = open(fname, "w")
  serialize(f, chkpointer)
  close(f)

  return nothing
end

"""
  Save to a specified checkpoint.  Throw error if checkpoint is not free.
  Users should not generally call this function directly.  Instead, they
  should prefer [`saveNextFreeCheckpoint`](@ref).

  This function automatically saves eqn.q_vec and the checkpointer to
  a file.  Any additional data should be in checkpoint_data.

  Note: mesh adaptation is not compatable with checkpointing
  #TODO; add a field to the mesh to record that it has been modified

  **Inputs**

   * checkpointer: the CheckPointer
   * checkpoint: the index of the checkpoint
   * mesh: an AbstractMesh object
   * sbp: SBP operator
   * eqn: an AbstractSolutionData
   * checkpoint_data: an AbstractCheckpointData.  This is the random bag
                      of data the user needs saved.

   **Inputs/Outputs**

    * opts: options dictionary


  Implementation Notes:

  Uses options dictionary keys described by [`Checkpointer`](@ref)
  Note that the checkpoint is eagerly marked as used, before finishing writing
  the checkpoint.  Upon restart the code needs to check if this checkpoint
  is really finished.
"""
function saveCheckpoint(checkpointer::Checkpointer, checkpoint::Int,
                         mesh::AbstractMesh,
                         sbp::AbstractSBP, eqn::AbstractSolutionData,
                         opts, checkpoint_data::AbstractCheckpointData)

  @assert checkpoint > 0
  @assert checkpoint <= checkpointer.ncheckpoints

  if checkpointer.status[checkpoint] != CheckpointFree
    throw(ErrorException("Cannot save to non-free checkpoint number $checkpoint"))
  end

  # wait for all processes to catch up before overwriting whatever data
  # might be in the checkpoint
  MPI.Barrier(mesh.comm)
  if mesh.myrank == 0
    deleteFlagFile(checkpointer, checkpoint)
  end

  # record the state change now, because we have started to overwrite things
  # the options dictionary handles the case of the program getting killed
  # while the checkpoint is being written
  markCheckpointUsed(checkpointer, checkpoint)

  # record which checkpoint to restart from
  opts["writing_checkpoint"] = checkpoint
  opts["writing_checkpoint_path"] = checkpointer.paths[checkpoint]

  # write the solution to the file
  # note: the residual is not saved
  solution_fname = joinpath_ascii(checkpointer.paths[checkpoint], Solution_fname)
  writeSolutionFiles(mesh, sbp, eqn, opts, solution_fname)

 
  # write the user data
  writeCheckpointData(checkpointer, checkpoint, checkpoint_data, mesh.myrank)


  # write input file
  if mesh.myrank == 0
    # create a file suitable for restarting
    opts_restart = copy(opts)
    opts_restart["is_restart"] = true  # general flag that all parts of
                                       # the code will have to look for

    opts_fname = joinpath_ascii(checkpointer.paths[checkpoint], Input_fname)
    opts_fname = make_input(opts_restart, opts_fname)

    # symlink it into the current directory so we know which checkpoint to
    # load on restart
    try 
      rm(Input_fname)
    end
    symlink(opts_fname, Input_fname)


    #TODO: deal with mesh adaptation
    # if only coordinates have changed, need to make sure that all the
    # solver-visible parts of the mesh object do not depend on the original
    # mesh coordinates (I think deciding elementL vs elementR for interfaces
    # does)
    # OR: write mesh coordinates to a file, load original mesh, then adapt
    # it.
    # it would be faster to have a flag in the mesh determining if it has been
    # adapted, to avoid needing to write the coordinates and read them back
    # in again

    # write checkpointer to file
    # this is needed in the case of restarting an unsteady adjoint computation

    writeCheckpointer(checkpointer, checkpoint)
  end

  # wait for all processes to finish writing
  MPI.Barrier(mesh.comm)

  # write flag file
  if mesh.myrank == 0
    writeFlagFile(checkpointer, checkpoint)
  end
  opts["most_recent_checkpoint"] = checkpoint
  opts["most_recent_checkpoint_path"] = checkpointer.paths[checkpoint]


  return nothing
end


"""
  Loads a specified checkpoint.  Only loads eqn.q_vec.  Does not load
  the checkpointer (because this is loading a checkpoint and not a restart).
  Also does not load the AbstractCheckpointData.
  See [`readCheckpointData`](@ref) for loading of the AbstractCheckpointData.
  Users should generally not call this function directly.  Users should
  prefer [`loadLastCheckpoint`](@ref) whenever possible.

  Note that loading a checkpoint does not mark it as free.  Users must
  explictly call [`freeCheckpoint`](@ref).

  **Inputs**

   * checkpointer: the Checkpointer object
   * checkpoint: the index of the checkpoint
   * mesh: an AbstractMesh object
   * sbp: the SBP operator
   * opts: the options dictionary

  **Inputs/Output**

   * eqn: an AbstractSolutionData, eqn.q_vec is overwritten
"""
function loadCheckpoint(checkpointer::Checkpointer, checkpoint::Int,
                        mesh::AbstractMesh,
                        sbp::AbstractSBP, eqn::AbstractSolutionData,
                        opts)

  @assert checkpoint > 0
  @assert checkpoint <= checkpointer.ncheckpoints

  if checkpointer.status[checkpoint] != CheckpointUsed
    throw(ErrorException("Cannot load unused checkpoint number $checkpoint"))
  end


  solution_fname = joinpath_ascii(checkpointer.paths[checkpoint], Solution_fname)
  readSolutionFiles(mesh, sbp, eqn, opts, solution_fname)
  
  # we can't load the AbstractCheckpointData here in a type stable way

  return nothing
end

# file name of flag file
global const FlagFile_fname = "flagfile.dat"

# value to write to flag file
global const FlagFile_val = Int64(42)

"""
  Write the flag file that signifies a checkpoint is complete
  Call on master process only!

  **Inputs**

   * checkpointer: the Checkpointer object
   * checkpoint: the checkpoint index
"""
function writeFlagFile(checkpointer::Checkpointer, checkpoint::Integer)

  fpath = joinpath_ascii(checkpointer.paths[checkpoint], FlagFile_fname)
  f = open(fpath, "w")
  write(f, FlagFile_val)
  close(f)

  return nothing
end

"""
  Returns true if the flag file exists and is consistent, returns false
  otherwise

  **Inputs**

   * checkpointer: the Checkpointer object
   * checkpoint: the checkpoint index
"""
function checkFlagFile(checkpointer::Checkpointer, checkpoint::Integer)

  fpath = checkpointer.paths[checkpoint]
  return checkFlagFile(fpath)
end

"""
  Sometimes need to check the flag file before the Checkpointer is available.
  See the other method ofr details.

  **Inputs**

   * fpath: path to the checkpoint directory
"""
function checkFlagFile(fpath::AbstractString)

  fpath = joinpath_ascii(fpath, FlagFile_fname)
  if !isfile(fpath)
    return false
  end

  f = open(fpath, "r")
  val = read(f, typeof(FlagFile_val))
  close(f)

  if val == FlagFile_val
    return true
  else
    return false
  end
end

"""
  Deletes the flag file if it exists.  Does not error if flag file does not
  exist.
  Call on master process only

  **Inputs**

   * checkpointer: the Checkpointer object
   * checkpoint: the checkpoint index
"""
function deleteFlagFile(checkpointer::Checkpointer, checkpoint::Integer)

  fpath = joinpath_ascii(checkpointer.paths[checkpoint], FlagFile_fname)
  if isfile(fpath)
    rm(fpath)
  end

  return nothing
end

"""
  Writes the AbstractCheckpointData to a file.  The file can be read by
  [`readCheckpointData`](@ref).

  **Inputs**

   * checkpoint: Checkpointer object
   * checkpoint: the checkpoint index
   * obj: the AbstractCheckpointData object
   * comm_rank: the MPI rank of this process
"""
function writeCheckpointData(chkpointer::Checkpointer, chkpoint::Integer,
                             obj::AbstractCheckpointData, comm_rank::Integer)

  fname = string(get_parallel_fname(CheckpointData_fname, comm_rank), ".dat")
  fpath = joinpath_ascii(chkpointer.paths[chkpoint], fname)

  f = open(fpath, "w")
  serialize(f, obj)
  close(f)

  return nothing
end

"""
  Marks a checkpoint as unused.

  **Inputs**

   * checkpointer: the Checkpointer object
   * checkpoint: the checkpoint index
"""
function markCheckpointUsed(chkpointer::Checkpointer, chkpoint::Integer)

  @assert chkpointer.status[chkpoint] == CheckpointFree

  chkpointer.status[chkpoint] = CheckpointUsed

  # shift everything right by one
  @assert chkpointer.history[end] == -1  # there has to be at least one unused
  last_idx = 0
  for i=1:chkpointer.ncheckpoints
    if chkpointer.history[i] == -1
      last_idx = i
      break
    end
  end
 


  for i=(last_idx-1):-1:1
    chkpointer.history[i+1] = chkpointer.history[i]
  end

  # add the checkpoint to the head of the line
  chkpointer.history[1] = chkpoint

  return nothing
end

#------------------------------------------------------------------------------
# API

"""
  This function reads an object of type AbstractSolutionData from a file
  and reuturns the object.  This is type unstable, so every AbstractCheckpointData
  implementation should create a constructor that calls this function and
  uses a type assertion to specify that the object must be of their concrete
  type.

  **Inputs**

   * chkpointer: a Checkpointer
   * chkpoint: the index of hte checkpoint
   * comm_rank: the MPI rank of this process

  **Outputs**

   * the [`AbstractCheckpointData`](@ref) loaded from the checkpoint
"""
function readCheckpointData(chkpointer::Checkpointer, chkpoint::Integer,
                            comm_rank::Integer)

  @assert chkpoint > 0
  @assert chkpoint <= chkpointer.ncheckpoints

  fname = string(get_parallel_fname(CheckpointData_fname, comm_rank), ".dat")
  fpath = joinpath_ascii(chkpointer.paths[chkpoint], fname)

  f = open(fpath, "r")
  obj = deserialize(f)
  close(f)

  return obj
end

"""
  Simple wrapper around [`readCheckpointData`](@ref) to loading the most
  recently saved checkpoint data

  **Inputs**

   * chkpointer: a Checkpointer
   * comm_rank: MPI rank of this process
"""
function readLastCheckpointData(chkpointer::Checkpointer, comm_rank::Integer)

  chkpoint = getLastCheckpoint(chkpointer)
  obj = readCheckpointData(chkpointer, chkpoint, comm_rank)

  return obj
end


"""
  This function saves to the next free checkpoint

  **Inputs**

   * checkpointer: the Checkpointer object
   * mesh: an AbstractMesh
   * sbp: an SBP operator
   * eqn: an AbstractSolutionData object
   * checkpoint_data: the users AbstractCheckpointData implementation

  **Inputs/Outputs**

   * opts: the options dictionary

  **Outputs**

   * checkpoint: the index of the checkpoint saved.
"""
function saveNextFreeCheckpoint(checkpointer::Checkpointer,
                                mesh::AbstractMesh,
                                sbp::AbstractSBP, eqn::AbstractSolutionData,
                                opts, checkpoint_data::AbstractCheckpointData)

  checkpoint = getNextFreeCheckpoint(checkpointer)
  saveCheckpoint(checkpointer, checkpoint, mesh, sbp, eqn, opts, checkpoint_data)

  return checkpoint

end

"""
  This function loads the most recently saved checkpoint.

  **Inputs**

   * checkpointer: the CheckPointer object
   * mesh: an AbstractMesh
   * sbp: an SBP operator
   * opts: the options dictionary

  **Inputs/Outputs**

   * eqn: an AbstractSolutionData object

  **Outputs**

   * the checkpoint index that was loaded
"""
function loadLastCheckpoint(checkpointer::Checkpointer,
                            mesh::AbstractMesh,
                            sbp::AbstractSBP, eqn::AbstractSolutionData,
                            opts)

  checkpoint = getLastCheckpoint(checkpointer)
  loadCheckpoint(checkpointer, checkpoint, mesh, sbp, eqn, opts)

  return checkpoint
end

"""
  This function returns the number of free checkpoints

  **Inputs**

   * checkpointer: the Checkpoint object
"""
function countFreeCheckpoints(checkpointer::Checkpointer)

  cnt = 0
  for i=1:checkpointer.ncheckpoints
    if checkpointer.status[i] == CheckpointFree
      cnt += 1
    end
  end

  return cnt
end

"""
  This function returns the index of the next free checkpoint,  or 0 if
  there are no free checkpoints

  **Input**

   * checkpointer: a Checkpointer
"""
function getNextFreeCheckpoint(checkpointer::Checkpointer)

  for i=1:checkpointer.ncheckpoints
    if checkpointer.status[i] == CheckpointFree
      return i
    end
  end

  return 0  # no free checkpoints
end


"""
  This function returns the index of the most recently written checkpoint

  **Inputs**
  
   * checkpointer: the Checkpointer object

  **Outputs**

   * the index of the most recently saved checkpoint
"""
function getLastCheckpoint(checkpointer::Checkpointer)

  return checkpointer.history[1]
end

"""
  Get the least recently written checkpoint

  **Inputs**

   * checkpointer: the Checkpointer

  **Outputs**

   * the index of the least recently written checkpoint
"""
function getOldestCheckpoint(checkpointer::Checkpointer)

  last_idx = checkpointer.ncheckpoints
  for i=checkpointer.ncheckpoints:-1:1
    if checkpointer.history[i] != -1
      last_idx = i
      break
    end
  end

  return checkpointer.history[last_idx]
end

"""
  Frees the least recently written checkpoint. Calling this function when
  all checkpoints are free is allowed.

  **Inputs**

   * checkpointer: the Checkpointer

  **Outputs**

   * returns the checkpoint freed (0 all checkpoints were free on entry)
"""
function freeOldestCheckpoint(checkpointer::Checkpointer)

  checkpoint = getOldestCheckpoint(checkpointer)
  if checkpoint != -1
    freeCheckpoint(checkpointer, checkpoint)
  end

  return checkpoint
end

"""
  This function frees a checkpoint (marks it as available to be overwritten)
  Unlike certain algorithms that use free lists (cough, cough, malloc)
  freeing an already free checkpoint is allowed.

  **Inputs**

   * checkpointer: the Checkpointer object
   * checkpoint: the index of the checkpoint to free

  The user must explictly free checkpoints (loading a checkpoint does not
  free it).
  Note that a free checkpoint can still be loaded for restarting if it has
  not been saved to yet.
"""
function freeCheckpoint(checkpointer::Checkpointer, checkpoint::Integer)

  @assert checkpoint > 0
  @assert checkpoint <= checkpointer.ncheckpoints

  if checkpointer.status[checkpoint] == CheckpointFree
    # nothing to do
    return nothing
  end

  # remove the checkpoint from the history
  checkpoint_idx = 0
  for i=1:checkpointer.ncheckpoints
    if checkpointer.history[i] == checkpoint
      checkpoint_idx = i
      break
    end
  end

  last_idx = 0
  for i=1:checkpointer.ncheckpoints
    if checkpointer.history[i] == -1
      last_idx = i - 1
      break
    end
  end

  # shift all entries after the checkpoint to be removed left by one
  for i=checkpoint_idx:last_idx
    checkpointer.history[i] = checkpointer.history[i+1]
  end

  # if we are freeing the last checkpoint, last_idx == 0
  if last_idx != 0
    checkpointer.history[last_idx] = -1  # the last one is now unused
  else
    checkpointer.history[end] = -1  
  end
  checkpointer.status[checkpoint] = CheckpointFree

  return nothing
end
