# io.jl: functions to do IO more effiently

@doc """
### Utils.BufferedIO

  This type provides a means to buffer IO in memory before writing it to a file.
  Data written to the object is stored in the IOBuffer until flush() is called, 
  when the buffer is (efficiently) dumped into the file
"""->
mutable struct BufferedIO{T <: IO}  <: IO
  fstream::T
  fbuf::IOBuffer
end

@doc """
### Utils.BufferedIO

  Constructor for BufferedIO type.  Takes an IOStream and creates an IOBuffer.
  If no IOStream is given, a dummy stream is created.  If a dummy stream is
  used, flush must never be called on the returned BufferedIO object

  Inputs:

    f: an IOStream object, defaults to a dummy stream

  Outputs:

    a BufferedIO object
"""->
function BufferedIO(f::IO=DevNull)
  buf = IOBuffer()
  fbuf = BufferedIO{typeof(f)}(f, buf)

  # register atexit hook to make sure any buffered data is flushed before
  # julia exits
  # this causes errors to be printed when julia exists, possible #10431
  atexit( () -> if isopen(fbuf.fstream)
                  close(fbuf)
                end
        )
  return fbuf
end


"""
  Alternative constructor for BufferedIO, emulating the open() function.
  This function creates the underlying file using open() and then creates
  a BufferedIO around it.

  Inputs:

    fname: AbstractString, name of file to open
    mode: file open mode, see documentation of open(), defaults to append

  Outputs:

    a BufferedIO object
"""
function BufferedIO(fname::AbstractString, mode::AbstractString="a")

  f = open(fname, mode)
  return BufferedIO(f)
end
# only provide write functionality, for now
import Base.write, Base.flush, Base.close, Base.isopen

# I think this is all that needs to be implemented, because it is the job of
# functions like println and print to convert things to arrays of UInt8s
write(io::BufferedIO, x::UInt8) = write(io.fbuf, x)

"""
  `Base` function `flush` extended for BufferedIO
"""
function flush(io::BufferedIO)
  write(io.fstream, take!(io.fbuf))
  flush(io.fstream)
end

"""
  `Base` function `close` extended for BufferedIO
"""
function close(io::BufferedIO)
  if isopen(io.fstream)
    flush(io)
    close(io.fstream)
  end
end

"""
  `Base` function `isopen` extended for BufferedIO
"""
function isopen(f::BufferedIO)
  return isopen(f.fstream)
end


#------------------------------------------------------------------------------
# create some useful streams

"""
  Buffered version of STDOUT.  This should *always* be used instead of STDOUT
"""
global const BSTDOUT = BufferedIO(STDOUT)

"""
  Buffered version of STDERR.  This should *always* be used instead of STDERR
"""
global const BSTDERR = BufferedIO(STDERR)

#=
#------------------------------------------------------------------------------
# create file io to numbered file with time stamp, unbuffered

"""
  This IO type writes the output to a numbered file (see the constructor),
  without buffering, prefixing the output with a time stamp using MPI_Wtime
  Note that only println() includes the time stamp prefix, print() does not.
"""
type DebugFileIO  <: IO
  fstream::IO  # TODO: make this an IOStream (concrete)
end

@doc """
### Utils.DebugFileIO

  Constructor for DebugFileIO type.  Creates a numbered file by appending
  the supplied rank to the file name (see get_parallel_fname in ODLCommonTools)

  Inputs:
    fname: file name, including extension
    mode: mode flags for opening the file (see docs for open() )
    rank: the rank of the process

  Outputs:
    a DebugFileIO object
"""->
function DebugFileIO(fname::AbstractString, mode::AbstractString, rank::Integer )

  f = open( get_parallel_fname(fname, rank), mode)
  fbuf = DebugFileIO(f)


  # register atexit hook to make sure any buffered data is flushed before
  # julia exits
  # this causes errors to be printed when julia exists, possible #10431
  atexit( () -> if isopen(fbuf.fstream)
                  close(fbuf)
                end
        )
  return fbuf
end

write(io::DebugFileIO, x::UInt8) = write(io.fstream, x)
flush(io::DebugFileIO) = flush(io.fstream)
close(io::DebugFileIO) = close(io.fstream)


import Base.println

function println(io::DebugFileIO, args...)
  tval = MPI.Wtime()
  println(io.fstream, @printf("%16e ", tval), args...)
end
=#

