# io.jl: functions to do IO more effiently

export BufferedIO

@doc """
### Utils.BufferedIO

  This type provides a means to buffer IO in memory before writing it to a file.
  Data written to the object is stored in the IOBuffer until flush() is called, 
  when the buffer is (efficiently) dumped into the file
"""->
type BufferedIO{T <: IO}  <: IO
  fstream::T
  fbuf::IOBuffer

end

@doc """
### Utils.BufferedIO

  Constructor for BufferedIO type.  Takes an IOStream and creates an IOBuffer.
  If no IOStream is given, a dummy stream is created.  If a dummy stream is
  used, flush must never be called on the returned BufferedIO object

  Inputs:
    f: and IOStream object, defaults to a dummy stream

  Outputs:
    a BufferedIO object
"""->
function BufferedIO(f::IO=DevNull)
  buf = IOBuffer()
  fbuf = BufferedIO{typeof(f)}(f, buf)

  # register atexit hook to make sure any buffered data is flushed before
  # julia exits
  # this causes errors to be printed when julia exists, possible #10431
#  atexit( () -> if isopen(fbuf.fstream)
#                  close(fbuf)
#                end
#        )
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
import Base.write, Base.flush, Base.close

# I think this is all that needs to be implemented, because it is the job of
# functions like println and print to convert things to arrays of UInt8s
write(io::BufferedIO, x::UInt8) = write(io.fbuf, x)


function flush(io::BufferedIO)
  write(io.fstream, takebuf_array(io.fbuf))
  flush(io.fstream)
end

function close(io::BufferedIO)
  flush(io)
  close(io.fstream)
end
