# io.jl: functions to do IO more effiently

export BufferedIO

@doc """
### Utils.BufferedIO

  This type provides a means to buffer IO in memory before writing it to a file.
  Data written to the object is stored in the IOBuffer until flush() is called, 
  when the buffer is (efficiently) dumped into the file
"""->
type BufferedIO  <: IO
  fstream::IOStream
  fbuf::IOBuffer

end

function BufferedIO(f::IOStream)
  fbuf = IOBuffer()
  return BufferedIO(f, fbuf)
end

# only provide write functionality, for now
import Base.print, Base.println, Base.write, Base.flush

#print(io::BufferedIO, x::UInt8) = print(io.fbuf, x)
#print(io::BufferedIO, x, y, args...) = print(io.fbuf, x, y, args...)

#println(io::BufferedIO, x) = println(io.fbuf, x)
#println(io::BufferedIO, x, y, args...) = println(io.fbuf, x, y, args...)

# I think this is all that needs to be implemented, because it is the job of
# functions like println and print to convert things to arrays of UInt8s
write(io::BufferedIO, x::UInt8) = write(io.fbuf, x)


function flush(io::BufferedIO)
  write(io.fstream, takebuf_array(io.fbuf))
  flush(io.fstream)
end


