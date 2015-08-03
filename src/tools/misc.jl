using ArrayViews


@doc """
 ### Tools rmfile

 Removes the file if it exists.  This is only necessary because it is an error
   to remove a file that does not exist.  Does not work on directories

   Arguments:
     fname::AbstractString  :  name of file

"""->
function rmfile(fname::AbstractString)

  if isfile(fname)
    rm(fname)
  end

  return nothing
end



function printbacktrace()
# print the current call stack
  bt = backtrace()
  s = sprint(io->Base.show_backtrace(io, bt))
  println("backtrace: ", s)

  return nothing
end


function checkZeroRows{T <: Number}(A::AbstractArray{T,2}, tol::FloatingPoint)
# checks each row of a matrix for zeros 
# 2d matrices only
# returns the integer number of zero rows, and a Bool
# array telling which rows have all zeros
(m,n) = size(A)

zero_mask = zeros(Bool, m)  # record whether each row has all zeros
zero_row = zeros(Bool, n)  # record results for a row
for i=1:m
  fill!(zero_row, false)
  for j=1:n
    if abs(A[i,j]) < tol  # if zero entry
      zero_row[j] = true
    end
  end  # end loop over column

  rowsum = sum(zero_row)
  zero_mask[i] = (rowsum == n)  # true if all zeros in row
end  # end loop over rows

return sum(zero_mask), zero_mask

end

function checkZeroColumns{T <: Number}(A::AbstractArray{T,2}, tol::FloatingPoint)
# checks each column of a matrix for zeros 
# 2d matrices only
# returns integer number of zero rows, and a Bool
# array telling which rows have all zeros
(m,n) = size(A)

zero_mask = zeros(Bool, n)  # record whether each column has all zeros
zero_row = zeros(Bool, m)  # record results for a column
for i=1:n
  fill!(zero_row, false)
  for j=1:m
    if abs(A[j, i]) < tol  # if zero entry
      zero_row[j] = true
    end
  end  # end loop over row

  rowsum = sum(zero_row)
  zero_mask[i] = (rowsum == n)  # true if all zeros in row
end  # end loop over columns

return sum(zero_mask), zero_mask

end

function checkIdenticalColumns{T <: Number}(A::AbstractArray{T,2}, colnum::Integer, tol::FloatingPoint)
# checks which columns are identical to column number col
# returns number of column identical and an array of bools telling which ones
# does not count column colnum as being the same as itself

(m,n) = size(A)
cnt = 0

col = view(A, :, colnum)
is_same = zeros(Bool, n)

for i=1:n  # loop over columns

  if i == colnum  # skip the specified column
    continue
  end
  col_i = view(A, :, i)
  diff_norm = norm(col - col_i)/m

  if diff_norm < tol
    is_same[i] = true
    cnt += 1
  end
end

return cnt, is_same

end

function findLarge{T <: Number}(A::AbstractArray{T,2}, tol::FloatingPoint)
  (m,n) = size(A)

  cnt = 0
  for i=1:m
    for j=1:n
      if A[i,j] > tol
	println(i, ", ", j)
	cnt += 1
      end
    end
  end

  return cnt
end
